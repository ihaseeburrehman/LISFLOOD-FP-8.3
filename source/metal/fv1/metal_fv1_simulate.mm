//
// LISFLOOD-FP Metal FV1 Solver - Implementation
// Objective-C++ for Metal API access
//

#ifdef __APPLE__

#include "metal_fv1_simulate.h"
#include "../../rain/rain.h"
#include "../../rain/rain.tpp"
#include "../../swe/fv1.h"              // CPU FV1 solver functions
#include "../../swe/fv1/modifiedvars.h" // initialise_Zstar
#include "../../swe/output.h"           // write_mass_stats, write_solution
#include "../../swe/stats.h"            // update_mass_stats, DomainVol
#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include <algorithm>
#include <cmath>

namespace lis {
namespace metal {
namespace fv1 {

// Threadgroup tile used by the flux-sharing sweeps (must match TG_X/TG_Y in
// metal_fv1_kernels.metal).
static const NSUInteger kTileX = 16;
static const NSUInteger kTileY = 16;

// Metal geometry structure (must match shader)
struct MetalGeometry {
  int32_t xsz;
  int32_t ysz;
  int32_t pitch;
  float dx;
  float dy;
  float tly;
  float blx;
  float bly;
};

// Metal physical params (must match shader)
struct MetalPhysicalParams {
  float g;
  float manning;
  float speed_thresh;
};

// Metal solver params (must match shader)
struct MetalSolverParams {
  float dt;
  float depth_thresh;
};

// Implementation class using PIMPL idiom
struct SimulatorImpl {
  // Metal objects
  id<MTLDevice> device;
  id<MTLCommandQueue> commandQueue;
  id<MTLLibrary> library;

  // Compute pipelines
  id<MTLComputePipelineState> applyFrictionPipeline;
  id<MTLComputePipelineState> computeZstarXPipeline;
  id<MTLComputePipelineState> computeZstarYPipeline;
  id<MTLComputePipelineState> updateFlowXPipeline;
  id<MTLComputePipelineState> updateFlowYPipeline;
  id<MTLComputePipelineState> computeDtPipeline;
  id<MTLComputePipelineState> zeroGhostCellsPipeline;
  id<MTLComputePipelineState> applyUniformRainPipeline;
  id<MTLComputePipelineState> applyRainPipeline; // For dynamic rain
  id<MTLComputePipelineState>
      updateGhostCellsPipeline; // Reflective ghost cells
  id<MTLComputePipelineState>
      drainNodataPipeline; // Drain water from NoData cells
  id<MTLComputePipelineState>
      accumulateOutflowPipeline; // Accumulate water leaving domain
  id<MTLComputePipelineState>
      updateMaxFieldPipeline; // Running per-cell maximum depth

  // Buffers for flow variables (double-buffered)
  id<MTLBuffer> H[2];
  id<MTLBuffer> HU[2];
  id<MTLBuffer> HV[2];
  id<MTLBuffer> H_int; // Intermediate state between x and y sweeps
  id<MTLBuffer> HU_int;
  id<MTLBuffer> HV_int;
  id<MTLBuffer> maxH; // Running maximum depth, for the .max output

  // Buffers for grid data
  id<MTLBuffer> DEM;
  id<MTLBuffer> Zstar_x;
  id<MTLBuffer> Zstar_y;
  id<MTLBuffer> manning;
  id<MTLBuffer> dt_field;
  id<MTLBuffer> minDtBuffer; // For atomic reduction

  // Parameter buffers
  id<MTLBuffer> geomBuffer;
  id<MTLBuffer> physicsBuffer;
  id<MTLBuffer> solverBuffer;

  // Boundary flux buffers
  id<MTLBuffer> flux_west;
  id<MTLBuffer> flux_east;
  id<MTLBuffer> flux_north;
  id<MTLBuffer> flux_south;

  // Rain buffers
  id<MTLBuffer> rainDataBuffer;
  id<MTLBuffer> rainGeometryBuffer;

  // Boundary condition type buffer, plus the per-edge-cell prescribed value
  // (stage or unit discharge) refreshed every step for the time-varying types.
  id<MTLBuffer> boundaryTypesBuffer;
  id<MTLBuffer> boundaryValuesBuffer;

  // Simulation state
  int current_buffer; // 0 or 1 for double buffering
  MetalGeometry geom;
  MetalPhysicalParams physics;
  MetalSolverParams solverParams;
  NUMERIC_TYPE time;
  NUMERIC_TYPE cfl;
  NUMERIC_TYPE max_dt;
  NUMERIC_TYPE nodata;
  // The GPU works in single precision, so the DEM is stored relative to the
  // lowest valid elevation and NoData cells are given a pit a fixed distance
  // below that, rather than the raw nodata_elevation (typically -9999).
  // The FV1 scheme only ever uses elevation *differences* (Z*, ETA - Z*,
  // Zdagger differences), so the shift is exact; what it buys is conditioning.
  // At -9999 a float32 ulp is ~9.8e-4 m, i.e. the same size as DepthThresh, so
  // (h + Z) - Z* inside the NoData region lost the depth entirely and the
  // scheme silently leaked mass. See dem_offset/gpu_nodata use in copyToGPU().
  float dem_offset;
  float gpu_nodata;
  bool has_manning_map;
  bool is_first_step;
  int verbose;
  id<MTLBuffer> volOutBuffer; // Buffer for tracking outflow volume
  double last_mass_sync_time = 0.0;
  double vol_out_total;       // Host-side cumulative outflow (m3)
  double last_rain_time;      // Simulation time of last rain upload
  double current_rain_sum;    // Sum of current rain grid for sub-basin (m3/s)
  double
      interval_out_total; // Outflow accumulated since last mass balance check

  size_t buffer_size;
  size_t num_cells;

  // CPU data pointers
  Arrays *arrays;
  BoundCs *boundCs;
  Pars *pars;
  ::Solver *solver;
  States *states;

  // Dynamic Rain object
  ::DynamicRain<std::allocator<NUMERIC_TYPE>> *dynamic_rain;

  SimulatorImpl()
      : device(nil), commandQueue(nil), library(nil), current_buffer(0),
        time(0), cfl(0.7), max_dt(10.0), nodata(-9999), dem_offset(0.0f),
        gpu_nodata(-1000.0f), has_manning_map(false),
        is_first_step(true), verbose(0), volOutBuffer(nil), vol_out_total(0.0), last_rain_time(-1.0),
        current_rain_sum(0.0), interval_out_total(0.0), arrays(nullptr),
        boundCs(nullptr), pars(nullptr), solver(nullptr), states(nullptr),
        dynamic_rain(nullptr) {}

  ~SimulatorImpl() {
    if (dynamic_rain) {
      delete dynamic_rain;
      dynamic_rain = nullptr;
    }
  }

  bool init(Fnames &filenames, Pars &pars, ::Solver &solver, int verb) {
    verbose = verb;

    // Initialize Dynamic Rain
    // Note: DynamicRain constructor reads file header.
    dynamic_rain = new ::DynamicRain<std::allocator<NUMERIC_TYPE>>(
        filenames.dynamicrainfilename, verbose);

    // Get default Metal device
    device = MTLCreateSystemDefaultDevice();
    if (!device) {
      if (verbose)
        printf("Metal: No GPU device found\n");
      return false;
    }

    if (verbose) {
      printf("Metal: Using GPU: %s\n", [[device name] UTF8String]);
      printf("Metal: Unified Memory: %s\n",
             device.hasUnifiedMemory ? "Yes" : "No");
    }

    // Create command queue
    commandQueue = [device newCommandQueue];
    if (!commandQueue) {
      if (verbose)
        printf("Metal: Failed to create command queue\n");
      return false;
    }

    // Load shader library
    NSError *error = nil;

    // 1. Try bundle (portable for apps)
    NSString *libraryPath =
        [[NSBundle mainBundle] pathForResource:@"metal_fv1_kernels"
                                        ofType:@"metallib"];

    // 2. Try current directory
    if (!libraryPath) {
      if ([[NSFileManager defaultManager]
              fileExistsAtPath:@"metal_fv1_kernels.metallib"]) {
        libraryPath = @"metal_fv1_kernels.metallib";
      }
    }

    // 3. Try directory of executable
    if (!libraryPath) {
      NSString *executablePath = [[NSProcessInfo processInfo] arguments][0];
      NSString *executableDir =
          [executablePath stringByDeletingLastPathComponent];
      NSString *path = [executableDir
          stringByAppendingPathComponent:@"metal_fv1_kernels.metallib"];
      if ([[NSFileManager defaultManager] fileExistsAtPath:path]) {
        libraryPath = path;
      }
    }

    if (libraryPath) {
      if (verbose) printf("Metal: Shader library: %s\n", [libraryPath UTF8String]);
      library = [device newLibraryWithFile:libraryPath error:&error];
    }

    // 4. Compile from source if the source is found relative to the executable
    if (!library) {
      NSString *executablePath = [[NSProcessInfo processInfo] arguments][0];
      NSString *executableDir =
          [executablePath stringByDeletingLastPathComponent];
      // Go up from build/ to find metal/fv1/
      NSString *shaderPath =
          [executableDir stringByAppendingPathComponent:
                             @"../metal/fv1/metal_fv1_kernels.metal"];
      NSString *shaderSource =
          [NSString stringWithContentsOfFile:shaderPath
                                    encoding:NSUTF8StringEncoding
                                       error:nil];

      if (shaderSource) {
        MTLCompileOptions *options = [[MTLCompileOptions alloc] init];
        options.fastMathEnabled = YES;
        library = [device newLibraryWithSource:shaderSource
                                       options:options
                                         error:&error];
      }
    }

    if (!library) {
      if (verbose)
        printf("Metal: Failed to load shader library: %s\n",
               [[error localizedDescription] UTF8String]);
      return false;
    }

    // Create compute pipelines
    if (!createPipelines()) {
      return false;
    }

    // Set up geometry
    geom.xsz = pars.xsz;
    geom.ysz = pars.ysz;
    geom.pitch = pars.xsz + 2; // Add ghost cells
    geom.dx = (float)pars.dx;
    geom.dy = (float)pars.dy;
    geom.tly = (float)pars.tly;
    geom.blx = (float)pars.blx;
    geom.bly = (float)pars.bly;

    // Physical params
    physics.g = (float)solver.g;
    physics.manning = (float)pars.FPn;
    physics.speed_thresh = (float)solver.SpeedThresh;

    // Solver params
    solverParams.dt = (float)solver.InitTstep;
    solverParams.depth_thresh = (float)solver.DepthThresh;

    cfl = solver.cfl;
    max_dt = solver.InitTstep;
    nodata = pars.nodata_elevation;

    // Calculate buffer sizes
    num_cells = (geom.xsz + 2) * (geom.ysz + 2);
    buffer_size = num_cells * sizeof(float);

    // Allocate buffers
    if (!allocateBuffers()) {
      return false;
    }

    // Allocate Rain Buffers
    if (dynamic_rain && dynamic_rain->enabled()) {
      const Geometry &rgeom = dynamic_rain->geometry();

      // Rain Data Buffer
      size_t rsize = rgeom.xsz * rgeom.ysz *
                     sizeof(NUMERIC_TYPE); // assuming standard float
      // Note: rain data is NUMERIC_TYPE (likely float or double).
      // If double, we might need conversion or kernels that handle double
      // (Metal supports float primarily, double on some). Ideally convert to
      // float if NUMERIC_TYPE is double. For now assume NUMERIC_TYPE is float
      // (as per user objective "single precision").

      rainDataBuffer =
          [device newBufferWithLength:rsize
                              options:MTLResourceStorageModeManaged];

      // Rain Geometry (MetalGeometry)
      MetalGeometry mRGeom;
      mRGeom.xsz = rgeom.xsz;
      mRGeom.ysz = rgeom.ysz;
      mRGeom.dx = (float)rgeom.dx;
      mRGeom.dy = (float)rgeom.dy;
      mRGeom.tly = (float)rgeom.tly;
      mRGeom.blx = (float)rgeom.blx;
      mRGeom.bly = (float)rgeom.bly;
      // Pitch not really used for rain grid usually linear, but let's set it
      mRGeom.pitch = rgeom.xsz;

      rainGeometryBuffer =
          [device newBufferWithBytes:&mRGeom
                              length:sizeof(MetalGeometry)
                             options:MTLResourceStorageModeShared];

      if (verbose) {
        printf("Metal: Dynamic Rain enabled. Grid: %d x %d\n", rgeom.xsz,
               rgeom.ysz);
      }
    }

    if (verbose) {
      printf("Metal: Grid size: %d x %d\n", geom.xsz, geom.ysz);
      printf("Metal: dx = %.2f m, dy = %.2f m\n", geom.dx, geom.dy);
      printf("Metal: g = %.2f m/s², CFL = %.2f, max_dt = %.2f s\n", physics.g,
             cfl, max_dt);
      printf("Metal: Buffer size: %.2f MB per field\n",
             buffer_size / (1024.0 * 1024.0));
    }

    return true;
  }

  bool createPipelines() {
    NSError *error = nil;

    auto createPipeline = [&](const char *name) -> id<MTLComputePipelineState> {
      id<MTLFunction> func =
          [library newFunctionWithName:[NSString stringWithUTF8String:name]];
      if (!func) {
        if (verbose)
          printf("Metal: Function '%s' not found\n", name);
        return nil;
      }
      id<MTLComputePipelineState> pipeline =
          [device newComputePipelineStateWithFunction:func error:&error];
      if (!pipeline && verbose) {
        printf("Metal: Failed to create pipeline for '%s'\n", name);
      }
      return pipeline;
    };

    applyFrictionPipeline = createPipeline("apply_friction");
    computeZstarXPipeline = createPipeline("compute_zstar_x");
    computeZstarYPipeline = createPipeline("compute_zstar_y");
    updateFlowXPipeline = createPipeline("update_flow_x");
    updateFlowYPipeline = createPipeline("update_flow_y");
    computeDtPipeline = createPipeline("compute_dt_per_cell");
    zeroGhostCellsPipeline = createPipeline("zero_ghost_cells");
    applyUniformRainPipeline = createPipeline("apply_uniform_rain");
    applyRainPipeline = createPipeline("apply_rain");
    updateGhostCellsPipeline = createPipeline("update_ghost_cells");
    drainNodataPipeline = createPipeline("drain_nodata_cells");
    accumulateOutflowPipeline = createPipeline("accumulate_outflow");
    updateMaxFieldPipeline = createPipeline("update_max_field");

    return updateFlowXPipeline && updateFlowYPipeline &&
           accumulateOutflowPipeline;
  }

  bool allocateBuffers() {
    MTLResourceOptions options = MTLResourceStorageModeShared; // Unified memory

    // Flow variables (double buffered)
    for (int i = 0; i < 2; i++) {
      H[i] = [device newBufferWithLength:buffer_size options:options];
      HU[i] = [device newBufferWithLength:buffer_size options:options];
      HV[i] = [device newBufferWithLength:buffer_size options:options];
      if (!H[i] || !HU[i] || !HV[i])
        return false;
    }

    // Intermediate buffers
    H_int = [device newBufferWithLength:buffer_size options:options];
    HU_int = [device newBufferWithLength:buffer_size options:options];
    HV_int = [device newBufferWithLength:buffer_size options:options];
    maxH = [device newBufferWithLength:buffer_size options:options];
    if (!maxH)
      return false;

    // Grid data
    DEM = [device newBufferWithLength:buffer_size options:options];
    Zstar_x = [device newBufferWithLength:buffer_size options:options];
    Zstar_y = [device newBufferWithLength:buffer_size options:options];
    dt_field = [device newBufferWithLength:buffer_size options:options];
    minDtBuffer = [device newBufferWithLength:sizeof(uint32_t)
                                      options:MTLResourceStorageModeShared];
    manning = [device newBufferWithLength:buffer_size options:options];

    // Zero out all buffers to ensure clean state
    for (int i = 0; i < 2; i++) {
      memset([H[i] contents], 0, buffer_size);
      memset([HU[i] contents], 0, buffer_size);
      memset([HV[i] contents], 0, buffer_size);
    }
    memset([H_int contents], 0, buffer_size);
    memset([HU_int contents], 0, buffer_size);
    memset([HV_int contents], 0, buffer_size);
    memset([DEM contents], 0, buffer_size);
    memset([Zstar_x contents], 0, buffer_size);
    memset([Zstar_y contents], 0, buffer_size);
    memset([dt_field contents], 0, buffer_size);
    memset([minDtBuffer contents], 0, sizeof(uint32_t));
    memset([manning contents], 0, buffer_size);

    // Parameter buffers
    geomBuffer = [device newBufferWithBytes:&geom
                                     length:sizeof(MetalGeometry)
                                    options:options];
    physicsBuffer = [device newBufferWithBytes:&physics
                                        length:sizeof(MetalPhysicalParams)
                                       options:options];
    solverBuffer = [device newBufferWithBytes:&solverParams
                                       length:sizeof(MetalSolverParams)
                                      options:options];

    // Allocate flux buffers
    flux_west = [device newBufferWithLength:geom.ysz * sizeof(float) * 3
                                    options:options];
    flux_east = [device newBufferWithLength:geom.ysz * sizeof(float) * 3
                                    options:options];
    flux_north = [device newBufferWithLength:geom.xsz * sizeof(float) * 3
                                     options:options];
    flux_south = [device newBufferWithLength:geom.xsz * sizeof(float) * 3
                                     options:options];

    // Allocate boundary condition types buffer
    int num_bc_cells = 2 * geom.xsz + 2 * geom.ysz;
    boundaryTypesBuffer =
        [device newBufferWithLength:num_bc_cells * sizeof(int32_t)
                            options:options];
    boundaryValuesBuffer =
        [device newBufferWithLength:num_bc_cells * sizeof(float)
                            options:options];

    // Allocate tracking buffers
    volOutBuffer = [device newBufferWithLength:8 * sizeof(float)
                                       options:options];

    memset([volOutBuffer contents], 0, 8 * sizeof(float));

    return DEM && Zstar_x && Zstar_y && geomBuffer && physicsBuffer &&
           solverBuffer && flux_west && flux_east && flux_north && flux_south &&
           boundaryTypesBuffer && volOutBuffer;
  }

  void copyToGPU(NUMERIC_TYPE *h_H, NUMERIC_TYPE *h_HU, NUMERIC_TYPE *h_HV,
                 NUMERIC_TYPE *h_DEM, NUMERIC_TYPE *h_Mann) {
    float *d_H = (float *)[H[current_buffer] contents];
    float *d_HU = (float *)[HU[current_buffer] contents];
    float *d_HV = (float *)[HV[current_buffer] contents];
    float *d_DEM = (float *)[DEM contents];
    float *d_Mann = (float *)[manning contents];

    // Zero all flow buffers first (including both double buffers)
    for (int i = 0; i < 2; i++) {
      memset([H[i] contents], 0, buffer_size);
      memset([HU[i] contents], 0, buffer_size);
      memset([HV[i] contents], 0, buffer_size);
    }
    memset([DEM contents], 0, buffer_size);
    memset([manning contents], 0, buffer_size);
    memset([maxH contents], 0, buffer_size);

    // Copy CPU data to the "active" part of the GPU buffers (1-indexed)
    // Handle NoData values properly
    float nodata_val = (float)nodata;

    // Pick the vertical datum for the GPU copy of the DEM. The FV1 scheme only
    // ever uses elevation *differences* (Z*, ETA - Z*, Zdagger differences), so
    // rebasing is exact; what it buys is float32 conditioning. At the common
    // nodata_elevation of -9999 a float32 ulp is ~9.8e-4 m, i.e. the same size
    // as DepthThresh, which wrecks (h + Z) - Z* inside the NoData region.
    //
    // The replacement value must keep the *sense* of nodata_elevation: LISFLOOD
    // defaults it to DEM_NO_DATA = +1e7, an impermeable wall, while clipped
    // catchments usually set it to -9999, a pit that drain_nodata empties.
    // Mapping both onto a pit would let water flood into cells the CPU treats
    // as solid, so the sign is preserved here.
    {
      double min_valid = 0.0, max_valid = 0.0;
      bool any_valid = false;
      for (int j = 0; j < geom.ysz; j++) {
        for (int i = 0; i < geom.xsz; i++) {
          double z = (double)h_DEM[j * geom.xsz + i];
          if (fabs(z - (double)nodata) < 1.0)
            continue;
          if (!any_valid) {
            min_valid = max_valid = z;
            any_valid = true;
          } else {
            if (z < min_valid) min_valid = z;
            if (z > max_valid) max_valid = z;
          }
        }
      }
      dem_offset = any_valid ? (float)min_valid : 0.0f;

      const float kMargin = 1000.0f; // far beyond any credible water depth
      if (!any_valid) {
        gpu_nodata = (float)nodata;
      } else if ((double)nodata > max_valid) {
        gpu_nodata = (float)(max_valid - min_valid) + kMargin;  // wall
      } else if ((double)nodata < min_valid) {
        gpu_nodata = -kMargin;                                   // pit
      } else {
        gpu_nodata = (float)((double)nodata - min_valid);        // as given
      }
      if (verbose)
        printf("Metal: DEM datum shifted by %.3f m; NoData (%.1f) mapped to "
               "%.1f m [%s]\n",
               dem_offset, (double)nodata, gpu_nodata,
               ((double)nodata > max_valid) ? "wall" : "pit");
    }

    for (int j = 0; j < geom.ysz; j++) {
      for (int i = 0; i < geom.xsz; i++) {
        int cpu_idx = j * geom.xsz + i;
        int gpu_idx = (j + 1) * geom.pitch + (i + 1);

        float dem = (float)h_DEM[cpu_idx];

        // Check if this is a NoData cell
        bool is_nodata = (fabs(dem - nodata_val) < 1.0f);

        if (is_nodata) {
          // NoData cell: set DEM to the GPU pit value, zero flow
          d_DEM[gpu_idx] = gpu_nodata;
          d_H[gpu_idx] = 0.0f;
          d_HU[gpu_idx] = 0.0f;
          d_HV[gpu_idx] = 0.0f;
          // Roughness still has to match the CPU here. NoData cells all sit at
          // nodata_elevation, so between two of them Z* is flat and water that
          // spills in spreads across them freely. LoadManningsn() gives those
          // cells Parptr->FPn, so zeroing n made the whole NoData region
          // frictionless and let water run far deeper into it than on the CPU.
          d_Mann[gpu_idx] =
              h_Mann ? (float)h_Mann[cpu_idx] : physics.manning;
        } else {
          // Valid data cell
          float h = (float)h_H[cpu_idx];

          d_H[gpu_idx] = h;
          d_HU[gpu_idx] = (float)h_HU[cpu_idx];
          d_HV[gpu_idx] = -(float)h_HV[cpu_idx];
          d_DEM[gpu_idx] = dem - dem_offset;

          if (h_Mann) {
            d_Mann[gpu_idx] = (float)h_Mann[cpu_idx];
          } else {
            d_Mann[gpu_idx] = physics.manning; // Use uniform value
          }
        }
      }
    }

    has_manning_map = (h_Mann != nullptr);

    // Fill DEM ghost cells (zero-gradient boundary conditions)
    // This is critical to avoid artificial "cliffs" at the boundary
    for (int j = 1; j <= geom.ysz; j++) {
      d_DEM[j * geom.pitch + 0] = d_DEM[j * geom.pitch + 1]; // West
      d_DEM[j * geom.pitch + geom.xsz + 1] =
          d_DEM[j * geom.pitch + geom.xsz]; // East
    }
    for (int i = 0; i < geom.pitch; i++) {
      d_DEM[0 * geom.pitch + i] = d_DEM[1 * geom.pitch + i]; // North
      d_DEM[(geom.ysz + 1) * geom.pitch + i] =
          d_DEM[geom.ysz * geom.pitch + i]; // South
    }

    // Populate boundary condition types buffer
    if (boundCs && boundCs->BC_Ident) {
      int32_t *d_bc = (int32_t *)[boundaryTypesBuffer contents];
      int num_bc_cells = 2 * geom.xsz + 2 * geom.ysz;
      // Pass the real BC_Ident through. Collapsing everything that was not
      // FREE into a reflective wall silently discarded QFIX/QVAR/HFIX/HVAR
      // boundaries, so prescribed inflows never entered the domain.
      for (int i = 0; i < num_bc_cells; i++) {
        d_bc[i] = (int32_t)boundCs->BC_Ident[i];
      }
    }

    // Compute Z-star fields (bed levels at interfaces)
    @autoreleasepool {
      id<MTLCommandBuffer> commandBuffer = [commandQueue commandBuffer];
      id<MTLComputeCommandEncoder> encoder =
          [commandBuffer computeCommandEncoder];

      MTLSize gridSize = MTLSizeMake(geom.xsz + 1, geom.ysz, 1);
      MTLSize threadGroupSize = MTLSizeMake(16, 16, 1);

      // Compute Zstar_x
      if (computeZstarXPipeline) {
        [encoder setComputePipelineState:computeZstarXPipeline];
        [encoder setBuffer:DEM offset:0 atIndex:0];
        [encoder setBuffer:Zstar_x offset:0 atIndex:1];
        [encoder setBuffer:geomBuffer offset:0 atIndex:2];
        [encoder dispatchThreads:gridSize
            threadsPerThreadgroup:threadGroupSize];
      }

      // Compute Zstar_y
      gridSize = MTLSizeMake(geom.xsz, geom.ysz + 1, 1);
      if (computeZstarYPipeline) {
        [encoder setComputePipelineState:computeZstarYPipeline];
        [encoder setBuffer:DEM offset:0 atIndex:0];
        [encoder setBuffer:Zstar_y offset:0 atIndex:1];
        [encoder setBuffer:geomBuffer offset:0 atIndex:2];
        [encoder dispatchThreads:gridSize
            threadsPerThreadgroup:threadGroupSize];
      }

      // Set initial ghost cells for H, HU, HV (reflective BCs)
      if (updateGhostCellsPipeline) {
        [encoder setComputePipelineState:updateGhostCellsPipeline];
        [encoder setBuffer:H[current_buffer] offset:0 atIndex:0];
        [encoder setBuffer:HU[current_buffer] offset:0 atIndex:1];
        [encoder setBuffer:HV[current_buffer] offset:0 atIndex:2];
        [encoder setBuffer:geomBuffer offset:0 atIndex:3];
        [encoder setBuffer:boundaryTypesBuffer offset:0 atIndex:4];

        int maxBoundarySize = std::max(geom.pitch, geom.ysz + 2);
        MTLSize ghostGridSize = MTLSizeMake(maxBoundarySize, 1, 1);
        MTLSize ghostThreadGroupSize = MTLSizeMake(64, 1, 1);
        [encoder dispatchThreads:ghostGridSize
            threadsPerThreadgroup:ghostThreadGroupSize];
      }

      [encoder endEncoding];
      [commandBuffer commit];
      [commandBuffer waitUntilCompleted];
    }

    if (verbose) {
      printf("Metal: Initialized Z-star fields\n");
    }
  }

  void syncMassTotals() {
    @autoreleasepool {
      id<MTLCommandBuffer> commandBuffer = [commandQueue commandBuffer];
      [commandBuffer commit];
      [commandBuffer waitUntilCompleted];

      float *totals = (float *)[volOutBuffer contents];
      const double out_boundary = (double)totals[0];
      const double out_drain = (double)totals[2];
      const double out = out_boundary + out_drain;
      const double in = totals[1];
      const double uniform_rain = totals[3];
      const double elapsed = (double)time - last_mass_sync_time;

      if (verbose && elapsed > 0.0) {
        printf("Metal mass: t=%.1f  boundary_out=%.4e  nodata_drain=%.4e  "
               "boundary_in=%.4e\n",
               (double)time, out_boundary, out_drain, in);
        fflush(stdout);
      }
      if (boundCs) {
        boundCs->VolOutMT += (NUMERIC_TYPE)out;
        boundCs->VolInMT += (NUMERIC_TYPE)(in + uniform_rain);
        pars->RainTotalLoss += (NUMERIC_TYPE)uniform_rain;
        if (elapsed > 0.0) {
          boundCs->Qout += (NUMERIC_TYPE)(out / elapsed);
          boundCs->Qin += (NUMERIC_TYPE)(in / elapsed);
        }
      }
      last_mass_sync_time = (double)time;

      // Zero the GPU counters for the next interval
      memset([volOutBuffer contents], 0, 8 * sizeof(float));
    }
  }

  void copyFromGPU(NUMERIC_TYPE *h_H, NUMERIC_TYPE *h_HU, NUMERIC_TYPE *h_HV) {
    // Sync mass balance first to ensure the final numbers are captured
    syncMassTotals();

    float *d_H = (float *)[H[current_buffer] contents];
    float *d_HU = (float *)[HU[current_buffer] contents];
    float *d_HV = (float *)[HV[current_buffer] contents];

    for (int j = 0; j < geom.ysz; j++) {
      for (int i = 0; i < geom.xsz; i++) {
        int cpu_idx = j * geom.xsz + i;
        int gpu_idx = (j + 1) * geom.pitch + (i + 1);

        h_H[cpu_idx] = (NUMERIC_TYPE)d_H[gpu_idx];
        h_HU[cpu_idx] = (NUMERIC_TYPE)d_HU[gpu_idx];
        h_HV[cpu_idx] = -(NUMERIC_TYPE)d_HV[gpu_idx];
      }
    }
  }

  void copyMaxFromGPU(NUMERIC_TYPE *h_maxH) {
    if (h_maxH == nullptr)
      return;

    @autoreleasepool {
      id<MTLCommandBuffer> commandBuffer = [commandQueue commandBuffer];
      [commandBuffer commit];
      [commandBuffer waitUntilCompleted];
    }

    const float *d_maxH = (const float *)[maxH contents];
    for (int j = 0; j < geom.ysz; j++) {
      for (int i = 0; i < geom.xsz; i++) {
        h_maxH[j * geom.xsz + i] =
            (NUMERIC_TYPE)d_maxH[(j + 1) * geom.pitch + (i + 1)];
      }
    }
  }

  void uploadState(NUMERIC_TYPE *h_H, NUMERIC_TYPE *h_HU, NUMERIC_TYPE *h_HV) {
    float *d_H = (float *)[H[current_buffer] contents];
    float *d_HU = (float *)[HU[current_buffer] contents];
    float *d_HV = (float *)[HV[current_buffer] contents];

    for (int j = 0; j < geom.ysz; j++) {
      for (int i = 0; i < geom.xsz; i++) {
        int cpu_idx = j * geom.xsz + i;
        int gpu_idx = (j + 1) * geom.pitch + (i + 1);

        d_H[gpu_idx] = (float)h_H[cpu_idx];
        d_HU[gpu_idx] = (float)h_HU[cpu_idx];
        d_HV[gpu_idx] = -(float)h_HV[cpu_idx];
      }
    }
  }

  void updateBoundaryFluxes() {
    // Extract fluxes from arrays->Qx, Qy (calculated by BCs on CPU)
    // Convert Total Q (m3/s) to Unit Flux (m2/s) and Momentum columns

    float *d_west = (float *)[flux_west contents];
    float *d_east = (float *)[flux_east contents];
    float *d_north = (float *)[flux_north contents];
    float *d_south = (float *)[flux_south contents];

    int xsz = geom.xsz;
    int ysz = geom.ysz;
    float dx = geom.dx;
    float dy = geom.dy;
    float g = physics.g;

    int pitch_x = xsz + 1; // Qx width (xsz+1)

    // West (i=0) and East (i=xsz) Interfaces
    for (int j = 0; j < ysz; j++) {
      // West
      int idx_west = j * pitch_x + 0;
      float Q_west = (float)arrays->Qx[idx_west];

      float h_west = (float)arrays->H[j * xsz + 0]; // Cell 0
      // BCs code for x-direction: *qptr = sign*...*Parptr->dy
      // So Qx is integrated over dy. So unit_q = Qx / dy.
      float unit_q_west = Q_west / dy;

      d_west[3 * j + 0] = unit_q_west; // F.H
      if (h_west > 0.001f) {
        // F.HU = u * F.H + 0.5*g*h^2 = (q/h)*q + ... = q^2/h + ...
        d_west[3 * j + 1] =
            (unit_q_west * unit_q_west / h_west) + 0.5f * g * h_west * h_west;
      } else {
        d_west[3 * j + 1] = 0.0f;
      }
      d_west[3 * j + 2] = 0.0f; // F.HV (approx 0 shear)

      // East
      int idx_east = j * pitch_x + xsz;
      float Q_east = (float)arrays->Qx[idx_east];
      float h_east = (float)arrays->H[j * xsz + (xsz - 1)]; // Last cell
      float unit_q_east = Q_east / dy;

      d_east[3 * j + 0] = unit_q_east;
      if (h_east > 0.001f) {
        d_east[3 * j + 1] =
            (unit_q_east * unit_q_east / h_east) + 0.5f * g * h_east * h_east;
      } else {
        d_east[3 * j + 1] = 0.0f;
      }
      d_east[3 * j + 2] = 0.0f;
    }

    // North (j=0) and South (j=ysz) Interfaces
    for (int i = 0; i < xsz; i++) {
      // North
      // Qy index: j * xsz + i. North i j=0. (Note: Qy is (xsz)*(ysz+1))
      // So index is just i.
      int idx_north = 0 * xsz + i;
      float Q_north = (float)arrays->Qy[idx_north];
      float h_north = (float)arrays->H[0 * xsz + i]; // Row 0

      // BCs for Qy multiplies by dx. So unit_q = Qy / dx.
      float unit_q_north = Q_north / dx;

      d_north[3 * i + 0] = unit_q_north; // F.H
      d_north[3 * i + 1] = 0.0f;         // F.HU (approx 0 shear)
      if (h_north > 0.001f) {
        // F.HV = v * F.H + ...
        d_north[3 * i + 2] = (unit_q_north * unit_q_north / h_north) +
                             0.5f * g * h_north * h_north;
      } else {
        d_north[3 * i + 2] = 0.0f;
      }

      // South
      int idx_south = ysz * xsz + i;
      float Q_south = (float)arrays->Qy[idx_south];
      float h_south = (float)arrays->H[(ysz - 1) * xsz + i]; // Last Row
      float unit_q_south = Q_south / dx;

      d_south[3 * i + 0] = unit_q_south;
      d_south[3 * i + 1] = 0.0f;
      if (h_south > 0.001f) {
        d_south[3 * i + 2] = (unit_q_south * unit_q_south / h_south) +
                             0.5f * g * h_south * h_south;
      } else {
        d_south[3 * i + 2] = 0.0f;
      }
    }
  }

  // Refresh the prescribed boundary values for this timestep. Stage boundaries
  // are shifted into the GPU's vertical datum so that (value - Z) is evaluated
  // consistently with the shifted DEM; discharge boundaries are in m2/s and
  // need no shift.
  void updateBoundaryValues() {
    if (!boundCs || !boundCs->BC_Ident || !boundaryValuesBuffer)
      return;
    float *d_val = (float *)[boundaryValuesBuffer contents];
    const int num_bc_cells = 2 * geom.xsz + 2 * geom.ysz;
    for (int i = 0; i < num_bc_cells; i++) {
      switch (boundCs->BC_Ident[i]) {
      case HFIX2:
        d_val[i] = (float)boundCs->BC_Val[i] - dem_offset;
        break;
      case HVAR3:
        d_val[i] = (float)InterpolateTimeSeries(boundCs->BC_TimeSeries[i], time) -
                   dem_offset;
        break;
      case QFIX4:
        d_val[i] = (float)boundCs->BC_Val[i];
        break;
      case QVAR5:
        d_val[i] =
            (float)InterpolateTimeSeries(boundCs->BC_TimeSeries[i], time);
        break;
      default:
        d_val[i] = 0.0f;
        break;
      }
    }
  }

  void applyPointSources() {
    if (!boundCs || boundCs->numPS <= 0)
      return;

    int xsz = geom.xsz;
    int ysz = geom.ysz;
    float dx = geom.dx;
    float dy = geom.dy;
    float dt = solverParams.dt;
    float cell_area = dx * dy;

    // Get direct pointer to GPU H buffer (unified memory)
    float *d_H = (float *)[H[current_buffer] contents];

    for (int i = 0; i < boundCs->numPS; i++) {
      int px = boundCs->xpi[i];
      int py = boundCs->ypi[i];

      if (px < 0 || px >= xsz || py < 0 || py >= ysz)
        continue;

      float val = 0.0f;
      bool is_q_type = false;
      if (boundCs->PS_Ident[i] == QVAR5) {
        val = (float)InterpolateTimeSeries(boundCs->PS_TimeSeries[i], time);
        is_q_type = true;
      } else if (boundCs->PS_Ident[i] == QFIX4) {
        val = (float)boundCs->PS_Val[i];
        is_q_type = true;
      } else if (boundCs->PS_Ident[i] == HVAR3) {
        float h_target =
            (float)InterpolateTimeSeries(boundCs->PS_TimeSeries[i], time);
        float z = arrays->DEM[py * xsz + px];
        float h_new = std::max(0.0f, h_target - z);
        float h_old = d_H[(py + 1) * geom.pitch + (px + 1)];
        val = (h_new - h_old) * cell_area / dt; // Calculate Q to match H
        is_q_type = false;                      // Already scaled to m3/s
      } else if (boundCs->PS_Ident[i] == HFIX2) {
        float h_target = (float)boundCs->PS_Val[i];
        float z = arrays->DEM[py * xsz + px];
        float h_new = std::max(0.0f, h_target - z);
        float h_old = d_H[(py + 1) * geom.pitch + (px + 1)];
        val = (h_new - h_old) * cell_area / dt; // Calculate Q to match H
        is_q_type = false;                      // Already scaled to m3/s
      }

      if (std::abs(val) > 1e-9f) {
        float Q_total = is_q_type ? (val * dx) : val;

        // GPU index (includes ghost cells)
        int gpu_idx = (py + 1) * geom.pitch + (px + 1);

        // Add volume: dH = Q_total * dt / Area
        float dH = (Q_total * dt) / cell_area;

        // Update GPU buffer directly (Unified memory)
        d_H[gpu_idx] += dH;

        // Update stats for mass balance (Qin is discharge m3/s)
        if (Q_total > 0) {
          boundCs->Qin += Q_total;
          boundCs->VolInMT += (NUMERIC_TYPE)(Q_total * dt);
        } else {
          boundCs->Qout -= Q_total;
          boundCs->VolOutMT -= (NUMERIC_TYPE)(Q_total * dt);
        }
      }
    }
  }

  NUMERIC_TYPE computeTimestep() {
    // Run compute_dt kernel
    @autoreleasepool {
      id<MTLCommandBuffer> commandBuffer = [commandQueue commandBuffer];
      id<MTLComputeCommandEncoder> encoder =
          [commandBuffer computeCommandEncoder];

      if (computeDtPipeline) {
        // Reset minDtBuffer to max value before kernel
        float max_dt_val = (float)max_dt;
        uint32_t init_uint;
        memcpy(&init_uint, &max_dt_val, sizeof(uint32_t));
        memcpy([minDtBuffer contents], &init_uint, sizeof(uint32_t));

        [encoder setComputePipelineState:computeDtPipeline];
        [encoder setBuffer:H[current_buffer] offset:0 atIndex:0];
        [encoder setBuffer:HU[current_buffer] offset:0 atIndex:1];
        [encoder setBuffer:HV[current_buffer] offset:0 atIndex:2];

        // Match kernel indices: 3:dt_field, 4:geom, 5:physics, 6:cfl, 7:max_dt
        [encoder setBuffer:dt_field offset:0 atIndex:3];
        [encoder setBuffer:geomBuffer offset:0 atIndex:4];
        [encoder setBuffer:physicsBuffer offset:0 atIndex:5];

        float cfl_val = (float)solver->cfl;
        [encoder setBytes:&cfl_val length:sizeof(float) atIndex:6];

        float max_dt_f = (float)max_dt;
        [encoder setBytes:&max_dt_f length:sizeof(float) atIndex:7];

        [encoder setBytes:&max_dt_f length:sizeof(float) atIndex:7];
        [encoder setBuffer:minDtBuffer
                    offset:0
                   atIndex:8]; // Bind minDt atomic buffer

        MTLSize gridSize = MTLSizeMake(geom.xsz, geom.ysz, 1);
        MTLSize threadGroupSize =
            MTLSizeMake(16, 16, 1); // Optimized for M1 Max
        [encoder dispatchThreads:gridSize
            threadsPerThreadgroup:threadGroupSize];
      }
      [encoder endEncoding];
      [commandBuffer commit];
      [commandBuffer waitUntilCompleted];
    }

    // Read back reduction result from global buffer (Instant O(1))
    uint32_t min_dt_uint;
    memcpy(&min_dt_uint, [minDtBuffer contents], sizeof(uint32_t));

    // Convert back to float
    // Note: If no thread updated it (e.g. all dry), it remains max_dt logic
    float min_dt = *(float *)&min_dt_uint;

    // Sanity limit
    if (min_dt <= 0.0f)
      min_dt = (float)max_dt;

    return (NUMERIC_TYPE)min_dt;
  }

  NUMERIC_TYPE step() {
    // 0. Compute Adaptive Timestep
    NUMERIC_TYPE adaptive_dt = computeTimestep();

    if (is_first_step) {
      adaptive_dt = states->adaptive_ts == ON ? (NUMERIC_TYPE)1e-6 : max_dt;
      is_first_step = false;
    }

    solverParams.dt = (float)adaptive_dt;

    // CRITICAL: Update solver->Tstep so BCs() uses the correct timestep
    solver->Tstep = adaptive_dt;

    // Update solver buffer
    memcpy([solverBuffer contents], &solverParams, sizeof(MetalSolverParams));

    int src = current_buffer;
    int dst = 1 - current_buffer;

    updateBoundaryValues();

    // --- APPLY POINT SOURCES DIRECTLY TO GPU ---
    // Point sources add water directly to current GPU buffer (unified memory)
    applyPointSources();

    // LoadRain has already converted the time series from mm/hour to m/s.
    // Match the legacy Rainfall routine, including depth-limited negative rain.
    if (states->rainfall == ON) {
      float rate = (float)InterpolateTimeSeries(arrays->rain, time);
      @autoreleasepool {
        id<MTLCommandBuffer> commandBuffer = [commandQueue commandBuffer];
        id<MTLComputeCommandEncoder> encoder = [commandBuffer computeCommandEncoder];
        [encoder setComputePipelineState:applyUniformRainPipeline];
        [encoder setBuffer:H[src] offset:0 atIndex:0];
        [encoder setBuffer:DEM offset:0 atIndex:1];
        [encoder setBuffer:geomBuffer offset:0 atIndex:2];
        [encoder setBytes:&rate length:sizeof(float) atIndex:3];
        float dt = solverParams.dt;
        [encoder setBytes:&dt length:sizeof(float) atIndex:4];
        [encoder setBytes:&gpu_nodata length:sizeof(float) atIndex:5];
        [encoder setBuffer:volOutBuffer offset:0 atIndex:6];
        [encoder dispatchThreads:MTLSizeMake(geom.xsz, geom.ysz, 1)
            threadsPerThreadgroup:MTLSizeMake(16, 16, 1)];
        [encoder endEncoding];
        [commandBuffer commit];
      }
    }

    // --- APPLY DYNAMIC RAIN ON GPU ---
    if (dynamic_rain && dynamic_rain->enabled() && applyRainPipeline) {
      dynamic_rain->update_time(time);

      if (time >= last_rain_time + 1.0 || last_rain_time < 0) {
        NUMERIC_TYPE *rdata = dynamic_rain->data();
        size_t num_rain_cells =
            dynamic_rain->geometry().xsz * dynamic_rain->geometry().ysz;
        float *d_rain = (float *)[rainDataBuffer contents];

        for (size_t k = 0; k < num_rain_cells; ++k)
          d_rain[k] = (float)rdata[k];
        const MetalGeometry &rg = *(const MetalGeometry *)[rainGeometryBuffer contents];
        double sum = 0.0;
        for (int j = 0; j < geom.ysz; ++j) {
          for (int i = 0; i < geom.xsz; ++i) {
            if (std::fabs(arrays->DEM[j * geom.xsz + i] - nodata) < 1e-6) continue;
            int ri = (int)std::floor((geom.blx - rg.blx + i * geom.dx) / rg.dx);
            int rj = (int)std::floor((j * geom.dy - (geom.tly - rg.tly)) / rg.dy);
            if (ri >= 0 && ri < rg.xsz && rj >= 0 && rj < rg.ysz)
              sum += std::max(0.0f, d_rain[rj * rg.xsz + ri]);
          }
        }
        current_rain_sum = sum;
        last_rain_time = time;

        [rainDataBuffer
            didModifyRange:NSMakeRange(0, num_rain_cells * sizeof(float))];
      }

      // Update mass balance (CPU side is much faster for global sums)
      double interval_rain_vol =
          current_rain_sum * (double)pars->dA * (double)adaptive_dt;
      boundCs->VolInMT += (NUMERIC_TYPE)interval_rain_vol;
      pars->RainTotalLoss += (NUMERIC_TYPE)interval_rain_vol;

      @autoreleasepool {
        id<MTLCommandBuffer> commandBuffer = [commandQueue commandBuffer];
        id<MTLComputeCommandEncoder> encoder =
            [commandBuffer computeCommandEncoder];
        [encoder setComputePipelineState:applyRainPipeline];
        [encoder setBuffer:H[src] offset:0 atIndex:0];
        [encoder setBuffer:DEM offset:0 atIndex:1];
        [encoder setBuffer:rainDataBuffer offset:0 atIndex:2];
        [encoder setBuffer:geomBuffer offset:0 atIndex:3];
        [encoder setBuffer:rainGeometryBuffer offset:0 atIndex:4];
        float dt_f = (float)adaptive_dt;
        [encoder setBytes:&dt_f length:sizeof(float) atIndex:5];
        float nd = gpu_nodata;
        [encoder setBytes:&nd length:sizeof(float) atIndex:6];

        [encoder dispatchThreads:MTLSizeMake(geom.xsz, geom.ysz, 1)
            threadsPerThreadgroup:MTLSizeMake(16, 16, 1)];
        [encoder endEncoding];
        [commandBuffer commit];
        // PERFORMANCE: No waitUntilCompleted here.
      }
    }

    // --- GPU KERNELS ---

    @autoreleasepool {
      id<MTLCommandBuffer> commandBuffer = [commandQueue commandBuffer];
      id<MTLComputeCommandEncoder> encoder =
          [commandBuffer computeCommandEncoder];

      MTLSize gridSize = MTLSizeMake(geom.xsz, geom.ysz, 1);
      MTLSize threadGroupSize = MTLSizeMake(16, 16, 1); // Optimized for M1 Max

      // Apply friction
      if (applyFrictionPipeline && (has_manning_map || physics.manning > 0)) {
        [encoder setComputePipelineState:applyFrictionPipeline];
        [encoder setBuffer:H[src] offset:0 atIndex:0];
        [encoder setBuffer:HU[src] offset:0 atIndex:1];
        [encoder setBuffer:HV[src] offset:0 atIndex:2];
        [encoder setBuffer:geomBuffer offset:0 atIndex:3];
        [encoder setBuffer:physicsBuffer offset:0 atIndex:4];
        [encoder setBuffer:solverBuffer offset:0 atIndex:5];
        [encoder setBuffer:manning offset:0 atIndex:6];
        [encoder dispatchThreads:gridSize
            threadsPerThreadgroup:threadGroupSize];
      }

      // Update ghost cells BEFORE flow calculations
      if (updateGhostCellsPipeline) {
        [encoder setComputePipelineState:updateGhostCellsPipeline];
        [encoder setBuffer:H[src] offset:0 atIndex:0];
        [encoder setBuffer:HU[src] offset:0 atIndex:1];
        [encoder setBuffer:HV[src] offset:0 atIndex:2];
        [encoder setBuffer:geomBuffer offset:0 atIndex:3];
        [encoder setBuffer:boundaryTypesBuffer offset:0 atIndex:4];

        // Dispatch enough threads for max(pitch, ysz+2)
        int maxBoundarySize = std::max(geom.pitch, geom.ysz + 2);
        MTLSize ghostGridSize = MTLSizeMake(maxBoundarySize, 1, 1);
        MTLSize ghostThreadGroupSize = MTLSizeMake(64, 1, 1);
        [encoder dispatchThreads:ghostGridSize
            threadsPerThreadgroup:ghostThreadGroupSize];
      }

      // X-sweep
      // NOTE: Pass nil for boundary flux buffers to use HLL wall BCs
      [encoder setComputePipelineState:updateFlowXPipeline];
      [encoder setBuffer:H[src] offset:0 atIndex:0];
      [encoder setBuffer:HU[src] offset:0 atIndex:1];
      [encoder setBuffer:HV[src] offset:0 atIndex:2];
      [encoder setBuffer:H_int offset:0 atIndex:3];
      [encoder setBuffer:HU_int offset:0 atIndex:4];
      [encoder setBuffer:HV_int offset:0 atIndex:5];
      [encoder setBuffer:DEM offset:0 atIndex:6];
      [encoder setBuffer:Zstar_x offset:0 atIndex:7];
      [encoder setBuffer:geomBuffer offset:0 atIndex:8];
      [encoder setBuffer:physicsBuffer offset:0 atIndex:9];
      [encoder setBuffer:solverBuffer offset:0 atIndex:10];
      [encoder setBuffer:boundaryTypesBuffer offset:0 atIndex:11];
      [encoder setBuffer:boundaryValuesBuffer offset:0 atIndex:12];
      [encoder setBuffer:volOutBuffer offset:0 atIndex:13];

      // Overlapping threadgroups: tiles share interface fluxes through
      // threadgroup memory, so they must overlap by one column in x (and one
      // row in y below), matching the CUDA kernel's block layout.
      {
        MTLSize tg = MTLSizeMake(kTileX, kTileY, 1);
        MTLSize groups = MTLSizeMake(
            (geom.xsz + 1 + (kTileX - 1) - 1) / (kTileX - 1),
            (geom.ysz + kTileY - 1) / kTileY, 1);
        [encoder dispatchThreadgroups:groups threadsPerThreadgroup:tg];
      }

      // Y-sweep
      // NOTE: Pass nil for boundary flux buffers to use HLL wall BCs
      [encoder setComputePipelineState:updateFlowYPipeline];
      [encoder setBuffer:H_int offset:0 atIndex:0];
      [encoder setBuffer:HU_int offset:0 atIndex:1];
      [encoder setBuffer:HV_int offset:0 atIndex:2];
      [encoder setBuffer:H[src] offset:0 atIndex:3];
      [encoder setBuffer:HU[src] offset:0 atIndex:4];
      [encoder setBuffer:HV[src] offset:0 atIndex:5];
      [encoder setBuffer:H[dst] offset:0 atIndex:6];
      [encoder setBuffer:HU[dst] offset:0 atIndex:7];
      [encoder setBuffer:HV[dst] offset:0 atIndex:8];
      [encoder setBuffer:DEM offset:0 atIndex:9];
      [encoder setBuffer:Zstar_y offset:0 atIndex:10];
      [encoder setBuffer:geomBuffer offset:0 atIndex:11];
      [encoder setBuffer:physicsBuffer offset:0 atIndex:12];
      [encoder setBuffer:solverBuffer offset:0 atIndex:13];
      [encoder setBuffer:boundaryTypesBuffer offset:0 atIndex:14];
      [encoder setBuffer:boundaryValuesBuffer offset:0 atIndex:15];
      [encoder setBuffer:volOutBuffer offset:0 atIndex:16];

      {
        MTLSize tg = MTLSizeMake(kTileX, kTileY, 1);
        MTLSize groups = MTLSizeMake(
            (geom.xsz + kTileX - 1) / kTileX,
            (geom.ysz + 1 + (kTileY - 1) - 1) / (kTileY - 1), 1);
        [encoder dispatchThreadgroups:groups threadsPerThreadgroup:tg];
      }

      [encoder endEncoding];
      [commandBuffer commit];
      [commandBuffer waitUntilCompleted];
    }

    // Drain water from NoData cells if enabled
    if (pars->drain_nodata == ON && drainNodataPipeline) {
      @autoreleasepool {
        id<MTLCommandBuffer> commandBuffer = [commandQueue commandBuffer];
        id<MTLComputeCommandEncoder> encoder =
            [commandBuffer computeCommandEncoder];
        [encoder setComputePipelineState:drainNodataPipeline];

        [encoder setBuffer:H[dst] offset:0 atIndex:0];
        [encoder setBuffer:HU[dst] offset:0 atIndex:1];
        [encoder setBuffer:HV[dst] offset:0 atIndex:2];
        [encoder setBuffer:DEM offset:0 atIndex:3];
        [encoder setBuffer:geomBuffer offset:0 atIndex:4];

        float nd = gpu_nodata;
        [encoder setBytes:&nd length:sizeof(float) atIndex:5];
        [encoder setBuffer:volOutBuffer offset:0 atIndex:6];

        MTLSize gridSize = MTLSizeMake(geom.xsz, geom.ysz, 1);
        MTLSize threadGroupSize =
            MTLSizeMake(16, 16, 1); // Optimized for M1 Max
        [encoder dispatchThreads:gridSize
            threadsPerThreadgroup:threadGroupSize];
        [encoder endEncoding];
        [commandBuffer commit];
        [commandBuffer waitUntilCompleted];
      }
    }

    // Track the running maximum depth (CPU equivalent: update_max_field() once
    // per timestep in fv1::solve)
    if (updateMaxFieldPipeline) {
      @autoreleasepool {
        id<MTLCommandBuffer> commandBuffer = [commandQueue commandBuffer];
        id<MTLComputeCommandEncoder> encoder =
            [commandBuffer computeCommandEncoder];
        [encoder setComputePipelineState:updateMaxFieldPipeline];
        [encoder setBuffer:H[dst] offset:0 atIndex:0];
        [encoder setBuffer:maxH offset:0 atIndex:1];
        [encoder setBuffer:geomBuffer offset:0 atIndex:2];
        [encoder dispatchThreads:MTLSizeMake(geom.xsz, geom.ysz, 1)
            threadsPerThreadgroup:MTLSizeMake(16, 16, 1)];
        [encoder endEncoding];
        [commandBuffer commit];
      }
    }

    current_buffer = dst;
    time += (NUMERIC_TYPE)adaptive_dt;

    // Report outflow to mass balance only when simulate() needs it (at
    // interval) We'll add a helper to sync this correctly.

    return (NUMERIC_TYPE)adaptive_dt;
  }
};

Simulator::Simulator(Fnames &filenames, Pars &pars, ::Solver &solver,
                     Arrays &arrays, States &states, BoundCs &boundCs,
                     int verbose) {
  impl = new SimulatorImpl();
  impl->arrays = &arrays;
  impl->boundCs = &boundCs;
  impl->pars = &pars;
  impl->solver = &solver;
  impl->states = &states;

  if (!impl->init(filenames, pars, solver, verbose)) {
    delete impl;
    impl = nullptr;
  }
}

Simulator::~Simulator() {
  if (impl)
    delete impl;
}

void Simulator::initialize(NUMERIC_TYPE *H, NUMERIC_TYPE *HU, NUMERIC_TYPE *HV,
                           NUMERIC_TYPE *DEM, NUMERIC_TYPE *Mann) {
  if (impl) {
    impl->copyToGPU(H, HU, HV, DEM, Mann);
  }
}

NUMERIC_TYPE Simulator::step() {
  if (impl) {
    return impl->step();
  }
  return 0;
}

bool Simulator::is_available() {
  id<MTLDevice> device = MTLCreateSystemDefaultDevice();
  return device != nil;
}

void Simulator::copy_to_host(NUMERIC_TYPE *H, NUMERIC_TYPE *HU,
                             NUMERIC_TYPE *HV) {
  if (impl) {
    impl->copyFromGPU(H, HU, HV);
  }
}

void Simulator::copy_max_to_host(NUMERIC_TYPE *maxH) {
  if (impl) {
    impl->copyMaxFromGPU(maxH);
  }
}

void simulate(Fnames &filenames, Files *files, States &states, Pars &pars,
              ::Solver &solver, BoundCs &boundCs, Stage &stage, Arrays &arrays,
              int verbose) {
  if (states.rainfallmask == ON || states.routing == ON) {
    fprintf(stderr, "ERROR: Metal FV1 does not support rainfallmask or rainfallrouting; "
                    "use unmasked rainfall or dynamicrainfile.\n");
    std::exit(EXIT_FAILURE);
  }
  if (states.rainfall == ON && !arrays.rain) {
    fprintf(stderr, "ERROR: Metal FV1 requires a loaded rainfall time series.\n");
    std::exit(EXIT_FAILURE);
  }
  Simulator sim(filenames, pars, solver, arrays, states, boundCs, verbose);

  // Check if Metal initialization succeeded
  if (!sim.is_valid()) {
    printf("ERROR: Metal GPU initialization failed! Check that:\n");
    printf("  1. Metal shader library (metal_fv1_kernels.metallib) exists\n");
    printf("  2. Metal GPU is available on this system\n");
    printf("  3. Single precision mode is enabled (NUMERIC_TYPE=float)\n");
    return;
  }

  sim.initialize(arrays.H, arrays.HU, arrays.HV, arrays.DEM, arrays.Manningsn);

  if (verbose) {
    printf("Metal GPU initialized successfully.\n");
    printf("Grid: %d x %d cells\n", pars.xsz, pars.ysz);
    printf("SaveInt: %.2f, MassInt: %.2f\n", pars.SaveInt, pars.MassInt);
    printf("Simulation time: %.2f seconds\n", solver.Sim_Time);
    printf("\nStarting GPU simulation...\n\n");
    fflush(stdout);
  }

  // Initialize Zstar fields (needed for bed slope source terms)
  ::fv1::initialise_Zstar(&pars, &arrays);

  // Main simulation loop
  time_t loop_start;
  time(&loop_start);

  solver.MinTstep = solver.InitTstep;
  solver.vol1 = DomainVol(&states, &pars, nullptr, &arrays, nullptr);

  while (solver.t < solver.Sim_Time) {
    // Reset Qin/Qout for this step (matches CPU behavior)
    zero_flux_stats(&boundCs);

    // Execute one timestep on GPU
    NUMERIC_TYPE dt = sim.step();

    if (dt <= 0) {
      if (verbose)
        printf("\nERROR: Metal simulator returned zero timestep.\n");
      break;
    }

    if (solver.t > 0) {
      solver.MinTstep = std::min(solver.MinTstep, dt);
    }
    solver.t += dt;
    solver.itCount += 1;

    // Progress output
    if (verbose && solver.itCount % 100 == 0) {
      printf("t=%.2f\tdt=%.6f\n", solver.t, dt);
      fflush(stdout);
    }

    // Mass balance output
    if (solver.t >= pars.MassTotal) {
      sim.copy_to_host(arrays.H, arrays.HU, arrays.HV);

      // Compute boundary fluxes for mass balance (Metal doesn't track these
      // internally) We need to call the CPU FV1 functions to populate FHx/FHy
      // arrays at boundaries
      ::fv1::update_Hstar(&pars, &arrays);
      ::fv1::update_fluxes(&pars, &solver, &arrays);
      ::fv1::update_fluxes_on_boundaries(&pars, &solver, &boundCs, &arrays);
      // GPU interface fluxes and point sources are already integrated per step.

      // Virtual gauge output
      if (states.gsection == ON && files->gau_fp) {
        fprintf(files->gau_fp, "%12.2" NUM_FMT "", solver.t); // print time
        for (int i = 0; i < stage.Ngauges; i++) {
          // Calculate discharge using unit fluxes (FHx, FHy)
          // CalcVirtualGauge sums up values in the unit flux arrays
          NUMERIC_TYPE unit_discharge =
              CalcVirtualGauge(i, pars.xsz + 1, arrays.FHx, arrays.FHy, &stage);

          // Convert unit discharge (m2/s) to discharge (m3/s)
          // Q = q * width
          NUMERIC_TYPE discharge = 0.0;
          if (stage.gauge_dir[i] == North || stage.gauge_dir[i] == South) {
            discharge = unit_discharge * pars.dx;
          } else {
            discharge = unit_discharge * pars.dy;
          }

          fprintf(files->gau_fp, " %10.3" NUM_FMT "", discharge);
        }
        fprintf(files->gau_fp, "\n");
        fflush(files->gau_fp);
      }

      update_mass_stats(&states, &pars, &solver, &boundCs, &arrays);
      write_mass_stats(files, &pars, &solver, &boundCs);
    }

    // Solution output
    if (solver.t >= pars.SaveTotal) {
      sim.copy_to_host(arrays.H, arrays.HU, arrays.HV);
      if (states.voutput == ON) {
        update_velocity(&pars, &solver, &arrays);
      }
      write_solution(&filenames, &states, &pars, &solver, &arrays);
    }
  }

  // Final copy and final solution output to ensure the last step is captured
  // correctly
  sim.copy_to_host(arrays.H, arrays.HU, arrays.HV);

  if (states.voutput == ON) {
    update_velocity(&pars, &solver, &arrays);
  }

  // Force a write of the final state (SaveNo will be correct)
  write_solution(&filenames, &states, &pars, &solver, &arrays);

  // The maximum depth field is tracked on the GPU every timestep; bring it
  // back and write the .max file, as fv1::solve does on the CPU.
  sim.copy_max_to_host(arrays.maxH);
  write_max_field(&filenames, &states, &pars, &arrays);

  time_t loop_end;
  time(&loop_end);
  double elapsed = difftime(loop_end, loop_start);

  if (verbose) {
    printf("\nloop time %lf\n", elapsed);
    printf("Finished.\n\n");
  }
}

} // namespace fv1
} // namespace metal
} // namespace lis

#endif // __APPLE__
