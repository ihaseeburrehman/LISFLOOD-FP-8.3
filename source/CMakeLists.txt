cmake_minimum_required(VERSION 3.13 FATAL_ERROR)
include(CheckLanguage)

if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES 61)
endif()

project(lisflood-fp)

set(_CONFIG config.default.cmake CACHE FILEPATH "CMake config file")
include(${_CONFIG})

if (NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()


enable_language(CXX)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)

add_library(lisflood-base INTERFACE)
target_sources(lisflood-base INTERFACE
	boundary.cpp
	ch_flow.cpp
	chkpnt.cpp
	fp_acc.cpp
	fp_flow.cpp
	fp_trent.cpp
	infevap.cpp
	input.cpp
	iterateq.cpp
	output.cpp
	pars.cpp
	por_flow.cpp
	sgc.cpp
	util.cpp
	utility.cpp
	weir_flow.cpp
	lisflood2/DataTypes.cpp
	lisflood2/file_tool.cpp
	lisflood2/lis2_output.cpp
	lisflood2/lisflood_processing.cpp
	lisflood2/sgm_fast.cpp
	lisflood2/windows_ext.cpp
	rain/rain.cpp
	swe/boundary.cpp
	swe/dg2.cpp
	swe/fields.cpp
	swe/flux.cpp
	swe/fv1.cpp
	swe/hll.cpp
	swe/input.cpp
	swe/output.cpp
	swe/stats.cpp
	swe/fv1/modifiedvars.cpp
	swe/dg2/dg2_output.cpp
	swe/dg2/fields.cpp
	swe/dg2/friction.cpp
	swe/dg2/modifiedvars.cpp
	swe/dg2/slope_limiter.cpp
)

add_executable(lisflood)
target_sources(lisflood PRIVATE lisflood.cpp)
target_link_libraries(lisflood PRIVATE lisflood-base)


check_language(CUDA)
if (CMAKE_CUDA_COMPILER)
	enable_language(CUDA)
	set_property(TARGET lisflood PROPERTY CUDA_SEPARABLE_COMPILATION ON)
	
	add_compile_definitions(CUDA THRUST_IGNORE_CUB_VERSION_CHECK)
    
	target_sources(lisflood PRIVATE
		cuda/cuda_boundary.cu
		cuda/cuda_dem.cu
		cuda/cuda_flow.cu
		cuda/cuda_geometry.cu
		cuda/cuda_hll.cu
		cuda/cuda_max_field.cu
		cuda/cuda_rain.cu
		cuda/cuda_sample.cu
		cuda/cuda_simulate.cu
		cuda/cuda_snapshot.cu
		cuda/cuda_solver.cu
		cuda/cuda_stats.cu
		cuda/cuda_unifiedallocator.cu
		cuda/cuda_util.cu
		cuda/ghostraster.cpp
		cuda/io.cpp
		cuda/sample.cpp
		cuda/stats.cpp
		cuda/fv1/cuda_fv1_flow.cu
		cuda/fv1/cuda_fv1_simulate.cu
		cuda/fv1/cuda_fv1_snapshot.cu
		cuda/fv1/cuda_fv1_solver.cu
		cuda/fv2/cuda_fv2_dem.cu
		cuda/fv2/cuda_fv2_flow.cu
		cuda/fv2/cuda_fv2_simulate.cu
		cuda/fv2/cuda_fv2_slope_limit.cu
		cuda/fv2/cuda_fv2_snapshot.cu
		cuda/fv2/cuda_fv2_solver.cu
		cuda/dg2/cuda_dg2_dem.cu
		cuda/dg2/cuda_dg2_flow.cu
		cuda/dg2/cuda_dg2_simulate.cu
		cuda/dg2/cuda_dg2_slope_limit.cu
		cuda/dg2/cuda_dg2_snapshot.cu
		cuda/dg2/cuda_dg2_solver.cu
		cuda/acc/cuda_acc_flow.cu
		cuda/acc/cuda_acc_simulate.cu
		cuda/acc/cuda_acc_snapshot.cu
		cuda/acc/cuda_acc_solver.cu
		cuda/acc/cuda_acc_max.cu
		cuda/acc_nugrid/cuda_acc_nugrid_simulate.cu
		cuda/acc_nugrid/get_num_blocks.cpp
		cuda/acc_nugrid/cuda_utils.cu
		cuda/acc_nugrid/read_gauge_points.cu
		cuda/acc_nugrid/read_point_srcs.cu
		cuda/acc_nugrid/read_bdy_conds.cu
		cuda/acc_nugrid/gen_bdy_morton_codes.cu
		cuda/acc_nugrid/preflag_details.cu
		cuda/acc_nugrid/get_modal_values.cu
		cuda/acc_nugrid/read_raster_file.cpp
		cuda/acc_nugrid/generate_all_morton_codes.cu
		cuda/acc_nugrid/sort_finest_scale_coefficients_z_order.cu
		cuda/acc_nugrid/copy_finest_coefficients.cu
		cuda/acc_nugrid/get_max_scale_coefficients.cu
		cuda/acc_nugrid/preflag_topo.cu
		cuda/acc_nugrid/encode_and_thresh_topo.cu
		cuda/acc_nugrid/get_reg_tree.cu
		cuda/acc_nugrid/reg_new.cu
		cuda/acc_nugrid/project_assem_sol.cu
		cuda/acc_nugrid/reinsert_assem_sol.cu
		cuda/acc_nugrid/traverse_tree_of_sig_details.cu
		cuda/acc_nugrid/reverse_z_order_assembled_solution.cu
		cuda/acc_nugrid/write_reals_to_file.cu
		cuda/acc_nugrid/write_raster_maps.cu
		cuda/acc_nugrid/write_soln_vtk.cu
		cuda/acc_nugrid/reverse_z_order_act_idcs.cu
		cuda/acc_nugrid/get_compaction_flags.cu
		cuda/acc_nugrid/compaction.cu
		cuda/acc_nugrid/count_neighbours.cu
		cuda/acc_nugrid/find_nonuniform_neighbours.cu
		cuda/acc_nugrid/count_interfaces_per_neighbours.cu
		cuda/acc_nugrid/sort_nghbr_assem_sol_row_major.cu
		cuda/acc_nugrid/init_neighbours.cu
		cuda/acc_nugrid/init_interfaces.cu
		cuda/acc_nugrid/find_interfaces.cu
		cuda/acc_nugrid/load_interface_q_vol.cu
		cuda/acc_nugrid/load_neighbour_h.cu
		cuda/acc_nugrid/compute_q.cu
		cuda/acc_nugrid/update_h.cu
		cuda/acc_nugrid/calculate_dt.cu
		cuda/acc_nugrid/get_dt_CFL.cu
		cuda/acc_nugrid/init_h.cu
		cuda/acc_nugrid/init_q.cu
		cuda/acc_nugrid/read_all_bdy_conds.cu
		cuda/acc_nugrid/read_num_stage_points.cpp
		cuda/acc_nugrid/read_num_bdy_cells.cpp
		cuda/acc_nugrid/read_num_point_srcs.cpp
		cuda/acc_nugrid/write_gauge_point_data.cu
		cuda/acc_nugrid/unifiedallocator.cu
		cuda/acc_nugrid/rain.cu
		cuda/acc_nugrid/write_mass_data.cu
		cuda/acc_nugrid/rain_uniform.cu
		cuda/acc_nugrid/update_max.cu
		cuda/acc_nugrid/write_max_maps.cu
		cuda/acc_nugrid/init_max.cu
		cuda/acc_nugrid/write_chkpnt.cu
		cuda/acc_nugrid/read_chkpnt.cu
	)
    
    add_subdirectory(cuda/adaptive)
    
    target_compile_options(lisflood PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--generate-line-info>)

    
	target_include_directories(lisflood PUBLIC cuda cuda/common)
endif()

 include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})



add_executable(generateDG2DEM preprocess/generateDG2DEM.cpp)
target_link_libraries(generateDG2DEM PRIVATE lisflood-base)

add_executable(generateDG2start preprocess/generateDG2start.cpp)
target_link_libraries(generateDG2start PRIVATE lisflood-base)

add_executable(DG2downscale postprocess/DG2downscale.cpp)
target_link_libraries(DG2downscale PRIVATE lisflood-base)

find_package(OpenMP)
if (OpenMP_CXX_FOUND)
	target_link_libraries(lisflood-base INTERFACE OpenMP::OpenMP_CXX)
endif()

if (UNIX)
	find_package(NUMA REQUIRED)
	target_include_directories(lisflood-base INTERFACE ${NUMA_INCLUDE_DIRS})
	target_link_libraries(lisflood-base INTERFACE ${NUMA_LIBRARIES})
endif()

if (_NETCDF)
	if (WIN32)
		target_include_directories(lisflood-base INTERFACE windep/netCDF4-64/include) # TODO: x64 only
		target_link_directories(lisflood-base INTERFACE windep/netCDF4-64/lib) # TODO: x64 only
		target_link_libraries(lisflood-base INTERFACE netcdf)
		file(GLOB DLLS "DLL's/*.dll")
		file(COPY ${DLLS} DESTINATION ${CMAKE_BINARY_DIR})
	else()
		find_package(NetCDF REQUIRED)
		target_include_directories(lisflood-base INTERFACE ${NetCDF_INCLUDE_DIRS})
		target_link_libraries(lisflood-base INTERFACE NetCDF::NetCDF)
	endif()
	target_compile_definitions(lisflood-base INTERFACE _NETCDF=1)
else()
	target_compile_definitions(lisflood-base INTERFACE _NETCDF=0)
endif()

if (CMAKE_BUILD_TYPE MATCHES Debug)
	target_compile_definitions(lisflood-base INTERFACE _DEBUG)
endif()
