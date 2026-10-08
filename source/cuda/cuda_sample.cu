#include "cuda_sample.cuh"
#include "cuda_util.cuh"
#include "fv1/cuda_fv1_flow.cuh"
#include "fv1/cuda_fv1_solver.cuh"
#include "acc/cuda_acc_flow.cuh" 
#include "acc/cuda_acc_solver.cuh" 

namespace lis
{
namespace cuda
{

__global__ void sample
(
	NUMERIC_TYPE* H,
	NUMERIC_TYPE* HU,
	NUMERIC_TYPE* HV,
	SampleBuffer sample_buf,
	NUMERIC_TYPE t
)
{
	sample_buf.time[cuda::sample_buf_idx] = t;

	for (int i=0; i<cuda::sample_points.count; i++)
	{
		int buf_offset = cuda::sample_buf_idx*cuda::sample_points.count;

		if (cuda::sample_points.inside_domain[i])
		{
			NUMERIC_TYPE Hval = H[cuda::sample_points.idx[i]];
			sample_buf.H[buf_offset + i] = Hval;

			if (Hval > cuda::solver_params.DepthThresh)
			{
				NUMERIC_TYPE HUval = HU[cuda::sample_points.idx[i]];
				NUMERIC_TYPE HVval = HV[cuda::sample_points.idx[i]];

				NUMERIC_TYPE speed = SQRT(HUval/Hval * HUval/Hval
						+ HVval/Hval * HVval/Hval);
				sample_buf.speed[buf_offset + i] = speed;
			}
		}
	}

	cuda::sample_buf_idx++;
}

__global__ void sample_ACC
(
	NUMERIC_TYPE* H,
	NUMERIC_TYPE* Vx,
	NUMERIC_TYPE* Vy,
	SampleBuffer sample_buf,
	NUMERIC_TYPE t
)
{
	sample_buf.time[cuda::sample_buf_idx] = t;

	for (int i = 0; i < cuda::sample_points.count; i++)
	{
		int buf_offset = cuda::sample_buf_idx * cuda::sample_points.count;

		if (cuda::sample_points.inside_domain[i])
		{
			NUMERIC_TYPE Hval = H[cuda::sample_points.idx[i]];
			sample_buf.H[buf_offset + i] = Hval;

			if (Hval > cuda::solver_params.DepthThresh)
			{
				NUMERIC_TYPE Vxval1 = Vx[cuda::sample_points.idx_Qx1[i]];
				NUMERIC_TYPE Vxval2 = Vx[cuda::sample_points.idx_Qx2[i]];
				NUMERIC_TYPE Vyval1 = Vy[cuda::sample_points.idx_Qy1[i]];
				NUMERIC_TYPE Vyval2 = Vy[cuda::sample_points.idx_Qy2[i]];

				NUMERIC_TYPE speed = SQRT(getmax(FABS(Vxval1), FABS(Vxval2)) * getmax(FABS(Vxval1), FABS(Vxval2))
					                    + getmax(FABS(Vyval1), FABS(Vyval2)) * getmax(FABS(Vyval1), FABS(Vyval2)));
				sample_buf.speed[buf_offset + i] = speed;
			}
		}
	}

	cuda::sample_buf_idx++;
}


}
}

lis::cuda::Sampler::Sampler
(
	SamplePoints& d_sample_points,
	int& sample_buf_idx,
	int verbose
)
:
d_sample_points(d_sample_points),
sample_buf_idx(sample_buf_idx),
stage_file(sample_points),
gauge_file(sample_points),
verbose(verbose)
{}

void lis::cuda::Sampler::load_sample_points
(
	const char* filename,
	Geometry& geometry,
	int pitch,
	int offset
)
{
	lis::Sample::initialise(sample_points, filename, geometry, pitch, offset,
			verbose);
	initialise_sample_points();

	allocate_pinned(sample_buf, sample_points.count);
	allocate_device(d_sample_buf, sample_points.count);

	active = true;
}

void lis::cuda::Sampler::open_stage_file
(
	const char* filename,
	int checkpoint,
	NUMERIC_TYPE t
)
{
	stage_file.open(filename, checkpoint, t);
}

void lis::cuda::Sampler::open_gauge_file
(
	const char* filename,
	int checkpoint,
	NUMERIC_TYPE t
)
{
	gauge_file.open(filename, checkpoint, t);

	write_speed = true;
}

void lis::cuda::Sampler::write_stage_header
(
	NUMERIC_TYPE* DEM,
	const char* sample_points_filename,
	int checkpoint,
	NUMERIC_TYPE t
)
{
	stage_file.write_header(DEM, sample_points_filename, checkpoint, t);
}

void lis::cuda::Sampler::write_gauge_header
(
	NUMERIC_TYPE* DEM,
	const char* sample_points_filename,
	NUMERIC_TYPE t
)
{
	gauge_file.write_header(DEM, sample_points_filename, t);
}

void lis::cuda::Sampler::sample
(
	NUMERIC_TYPE* H,
	NUMERIC_TYPE* HU,
	NUMERIC_TYPE* HV,
	NUMERIC_TYPE t
)
{
	if (active) lis::cuda::sample<<<1, 1>>>(H, HU, HV, d_sample_buf, t);
}

void lis::cuda::Sampler::sample_ACC
(
	NUMERIC_TYPE* H,
	NUMERIC_TYPE* Vx,
	NUMERIC_TYPE* Vy,
	NUMERIC_TYPE t
)
{
	if (active) lis::cuda::sample_ACC << <1, 1 >> > (H, Vx, Vy, d_sample_buf, t);
}

void lis::cuda::Sampler::write_if_buffer_full()
{
	if (active && buffer_full())
	{
		if (verbose == ON) printf("SampleBuffer full: flushing to disk\n");
		copy_buffer();
		stage_file.write(sample_buf, sample_buf_idx);
		if (write_speed) gauge_file.write(sample_buf, sample_buf_idx);
		sample_buf_idx = 0;
	}
}

void lis::cuda::Sampler::write()
{
	if (active)
	{
		if (verbose == ON) printf("Flushing SampleBuffer to disk\n");
		copy_buffer();
		stage_file.write(sample_buf, sample_buf_idx);
		if (write_speed) gauge_file.write(sample_buf, sample_buf_idx);
	}
}

bool lis::cuda::Sampler::buffer_full()
{
	return sample_buf_idx == sample_buf.size;
}

void lis::cuda::Sampler::copy_buffer()
{
	SampleBuffer& dst = sample_buf;
	SampleBuffer& src = d_sample_buf;

	cuda::copy(dst.time, src.time, sample_buf_idx*sizeof(NUMERIC_TYPE));
	cuda::copy(dst.H, src.H,
			sample_buf_idx*sample_points.count*sizeof(NUMERIC_TYPE));
	cuda::copy(dst.speed, src.speed,
			sample_buf_idx*sample_points.count*sizeof(NUMERIC_TYPE));
}

void lis::cuda::Sampler::initialise_sample_points()
{
	SamplePoints& d_dst = d_sample_points;
	SamplePoints& h_src = sample_points;

	SamplePoints temp;
	int count = temp.count = h_src.count;

	temp.x = static_cast<NUMERIC_TYPE*>(
			malloc_device(count*sizeof(NUMERIC_TYPE)));
	temp.y = static_cast<NUMERIC_TYPE*>(
			malloc_device(count*sizeof(NUMERIC_TYPE)));
	temp.idx = static_cast<int*>(malloc_device(count*sizeof(int)));
	temp.idx_Qx1 = static_cast<int*>(malloc_device(count * sizeof(int))); 
	temp.idx_Qx2 = static_cast<int*>(malloc_device(count * sizeof(int))); 
	temp.idx_Qy1 = static_cast<int*>(malloc_device(count * sizeof(int))); 
	temp.idx_Qy2 = static_cast<int*>(malloc_device(count * sizeof(int))); 
	temp.inside_domain = static_cast<bool*>(malloc_device(count*sizeof(bool)));

	cuda::copy(temp.x, h_src.x, count*sizeof(NUMERIC_TYPE));
	cuda::copy(temp.y, h_src.y, count*sizeof(NUMERIC_TYPE));
	cuda::copy(temp.idx, h_src.idx, count*sizeof(int));
	cuda::copy(temp.idx_Qx1, h_src.idx_Qx1, count * sizeof(int)); 
	cuda::copy(temp.idx_Qx2, h_src.idx_Qx2, count * sizeof(int)); 
	cuda::copy(temp.idx_Qy1, h_src.idx_Qy1, count * sizeof(int)); 
	cuda::copy(temp.idx_Qy2, h_src.idx_Qy2, count * sizeof(int)); 
	cuda::copy(temp.inside_domain, h_src.inside_domain, count*sizeof(bool));

	copy_to_symbol(d_dst, &temp, sizeof(SamplePoints));
}

void lis::cuda::Sampler::free
(
	SamplePoints& d_sample_points
)
{
	cuda::free_device(d_sample_points.x);
	cuda::free_device(d_sample_points.y);
	cuda::free_device(d_sample_points.idx);
	cuda::free_device(d_sample_points.idx_Qx1); 
	cuda::free_device(d_sample_points.idx_Qx2); 
	cuda::free_device(d_sample_points.idx_Qy1); 
	cuda::free_device(d_sample_points.idx_Qy2); 
	cuda::free_device(d_sample_points.inside_domain);
}

void lis::cuda::Sampler::allocate_pinned
(
	SampleBuffer& buf,
	int points,
	int size
)
{
	buf.size = size;
	buf.time = static_cast<NUMERIC_TYPE*>(
			malloc_pinned(size*sizeof(NUMERIC_TYPE)));
	buf.H = static_cast<NUMERIC_TYPE*>(
			malloc_pinned(size*points*sizeof(NUMERIC_TYPE)));
	buf.speed = static_cast<NUMERIC_TYPE*>(
			malloc_pinned(size*points*sizeof(NUMERIC_TYPE)));
}

void lis::cuda::Sampler::allocate_device
(
	SampleBuffer& buf,
	int points,
	int size
)
{
	buf.size = size;
	buf.time = static_cast<NUMERIC_TYPE*>(
			malloc_device(size*sizeof(NUMERIC_TYPE)));
	buf.H = static_cast<NUMERIC_TYPE*>(
			malloc_device(size*points*sizeof(NUMERIC_TYPE)));
	buf.speed = static_cast<NUMERIC_TYPE*>(
			malloc_device(size*points*sizeof(NUMERIC_TYPE)));
}

void lis::cuda::Sampler::free_pinned
(
	SampleBuffer& buf
)
{
	cuda::free_pinned(buf.time);
	cuda::free_pinned(buf.H);
	cuda::free_pinned(buf.speed);
}

void lis::cuda::Sampler::free_device
(
	SampleBuffer& buf
)
{
	cuda::free_device(buf.time);
	cuda::free_device(buf.H);
	cuda::free_device(buf.speed);
}

lis::cuda::Sampler::~Sampler()
{
	if (active)
	{
		lis::Sample::free(sample_points);
		free(d_sample_points);
		free_pinned(sample_buf);
		free_device(d_sample_buf);
	}
}
