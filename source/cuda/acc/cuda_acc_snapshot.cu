#include "cuda_acc_snapshot.cuh"
#include "io.h"

lis::cuda::acc::Snapshot::Snapshot
(
	const char* resrootname,
	NUMERIC_TYPE interval,
	NUMERIC_TYPE next_save,
	int counter,
	NUMERIC_TYPE* DEM,
	Flow& U,
	Geometry& geometry,
	int pitch,
	int offset,
	int verbose,
	int precision
)
:
cuda::Snapshot<Flow>(resrootname, interval, next_save, counter, U, geometry,
		pitch, offset, verbose, precision),
DEM(DEM),
writer_elev(std::async(std::launch::async, []{})),
writer_H(std::async(std::launch::async, []{})),
writer_HU(std::async(std::launch::async, []{})),
writer_HV(std::async(std::launch::async, []{})),
writer_U(std::async(std::launch::async, []{})),
writer_V(std::async(std::launch::async, []{}))
{}

void lis::cuda::acc::Snapshot::write(const int& call_gzip)
{
	writer_H = write_async(U.H, ".wd", 0, call_gzip);

	if (write_elevation)
	{
		writer_elev = write_async(ElevationWriter(U.H, DEM), ".elev", 0, call_gzip);
	}

	if (write_discharge)
	{
		writer_HU = write_async(U.Qxold, ".Qx", 1, call_gzip);
		writer_HV = write_async(U.Qyold, ".Qy", 2, call_gzip);
	}

	if (write_velocity)
	{
		writer_U = write_async(U.Vx, ".Vx", 1, call_gzip);
		writer_V = write_async(U.Vy, ".Vy", 2, call_gzip);
	}
}

void lis::cuda::acc::Snapshot::write_maxes()
{
	write_max(U.initHtm, ".inittm", 0);
	write_max(U.totalHtm, ".totaltm", 0);
	write_max(U.maxH, ".max", 0);
	write_max(U.maxHtm, ".maxtm", 0);


//	writer_elev = write_async_maxes(ElevationWriter(U.H, DEM), ".elev", 0);

}

void lis::cuda::acc::Snapshot::wait()
{
	writer_elev.wait();
	writer_H.wait(); 
	writer_HU.wait(); 
	writer_HV.wait();
    writer_U.wait(); 
	writer_V.wait(); 
}
