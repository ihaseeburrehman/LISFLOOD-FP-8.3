template<typename F>
lis::cuda::Snapshot<F>::Snapshot
(
	const char* resrootname,
	NUMERIC_TYPE interval,
	NUMERIC_TYPE next_save,
	int counter,
	F& U,
	Geometry& geometry,
	int pitch,
	int offset,
	int verbose,
	int precision
)
:
resrootname(resrootname),
interval(interval),
next_save(next_save),
counter(counter),
U(U),
geometry(geometry),
pitch(pitch),
offset(offset),
verbose(verbose),
precision(precision)
{}

template<typename F>
void lis::cuda::Snapshot<F>::write_if_needed
(
	F& d_U,
	NUMERIC_TYPE t,
	const int& call_gzip
)
{
	if (t >= next_save) {
		wait();
		F::copy(U, d_U, geometry);
		if (verbose) printf("Writing snapshot at t=%f\n", t);

		write(call_gzip);

		next_save += interval;
		counter++;
	}
}


template<typename F>
void lis::cuda::Snapshot<F>::write_max_files
(
	F& d_U
)
{

		F::copy(U, d_U, geometry);
		if (verbose) printf("Writing max files\n");

		write_maxes();

}


template<typename F>
void lis::cuda::Snapshot<F>::enable_elevation_writer()
{
	write_elevation = true;
}

template<typename F>
void lis::cuda::Snapshot<F>::enable_discharge_writer()
{
	write_discharge = true;
}

template<typename F>
void lis::cuda::Snapshot<F>::enable_velocity_writer
(
	NUMERIC_TYPE DepthThresh
)
{
	write_velocity = true;
	this->DepthThresh = DepthThresh;
}

template<typename F>
template<typename T>
std::future<void> lis::cuda::Snapshot<F>::write_async
(
	T field,
	const char* suffix,
	int outflag,
	const int& call_gzip
)
{	
	return lis::Snapshot::write_async(field, geometry, pitch, offset,
			resrootname, counter, suffix, verbose, outflag, call_gzip, precision);
}

template<typename F>
template<typename T>
void lis::cuda::Snapshot<F>::write_max
(
	T field,
	const char* suffix,
	int outflag
)
{

	return lis::Snapshot::write_max(field, geometry, pitch, offset, resrootname, suffix,
		verbose, outflag, precision);
}
