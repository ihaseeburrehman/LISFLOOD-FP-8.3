template<typename F>
void lis::AsciiRaster::write
(
	FILE* file,
	F array,
	Geometry& geometry,
	int pitch,
	int offset,
	NUMERIC_TYPE no_data_value,
	int outflag, 
	int precision 
	
)
{
	if (outflag == 0) {
		fprintf(file, "ncols         %i\n", geometry.xsz);
		fprintf(file, "nrows         %i\n", geometry.ysz);
		fprintf(file, "xllcorner     %.*" NUM_FMT"\n", precision, geometry.blx);
		fprintf(file, "yllcorner     %.*" NUM_FMT"\n", precision, geometry.bly);
		fprintf(file, "cellsize      %.*" NUM_FMT"\n", precision, geometry.dx);
		fprintf(file, "NODATA_value  %.*" NUM_FMT"\n", precision, no_data_value);

		for (int j = 0; j < geometry.ysz; j++)
		{
			for (int i = 0; i < geometry.xsz; i++)
			{
				fprintf(file, "%.*" NUM_FMT "\t", precision,
					array[j * pitch + i + offset]);
			}
			fprintf(file, "\n");
		}
	} 
	else if (outflag == 1) {
		fprintf(file, "ncols         %i\n", geometry.xsz + 1);
		fprintf(file, "nrows         %i\n", geometry.ysz);
		fprintf(file, "xllcorner     %.*" NUM_FMT"\n", precision, geometry.blx - (geometry.dx / C(2.0)));
		fprintf(file, "yllcorner     %.*" NUM_FMT"\n", precision, geometry.bly);
		fprintf(file, "cellsize      %.*" NUM_FMT"\n", precision, geometry.dx);
		fprintf(file, "NODATA_value  %.*" NUM_FMT"\n", precision, no_data_value);

		for (int j = 0; j < geometry.ysz; j++)
		{
			for (int i = 0; i < geometry.xsz + 1; i++)
			{
				fprintf(file, "%.*" NUM_FMT "\t", precision,
					array[j * (pitch + 1) + i + offset]);
			}
			fprintf(file, "\n");
		}
	}
	else if (outflag == 2) {
		fprintf(file, "ncols         %i\n", geometry.xsz);
		fprintf(file, "nrows         %i\n", geometry.ysz + 1);
		fprintf(file, "xllcorner     %.*" NUM_FMT"\n", precision, geometry.blx);
		fprintf(file, "yllcorner     %.*" NUM_FMT"\n", precision, geometry.bly - (geometry.dx / C(2.0)));
		fprintf(file, "cellsize      %.*" NUM_FMT"\n", precision, geometry.dx);
		fprintf(file, "NODATA_value  %.*" NUM_FMT"\n", precision, no_data_value);

		for (int j = 0; j < geometry.ysz + 1; j++)
		{
			for (int i = 0; i < geometry.xsz; i++)
			{
				fprintf(file, "%.*" NUM_FMT "\t", precision,
					array[j * (pitch + 1) + i + offset]);
			}
			fprintf(file, "\n");
		}
	}
}

template<typename F>
std::future<void> lis::Snapshot::write_async
(
	F array,
	Geometry& geometry,
	int pitch,
	int offset,
	const char* prefix,
	int counter,
	const char* suffix,
	int verbose,
	int outflag,
	const int& call_gzip,
	int precision,
	NUMERIC_TYPE no_data_value
)
{
	return std::async(std::launch::async, [=, &geometry]{
		Snapshot::write(array, geometry, pitch, offset, prefix, counter, suffix,
				verbose, outflag, call_gzip, precision, no_data_value);
	});
}


template<typename F>
void lis::Snapshot::write
(
	F array,
	Geometry& geometry,
	int pitch,
	int offset,
	const char* prefix,
	int counter,
	const char* suffix,
	int verbose,
	int outflag,
	const int& call_gzip,
	int precision,
	NUMERIC_TYPE no_data_value
	
)
{
	char filename[800];
	char tmp_sys_com[255];

	snprintf(filename, 800*sizeof(char), "%s-%.4d%s", prefix, counter, suffix);

	FILE* file = fopen_or_die(filename, "wb");
	
	AsciiRaster::write(file, array, geometry, pitch, offset, no_data_value,
		     outflag, precision); 

	fclose(file);

	// check if we need to zip the file up
	if (call_gzip == ON)
	{
		sprintf(tmp_sys_com, "%s%s", "gzip -9 -f ", filename);
		system(tmp_sys_com);
	}
}

template<typename F>
void lis::Snapshot::write_max
(
	F array,
	Geometry& geometry,
	int pitch,
	int offset,
	const char* prefix,
	const char* suffix,
	int verbose,
	int outflag,
	int precision,
	NUMERIC_TYPE no_data_value

)
{
	char filename[800];
	snprintf(filename, 800 * sizeof(char), "%s%s", prefix, suffix);

	FILE* file = fopen_or_die(filename, "wb");

	AsciiRaster::write(file, array, geometry, pitch, offset, no_data_value, outflag, 
		precision);

	fclose(file);

}

