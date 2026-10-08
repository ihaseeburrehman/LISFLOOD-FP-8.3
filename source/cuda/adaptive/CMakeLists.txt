add_subdirectory(input)
add_subdirectory(mra)
add_subdirectory(neighbours)
add_subdirectory(operators)
add_subdirectory(output)
add_subdirectory(utilities)
add_subdirectory(zorder)

target_sources(
	    lisflood
		PRIVATE
		cuda_adaptive_simulate.cu
	)