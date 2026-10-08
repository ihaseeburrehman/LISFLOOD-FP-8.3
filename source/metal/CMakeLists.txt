# Metal GPU back end for the FV1 solver (Apple Silicon)

find_library(METAL_FRAMEWORK Metal)
find_library(FOUNDATION_FRAMEWORK Foundation)

if(METAL_FRAMEWORK AND FOUNDATION_FRAMEWORK)
    # The kernels are compiled to a metallib next to the executable, which loads
    # it at run time. -fno-fast-math keeps the GPU arithmetic identical to the
    # CPU solver's, so the verification cases are bit-for-bit reproducible.
    set(METALLIB ${CMAKE_BINARY_DIR}/metal_fv1_kernels.metallib)
    add_custom_command(
        OUTPUT ${METALLIB}
        COMMAND xcrun -sdk macosx metal -fno-fast-math
                -c ${CMAKE_CURRENT_SOURCE_DIR}/fv1/metal_fv1_kernels.metal
                -o ${CMAKE_BINARY_DIR}/metal_fv1_kernels.air
        COMMAND xcrun -sdk macosx metallib ${CMAKE_BINARY_DIR}/metal_fv1_kernels.air
                -o ${METALLIB}
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/fv1/metal_fv1_kernels.metal
        COMMENT "Compiling Metal FV1 kernels")
    add_custom_target(metal_shaders DEPENDS ${METALLIB})

    add_library(lisflood-metal STATIC
        ${CMAKE_CURRENT_SOURCE_DIR}/metal_integrate.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/fv1/metal_fv1_simulate.mm)
    add_dependencies(lisflood-metal metal_shaders)
    target_include_directories(lisflood-metal PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR})
    target_link_libraries(lisflood-metal PUBLIC ${METAL_FRAMEWORK} ${FOUNDATION_FRAMEWORK} lisflood-base)
    target_compile_options(lisflood-metal PRIVATE -fobjc-arc)
    target_compile_definitions(lisflood-metal PUBLIC METAL_ENABLED)

    set(METAL_ENABLED TRUE PARENT_SCOPE)
    set(METAL_LIBRARY lisflood-metal PARENT_SCOPE)
else()
    message(WARNING "Metal framework not found: building without the Metal back end")
endif()
