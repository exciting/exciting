# Build all the necessary executables
set(MOD_DIR "${CMAKE_BINARY_DIR}/include")
set(linkdir ${MOD_DIR})
set(LINKS deviceacc)

add_executable(test_lu_inversion test_lu_inversion.f90)
target_include_directories(test_lu_inversion PUBLIC ${linkdir})
target_link_libraries(test_lu_inversion ${LINKS} ${devacc_link_libs})
set_property(TARGET test_lu_inversion PROPERTY LINKER_LANGUAGE Fortran)

add_executable(test_fft test_fft.f90)
target_include_directories(test_fft PUBLIC ${linkdir})
target_link_libraries(test_fft ${LINKS} ${devacc_link_libs})
set_property(TARGET test_fft PROPERTY LINKER_LANGUAGE Fortran)

add_executable(test_zgemm test_zgemm.f90)
target_include_directories(test_zgemm PUBLIC ${linkdir})
target_link_libraries(test_zgemm ${LINKS} ${devacc_link_libs})
set_property(TARGET test_zgemm PROPERTY LINKER_LANGUAGE Fortran)

add_executable(test_zgemm_filled_in_gpu test_zgemm_filled_in_gpu.f90)
target_include_directories(test_zgemm_filled_in_gpu PUBLIC ${linkdir})
target_link_libraries(test_zgemm_filled_in_gpu ${LINKS} ${devacc_link_libs})
set_property(TARGET test_zgemm_filled_in_gpu PROPERTY LINKER_LANGUAGE Fortran)

add_executable(test_do_loop_gpu test_do_loop_gpu.f90)
target_include_directories(test_do_loop_gpu PUBLIC ${linkdir})
target_link_libraries(test_do_loop_gpu ${LINKS} ${devacc_link_libs})
set_property(TARGET test_do_loop_gpu PROPERTY LINKER_LANGUAGE Fortran)

# Enable testing
include(CTest)
add_test(NAME test_lu_inversion COMMAND test_lu_inversion)
add_test(NAME test_fft COMMAND test_fft)
add_test(NAME test_zgemm COMMAND test_zgemm)
add_test(NAME test_zgemm_filled_in_gpu COMMAND test_zgemm_filled_in_gpu)
add_test(NAME test_do_loop_gpu COMMAND test_do_loop_gpu)
