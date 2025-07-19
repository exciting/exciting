# Build all the necessary executables
set(MOD_DIR "${CMAKE_BINARY_DIR}/include")
set(linkdir ${MOD_DIR})

set(LINKS IDieL ${LINKS_DEVICEACC})

add_executable(test_crystal test_crystal.f90)
target_include_directories(test_crystal PUBLIC ${linkdir})
target_link_libraries(test_crystal ${LINKS})
set_target_properties(test_crystal PROPERTIES LINKER_LANGUAGE Fortran)

add_executable(test_quadrature_lebedev test_quadrature_lebedev.f90)
target_include_directories(test_quadrature_lebedev PUBLIC ${linkdir})
target_link_libraries(test_quadrature_lebedev ${LINKS})
set_target_properties(test_quadrature_lebedev PROPERTIES LINKER_LANGUAGE Fortran)

add_executable(test_dielectric_average_iso test_dielectric_average_iso.f90)
execute_process( COMMAND ${CMAKE_COMMAND} -E tar xzf ${PROJECT_SOURCE_DIR}/test_resources/dielectric/Binv.dat.tar.gz WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
execute_process( COMMAND ${CMAKE_COMMAND} -E tar xzf ${PROJECT_SOURCE_DIR}/test_resources/dielectric/head.dat.tar.gz WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
execute_process( COMMAND ${CMAKE_COMMAND} -E tar xzf ${PROJECT_SOURCE_DIR}/test_resources/dielectric/wingL.dat.tar.gz WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
execute_process( COMMAND ${CMAKE_COMMAND} -E tar xzf ${PROJECT_SOURCE_DIR}/test_resources/dielectric/wingU.dat.tar.gz WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(test_dielectric_average_iso PUBLIC ${linkdir})
target_link_libraries(test_dielectric_average_iso ${LINKS})
set_target_properties(test_dielectric_average_iso PROPERTIES LINKER_LANGUAGE Fortran)

add_executable(test_inverse test_inverse.f90)
target_include_directories(test_inverse PUBLIC ${linkdir})
target_link_libraries(test_inverse ${LINKS})
set_target_properties(test_inverse PROPERTIES LINKER_LANGUAGE Fortran)

add_executable(test_sph_harm test_sph_harm.f90)
file(COPY ${PROJECT_SOURCE_DIR}/test_resources/spherical_harmonics/real_part.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${PROJECT_SOURCE_DIR}/test_resources/spherical_harmonics/imag_part.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(test_sph_harm PUBLIC ${linkdir})
target_link_libraries(test_sph_harm ${LINKS})
set_target_properties(test_sph_harm PROPERTIES LINKER_LANGUAGE Fortran)

add_executable(test_spline test_spline.f90)
target_include_directories(test_spline PUBLIC ${linkdir})
target_link_libraries(test_spline ${LINKS})
set_target_properties(test_spline PROPERTIES LINKER_LANGUAGE Fortran)

add_executable(test_circular_expansion test_circular_expansion.f90)
target_include_directories(test_circular_expansion PUBLIC ${linkdir})
target_link_libraries(test_circular_expansion ${LINKS})
set_target_properties(test_circular_expansion PROPERTIES LINKER_LANGUAGE Fortran)

add_executable(test_2d test_2d.f90)
execute_process( COMMAND ${CMAKE_COMMAND} -E tar xzf ${PROJECT_SOURCE_DIR}/test_resources/dielectric/2d.tar.gz  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(test_2d PUBLIC ${linkdir})
target_link_libraries(test_2d ${LINKS})
set_target_properties(test_2d PROPERTIES LINKER_LANGUAGE Fortran)

# Enable testing
include(CTest)
add_test(NAME test_crystal COMMAND test_crystal)
add_test(NAME test_quadrature_lebedev COMMAND test_quadrature_lebedev)
add_test(NAME test_dielectric_average_iso COMMAND test_dielectric_average_iso)
add_test(NAME test_inverse COMMAND test_inverse)
add_test(NAME test_sph_harm COMMAND test_sph_harm)
add_test(NAME test_spline COMMAND test_spline)
add_test(NAME test_circular_expansion COMMAND test_circular_expansion)
add_test(NAME test_2d COMMAND test_2d)

