# Build all the necessary executables
set(MOD_DIR "${CMAKE_BINARY_DIR}/include")
set(linkdir ${MOD_DIR})

set(LINKS IDieL)

execute_process(COMMAND wget -q -O diel2d.tar.gz https://box.hu-berlin.de/f/57340f6837ba477b9b15/?dl=1 WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
execute_process(COMMAND wget -q -O diel3d.tar.gz https://box.hu-berlin.de/f/e904f7489e5e4ef981c7/?dl=1 WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
execute_process( COMMAND ${CMAKE_COMMAND} -E tar xzf ${CMAKE_CURRENT_BINARY_DIR}/diel2d.tar.gz WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
execute_process( COMMAND ${CMAKE_COMMAND} -E tar xzf ${CMAKE_CURRENT_BINARY_DIR}/diel3d.tar.gz WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_executable(procedures_3d procedures_3d.f90)
target_include_directories(procedures_3d PUBLIC ${linkdir})
target_link_libraries(procedures_3d ${LINKS})

add_executable(procedures_2d procedures_2d.f90)
target_include_directories(procedures_2d PUBLIC ${linkdir})
target_link_libraries(procedures_2d ${LINKS})

