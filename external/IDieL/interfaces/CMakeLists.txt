# This first creates a CXX library and then uses it to create python bindings
# to IDieL

set(pybind11_DIR ${PYBIND_CMAKE_DIR})
find_package(pybind11 REQUIRED CONFIG)

set(SRC_INTERFACES interface.f90
                   interface.cpp)

set(MOD_DIR "${CMAKE_BINARY_DIR}/include")

add_library(IDieL_interfaces SHARED ${SRC_INTERFACES})
target_include_directories(IDieL_interfaces PUBLIC ${MOD_DIR} ${pybind11_INCLUDE_DIRS})
if (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    target_link_libraries(IDieL_interfaces IDieL ifcore)
endif()

if (NVIDIA OR AMD)
    set_target_properties(IDieL_interfaces PROPERTIES
        INSTALL_RPATH "${MAGMA_DIR}/lib/libmagma.so:${INSTALL_DIR}/lib/"
        INSTALL_RPATH_USE_LINK_PATH TRUE
    )
else()
    set_target_properties(IDieL_interfaces PROPERTIES
        INSTALL_RPATH "${INSTALL_DIR}/lib/"
        INSTALL_RPATH_USE_LINK_PATH TRUE
    )
endif()

pybind11_add_module(IDieLPy MODULE interface.cpp)
target_include_directories(IDieLPy PRIVATE ${MOD_DIR} ${pybind11_INCLUDE_DIRS} ${Python3_INCLUDE_DIRS})

if (MKL)
    target_link_libraries(IDieLPy PRIVATE IDieL_interfaces ${Python3_LIBRARIES} mkl_rt)
endif()

if (NVIDIA OR AMD)
    set_target_properties(IDieLPy PROPERTIES
        INSTALL_RPATH "${INSTALL_DIR}/lib:${INSTALL_DIR}/interfaces:${SPGLIBDIR}" 
        INSTALL_RPATH_USE_LINK_PATH TRUE
    )
else()
    set_target_properties(IDieLPy PROPERTIES
        INSTALL_RPATH "${INSTALL_DIR}/lib:${INSTALL_DIR}/interfaces:${SPGLIBDIR}"
        INSTALL_RPATH_USE_LINK_PATH TRUE
    )
endif(GPU)

# Install
install(TARGETS IDieL_interfaces DESTINATION ${INSTALL_DIR}/interfaces)
install(TARGETS IDieLPy DESTINATION ${INSTALL_DIR}/interfaces)
