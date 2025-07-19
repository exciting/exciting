# This CMake script is responsible for creating the documentation for IDieL

# Option to enable documentation generation.
option(DOCUMENTATION "Generate documentation for the exciting project" OFF)

if (DOCUMENTATION)
  find_package(Doxygen REQUIRED)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile @ONLY)
  file(INSTALL ${CMAKE_CURRENT_SOURCE_DIR}/doc/logo.png DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/doc)
  add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc
    COMMENT "Generating documentation with Doxygen" VERBATIM)
endif()
