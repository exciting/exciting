# This CMake script is responsible for creating the documentation for exciting.

# Option to enable documentation generation.
option(DOCUMENTATION "Generate documentation for the exciting project" OFF)

if (DOCUMENTATION)
  find_program(FORD_EXECUTABLE ford)
  if(FORD_EXECUTABLE)
    message(STATUS "Found ford: ${FORD_EXECUTABLE}")
  else()
    message(FATAL_ERROR "FORD executable not found. It is required for documentation.")
  endif()
  # Define the custom target to generate documentation
  add_custom_target(GenerateExcitingDocs ALL
      COMMAND ${CMAKE_COMMAND} -E echo "Calling ford to generate exciting documentation..."
      COMMAND ${CMAKE_COMMAND} -E env ford -o exciting_ford/ docs/ford_settings.md
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
      COMMENT "Generating documentation (1) for exciting subroutines"
  )
  
  # Add a post-build message to indicate completion
  add_custom_command(TARGET GenerateExcitingDocs
      POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E echo "Documentation generation complete."
  )

endif()

