cmake_minimum_required (VERSION 2.6) 

if(NOT CMAKE_BUILD_TYPE)
  message("Type of build not specified. Use 'Release' by default.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra -std=gnu++0x")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}") # -Wconversion")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -march=native")


add_executable (hapchat_core
  entry.h
  entry.cpp
  new_columnreader.h
  new_columnreader.cpp
  blockreader.h
  blockreader.cpp
  combinations.h
  balanced_combinations.h
  balanced_combinations.cpp
  binomial.h
  binomial.cpp
  basic_types.h
  basic_types.cpp
  log.h
  hapchat_core.cpp
)


install(TARGETS hapchat_core RUNTIME DESTINATION bin)



if(NOT DEFINED DEFAULT_ENABLE_DISPLAY_REVISION)
  # Test whether the code is used in a repository if not autorevision will
  # fail and should be disabled by default. If inside a repository enable
  # the display of revision numbers by default.
  execute_process(
    COMMAND  ${CMAKE_SOURCE_DIR}/../scripts/autorevision.sh -th > ${CMAKE_CURRENT_BINARY_DIR}/revision.dummy
    WORKING_DIRECTORY  ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE  ENABLE_DISPLAY_REVISION_TEST_OUTPUT
    ERROR_VARIABLE   ENABLE_DISPLAY_REVISION_TEST_ERRNO
    )
  if("${ENABLE_DISPLAY_REVISION_TEST_ERRNO}" STREQUAL "")
    set(DEFAULT_ENABLE_DISPLAY_REVISION true)
  else()
    set(DEFAULT_ENABLE_DISPLAY_REVISION false)
  endif()
  unset(ENABLE_DISPLAY_REVISION_TEST_OUTPUT)
  unset(ENABLE_DISPLAY_REVISION_TEST_ERRNO)
endif()

option(
  ENABLE_DISPLAY_REVISION
  "Enable the display of the revision number while starting the program, only enable it when in a checkout"
  ${DEFAULT_ENABLE_DISPLAY_REVISION}
  )

# Depending on the flag create a real file or an empty dummy.
#
# Since the code is executed every build run of Wesnoth we need to make sure
# that no modifications don't turn into a rebuild. Therefore a the dummy is
# created and the second target only copies the file if modified.
if(ENABLE_DISPLAY_REVISION)
  add_custom_target(revision_dummy
    COMMAND  ${CMAKE_SOURCE_DIR}/../scripts/autorevision.sh -th > ${CMAKE_CURRENT_BINARY_DIR}/revision.dummy
    WORKING_DIRECTORY  ${CMAKE_SOURCE_DIR}
    )
  add_custom_target(revision
    COMMAND  ${CMAKE_COMMAND} -E
      copy_if_different
        ${CMAKE_CURRENT_BINARY_DIR}/revision.dummy
        ${CMAKE_CURRENT_BINARY_DIR}/revision.h
    )
  add_dependencies(revision revision_dummy)
  add_dependencies(hapchat_core revision)
  set_source_files_properties(
    hapchat_core.cpp
    PROPERTIES
    COMPILE_DEFINITIONS "LOAD_REVISION"
    )
  include_directories( ${CMAKE_CURRENT_BINARY_DIR} )
endif()
