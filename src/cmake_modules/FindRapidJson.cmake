# Try to find rapidjson
#
# Sets the following variables:
# RAPIDJSON_FOUND
# RAPIDJSON_INCLUDE_DIRS - directories with rapidjson headers
# RAPIDJSON_DEFINITIONS - rapidjson compiler flags

find_path(rapidjson_header_paths_tmp
  NAMES
    document.h
  PATH_SUFFIXES
  include
  rapidjson/include
	PATHS
    ${RAPIDJSON_ROOT_DIR}
    ${RAPIDJSON_ROOT_DIR}/include
    ${RAPIDJSON_ROOT_DIR}/rapidjson/include
    ${RAPIDJSON_ROOT_DIR}/include/rapidjson
    $ENV{RAPIDJSON_ROOT_DIR}
    $ENV{RAPIDJSON_ROOT_DIR}/include
    $ENV{RAPIDJSON_ROOT_DIR}/rapidjson
    $ENV{RAPIDJSON_ROOT_DIR}/include/rapidjson
    )

get_filename_component(RAPIDJSON_INCLUDE_DIRS ${rapidjson_header_paths_tmp} PATH)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(rapidjson
  REQUIRED_VARS RAPIDJSON_INCLUDE_DIRS
  )

mark_as_advanced(RAPIDJSON_FOUND)

