#########################################################################################
# By Dave Favis-Mortlock, May 2018
#
# FindCShoreLib
# --------
#
# This module defines the following variables:
#
# CSHORELIB_INCLUDE_DIR - The libcshore include directories
# CSHORELIB_DIR - The libraries needed to use libcshore
#
# If CSHORE_LIBRARY_LC is "static", then this module also finds GFortran and defines:
# LIBGFORTRAN_LIBRARIES
# LIBQUADMATH_LIBRARIES
# GFORTRAN_INCLUDE_DIR
#
#########################################################################################
# Control verbosity
set(CMAKE_REQUIRED_QUIET true)
#set(CMAKE_REQUIRED_QUIET false)

if (CSHORE_LIBRARY_LC STREQUAL "static")
   set (CSHORE_LIB_NAME "libcshore.a")

   if (NOT CMAKE_REQUIRED_QUIET)
      message(STATUS "Looking for CShore library with static linkage (${CSHORE_LIB_NAME})")
   endif ()

else ()
   set (CSHORE_LIB_NAME "libcshore.so")

   if (NOT CMAKE_REQUIRED_QUIET)
      message(STATUS "Looking for CShore library with dynamic linkage (${CSHORE_LIB_NAME})")
   endif ()
endif ()

set (CSHORE_HEADER_NAME "cshore.h")

find_path(CSHORELIB_INCLUDE_DIR ${CSHORE_HEADER_NAME} HINTS ${CMAKE_SOURCE_DIR}/inc ${CMAKE_INSTALL_DIR}/inc PATH_SUFFIXES inc)
message (STATUS, "CSHORELIB_INCLUDE_DIR=${CSHORELIB_INCLUDE_DIR}")

find_library(CSHORELIB_DIR NAMES ${CSHORE_LIB_NAME} HINTS ${CMAKE_SOURCE_DIR}/lib ${CMAKE_INSTALL_DIR}/lib $ENV{LD_LIBRARY_PATH})
message (STATUS, "CSHORELIB_DIR=${CSHORELIB_DIR}")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CShoreLib DEFAULT_MSG CSHORELIB_DIR CSHORELIB_INCLUDE_DIR)

#message (STATUS, "CShoreLib_FOUND=${CShoreLib_FOUND}")

# if (CSHORELIB_FOUND)
#    message(STATUS "Found ${CSHORE_LIB_NAME} at ${CSHORELIB_DIR}")
# else ()
#   message(FATAL_ERROR "Could not find ${CSHORE_LIB_NAME}")
# endif ()

if (NOT CSHORELIB_FOUND)
  message(FATAL_ERROR "Could not find ${CSHORE_LIB_NAME}")
endif ()

# OK, all found. For neatness and consistency, also specify plural version of variables
set(CSHORELIB_DIRS ${CSHORELIB_DIR})
set(CSHORELIB_INCLUDE_DIRS ${CSHORELIB_INCLUDE_DIR})

mark_as_advanced(CSHORELIB_INCLUDE_DIR CSHORELIB_DIR)

if (CSHORE_LIBRARY_LC STREQUAL "static")
   # For static linkage, the CShore library needs GFortran on Linux TODO Windows, OSX?
   find_package (GFortranLibs REQUIRED)

#    message (STATUS, "GFORTRANLIBS_FOUND=${GFORTRANLIBS_FOUND}")
#    message (STATUS, "LIBGFORTRAN_LIBRARIES=${LIBGFORTRAN_LIBRARIES}")
#    message (STATUS, "LIBQUADMATH_LIBRARIES=${LIBQUADMATH_LIBRARIES}")
#    message (STATUS, "GFORTRAN_INCLUDE_DIR=${GFORTRAN_INCLUDE_DIR}")
#
#    if (GFORTRANLIBS_FOUND)
#       message(STATUS "Found GNU Fortran")
#    else ()
#       message(FATAL_ERROR "Could not find GNU Fortran")
#    endif ()
endif ()

#########################################################################################
