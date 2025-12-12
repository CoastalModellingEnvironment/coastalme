################################################################################
#
# CMakeLists.txt for CoastalME(Coastal Modelling Environment)
#
################################################################################
# Set CMake required version
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

# This is needed if we want to use relative paths for e.g. the install dir.
# Such as ${CMAKE_SOURCE_DIR}/.. This was deprecated in CMake 3.31
if(${CMAKE_VERSION} VERSION_GREATER "3.31")
  cmake_policy(SET CMP0177 OLD)
endif()

# Create a .json file for use with clang-tidy
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

################################################################################
# Detect the build system
if(UNIX)
  if(NOT APPLE AND NOT CYGWIN)
    set(LINUX TRUE)
    set(SYSTEM "LINUX")
    message("Creating Linux-Unix Makefile for CoastalME")

    # Choose the Linux compiler
    string(TOLOWER "${COMPILER}" COMPILER_LC)
    if(COMPILER_LC STREQUAL "gnu")
      set(CMAKE_C_COMPILER gcc)
      set(CMAKE_CXX_COMPILER g++)
    elseif(COMPILER_LC STREQUAL "clang")
      set(CMAKE_C_COMPILER clang)
      set(CMAKE_CXX_COMPILER clang++)
    else()
      message(WARNING "Linux compiler not specified, falling back to gcc")
      set(CMAKE_C_COMPILER gcc)
      set(CMAKE_CXX_COMPILER g++)
    endif()

  elseif(CYGWIN)
    message("Creating Cygwin Makefile for CoastalME")

  elseif(APPLE)
    set(MACOS TRUE)
    set(SYSTEM "MACOS")
    message("Creating macOS-Unix Makefile for CoastalME")
  endif()

elseif(WIN32)
  message(
  FATAL_ERROR
  "CoastalME cannot be built using the Visual Studio compiler"
  )

elseif(MINGW)
  message(
  FATAL_ERROR
  "CoastalME cannot be builtusing the MinGW compiler"
  )
endif()

message("")

################################################################################
# First set up project-specific defaults
project(CME)

# Assume that we are in the src folder
set(CMAKE_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
set(CMAKE_INSTALL_DIR "${CMAKE_SOURCE_DIR}/..")

# May be changed by the user
# set(CMAKE_SOURCE_DIR "$ENV{HOME}/coast/CoastalME/src")
# set(CMAKE_INSTALL_DIR "$ENV{HOME}/coast/CoastalME")

# We are using C++
enable_language(CXX)

# find all of the cpp files
file(GLOB CME_SOURCE_FILES *.cpp)
set(CME_EXECUTABLE cme)

# Set the path for CoastalME's own cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

################################################################################
# Detect the compiler and set compiler-specific flags etc.
message("Checking compiler")

# We need c++17(for filesystem stuff, not that important howerver),
# so first set the flag
# (note that later on we check the compiler for c++17 capability)
add_compile_options(-std=c++17)

# Use 64-bit libraries
set(FIND_LIBRARY_USE_LIB64_PATHS true)

# Set the CoastalME CPU flag
# TODO Need to change this for different CPU architectures
if(APPLE)
  # assume only apple silicon exitists
  add_compile_options(-DCPU -Dxarm64)
else()
  add_compile_options(-DCPU -Dx86)
endif()

# Convert build type string to lower case
string(TOLOWER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_LC)
# Alert the user to their current compiler
message(
  "C++ compiler is ${CMAKE_CXX_COMPILER_ID} version ${CMAKE_CXX_COMPILER_VERSION}"
)

if(CMAKE_BUILD_TYPE_LC STREQUAL "debug")
  message("Applying universal debug flags")
  add_compile_options(
      -march=native
      -O0
      -fno-inline
      -W
      -Wall
      -Wextra
      -Wpointer-arith
      -Wconversion
      -Wcast-qual
      -Wcast-align
      -Wwrite-strings
      -Wredundant-decls
      -Wno-strict-overflow
      -Wshadow
      -Wuninitialized
      -Wnull-dereference
      -Wformat
      -Wformat-overflow
      -Wformat-signedness
      -Wempty-body
      -Wfloat-equal
  )
  if(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|LCC|QCC)$")
    message("Applying additional ${CMAKE_BUILD_TYPE_LC} flags for gcc")
    # ==========================================================================
    # Also use these with debug builds to check for vulnerabilities at compile time, see https://best.openssf.org/Compiler-Hardening-Guides/Compiler-Options-Hardening-Guide-for-C-and-C++.html
    # ==========================================================================
    add_compile_options(-Wuseless-cast)
    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "13.0")
      add_compile_options(
        # -O2
        -Wformat
        -Wformat=2
        -Wimplicit-fallthrough
        -Werror=format-security
        -U_FORTIFY_SOURCE
        -D_FORTIFY_SOURCE=3
        -D_GLIBCXX_ASSERTIONS
        -fstrict-flex-arrays=3
        -fstack-clash-protection
        -fstack-protector-strong
      )
      add_link_options(
        -Wl,-z,nodlopen
        -Wl,-z,noexecstack
        -Wl,-z,relro
        -Wl,-z,now
        -Wl,--as-needed
        -Wl,--no-copy-dt-needed-entries
      )
      # -Wtrampolines
      # -Wbidi-chars=any
      # -fPIE -pie
      #(for x86_64) -fcf-protection=full
    endif()

    #      add_compile_options(-Weffc++)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(Apple)?Clang$")
    message("Applying additional ${CMAKE_BUILD_TYPE_LC} flags for Clang")
  endif()

elseif(CMAKE_BUILD_TYPE_LC STREQUAL "relwithdebinfo")
  # Production build
  add_compile_options(
      -march=native
      -O3
      -g3
      -gdwarf-3
      -fopenmp
      # -foffload=nvidia-ptx
      -ffast-math
      -fno-omit-frame-pointer
      # -Wall
      -fno-delete-null-pointer-checks
      -fno-strict-overflow
      -fno-strict-aliasing
      -ftrivial-auto-var-init=zero
      -flto=auto
  )
  add_link_options(
      -flto=auto
      -rdynamic
  )
  if(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|LCC|QCC)$")
    # Production build using gcc
    message("Applying additional ${CMAKE_BUILD_TYPE_LC} flags for gcc")

  elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(Apple)?Clang$")
    message("Applying additional ${CMAKE_BUILD_TYPE_LC} flags for Clang")
  endif()

elseif(CMAKE_BUILD_TYPE_LC STREQUAL "prerelease")
  if(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|LCC|QCC)$")
    message("Applying additional ${CMAKE_BUILD_TYPE_LC} flags for gcc")
    # ========================================================================
    # Pre-release build using gcc, do this before a release to check for  memory access problems. Note that -fsanitize=address(and the other sanitize options?) DO NOT work if run under gdb
    # ========================================================================
    add_compile_options(-march=native)
    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "5.0")
      add_compile_options(
        -g
        -O0
        -fno-omit-frame-pointer
        -fsanitize=undefined
        -fsanitize=shift
        -fsanitize=integer-divide-by-zero
        -fsanitize=unreachable
        -fsanitize=null
        -fsanitize=return
        -fsanitize=vla-bound
        -fsanitize=float-divide-by-zero
        -fsanitize=float-cast-overflow
        -fsanitize=bounds
        -fsanitize=alignment
        -fsanitize=object-size
        -fsanitize=vptr
        -fsanitize=address
      )
      # Note  -fsanitize-ignorelist=sanitize_ignorelist.txt not yet working for gcc, need this to hide GDAL errors
    endif()

    # =======================================================================
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(Apple)?Clang$")
    message("Applying additional ${CMAKE_BUILD_TYPE_LC} flags for Clang")
  endif()

elseif(CMAKE_BUILD_TYPE_LC STREQUAL "release")
  # Production build
  add_compile_options(
      -march=native
      -O3
      -Wall
      -fno-delete-null-pointer-checks
      -fno-strict-overflow
      -fno-strict-aliasing
      -ftrivial-auto-var-init=zero
      -flto
  )
  add_link_options(
      -flto=auto
      -rdynamic
  )
  if(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|LCC|QCC)$")
    # Production build using gcc
    message("Applying additional ${CMAKE_BUILD_TYPE_LC} flags for gcc")

  elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(Apple)?Clang$")
    message("Applying additional ${CMAKE_BUILD_TYPE_LC} flags for Clang")
  endif()
endif()

# Finally set some gcc-specific linker flags, for all build types
if(UNIX)
  # add_link_options(-Wl)
  if(LINUX)
    # add_link_options(--as-needed -Wl,--no-undefined)
  elseif(APPLE)
  endif()
endif()

if(MSVC)
  message("Compiler is Visual Studio")
  add_compile_options(/W3)
endif()

message("")

################################################################################
# If not specified, assume a Debug build
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE_LC debug)
endif()

# Sort out compiler flags for various build types
if(CMAKE_BUILD_TYPE_LC STREQUAL "debug")
  set(CSHORE_BUILD "DEBUG")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS} ${COMPILE_OPTIONS}")
  add_link_options(-g)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CME_BINARY_DIR}/Debug CACHE PATH "Directory for Debug builds" FORCE)
  # For Debug builds on any platform, add the model's "internal" debug flag(used for extra logging) etc.
  add_compile_options(-D_DEBUG)

elseif(CMAKE_BUILD_TYPE_LC STREQUAL "relwithdebinfo")
  set(CSHORE_BUILD "DEBUG")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS}")
  add_link_options(-g)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CME_BINARY_DIR}/Debug CACHE PATH "Directory for Debug builds" FORCE)
  # For Debug builds on any platform, add the model's "internal" debug flag(used for extra logging) etc.
  add_compile_options(-D_DEBUG)

elseif(CMAKE_BUILD_TYPE_LC STREQUAL "release")
  set(CSHORE_BUILD "RELEASE")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS}")
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CME_BINARY_DIR}/Release CACHE PATH "Directory for Release builds" FORCE)

elseif(CMAKE_BUILD_TYPE_LC STREQUAL "gcov")
  set(CSHORE_BUILD "DEBUG")
  #   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS}")
  add_compile_options(-O0 --coverage)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CME_BINARY_DIR}/gcov CACHE PATH "Directory for gcov builds" FORCE)

elseif(CMAKE_BUILD_TYPE_LC STREQUAL "valgrind")
  set(CSHORE_BUILD "DEBUG")
  #   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS}")
  add_compile_options(-pg)
  add_link_options(-pg -v)
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg -v")
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CME_BINARY_DIR}/Callgrind CACHE PATH "Directory for Callgrind builds" FORCE)
endif()

################################################################################
# Next, find libraries
message("Finding libraries")

################################################################################
# Find OpenMP for parallelization. First try to find OpenMP using the standard CMake module
find_package(OpenMP)

if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP found via find_package: ${OpenMP_CXX_VERSION}")
  set(LIBS ${LIBS} OpenMP::OpenMP_CXX)
  message(STATUS "OpenMP enabled for parallel grid operations")
else()

  # Reset CMAKE_REQUIRED variables
  unset(CMAKE_REQUIRED_FLAGS)
  unset(CMAKE_REQUIRED_LIBRARIES)
  message(WARNING "OpenMP not found - grid operations will run sequentially")
endif()

################################################################################
# CoastalME cannot be built without GDAL, which is linked dynamically
find_package(GDAL REQUIRED)
if(GDAL_FOUND)
  # OK, we found GDAL
  message(STATUS "GDAL_INCLUDE_DIRS=${GDAL_INCLUDE_DIRS}")
  message(STATUS "GDAL_LIBRARIES=${GDAL_LIBRARIES}")
  message(STATUS "GDAL_VERSION=${GDAL_VERSION}")

  # Now check the GDAL version: CoastalME can only be built with GDAL 2.1 or above
  if(GDAL_VERSION VERSION_LESS "2.1")
    message(FATAL_ERROR "CoastalME requires GDAL 2.1 or later. You have version ${GDAL_VERSION}. Please upgrade")
  endif()

  set(LIBS ${LIBS} ${GDAL_LIBRARIES})
  set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${GDAL_INCLUDE_DIRS})
  message(STATUS "LIBS=${LIBS}")
  message(STATUS "CMAKE_INCLUDE_PATH=${CMAKE_INCLUDE_PATH}")
endif()

################################################################################
# CoastalME also requires the CShore library. This can be linked statically or dynamically(i.e. shared). If not specified, assume a shared library
if(NOT CSHORE_LIBRARY)
  set(CSHORE_LIBRARY_LC "shared")
else()
  # The library link type was specified, so create a lower-case version
  string(TOLOWER "${CSHORE_LIBRARY}" CSHORE_LIBRARY_LC )
endif()

# If not specified, assume that the CShore library communicates with CoastalME via argument passing
if(NOT CSHORE_INOUT)
  set(CSHORE_INOUT_LC "arg")
else()
  # The communication method was specified, so make a lower-case version
  string(TOLOWER "${CSHORE_INOUT}" CSHORE_INOUT_LC )
endif()

if(CSHORE_INOUT_LC STREQUAL "file")
  set(CSHORE_LINKTYPE "FILEINOUT")
  add_compile_options(-DCSHORE_FILE_INOUT)
elseif(CSHORE_INOUT_LC STREQUAL "arg")
  set(CSHORE_LINKTYPE "ARGINOUT")
  add_compile_options(-DCSHORE_ARG_INOUT)
elseif(CSHORE_INOUT_LC STREQUAL "both")
  set(CSHORE_LINKTYPE "ARGINBOTHOUT")
  add_compile_options(-DCSHORE_BOTH)
else()
  message(FATAL_ERROR "Invalid value specified for communication with the CShore library: ${CSHORE_INOUT}")
endif()

# set link type
# set buildtype
if(CSHORE_LIBRARY_LC STREQUAL "shared")
  set(EXT "so")
else()
  set(EXT "a")
endif()

set(CSHORE_FILE "${CMAKE_SOURCE_DIR}/lib/libcshore.${EXT}.${CSHORE_BUILD}.${SYSTEM}.${CSHORE_LINKTYPE}")
message(STATUS "CSHORE_FILE=${CSHORE_FILE}")
execute_process(COMMAND ln -s -f "${CSHORE_FILE}" "${CMAKE_SOURCE_DIR}/lib/libcshore.${EXT}")

if(WIN32)
  # TODO
endif()

# OK, now find the CShore library (and GFortran, if the CShore linkage is static)
find_package(CShoreLib REQUIRED)

#message(STATUS "CSHORELIB_DIR=${CSHORELIB_DIR}")
set(LIBS ${LIBS} ${CSHORELIB_DIR})
message(STATUS "CSHORELIB_INCLUDE_DIR=${CSHORELIB_INCLUDE_DIR}")
set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${CSHORELIB_INCLUDE_DIR})

if(CSHORE_LIBRARY_LC STREQUAL "static")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
  set(LIBS ${LIBS} ${LIBGFORTRAN_LIBRARIES} ${LIBQUADMATH_LIBRARIES})
  set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${GFORTRAN_INCLUDE_DIR})
     # include_directories(${GFORTRAN_INCLUDE_DIR})
     # target_link_libraries(${CME_EXECUTABLE} ${LIBS} ${LIBGFORTRAN_LIBRARIES})
     # target_link_libraries(${CME_EXECUTABLE} ${LIBS} ${LIBQUADMATH_LIBRARIES})
else()
  # Assume that libcshore.so is either in one of the system library dirs, or is in a location specified by export LD_LIBRARY_PATH
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

  # We have the absolute path to the CShore library, now get the relative path from the install dir
  file(RELATIVE_PATH CSHORELIB_REL_FILEPATH ${CMAKE_INSTALL_DIR} ${CSHORELIB_DIR})

  # Remove the file name
  get_filename_component(CSHORELIB_REL_PATH ${CSHORELIB_REL_FILEPATH} DIRECTORY)

  # And get the absolute path(without the file name)
  get_filename_component(CSHORELIB_ABS_PATH ${CSHORELIB_DIR} DIRECTORY)

  # Get GDAL library directory
  get_filename_component(GDAL_LIB_DIR ${GDAL_LIBRARIES} DIRECTORY)

  # Set platform-specific RPATH configuration
  if(APPLE)
    # macOS uses @loader_path instead of $ORIGIN
    set(CMAKE_INSTALL_RPATH "@loader_path/${CSHORELIB_REL_PATH}:${CSHORELIB_ABS_PATH}:${GDAL_LIB_DIR}:/opt/homebrew/lib:/usr/local/lib")
  else()
    # Linux and other Unix systems use $ORIGIN
    set(CMAKE_INSTALL_RPATH "$ORIGIN/${CSHORELIB_REL_PATH}:$$ORIGIN/${CSHORELIB_REL_PATH}:${CSHORELIB_ABS_PATH}:${GDAL_LIB_DIR}:/usr/local/lib")
  endif()
endif()

################################################################################
# The important bits
include_directories(SYSTEM ${CMAKE_INCLUDE_PATH})
include_directories(${CME_SOURCE_DIR})
add_executable(${CME_EXECUTABLE} ${CME_SOURCE_FILES})
target_link_libraries(${CME_EXECUTABLE} ${LIBS})

install(TARGETS ${CME_EXECUTABLE} RUNTIME DESTINATION ${CMAKE_INSTALL_DIR})

# Check to see if the compiler supports c++17
set_property(TARGET ${CME_EXECUTABLE} PROPERTY CXX_STANDARD 17)
set_property(TARGET ${CME_EXECUTABLE} PROPERTY CXX_STANDARD_REQUIRED ON)

################################################################################

################################################################################
# Tell the user what we have found
message("")
message("Writing Makefile")
message(STATUS "CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
message(STATUS "CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}")
message(STATUS "CMAKE_EXE_LINKER_FLAGS=${CMAKE_EXE_LINKER_FLAGS}")
message(STATUS "CMAKE_SOURCE_DIR = ${CMAKE_SOURCE_DIR}")
message(STATUS "CMAKE_INSTALL_DIR = ${CMAKE_INSTALL_DIR}")
# message(STATUS "CMAKE_RUNTIME_OUTPUT_DIRECTORY=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
message(STATUS "CMAKE_INCLUDE_PATH=${CMAKE_INCLUDE_PATH}")
message(STATUS "LIBS=${LIBS}")
# message(STATUS "CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}")
# message(STATUS "CMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}")
message(STATUS "CMAKE_INSTALL_RPATH=${CMAKE_INSTALL_RPATH}")
################################################################################

