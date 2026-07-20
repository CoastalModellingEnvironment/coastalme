#!/bin/bash

clang-tidy *.cpp -header-filter=.* -checks=-*,clang-analyzer-*,performance-*,-performance-avoid-endl,-clang-analyzer-core.NullDereference,misc-*  > 000_clang-analyze_advice.txt
clang-tidy *.h -header-filter=.* -checks=-*,clang-analyzer-*,performance-*,-performance-avoid-endl,-clang-analyzer-core.NullDereference,misc-*  >> 000_clang-analyze_advice.txt

# If get "error: 'omp.h' file not found [clang-diagnostic-error]" then need to install libomp-dev
# Some warnings seem crazy with openmp-*

# Run with misc-* occasionally, but note spurious errors especially re. headers (e.g. DO NOT remove "coast.h" #includes) and OpenMP

# Sometime, work on modernize-* results

# Not sure if readability-* results are useful
# Not sure if cppcoreguidelines-* results are useful
# Got no results from portability-*
