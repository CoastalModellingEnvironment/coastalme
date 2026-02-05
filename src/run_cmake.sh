#!/bin/bash
verbose='false'
cflag="false"
iflag="false"
while getopts 'civh' flag; do
	case "${flag}" in
	c) cflag="true" ;;
	i) iflag="true" ;;
	v) verbose='true' ;;
	h)
		printf "Usage: %s: [-c: clean] [-i: install] args\n" $0
		exit 2
		;;
	*)
		printf "Usage: %s [-c: clean] [-i: install] args\n" $0
		exit 2
		;;
	esac
done

# Change this to change build type
buildtype=DEBUG
#buildtype=RELEASE
#buildtype=PRERELEASE
#buildtype=RELWITHDEBINFO
#buildtype=MINSIZEREL            # Not yet implemented in CMakeLists.txt
#buildtype=GCOV
#buildtype=CALLGRIND

# Change this to select the Linux compiler
compiler=GNU
# compiler=CLANG

# Change this to select the CShore library type
#cshorelibrary=STATIC
cshorelibrary=SHARED

# Change this to select CShore input/output method
#cshoreinout=FILE
cshoreinout=ARG
#cshoreinout=BOTH

# Always build CShore
echo ""
cd cshore || exit
if [ "$cflag" = "true" ]; then
	rm -f ../lib/*
	make clean
	./make_cshore_lib.sh
fi

cd ..
# Note: The cshore Makefile now correctly names libraries for MacOS automatically
echo ""

# Now run CMake for CoastalME
echo ""
echo "================================================================="
echo ""

# On Mac, switch to gcc15 to ensure openMP compatability
# if [ "$OSTYPE" = "darwin"* ]; then
# 	# export CC=gcc-15
# 	# export CXX=g++-15
# 	# CMAKE_COMPILER_ARGS="-DCMAKE_C_COMPILER=gcc-15 -DCMAKE_CXX_COMPILER=g++-15"
# else
# 	#    echo "Using $compiler compiler, OpenMP may not be available"
# 	#    echo "For OpenMP support, install GCC or install libomp for Clang"
# 	CMAKE_COMPILER_ARGS=""
# fi
echo ""

if [ "$cflag" = "true" ]; then
	make clean
	rm -f CMakeCache.txt
fi
# CMAKE_COMPILER_ARGS="PKG_CPPFLAGS='-Xclang -fopenmp' PKG_LIBS=-lomp"
# cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=$buildtype -DCOMPILER=$compiler -DCSHORE_LIBRARY=$cshorelibrary -DCSHORE_INOUT=$cshoreinout $CMAKE_COMPILER_ARGS .
cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=$buildtype -DCOMPILER=$compiler -DCSHORE_LIBRARY=$cshorelibrary -DCSHORE_INOUT=$cshoreinout -DCMAKE_VERBOSE_MAKEFILE=OFF $CMAKE_COMPILER_ARGS .
#cmake -DCMAKE_BUILD_TYPE=$buildtype -DCOMPILER=$compiler -DCSHORE_LIBRARY=$cshorelibrary -DCSHORE_INOUT=$cshoreinout $CMAKE_COMPILER_ARGS . -G"CodeBlocks - Unix Makefiles"
# Or Ninja?
#cmake -G Ninja -DCMAKE_BUILD_TYPE=$buildtype -DCOMPILER=$compiler -DCSHORE_LIBRARY=$cshorelibrary -DCSHORE_INOUT=$cshoreinout $CMAKE_COMPILER_ARGS .

if [ "$iflag" = "true" ]; then
	make install
# 	ninja install #>output.txt
	if [[ $OSTYPE == 'darwin'* ]]; then
		# Let's sign to enable profiling
		codesign -s - -f --entitlements ../debug.plist ../cme
		echo "signed to enable profiling"
	fi
fi

echo ""
echo "================================================================="
echo ""

# Convert to upper case.
build_uc=$(echo $buildtype | tr 'a-z' 'A-Z')
compiler_uc=$(echo $compiler | tr 'a-z' 'A-Z')
cshorelibrary_uc=$(echo $cshorelibrary | tr 'a-z' 'A-Z')
cshoreinout_uc=$(echo $cshoreinout | tr 'a-z' 'A-Z')

echo "Finished CMake (Compiler=$compiler_uc, Build type=$build_uc, CShore library=$cshorelibrary_uc, CShore input/output method=$cshoreinout_uc)"
echo ""
echo "================================================================="

# Some extra messages
if [ "$buildtype" = "CALLGRIND" ]; then
	echo "When the build has finished, use valgrind/callgrind as follows:"
	echo ""
	echo "To check for memory leaks:"
	echo ""
	echo "   valgrind --leak-check=yes --suppressions=system-libs.supp --track-origins=yes ./cme &> valgrind.txt"
	echo ""
	echo "Then look at valgrind.txt"
	echo ""
	echo "Or to check coverage:"
	echo ""
	echo "   valgrind --tool=callgrind ./cme"
	echo ""
	echo "Then run:"
	echo ""
	echo "   callgrind_annotate --auto=yes callgrind.out.XXXXX > ./profile/callgrind/callgrind.txt"
	echo ""
	echo "where XXXXX is the number of the callgrind.out.XXXXX that was produced by valgrind. Then look at ./profile/callgrind.txt"
	echo ""
fi

if [ "$buildtype" = "GCOV" ]; then
	echo "When the build has finished, use gcov/lcov as follows:"
	echo ""
	echo "   ./cme"
	echo "   lcov --capture --directory ./src/CMakeFiles/cme.dir/ --output-file ./profile/lcov_output/coverage.info"
	echo "   cd ./profile/lcov_output"
	echo "   genhtml coverage.info"
	echo ""
	echo "Then look at index.html in your browser"
	echo ""
fi

if [ "$buildtype" = "PRERELEASE" ]; then
	echo "When the build has finished:"
	echo ""
	echo "   ./cme 2> sanitize.log"
	echo ""
	echo "Then look at sanitize.log"
	echo "DO NOT run using debugger e.g. gdb"
	echo "Note that we get a large number of GDAL leaks since sanitize-blacklist is not supported by gcc"
	echo ""
fi

#########################################################################################
