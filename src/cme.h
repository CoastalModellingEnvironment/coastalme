/*!
   \mainpage
   \section intro_sec Introduction

   <b>CoastalME</b> (Coastal Modelling Environment) simulates the long-term behaviour of a coast. This initial version considers only simple soft cliff cross-shore effects. However, development of CoastalME is ongoing. Watch this space!\n\n

   CoastalME was devised and constructed by Andres Payo Garcia (British Geological Survey: agarcia@bgs.ac.uk) and David Favis-Mortlock (British Geological Survey: dfm1@bgs.ac.uk). We are very grateful to the following for support, assistance, and inspiration: Tom Ashby, Manuel Cobos Budia, Wilf Chun, Mark Dickson, Jim W. Hall, Martin D. Hurst, Matthew Ives, Robert J. Nicholls, Ian Townend, and Mike
   J.A. Walkden.\n\n

   See <a href="https://github.com/coastalme/coastalme" target="_blank">https://github.com/coastalme/coastalme</a> for the stable release version, and the unstable development version, of the source code.\n
   \n
   From Shingle Street\n
   To Orford Ness\n
   The waves maraud,\n
   The winds oppress,\n
   The earth can’t help\n
   But acquiesce\n
   For this is east\n
   And east means loss,\n
   A lessening shore, receding ground,\n
   Three feet gone last year, four feet this\n
   Where land runs out and nothing’s sound.\n
   Nothing lasts long on Shingle Street.\n
   \n
   By Blake Morrison (2018). See <a href="https://www.penguin.co.uk/books/419911/shingle-street-by-morrison-blake/9780701188771" target="_blank">https://www.penguin.co.uk/books/419911/shingle-street-by-morrison-blake/9780701188771</a>\n

   \section install_sec Installing CoastalME

   \subsection install_step1 Obtaining the source code

   CoastalME builds easily using Linux. If you wish to run CoastalME on Windows, then we currently recommend using the Windows Subsystem Linux (WSL) software to do this.

   Create a local copy of the github repository, for example by downloading a zipfile, then unpacking it or cloning. We suggest unpacking it to something like "/home/YOUR NAME/Projects/CoastalME/", this is then your CoastalME folder.

   git clone https://github.com/coastalme/coastalme

   \subsection install_step2 Building CoastalME

   In a terminal window (i.e. at a command-line prompt) move to the CoastalME folder.

   Then move to the the src folder

   cd CoastalME/src

   and run run_cmake.sh

   ./run_cmake.sh

   If you get a "Permission denied" message: -bash: ./run_cmake.sh: Permission denied you will have to grant executable permission using chmod a+x run_cmake.sh, chmod a+x cshore/make_cshore.sh and then
   ./run_cmake.sh

   This will build CShore, look for GDAL, and write the CMake files. If you see error messages about missing software (for example, telling you that CMake cannot be found or is too old, or GDAL cannot be found or is too old) then you need to install or update the software that is causing the problem.

   Next, run

   make install

   This will create an executable file called cme in the CoastalME folder.

   \section run_sec Running CoastalME

   \subsection run_step1 Specifying the input file

   Edit cme.ini to tell CoastalME which input file to read (for example, in/test_suite/minimal_wave_angle_230/minimal.dat).

   \subsection run_step2 Running CoastalME

   Leave the src folder, and run cme

   cd ..
   ./cme

   Output will appear in the "Path for output" folder.

   \subsection run_step3 Running CoastalME's test suite

   To check that your installation is running correctly, you can run a suite of pre-defined tests by running the following commands:

   chmod a+x run_test_suite.sh
   ./run_test_suite.sh

   The `chmod` comand ensures that you have permission to execute the run_test_suite.sh file.

   \subsection run_step4 Managing CoastalME's output

   Once you have CoastalME (CME) up and running, you can reduce the quantity of output (it can be overwhelming!) in several ways.

   Change "Content of log file" in the CME input file for any of the test suite runs (the name of this input file is listed in cme.ini, both are simple text files). If you set "Content of log file" to zero, then CME won't output a log file; setting it to 4 (all output) is really only useful to developers.

   Change "GIS vector files to output" and "GIS vector files to output" in the CME input file. These are both set to "all" in the test suite files on GitHub. Instead of "all" you can list the space-separated codes for only the GIS output that you want to see. A list of CME GIS output codes is in codes.txt.

   Enjoy!

   \file cme.h
   \brief This file contains global definitions for CoastalME
*/

/*
   NOTE Before releasing a new version, do a pre-release build to check for memory leaks with -fsanitize options enabled (see CMakeLists.txt) then run ./cme 2> sanitize.txt NOT UNDER DEBUG (i.e. not using gdb)

   TODOLIST
 ***********************************************************************************************************
   DOCUMENTATION
   TODO 001 Add more Doxygen information about all classes
   TODO 007 We now have setup and surge info from CShore (thanks to Manuel). But what shall we do with this info? "The variable VdWaveSetupSurge() represents the sea level rise due to wave effects (setup) and storm surge. CSHORE calculates them together and they can’t be separated. That’s what the VdWaveSetupSurge variable is. That’s why you saw my initial efforts to try to separate both variables from CSHORE commented out, which is impossible. What is possible is to get the RunUp from CSHORE, but since it uses an empirical formula for that, I finally decided to calculate it separately. To your question about whether you should remove VdStormSurge, the answer is yes. I left it because I still intend at some point to extract the cross-shore transport from CSHORE and balance it in CME with the longshore and cross-shore transports without needing the Dean profile. From my point of view, this would be even more realistic, though at first it will surely drive us crazy."

   USER INPUT
   TODO 000 Should user input be split in two main files: one for frequently-changed things, one for rarely-changed things? If so, what should go into each file ('testing only' OK, but what else?)
   TODO 011 Should this constant be a user input? If so, TODO 071
   TODO 036 Read in changed deep water wave values (need TODO 071)
   TODO 022 Get intervention update working (need TODO 071)
   TODO 042 Should we have a smallest valid input for KLS in the CERC equation?
   TODO 045 Method of getting depth of closure value needs to be a user input (need TODO 071)
   TODO 049 Handle other command line parameters e.g. path to .ini file, path to datafile
   TODO 035 Also handle other EPSG for vector spatial reference systems
   TODO 054 Choose more files to omit from "usual" raster output
   TODO 069 Enable ability to represent intervention structures which have their foundation embedded in consolidated sediment. In other words, with the elevation of the base of the intervention structure *below* the top of all consolidated sediment layers. Will need some sanity checking of elevations
   TODO 071 If the user input file format is changed, write a Python script to convert from the old file format to the new
   TODO 083 Get all three kinds of sediment input events working correctly

   ERROR HANDLING
   TODO 038 Do better error handling if insufficient memory
   TODO 004 Improve error handling of situation where we have a valid shadow zone but cannot find a neighbouring cell which is 'under' the coastline
   TODO 006 Check GDALGridCreate() with only start-of-coast or an end-of-coast profiles
   TODO 009 Decide what to do when we have eroded down to basement
   TODO 017 Extra safety check needed, make sure that each point is within valid grid
   TODO 018 Improve situation where new landwards point on parallel profile is not within the raster grid
   TODO 019 Improve situation where Dean profile has a near-zero elevation difference
   TODO 020 Check calculation of elevation of coast point of Dean parallel profile
   TODO 021 Improve situation where all layers have zero thickness
   TODO 025 Improve situation where this point has only zero thickness layers
   TODO 026 Check situation where cell in parallel profile is not in a polygon
   TODO 028 Give a warning if raster input layer has several bands
   TODO 053 Improve handling of situation where landward elevation of profile is -ve
   TODO 055 Maybe add a safety check here?
   TODO 080 Do we get -ve breaking wave heights here?
   TODO 084 Improve handling of situation where consecutive profile points are same distance from shoreline

   THEORY/EFFICIENCY
   TODO 002 Do we really need D50 for drift landform class? What do we need for drift?
   TODO 005 Maybe give every coast point a value for end-of-profile wave height and direction instead of for deep water wave height and direction
   TODO 010 Do we also need to update the active zone cells?
   TODO 012 Change finding of adjacent polygons, and calculation of the length of shared normals, when we make polygon seaward length determined by depth of closure
   TODO 013 Change calculation (need user input?) of coastline smoothing convexity threshold
   TODO 014 Profile spacing, could try gradually increasing the profile spacing with increasing concavity, and decreasing the profile spacing with increasing convexity
   TODO 016 Check mass balance for recirculating unconsolidated sediment option
   TODO 024 Should we calculate platform erosion on a profile that has hit dry land?
   TODO 044 Implement estuaries
   TODO 051 Implement other ways of calculating depth of closure, see TODO 045
   TODO 056 Check this please Andres
   TODO 059 Implement dune landform class
   TODO 060 Remove 'magic numbers' from code here
   TODO 061 Is this safety check to depth of breaking a reasonable thing to do?
   TODO 066 Should this be for all layers? Check
   TODO 067 Suspended fine sediment never decreases i.e. no suspended fine sediment ever leaves the grid. Is this OK?
   TODO 070 Change CShore to use allocatable arrays (https://fortran-lang.org/en/learn/best_practices/allocatable_arrays/) so that the number of points in the CShore output profiles can either be a user input, or determined by e.g. the physical length of the profile. At present, max is NN = 1000 in cshore_wrapper.f03
   TODO 075 What if bedrock sticks above Dean profile?
   TODO 076 When doing parallel profiles, start from the profile which is closest to a right angle with the coast
   TODO 077 As traverse between the bounding profiles creating parallel profiles, gradually change the parallel profile orientation based on distance weighting of two bounding profiles
   TODO 078 At present, we don't allow cliff collapse onto interventions. Is this realistic? Should it be different for different types on intervention?
   TODO 089 Why do we get patches of sediment in the sea?
   TODO 086 Try these as a more efficient replacement for GDALGridCreate(): https://github.com/delfrrr/delaunator-cpp https://www.cs.cmu.edu/~quake/triangle.html https://github.com/greenm01/poly2tri https://gts.sourceforge.net/index.html
   TODO 088 In (almost) all whole-grid loops, immediately continue if cell is hinterland (but not when calculating cliff collapse)
   TODO 090 At present, sediment cannot move from a given coastline polygon to a polygon belonging to another coastline. Is this always true?
   TODO 092 If we have only fine sediment, the surface formed as the coast recedes inland is dead level (because fine sediment goes to suspension, and hence the Dean profile stuff does not operate). This causes problems with profile creation and CShore. Need to impose a small slope here somehow
   TODO 093 There are a number of cell attributes that are really only useful for debugging. To keep memory usage down on release versions, need to flag these attributes and their methods so that they are included only in debug versions
   TODO 094 Problems with sediment recirculation when large volumes of sedimentare introduced at the input end of the coast: need to spread this input sediment over the whole of the grid-end polygon
   TODO 095 Parallel profiles seem to create "streaks" of above-water sediment, these streaks interfere with subsequent profile creation (profiles which hit steaks arte marked as invalid). Investigate this
   TODO 096 It is OK to specify OGRFieldDefn objects as const in recent versions (e.g. 14.2) of g++, but in older versions of g++ (e.g. 12.2) have problems with const here. Maybe make the const a pre-processor condion?

   OUTPUT
   TODO 065 Get GPKG output working: GDAL 3.9.1 does not yet implement this correctly. Currently is OK for vector output (but is very slow), not yet working for raster output
   TODO 063 Add NetCDF support, see https://trac.osgeo.org/gdal/wiki/NetCDF
   TODO 064 Add support for grids that are not oriented N-S and W-E, but which are still rectangular. See https://gdal.org/en/stable/tutorials/geotransforms_tut.html
   TODO 031 Get raster slice output working with multiple slices
   TODO 032 Improve output scaling for DBL_NODATA situation
   TODO 033 Also test and configure (e.g. by passing open() options) other vector output file formats
   TODO 034 Also test and configure (e.g. by passing open() options) other raster output file formats
   TODO 043 When outputting profiles, how do we deal with randomness of profile spacing (since profile location is determined by curvature)?
   TODO 052 Improve saving of profiles and parallel profiles
   TODO 062 Show end-of-iteration number of cells with sediment somewhere
   TODO 068 Only show output in log file that is relevant to processes being simulated
   TODO 074 Output history of what landforms are on a particular cell or cells. User inputs cell(s), how?
   TODO 082 Also show m_dStartIterUnconsFineAllCells etc. in log file

   096 is max

   COMPLETED
   TODO 003 Make coastline curvature moving window size a user input DONE in 1.1.22
   TODO 046 Why is cliff collapse eroded during deposition (three size classes) no longer calculated? DONE IN 1.1.22
   TODO 058 Dave to check this DONE in 1.1.22
   TODO 039 Rewrite reading of multiple random number seeds DONE in 1.2.1, 8 Nov 2024
   TODO 041 Read in SWL per-timestep
   BUG 002 Useless output e.g. clay layers even if no clay input DONE in 1.1.21
   BUG 003 Use mean SWL for elevations of Dean profiles DONE in 1.2.1, 27 Nov 2024
   BUG 004 Don't smooth intervention coastline DONE 1.2.1, 27 Nov 2024
   TODO 073 If output dir does not exist, then create it (ask user first) DONE 1.2.2, 28 Nov 2024
   TODO 047 Where is the GDAL description for the deep water wave stations vector file? DONE 1.2.3, 2 Dec 2024
   TODO 048 Where is the GDAL description for the flood input locations point or vector file? DONE 1.2.3, 2 Dec 2024
   TODO 027 Sort out GDAL problem with raster reference units DONE 1.2.3, 2 Dec 2024
   TODO 079 Do sanity checking on wave and tide input DONE 1.2.3, 2 Dec 2024
   TODO 072 CShore crashes occasionally, is it because of -ve Z values here? DONE 1.2.3, 2 Dec 2024
   TODO 050 Update for recent versions of Windows DONE 1.2.3, 2 Dec 2024
   TODO 037 Need more info on nFindIndex() DONE 1.2.3, 2 Dec 2024 Improve coast normals DONE 1.2.3, 20 Dec 2024
   TODO 057 Check this please Manuel DONE 1.2.4, 4 Jan 2025
   TODO 087 Is there a problem if profile is not long enough for user-input depth of closure? DONE 1.3.0 2 Feb 2025
   TODO 023 Only calculate shore platform erosion if cell is in a polygon
*/

#ifndef CME_H
#define CME_H
/* ===============================================================================================================================
   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
===============================================================================================================================*/
#include <csignal>      // Only for programmatic debug breakpoints using gdb i.e. std::raise(SIGINT);

#include <climits>

#include <string>
using std::string;

#include <sstream>
using std::ostringstream;

#include <ostream>
using std::ostream;

#include <vector>

//===================================================== platform-specific stuff =================================================
#ifdef _WIN32
   #define access _access
   #define F_OK 0 // Test for file existence
#endif

#ifdef _MSC_VER
   // MS Visual C++ compiler, byte order is IEEE little-endian
   #ifdef _DEBUG
      #include <crtdbg.h> // useful
   #endif

   // clock_t is a signed long: see <time.h>
   long const CLOCK_T_MIN = LONG_MIN;
   double const CLOCK_T_RANGE = static_cast<double>(LONG_MAX) - static_cast<double>(CLOCK_T_MIN);

   #ifdef _M_ALPHA
      string const PLATFORM = "MS Visual C++ for Alpha";
   #elif defined _M_IX86
      string const PLATFORM = "MS Visual C++ for Intel x86";
   #elif defined _M_MPPC
      string const PLATFORM = "MS Visual C++ for Power PC";
   #elif defined _M_MRX000
      string const PLATFORM = "MS Visual C++ for MIPS";
   #else
      string const PLATFORM = "MS Visual C++ for unknown CPU";
   #endif

#elif defined __GNUG__
   // GNU compiler
   #ifndef CPU
      #error "CPU not defined"
   #else
      #ifdef x86
         // Intel x86, byte order is little-endian
         string const PLATFORM = "GNU Compiler for Intel x86";
         // clock_t is an unsigned long: see <time.h>
         unsigned long const CLOCK_T_MIN = 0;
         double const CLOCK_T_RANGE = static_cast<double>(ULONG_MAX);
      #elif defined rs6000
         // IBM RS-6000, byte order is big-endian
         string const PLATFORM = "GNU complier for IBM RS-6000";
         // clock_t is a signed long: see <time.h> NEED TO CHECK
         long const CLOCK_T_MIN = LONG_MIN;
         double const CLOCK_T_RANGE = static_cast<double>(LONG_MAX) - static_cast<double>(CLOCK_T_MIN);
      #elif defined ultrasparc
         // Sun UltraSparc, byte order is big-endian
         string const PLATFORM = "GNU compiler for Sun UltraSPARC";
         // clock_t is a signed long: see <time.h>
         long const CLOCK_T_MIN = LONG_MIN;
         double const CLOCK_T_RANGE = static_cast<double>(LONG_MAX) - static_cast<double>(CLOCK_T_MIN);
      #else
         // Something else
         string const PLATFORM = "GNU compiler for unknown CPU";
         // clock_t is a signed long: NEED TO CHECK <time.h>
         long const CLOCK_T_MIN = LONG_MIN;
         double const CLOCK_T_RANGE = static_cast<double>(LONG_MAX) - static_cast<double>(CLOCK_T_MIN);
      #endif
   #endif

#elif defined __clang__
   // Clang compiler
   #ifndef CPU
      #error "CPU not defined"
      #else
         #ifdef x86
            // Intel x86, byte order is little-endian
            string const PLATFORM = "Clang compiler for Intel x86";
            // clock_t is an unsigned long: see <time.h>
            unsigned long const CLOCK_T_MIN = 0;
            double const CLOCK_T_RANGE = static_cast<double>(ULONG_MAX);
         #else
            // Something else
            string const PLATFORM = "Clang compiler for unknown CPU";
            // clock_t is a signed long: NEED TO CHECK <time.h>
            long const CLOCK_T_MIN = LONG_MIN;
            double const CLOCK_T_RANGE = static_cast<double>(LONG_MAX) - static_cast<double>(CLOCK_T_MIN);
         #endif
      #endif

#elif defined __MINGW32__
   // Minimalist GNU for Windows
   //   #define __USE_MINGW_ANSI_STDIO 1        // Fix long doubles output problem,
   //   see http://stackoverflow.com/questions/7134547/gcc-printf-and-long-double-leads-to-wrong-output-c-type-conversion-messes-u
   #define WEXITSTATUS(x) ((x) & 0xff)

#elif defined __HP_aCC
   // HP-UX aCC, byte order is big-endian, can be either 32-bit or 64-bit
   string const PLATFORM = "HP-UX aC++";
   // clock_t is an unsigned long: see <time.h>
   unsigned long const CLOCK_T_MIN = 0;
   #ifdef __ia64
      // However, clock_t is a 32-bit unsigned long and we are using 64-bit unsigned longs here
      double const CLOCK_T_RANGE = 4294967295UL; // crude, improve
   #else
      double const CLOCK_T_RANGE = static_cast<double>(ULONG_MAX);
   #endif
#endif

// TODO: Check
// #if defined(WIN32)
// #define STRCASECMP(a, b) (_stricmp(a, b))
// #define STRNCASECMP(a, b, n) (_strnicmp(a, b, n))
// #else
// /** Alias for strcasecmp() */
// #define STRCASECMP(a, b) (strcasecmp(a, b))
// /** Alias for strncasecmp() */
// #define STRNCASECMP(a, b, n) (strncasecmp(a, b, n))
// //#  endif
// /** Alias for strncasecmp() == 0 */
// #define EQUALN(a, b, n) (STRNCASECMP(a, b, n) == 0)
// /** Alias for strcasecmp() == 0 */
// #define EQUAL(a, b) (STRCASECMP(a, b) == 0)
// #endif

//===================================================== hard-wired constants ====================================================
char const COLON = ':';
char const COMMA = ',';
char const DASH = '-';
char const PATH_SEPARATOR = '/'; // Works for Windows too!
char const QUOTE1 = ';';
char const QUOTE2 = '#';
char const SLASH = '/';
char const SPACE = ' ';
char const TILDE = '~';

// TESTING options
bool const ACCEPT_TRUNCATED_PROFILES = true;
bool const CREATE_SHADOW_ZONE_IF_HITS_GRID_EDGE = true;  // If shadow line tracing hits grid edge, create shadow zone?
bool const SAVE_CSHORE_OUTPUT = true;                    // If CSHORE_FILE_INOUT or CSHORE_BOTH, append all CShore output files to a whole-run master
bool const USE_DEEP_WATER_FOR_SHADOW_LINE = true;        // Use deep water wave orientation in determining shadow line orientation?

// Not likely that user will need to change these
int const NUMBER_OF_RNGS = 2;                            // Number of random number generators
int const SAVEMAX = 100000;                              // Maximum number of saves of spatial output
int const BUF_SIZE = 2048;                               // Max length (inc. terminating NULL) of any C-type string
int const CAPE_POINT_MIN_SPACING = 10;                   // In cells: for shadow zone stuff, cape points must not be closer than this
int const CLOCK_CHECK_ITERATION = 5000;                  // If have done this many timesteps then reset the CPU time running total
int const COAST_LENGTH_MAX_CONST = 3;                    // When tracing coast, quit if coast length > this * smallest grid edge size
int const COAST_LENGTH_MIN_CONST = 1;                    // Ignore very short coasts less than this * profile spacing

//! The size of the arrays output by CShore. If this is changed, then must also set the same value on line 12 of cshore_wrapper.f03 (integer, parameter :: NN = 1000, NL = 1) and recompile CShore. Eventually we should move to dynamically allocated arrays TODO 070
int const CSHORE_ARRAY_OUT_SIZE = 1000;

int const FLOOD_FILL_START_OFFSET = 2;                   // In cells: cell-by-cell fill starts this distance inside polygon
int const GRID_MARGIN = 10;                              // Ignore this many along-coast grid-edge points re. shadow zone calcs
int const INT_NODATA = -9999;                            // CME's internal NODATA value for ints
int const MAX_LEN_SHADOW_LINE_TO_IGNORE = 200;           // In cells: if can't find cell-by-cell fill start point, continue if short shadow line
int const MAX_NUM_PREV_ORIENTATION_VALUES = 10;          // Max length of deque used in tracing shadow boundary
int const MAX_NUM_SHADOW_ZONES = 10;                     // Consider at most this number of shadow zones
int const MIN_INLAND_OFFSET_UNCONS_EROSION = 5;          // Used in estimation of beach erosion
int const MIN_PARALLEL_PROFILE_SIZE = 3;                 // In cells: min size for valid unconsolidated sediment parallel profile
int const MIN_PROFILE_SIZE = 3;                          // In cells: min size for valid unconsolidated sediment profile
int const DEFAULT_PROFILE_SPACING = 20;                  // In cells: profile creation does not work well if profiles are too closely spaced
int const SAVGOL_POLYNOMIAL_MAX_ORDER = 6;               // Maximum order of Savitzky-Golay smoothing polynomial

// Log file detail level
int const NO_LOG_FILE = 0;
int const LOG_FILE_LOW_DETAIL = 1;
int const LOG_FILE_MIDDLE_DETAIL = 2;
int const LOG_FILE_HIGH_DETAIL = 3;
int const LOG_FILE_ALL = 4;

// Direction codes
int const NO_DIRECTION = 0;
int const NORTH = 1;
int const NORTH_EAST = 2;
int const EAST = 3;
int const SOUTH_EAST = 4;
int const SOUTH = 5;
int const SOUTH_WEST = 6;
int const WEST = 7;
int const NORTH_WEST = 8;

int const DIRECTION_DOWNCOAST = 0;     // Down-coast, i.e. along the coast so that the index of coastline points INCREASES
int const DIRECTION_UPCOAST = 1;       // Up-coast, i.e. along the coast so that the index of coastline points DECREASES

// Handedness codes, these show which side the sea is on when travelling down-coast (i.e. in the direction in which coastline point numbers INCREASE)
int const NULL_HANDED = -1;
int const RIGHT_HANDED = 0;
int const LEFT_HANDED = 1;

// Sediment texture codes
int const TEXTURE_FINE = 0;
int const TEXTURE_SAND = 1;
int const TEXTURE_COARSE = 2;

// Time unit codes
int const TIME_UNKNOWN = -1;
int const TIME_HOURS = 0;
int const TIME_DAYS = 1;
int const TIME_WEEKS = 2;
int const TIME_MONTHS = 3;
int const TIME_YEARS = 4;

// Intervention input and output codes
int const IO_INTERVENTION_NONE = 0;
int const IO_INTERVENTION_STRUCT = 1;
int const IO_INTERVENTION_NON_STRUCT = 2;

// Landform category codes
int const LF_UNKNOWN = 0;
int const LF_HINTERLAND = 1;
int const LF_SEA = 2;
int const LF_CLIFF = 6;
int const LF_DRIFT_TALUS = 9;
int const LF_DRIFT_BEACH = 10;
int const LF_DRIFT_DUNES = 11;                        // TODO 059 Implement dune landform class
int const LF_INTERVENTION_STRUCT = 12;
int const LF_INTERVENTION_NON_STRUCT = 13;
int const LF_ISLAND = 14;                             // Not yet implemented
int const LF_SEDIMENT_INPUT_UNCONSOLIDATED = 18;
int const LF_SEDIMENT_INPUT_CONSOLIDATED = 19;

// GIS raster input codes
int const FINE_CONS_RASTER = 1;
int const SAND_CONS_RASTER = 2;
int const COARSE_CONS_RASTER = 3;
int const FINE_UNCONS_RASTER = 4;
int const SAND_UNCONS_RASTER = 5;
int const COARSE_UNCONS_RASTER = 6;
int const SUSP_SED_RASTER = 7;
int const LANDFORM_RASTER = 8;
int const INTERVENTION_CLASS_RASTER = 9;
int const INTERVENTION_HEIGHT_RASTER = 10;

// GIS vector data type codes
int const VEC_FIELD_DATA_ANY = 0;
int const VEC_FIELD_DATA_INT = 1;
int const VEC_FIELD_DATA_REAL = 2;
int const VEC_FIELD_DATA_STRING = 3;
int const VEC_FIELD_DATA_OTHER = 4;

// GIS vector geometry codes
int const VEC_GEOMETRY_POINT = 1;
int const VEC_GEOMETRY_LINE = 2;
int const VEC_GEOMETRY_POLYGON = 3;
int const VEC_GEOMETRY_OTHER = 4;

// GIS vector input codes and constraints
int const DEEP_WATER_WAVE_STATIONS_VEC = 1;
int const DEEP_WATER_WAVE_STATIONS_MAX_LAYER = 1;
int const DEEP_WATER_WAVE_STATIONS_POINT_GEOMETRY = VEC_GEOMETRY_POINT;
int const SEDIMENT_INPUT_EVENT_LOCATION_VEC = 2;
int const SEDIMENT_INPUT_EVENT_LOCATION_MAX_LAYER = 1;
int const SEDIMENT_INPUT_EVENT_LOCATION_POINT_GEOMETRY = VEC_GEOMETRY_POINT;
int const FLOOD_LOCATION_POINT_GEOMETRY = VEC_GEOMETRY_POINT;
int const SEDIMENT_INPUT_EVENT_LOCATION_LINE_GEOMETRY = VEC_GEOMETRY_LINE;
int const FLOOD_LOCATION_VEC = 3;
int const FLOOD_LOCATION_MAX_LAYER = 1;

// GIS raster output codes
int const RASTER_PLOT_ACTIVE_ZONE = 1;
int const RASTER_PLOT_ACTUAL_BEACH_EROSION = 2;
int const RASTER_PLOT_ACTUAL_PLATFORM_EROSION = 3;
int const RASTER_PLOT_AVG_SEA_DEPTH = 4;
int const RASTER_PLOT_AVG_SUSPENDED_SEDIMENT = 5;
int const RASTER_PLOT_AVG_WAVE_HEIGHT = 6;
int const RASTER_PLOT_AVG_WAVE_ANGLE = 7;
int const RASTER_PLOT_BASEMENT_ELEVATION = 8;
int const RASTER_PLOT_BEACH_DEPOSITION = 9;
int const RASTER_PLOT_BEACH_MASK = 10;
int const RASTER_PLOT_BEACH_PROTECTION = 11;
int const RASTER_PLOT_CLIFF_COLLAPSE_DEPOSITION_COARSE = 12;
int const RASTER_PLOT_CLIFF_COLLAPSE_DEPOSITION_SAND = 13;
int const RASTER_PLOT_CLIFF_COLLAPSE_EROSION_COARSE = 14;
int const RASTER_PLOT_CLIFF_COLLAPSE_EROSION_FINE = 15;
int const RASTER_PLOT_CLIFF_COLLAPSE_EROSION_SAND = 16;
int const RASTER_PLOT_CLIFF_COLLAPSE_TIMESTEP = 17;
int const RASTER_PLOT_CLIFF_NOTCH_ALL = 18;
int const RASTER_PLOT_CLIFF_TOE = 19;
int const RASTER_PLOT_COARSE_CONSOLIDATED_SEDIMENT = 20;
int const RASTER_PLOT_COARSE_UNCONSOLIDATED_SEDIMENT = 21;
int const RASTER_PLOT_COAST = 22;
int const RASTER_PLOT_DEEP_WATER_WAVE_HEIGHT = 23;
int const RASTER_PLOT_DEEP_WATER_WAVE_ANGLE = 24;
int const RASTER_PLOT_DEEP_WATER_WAVE_PERIOD = 25;
int const RASTER_PLOT_FINE_CONSOLIDATED_SEDIMENT = 26;
int const RASTER_PLOT_FINE_UNCONSOLIDATED_SEDIMENT = 27;
int const RASTER_PLOT_INTERVENTION_CLASS = 28;
int const RASTER_PLOT_INTERVENTION_HEIGHT = 29;
int const RASTER_PLOT_INUNDATION_MASK = 30;
int const RASTER_PLOT_LANDFORM = 31;
int const RASTER_PLOT_NORMAL_PROFILE = 32;
int const RASTER_PLOT_TOP_ELEV_INC_SEA = 33;
int const RASTER_PLOT_POLYGON = 34;
int const RASTER_PLOT_POLYGON_GAIN_OR_LOSS = 35;
int const RASTER_PLOT_POLYGON_UPDRIFT_OR_DOWNDRIFT = 36;
int const RASTER_PLOT_POTENTIAL_BEACH_EROSION = 37;
int const RASTER_PLOT_POTENTIAL_PLATFORM_EROSION = 38;
int const RASTER_PLOT_POTENTIAL_PLATFORM_EROSION_MASK = 39;
int const RASTER_PLOT_SAND_CONSOLIDATED_SEDIMENT = 40;
int const RASTER_PLOT_SAND_UNCONSOLIDATED_SEDIMENT = 41;
int const RASTER_PLOT_SEA_DEPTH = 42;
int const RASTER_PLOT_SEDIMENT_INPUT = 43;
int const RASTER_PLOT_SED_TOP_INC_TALUS_ELEV = 44;
// int const RASTER_PLOT_SETUP_SURGE_FLOOD_MASK = 45;
// int const RASTER_PLOT_SETUP_SURGE_RUNUP_FLOOD_MASK = 46;
int const RASTER_PLOT_SHADOW_DOWNDRIFT_ZONE = 47;
int const RASTER_PLOT_SHADOW_ZONE = 48;
int const RASTER_PLOT_SLICE = 49;
int const RASTER_PLOT_SLOPE_FOR_CLIFF_TOE = 50;
int const RASTER_PLOT_CONS_SED_SLOPE = 51;
int const RASTER_PLOT_SUSPENDED_SEDIMENT = 52;
int const RASTER_PLOT_TALUS = 53;
int const RASTER_PLOT_TOTAL_ACTUAL_BEACH_EROSION = 54;
int const RASTER_PLOT_TOTAL_ACTUAL_PLATFORM_EROSION = 55;
int const RASTER_PLOT_TOTAL_BEACH_DEPOSITION = 56;
int const RASTER_PLOT_TOTAL_CLIFF_COLLAPSE_DEPOSITION_COARSE = 57;
int const RASTER_PLOT_TOTAL_CLIFF_COLLAPSE_DEPOSITION_SAND = 58;
int const RASTER_PLOT_TOTAL_CLIFF_COLLAPSE_EROSION_COARSE = 59;
int const RASTER_PLOT_TOTAL_CLIFF_COLLAPSE_EROSION_FINE = 60;
int const RASTER_PLOT_TOTAL_CLIFF_COLLAPSE_EROSION_SAND = 61;
int const RASTER_PLOT_TOTAL_POTENTIAL_BEACH_EROSION = 62;
int const RASTER_PLOT_TOTAL_POTENTIAL_PLATFORM_EROSION = 63;
int const RASTER_PLOT_WAVE_FLOOD_LINE = 64;
int const RASTER_PLOT_WAVE_HEIGHT = 65;
int const RASTER_PLOT_WAVE_ANGLE = 66;

// GIS vector output codes
// int const VECTOR_PLOT_FLOOD_SWL_SETUP_LINE = 19;
int const VECTOR_PLOT_AVG_WAVE_ANGLE_AND_HEIGHT = 1;
int const VECTOR_PLOT_BREAKING_WAVE_HEIGHT = 2;
// int const VECTOR_PLOT_CLIFF_EDGE = 3;
int const VECTOR_PLOT_CLIFF_NOTCH_ACTIVE = 4;
int const VECTOR_PLOT_COAST = 5;
int const VECTOR_PLOT_COAST_CURVATURE = 6;
int const VECTOR_PLOT_COAST_SWL_HIGHEST = 7;
int const VECTOR_PLOT_COAST_SWL_LOWEST = 8;
int const VECTOR_PLOT_DEEP_WATER_WAVE_ANGLE_AND_HEIGHT = 9;
int const VECTOR_PLOT_DOWNDRIFT_ZONE_BOUNDARY = 10;
int const VECTOR_PLOT_FLOOD_LINE = 11;
int const VECTOR_PLOT_INVALID_NORMALS = 12;
int const VECTOR_PLOT_MEAN_WAVE_ENERGY = 13;
int const VECTOR_PLOT_NORMALS = 14;
int const VECTOR_PLOT_POLYGON_BOUNDARY = 15;
int const VECTOR_PLOT_POLYGON_NODES = 16;
int const VECTOR_PLOT_RUN_UP = 17;
int const VECTOR_PLOT_SHADOW_ZONE_BOUNDARY = 18;
int const VECTOR_PLOT_STORM_SURGE = 19;
int const VECTOR_PLOT_WAVE_ANGLE_AND_HEIGHT = 20;
int const VECTOR_PLOT_WAVE_ENERGY_SINCE_COLLAPSE = 21;
int const VECTOR_PLOT_WAVE_SETUP = 22;
int const VECTOR_PLOT_WAVE_TRANSECT_POINTS = 23;
// int const VECTOR_PLOT_FLOOD_SWL_SETUP_SURGE_LINE = 24;
// int const VECTOR_PLOT_FLOOD_SWL_SETUP_SURGE_RUNUP_LINE = 25;

// Return codes
int const RTN_OK = 0;
int const RTN_HELP_ONLY = 1;
int const RTN_CHECK_ONLY = 2;
int const RTN_USER_ABORT = 3;
int const RTN_ERR_BADPARAM = 4;
int const RTN_ERR_INI = 5;
int const RTN_ERR_CMEDIR = 6;
int const RTN_ERR_RUNDATA = 7;
int const RTN_ERR_SCAPE_SHAPE_FUNCTION_FILE = 8;
int const RTN_ERR_TIDEDATAFILE = 9;
int const RTN_ERR_LOGFILE = 10;
int const RTN_ERR_OUTFILE = 11;
int const RTN_ERR_TSFILE = 12;
int const RTN_ERR_DEMFILE = 13;
int const RTN_ERR_RASTER_FILE_READ = 14;
int const RTN_ERR_VECTOR_FILE_READ = 15;
int const RTN_ERR_MEMALLOC = 16;
int const RTN_ERR_RASTER_GIS_OUT_FORMAT = 17;
int const RTN_ERR_VECTOR_GIS_OUT_FORMAT = 18;
int const RTN_ERR_TEXT_FILE_WRITE = 19;
int const RTN_ERR_RASTER_FILE_WRITE = 20;
int const RTN_ERR_VECTOR_FILE_WRITE = 21;
int const RTN_ERR_TIMESERIES_FILE_WRITE = 22;
int const RTN_ERR_LINETOGRID = 23;
int const RTN_ERR_PROFILE_SPACING = 24;
int const RTN_ERR_CHECK_INTERSECTION = 25;
int const RTN_ERR_PROFILE_ENDPOINT_IS_INLAND = 26;
int const RTN_ERR_NO_SOLUTION_FOR_ENDPOINT = 27;
int const RTN_ERR_PROFILE_END_INSUFFICIENT_DEPTH = 28;
int const RTN_ERR_NO_PROFILES_1 = 29;
int const RTN_ERR_NO_PROFILES_2 = 30;
int const RTN_ERR_NOSEACELLS = 31;
int const RTN_ERR_GRID_TO_LINE = 32;
int const RTN_ERR_NO_COAST = 34;
int const RTN_ERR_PROFILE_WRITE = 35;
int const RTN_ERR_TIME_UNITS = 36;
int const RTN_ERR_CLIFF_NOTCH = 37;
int const RTN_ERR_CLIFF_CANNOT_DEPOSIT_ALL = 38;
int const RTN_ERR_BAD_INDEX = 39;
int const RTN_ERR_EDGE_OF_GRID = 40;
int const RTN_ERR_NO_SEAWARD_END_OF_PROFILE_BEACH_EROSION = 42;
int const RTN_ERR_NO_SEAWARD_END_OF_PROFILE_UPCOAST_BEACH_DEPOSITION = 43;
int const RTN_ERR_NO_SEAWARD_END_OF_PROFILE_DOWNCOAST_BEACH_DEPOSITION = 44;
int const RTN_ERR_LANDFORM_TO_GRID = 45;
int const RTN_ERR_NO_TOP_LAYER_DURING_WAVE_CALC = 46;
int const RTN_ERR_NO_TOP_LAYER_DURING_BEACH_CALC = 47;
int const RTN_ERR_NO_TOP_LAYER_DURING_PLATFORM_CALC = 48;
int const RTN_ERR_NO_TOP_LAYER_DURING_CLIFF_COLLAPSE_CALC = 49;
int const RTN_ERR_NO_ADJACENT_POLYGON = 50;
int const RTN_ERR_BAD_MULTILINE = 51;
int const RTN_ERR_CANNOT_INSERT_POINT = 52;
int const RTN_ERR_CANNOT_ASSIGN_COASTAL_LANDFORM = 53;
int const RTN_ERR_SHADOW_ZONE_FLOOD_FILL_NOGRID = 54;
int const RTN_ERR_SHADOW_ZONE_FLOOD_START_POINT = 55;
int const RTN_ERR_CSHORE_EMPTY_PROFILE = 56;
int const RTN_ERR_CSHORE_FILE_INPUT = 57;
int const RTN_ERR_READING_CSHORE_FILE_OUTPUT = 58;
int const RTN_ERR_WAVE_INTERPOLATION_LOOKUP = 59;
int const RTN_ERR_GRIDCREATE = 60;
int const RTN_ERR_COAST_CANT_FIND_EDGE_CELL = 61;
int const RTN_ERR_CSHORE_ERROR = 62;
int const RTN_ERR_NO_CELL_UNDER_COASTLINE = 63;
int const RTN_ERR_OPEN_DEEP_WATER_WAVE_DATA = 64;
int const RTN_ERR_READING_DEEP_WATER_WAVE_DATA = 65;
int const RTN_ERR_BOUNDING_BOX = 66;
int const RTN_ERR_READING_SEDIMENT_INPUT_EVENT = 67;
int const RTN_ERR_SEDIMENT_INPUT_EVENT = 68;
int const RTN_ERR_SEDIMENT_INPUT_EVENT_LOCATION = 69;
int const RTN_ERR_WAVESTATION_LOCATION = 70;
int const RTN_ERR_FLOOD_LOCATION = 71;
int const RTN_ERR_CLIFF_NOT_IN_POLYGON = 72;
int const RTN_ERR_CELL_MARKED_PROFILE_COAST_BUT_NOT_PROFILE = 73;
int const RTN_ERR_TRACING_FLOOD = 74;
int const RTN_ERR_NO_START_FINISH_POINTS_TRACING_COAST = 75;
int const RTN_ERR_NO_VALID_COAST = 76;
int const RTN_ERR_COAST_TRACING_REPEATING = 77;
int const RTN_ERR_COAST_TRACING_SAME_START_FINISH = 78;
int const RTN_ERR_COAST_TRACING_ZERO_LENGTH = 79;
int const RTN_ERR_COAST_TRACING_TOO_SHORT = 80;
int const RTN_ERR_COAST_TRACING_TOO_LONG = 81;
int const RTN_ERR_COAST_TRACING_OFF_EDGE = 82;
int const RTN_ERR_CELL_NOT_FOUND_IN_HIT_PROFILE_DIFFERENT_COASTS = 83;
int const RTN_ERR_POINT_NOT_FOUND_IN_MULTILINE_DIFFERENT_COASTS = 84;
int const RTN_ERR_CELL_NOT_FOUND_IN_HIT_PROFILE = 85;
int const RTN_ERR_CELL_IN_POLY_BUT_NO_POLY_COAST = 86;
int const RTN_ERR_CLIFF_TALUS_TO_UNCONS = 87;
int const RTN_ERR_COAST_TRACING_SAME_EDGE_START_FINISH = 88;
int const RTN_ERR_POLYGON_FILL_START_POINT = 89;
int const RTN_ERR_UNKNOWN = 999;

// Coast-normal profile status codes
int const PROFILE_STATUS_OK = 0;                         // No problem, profile is valid
int const PROFILE_STATUS_TOO_SHORT = 1;                  // Is this profile invalid because it is too short?
int const PROFILE_STATUS_HIT_LAND = 2;                   // Has this profile invalid because it hit land?
int const PROFILE_STATUS_HIT_COAST = 3;                  // Has this profile invalid because it hit a coast?
int const PROFILE_STATUS_HIT_PROFILE = 4;                // Has this profile invalid because it hit another profile?
int const PROFILE_STATUS_HIT_INTERVENTION = 5;           // Has this profile invalid because it hit an intervention?
int const PROFILE_STATUS_TOO_MANY_SEGMENTS = 6;          // Too many line segments for this profile?

// Elevation and 'slice' codes
int const ELEV_IN_BASEMENT = -1;
int const ELEV_ABOVE_SEDIMENT_TOP = -2;
int const NO_NONZERO_THICKNESS_LAYERS = -3;

// Vector smoothing codes
int const SMOOTH_NONE = 0;
int const SMOOTH_RUNNING_MEAN = 1;
int const SMOOTH_SAVITZKY_GOLAY = 2;

// Grid-edge boundary treatment for unconsolidated sediment movement
int const GRID_EDGE_CLOSED = 0;
int const GRID_EDGE_OPEN = 1;
int const GRID_EDGE_RECIRCULATE = 2;

// Model for wave propagation
int const WAVE_MODEL_COVE = 0;
int const WAVE_MODEL_CSHORE = 1;

// Equation for estimating erosion of unconsolidated sediment
int const UNCONS_SEDIMENT_EQUATION_CERC = 0;
int const UNCONS_SEDIMENT_EQUATION_KAMPHUIS = 1;

int const PROFILE_CHECK_DIST_FROM_COAST = 5;                // Used in checking shoreline-normal profiles for intersection
int const GAP_BETWEEN_DIFFERENT_COAST_PROFILES = 30;        // In cells, is the gap between profile ends belonging to different coasts

int const NUM_DAYS_FOR_MEAN_HIGH_WATER_CALC = 30;           // Number of days to average daily high water elevation

int const NO_RUNUP_EQUATION = 0;                            // Runup not considered
int const RUNUP_EQUATION_NIELSEN_HANSLOW = 1;               // Runup equation is Nielsen, P. & Hanslow, D. J. 1991. Wave Runup Distributions on Natural Beaches. Journal of Coastal Research, 7, 1139-1152.  *** & DHI (2004) ???
int const RUNUP_EQUATION_MASE = 2;                          // Runup equation is Mase, H. 1989. Random Wave Runup Height on Gentle Slope. Journal of Waterway, Port, Coastal, and Ocean Engineering, 115, 649-661.
int const RUNUP_EQUATION_STOCKDON = 3;                      // Runup equation is Stockdon, H. F., Holman, R. A., Howd, P. A. & Sallenger JR, A. H. 2006. Empirical parameterization of setup, swash, and runup. Coastal Engineering, 53, 573-588.

//! Maximum number of avalanche iterations per timestep. Safety limit to prevent infinite loops in case of numerical issues.
int const MAX_SLUMP_ITERATIONS = 100;

//! Maximum number of along-coast passes to identify and correct coast-normal profile intersectionms
int const MAX_ALONG_COAST_PASSES = 20;

//! Maximum number of line segments into which a profile can be divided
int const MAX_SEGMENTS = 10;

unsigned long const MASK = 0xfffffffful;
unsigned long const SEDIMENT_INPUT_EVENT_ERROR = -1;
unsigned long const UNSIGNED_LONG_NODATA = 9999;

double const PI = 3.141592653589793238462643;
double const SQRT2 = 1.414213562;

double const D50_FINE_DEFAULT = 0.0625;                     // In mm
double const D50_SAND_DEFAULT = 0.42;                       // In mm
double const D50_COARSE_DEFAULT = 19.0;                     // In mm

double const BEACH_PROTECTION_HB_RATIO = 0.23;              // The beach protection factor is this times breaking depth
double const WALKDEN_HALL_PARAM_1 = 3.25;                   // First parameter in Equation 4 from Walkden & Hall, 2005
double const WALKDEN_HALL_PARAM_2 = 1.50;                   // Second parameter in Equation 4 from Walkden & Hall, 2005
double const TAN_BERM_SLOPE = 0.176;                        // Seaward slope of top of berm (that part of beach coastward of the start of the Dean profile)

double const DEPTH_OVER_DB_INCREMENT = 0.001;               // Depth over DB increment for erosion potential look-up function
double const INVERSE_DEPTH_OVER_DB_INCREMENT = 1000;        // Inverse of the above
double const DEAN_POWER = 2.0 / 3.0;                        // Dean profile exponent

// TODO 011 Let the user define these CShore input parameters
double const CSHORE_FRICTION_FACTOR = 0.015;                // Friction factor for CShore model
double const CSHORE_SURGE_LEVEL = 0.0;                      // TODO 007

double const TOLERANCE = 1e-7;                              // For bFPIsEqual, if too small (e.g. 1e-10), get spurious "rounding" errors
double const SED_ELEV_TOLERANCE = 1e-5;                     // Used to compare depth-equivalent sediment amounts in bFPIsEqual and elsewhere
double const MASS_BALANCE_TOLERANCE = 1e-5;                 // For bFPIsEqual, used to compare for mass balance checks
double const STRAIGHT_COAST_MAX_DETAILED_CURVATURE = -5;
double const STRAIGHT_COAST_MAX_SMOOTH_CURVATURE = -1;
double const MIN_LENGTH_OF_SHADOW_ZONE_LINE = 10;           // Used in shadow line tracing
double const MAX_LAND_LENGTH_OF_SHADOW_ZONE_LINE = 5;       // Used in shadow line tracing
double const CLIFF_COLLAPSE_HEIGHT_INCREMENT = 0.1;         // Increment the fractional height of the cliff talus Dean profile, if we have not been able to deposit enough
double const INTERVENTION_PROFILE_SPACING_FACTOR = 0.2;     // Profile spacing on interventions works better if it is smaller than profile spacing on coastline

double const CLIFF_NOTCH_CUTOFF_DISTANCE = 2;               // Cut-off SWL distance (m), measured downwards from the cliff notch apex: below this there is no notch incision
double const DBL_NODATA = -9999;

//! Angle of repose for sediment (degrees). This is a typical value for sand. Finer sediments may have lower angles (~28°), coarser sediments higher (~37°).
double const ANGLE_OF_REPOSE_DEG = 33.0;

//! Angle of repose in radians (pre-calculated for efficiency)
double const ANGLE_OF_REPOSE_RAD = 0.5759586531;  // 33° * π/180

//! Tangent of angle of repose (pre-calculated for direct slope comparison)
double const TAN_ANGLE_OF_REPOSE = 0.6494075931;  // tan(33°)

//! Minimum sediment volume (m³) to trigger avalanche redistribution. Prevents processing of trivially small amounts that don't affect morphology.
double const MIN_SLUMP_VOLUME = 0.001;  // 1 mm average depth over 1 m² cell

//! Fraction of excess sediment to redistribute per iteration. 0.5 = move 50% of unstable sediment each iteration. Lower values (0.2-0.3) are more stable but slower to converge. Higher values (0.6-0.8) converge faster but may overshoot.
double const SLUMP_REDISTRIBUTION_FRACTION = 0.5;

//! If there is lees than this depth (m) of cliff-collapse talus on a cell, then remove all talus from the cell
double const MIN_TALUS_DEPTH = 0.01;

string const PROGRAM_NAME = "Coastal Modelling Environment (CoastalME) version 1.4.1 (07 Jun 2026)";
string const PROGRAM_NAME_SHORT = "CME";
string const CME_INI = "cme.ini";
string const CME_YAML = "cme.yaml";

string const COPYRIGHT = "(C) 2026 Andres Payo, David Favis-Mortlock, and Wilf Chun";
string const LINE = "-------------------------------------------------------------------------------";
string const DISCLAIMER1 = "This program is distributed in the hope that it will be useful, but WITHOUT ANY";
string const DISCLAIMER2 = "WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A";
string const DISCLAIMER3 = "PARTICULAR PURPOSE. See the GNU General Public License for more details. You";
string const DISCLAIMER4 = "should have received a copy of the GNU General Public License along with this";
string const DISCLAIMER5 = "program; if not, contact the Free Software Foundation, Inc., 675 Mass Ave,";
string const DISCLAIMER6 = "Cambridge, MA 02139, USA.";

string const ABOUT = "simulates the long-term behaviour of a coast. This initial version considers only simple soft cliff cross-shore effects";
string const THANKS = "Many thanks to:\n\tTom Ashby\n\tManuel Cobos Budia\n\tMark Dickson\n\tJim W. Hall\n\tMartin D. Hurst\n\tMatthew Ives\n\tRobert J. Nicholls\n\tIan Townend\n\tMike J.A. Walkden";
string const GDAL_DRIVERS = "GDAL drivers";

string const USAGE = "Usage: cme [OPTION]...";
string const USAGE1 = "  --gdal             List GDAL drivers";
string const USAGE2 = "  --about            Information about this program";
string const USAGE3 = "  --help             Display this text";
string const USAGE4 = "  --home=DIRECTORY   Specify the location of the .ini file etc.";
string const USAGE5 = "  --datafile=FILE    Specify the location and name of the main datafile";
string const USAGE6 = "  --yaml             Use YAML format for the main datafile";

string const START_NOTICE = "- Started on ";
string const INITIALIZING_NOTICE = "- Initializing";
string const READING_FILE_LOCATIONS = "  - Reading file locations: ";
string const READING_RUN_DATA = "  - Reading run data file: ";
string const READING_BASEMENT = "  - Reading basement DEM: ";
string const READING_RASTER_FILES = "  - Reading raster GIS files";
string const READING_LANDFORM_FILE = "    - Landform class: ";
string const READING_INTERVENTION_CLASS_FILE = "    - Intervention class: ";
string const READING_INTERVENTION_HEIGHT_FILE = "    - Intervention height: ";
string const READING_SUSPENDED_SEDIMENT_FILE = "    - Suspended sediment: ";
string const READING_UNCONS_FINE_SEDIMENT_FILE = "    - Unconsolidated fine sediment (layer ";
string const READING_UNCONS_SAND_SEDIMENT_FILE = "    - Unconsolidated sand sediment (layer ";
string const READING_UNCONS_COARSE_SEDIMENT_FILE = "    - Unconsolidated coarse sediment (layer ";
string const READING_CONS_FINE_SEDIMENT_FILE = "    - Consolidated fine sediment (layer ";
string const READING_CONS_SAND_SEDIMENT_FILE = "    - Consolidated sand sediment (layer ";
string const READING_CONS_COARSE_SEDIMENT_FILE = "    - Consolidated coarse sediment (layer ";
string const READING_VECTOR_FILES = "  - Reading vector GIS files";
string const READING_DEEP_WATER_WAVE_FILE = "    - Deep water wave values: ";
string const READING_SED_INPUT_EVENT_FILE = "    - Sediment input event values: ";
string const READING_FLOOD_LOCATION = "    - Characteristic locations for flood: ";
string const READING_SCAPE_SHAPE_FUNCTION_FILE = "  - Reading SCAPE shape function file";
string const READING_TIDE_DATA_FILE = "  - Reading tide data file: ";
string const ALLOCATE_MEMORY = "  - Allocating memory for raster grid";
string const ADD_LAYERS = "  - Adding sediment layers to raster grid";
string const INITIALIZING_FINAL = "  - Finishing initialization";
string const RUN_NOTICE = "- Running simulation";
string const SIMULATING = "\r  - Simulating ";
string const FINAL_OUTPUT = "- Writing final output";
string const SEND_EMAIL = "  - Sending email to ";
string const RUN_END_NOTICE = "Run ended at ";
string const PRESS_KEY = "Press any key to continue...";

string const ERROR_NOTICE = " with error code ";
string const EMAIL_ERROR = "Could not send email";

string const SCAPE_DIR = "scape/";
string const SCAPE_SHAPE_FUNCTION_FILE = "ShapeFunction.dat";
string const EROSION_POTENTIAL_LOOKUP_FILE = "ErosionPotential.csv";

string const CSHORE_DIR = "cshore/";
string const CSHORE_INFILE = "infile";

string const ERR = "*** ERROR ";
string const WARN = "WARNING ";
string const NOTE = "      Note ";

string const MASS_BALANCE_ERROR = "MASS BALANCE ERROR";

string const PER_ITER_HEAD1 = "<-----ELAPSED----><--SEA-><----POTENTIAL---><-----------ACTUAL-----------><-----POTENTIAL-----><------------ACTUAL-------------><-----------ACTUAL------------><--SEDIMENT--><---CLIFF COLLAPSE---><SUSP>";

string const PER_ITER_HEAD2 = "       TIME         DEPTH  PLATFORM EROSION        PLATFORM EROSION          BEACH EROSION              BEACH EROSION                  BEACH DEPOSITION         INPUT EVENT     EROSION  DEPOSITION  SED";

string const PER_ITER_HEAD3 = "Time  Hours  Years    Avg  % Sea   All  Erod % Sea   All  Erod <--sea avg->  % Sea   All   Erod  % Sea   All   Erod <--sea avg->  % Sea   All  Depos  <--sea-->                <-coast avg-><--sea-->";

string const PER_ITER_HEAD4 = "Step                        Area   Sea  Area  Area   Sea  Area   F   S   C    Area   Sea   Area   Area   Sea   Area   F   S   C    Area   Sea   Area    S   C     F   S   C     F   S   C    S   C    F";

string const PER_ITER_HEAD5 = "                                   Avg   Avg         Avg   Avg                              Avg          Avg    Avg                       Avg    Avg";

string const PER_ITER_HEAD = "PER-ITERATION RESULTS =============================================================================================================================================================================================";

string const PER_ITER_CSV_HEAD = "Timestep,Hours, Years, AvgSeaDepth_m, PotPlatformErosion_PctSeaArea, PotPlatformErosion_AllAvg_mm, PotPlatformErosion_ErodAvg_mm, ActPlatformErosion_PctSeaArea, ActPlatformErosion_AllAvg_mm, ActPlatformErosion_ErodAvg_mm, ActPlatformErosion_Fine_mm, ActPlatformErosion_Sand_mm, ActPlatformErosion_Coarse_mm, PotBeachErosion_PctSeaArea, PotBeachErosion_AllAvg_mm, PotBeachErosion_ErodAvg_mm, ActBeachErosion_PctSeaArea, ActBeachErosion_AllAvg_mm, ActBeachErosion_ErodAvg_mm, ActBeachErosion_Fine_mm, ActBeachErosion_Sand_mm, ActBeachErosion_Coarse_mm, BeachDeposition_PctSeaArea, BeachDeposition_AllAvg_mm, BeachDeposition_DepAvg_mm, BeachDeposition_Sand_mm, BeachDeposition_Coarse_mm, SedimentInput_Fine, SedimentInput_Sand, SedimentInput_Coarse, CliffCollapse_Fine_mm, CliffCollapse_Sand_mm, CliffCollapse_Coarse_mm, CliffDeposition_Sand_mm, CliffDeposition_Coarse_mm, SuspendedSediment_mm, GISEvents";

string const ENDHYDROLOGYHEAD = "END OF SIMULATION: HYDROLOGY ======================================================================================================================================================================================";
string const ENDSEDIMENTHEAD = "END OF SIMULATION: SEDIMENT MOVEMENT ==============================================================================================================================================================================";
string const PERFORMHEAD = "END OF SIMULATION: PERFORMANCE ====================================================================================================================================================================================";

string const OUTEXT = ".out";
string const LOGEXT = ".log";
string const CSVEXT = ".csv";

string const DEEP_WATER_WAVE_STATION_ID = "id";
string const SEDIMENT_INPUT_EVENT_LOCATION_ID = "id";
string const FLOOD_LOCATION_ID = "id";

// GIS raster output user codes
string const RASTER_ACTIVE_ZONE_CODE = "active_zone";
string const RASTER_ACTIVE_ZONE_NAME = "active_zone";
string const RASTER_ACTUAL_BEACH_EROSION_CODE = "actual_beach_erosion";
string const RASTER_ACTUAL_BEACH_EROSION_NAME = "actual_beach_erosion";
string const RASTER_ACTUAL_PLATFORM_EROSION_CODE = "actual_platform_erosion";
string const RASTER_ACTUAL_PLATFORM_EROSION_NAME = "actual_platform_erosion";
string const RASTER_ALL_OUTPUT_CODE = "all";
string const RASTER_AVG_SEA_DEPTH_CODE = "avg_sea_depth";
string const RASTER_AVG_SEA_DEPTH_NAME = "avg_sea_depth";
string const RASTER_AVG_SUSP_SED_CODE = "avg_susp_sed";
string const RASTER_AVG_SUSP_SED_NAME = "avg_susp_sed";
string const RASTER_AVG_WAVE_HEIGHT_CODE = "avg_wave_height";
string const RASTER_AVG_WAVE_HEIGHT_NAME = "avg_wave_height";
string const RASTER_AVG_WAVE_ANGLE_CODE = "avg_wave_angle";
string const RASTER_AVG_WAVE_ANGLE_NAME = "avg_wave_angle";
string const RASTER_BASEMENT_ELEVATION_CODE = "basement_elevation";
string const RASTER_BASEMENT_ELEVATION_NAME = "basement_elevation";
string const RASTER_BEACH_DEPOSITION_CODE = "beach_deposition";
string const RASTER_BEACH_DEPOSITION_NAME = "beach_deposition";
string const RASTER_BEACH_MASK_CODE = "beach_mask";
string const RASTER_BEACH_MASK_NAME = "beach_mask";
string const RASTER_BEACH_PROTECTION_CODE = "beach_protection";
string const RASTER_BEACH_PROTECTION_NAME = "beach_protection";
string const RASTER_CLIFF_COLLAPSE_DEPOSITION_COARSE_CODE = "cliff_collapse_talus_deposition_coarse";
string const RASTER_CLIFF_COLLAPSE_DEPOSITION_COARSE_NAME = "cliff_collapse_talus_deposition_coarse";
string const RASTER_CLIFF_COLLAPSE_DEPOSITION_SAND_CODE = "cliff_collapse_talus_deposition_sand";
string const RASTER_CLIFF_COLLAPSE_DEPOSITION_SAND_NAME = "cliff_collapse_talus_deposition_sand";
string const RASTER_CLIFF_COLLAPSE_EROSION_COARSE_CODE = "cliff_collapse_erosion_coarse";
string const RASTER_CLIFF_COLLAPSE_EROSION_COARSE_NAME = "cliff_collapse_erosion_coarse";
string const RASTER_CLIFF_COLLAPSE_EROSION_FINE_CODE = "cliff_collapse_erosion_fine";
string const RASTER_CLIFF_COLLAPSE_EROSION_FINE_NAME = "cliff_collapse_erosion_fine";
string const RASTER_CLIFF_COLLAPSE_EROSION_SAND_CODE = "cliff_collapse_erosion_sand";
string const RASTER_CLIFF_COLLAPSE_EROSION_SAND_NAME = "cliff_collapse_erosion_sand";
string const RASTER_CLIFF_COLLAPSE_TIMESTEP_CODE = "cliff_collapse_timestep";
string const RASTER_CLIFF_COLLAPSE_TIMESTEP_NAME = "cliff_collapse_timestep_all";
string const RASTER_CLIFF_NOTCH_ALL_CODE = "cliff_notch_all";
string const RASTER_CLIFF_NOTCH_ALL_NAME = "cliff_notch_all";
string const RASTER_COARSE_CONS_CODE = "cons_sed_coarse";
string const RASTER_COARSE_CONS_NAME = "cons_sed_coarse";
string const RASTER_COARSE_UNCONS_CODE = "uncons_sed_coarse";
string const RASTER_COARSE_UNCONS_NAME = "uncons_sed_coarse";
string const RASTER_COAST_CODE = "rcoast";
string const RASTER_COAST_NAME = "rcoast";
string const RASTER_COAST_NORMAL_CODE = "rcoast_normal";
string const RASTER_COAST_NORMAL_NAME = "rcoast_normal";
string const RASTER_DEEP_WATER_WAVE_HEIGHT_CODE = "deep_water_wave_height";
string const RASTER_DEEP_WATER_WAVE_HEIGHT_NAME = "deep_water_wave_height";
string const RASTER_DEEP_WATER_WAVE_ANGLE_CODE = "deep_water_wave_angle";
string const RASTER_DEEP_WATER_WAVE_ANGLE_NAME = "deep_water_wave_angle";
string const RASTER_DEEP_WATER_WAVE_PERIOD_CODE = "deep_water_wave_period";
string const RASTER_DEEP_WATER_WAVE_PERIOD_NAME = "deep_water_wave_period";
string const RASTER_FINE_CONS_CODE = "cons_sed_fine";
string const RASTER_FINE_CONS_NAME = "cons_sed_fine";
string const RASTER_FINE_UNCONS_CODE = "uncons_sed_fine";
string const RASTER_FINE_UNCONS_NAME = "uncons_sed_fine";
string const RASTER_INTERVENTION_CLASS_CODE = "intervention_class";
string const RASTER_INTERVENTION_CLASS_NAME = "intervention_class";
string const RASTER_INTERVENTION_HEIGHT_CODE = "intervention_height";
string const RASTER_INTERVENTION_HEIGHT_NAME = "intervention_height";
string const RASTER_INUNDATION_MASK_CODE = "inundation_mask";
string const RASTER_INUNDATION_MASK_NAME = "inundation_mask";
string const RASTER_LANDFORM_CODE = "landform_class";
string const RASTER_LANDFORM_NAME = "landform_class";
string const RASTER_POLYGON_CODE = "polygon_raster";
string const RASTER_POLYGON_GAIN_OR_LOSS_CODE = "polygon_gain_or_loss";
string const RASTER_POLYGON_GAIN_OR_LOSS_NAME = "polygon_gain_or_loss";
string const RASTER_POLYGON_NAME = "polygon_raster";
string const RASTER_POLYGON_UPDRIFT_OR_DOWNDRIFT_CODE = "polygon_updrift_or_downdrift";
string const RASTER_POLYGON_UPDRIFT_OR_DOWNDRIFT_NAME = "polygon_updrift_or_downdrift";
string const RASTER_POTENTIAL_BEACH_EROSION_CODE = "potential_beach_erosion";
string const RASTER_POTENTIAL_BEACH_EROSION_NAME = "potential_beach_erosion";
string const RASTER_POTENTIAL_PLATFORM_EROSION_CODE = "potential_platform_erosion";
string const RASTER_POTENTIAL_PLATFORM_EROSION_MASK_CODE = "potential_platform_erosion_mask";
string const RASTER_POTENTIAL_PLATFORM_EROSION_MASK_NAME = "potential_platform_erosion_mask";
string const RASTER_POTENTIAL_PLATFORM_EROSION_NAME = "potential_platform_erosion";
string const RASTER_SAND_CONS_CODE = "cons_sed_sand";
string const RASTER_SAND_CONS_NAME = "cons_sed_sand";
string const RASTER_SAND_UNCONS_CODE = "uncons_sed_sand";
string const RASTER_SAND_UNCONS_NAME = "uncons_sed_sand";
string const RASTER_SEA_DEPTH_CODE = "sea_depth";
string const RASTER_SEA_DEPTH_NAME = "sea_depth";
string const RASTER_SEDIMENT_INPUT_EVENT_CODE = "sediment_input_total";
string const RASTER_SEDIMENT_INPUT_EVENT_NAME = "sediment_input_total";
string const RASTER_SEDIMENT_TOP_CODE = "sediment_top_elevation";
string const RASTER_SEDIMENT_TOP_ELEVATION_NAME = "sediment_top_elevation";
// string const RASTER_SETUP_SURGE_FLOOD_MASK_CODE = "flood_setup_surge_mask";
// string const RASTER_SETUP_SURGE_FLOOD_MASK_NAME = "flood_setup_surge_mask";
// string const RASTER_SETUP_SURGE_RUNUP_FLOOD_MASK_CODE = "flood_setup_surge_runup_mask";
// string const RASTER_SETUP_SURGE_RUNUP_FLOOD_MASK_NAME = "flood_setup_surge_runup_mask";
string const RASTER_SHADOW_DOWNDRIFT_ZONE_CODE = "shadow_downdrift_zones";
string const RASTER_SHADOW_DOWNDRIFT_ZONE_NAME = "shadow_downdrift_zones";
string const RASTER_SHADOW_ZONE_CODE = "shadow_zones";
string const RASTER_SHADOW_ZONE_NAME = "shadow_zones";
string const RASTER_SLICE_CODE = "slice";
string const RASTER_SLICE_NAME = "slice";
string const RASTER_SLOPE_FOR_CLIFF_TOE_NAME = "toe_slope";       // Note no code for this, because is chosen by m_bCliffToeLocate in input file
string const RASTER_SLOPE_OF_CONSOLIDATED_SEDIMENT_CODE = "cons_sediment_slope";
string const RASTER_SLOPE_OF_CONSOLIDATED_SEDIMENT_NAME = "cons_sediment_slope";
string const RASTER_SUSP_SED_CODE = "susp_sed";
string const RASTER_SUSP_SED_NAME = "susp_sed";
string const RASTER_TALUS_CODE = "talus";
string const RASTER_TALUS_NAME = "talus";
string const RASTER_TOP_ELEVATION_INC_SEA_CODE = "top_elevation";
string const RASTER_TOP_ELEVATION_INC_SEA_NAME = "top_elevation_inc_sea";
string const RASTER_TOTAL_ACTUAL_BEACH_EROSION_CODE = "total_actual_beach_erosion";
string const RASTER_TOTAL_ACTUAL_BEACH_EROSION_NAME = "total_actual_beach_erosion";
string const RASTER_TOTAL_ACTUAL_PLATFORM_EROSION_CODE = "total_actual_platform_erosion";
string const RASTER_TOTAL_ACTUAL_PLATFORM_EROSION_NAME = "total_actual_platform_erosion";
string const RASTER_TOTAL_BEACH_DEPOSITION_CODE = "total_beach_deposition";
string const RASTER_TOTAL_BEACH_DEPOSITION_NAME = "total_beach_deposition";
string const RASTER_TOTAL_CLIFF_COLLAPSE_DEPOSITION_COARSE_CODE = "total_cliff_collapse_talus_deposition_coarse";
string const RASTER_TOTAL_CLIFF_COLLAPSE_DEPOSITION_COARSE_NAME = "total_cliff_collapse_talus_deposition_coarse";
string const RASTER_TOTAL_CLIFF_COLLAPSE_DEPOSITION_SAND_CODE = "total_cliff_collapse_talus_deposition_sand";
string const RASTER_TOTAL_CLIFF_COLLAPSE_DEPOSITION_SAND_NAME = "total_cliff_collapse_talus_deposition_sand";
string const RASTER_TOTAL_CLIFF_COLLAPSE_EROSION_COARSE_CODE = "total_cliff_collapse_erosion_coarse";
string const RASTER_TOTAL_CLIFF_COLLAPSE_EROSION_COARSE_NAME = "total_cliff_collapse_erosion_coarse";
string const RASTER_TOTAL_CLIFF_COLLAPSE_EROSION_FINE_CODE = "total_cliff_collapse_erosion_fine";
string const RASTER_TOTAL_CLIFF_COLLAPSE_EROSION_FINE_NAME = "total_cliff_collapse_erosion_fine";
string const RASTER_TOTAL_CLIFF_COLLAPSE_EROSION_SAND_CODE = "total_cliff_collapse_erosion_sand";
string const RASTER_TOTAL_CLIFF_COLLAPSE_EROSION_SAND_NAME = "total_cliff_collapse_erosion_sand";
string const RASTER_TOTAL_POTENTIAL_BEACH_EROSION_CODE = "total_potential_beach_erosion";
string const RASTER_TOTAL_POTENTIAL_BEACH_EROSION_NAME = "total_potential_beach_erosion";
string const RASTER_TOTAL_POTENTIAL_PLATFORM_EROSION_CODE = "total_potential_platform_erosion";
string const RASTER_TOTAL_POTENTIAL_PLATFORM_EROSION_NAME = "total_potential_platform_erosion";
string const RASTER_USUAL_OUTPUT_CODE = "usual";
string const RASTER_WAVE_FLOOD_LINE_CODE = "wave_flood_line";
string const RASTER_WAVE_FLOOD_LINE_NAME = "wave_flood_line";
string const RASTER_WAVE_HEIGHT_CODE = "wave_height";
string const RASTER_WAVE_HEIGHT_NAME = "wave_height";
string const RASTER_WAVE_ANGLE_CODE = "wave_angle";
string const RASTER_WAVE_ANGLE_NAME = "wave_angle";
string const RASTER_WAVE_PERIOD_CODE = "wave_period";
string const RASTER_WAVE_PERIOD_NAME = "wave_period";

// GIS raster output titles
string const RASTER_PLOT_ACTIVE_ZONE_TITLE = "Active zone";
string const RASTER_PLOT_ACTUAL_BEACH_EROSION_TITLE = "Actual (constrained) beach erosion depth";
string const RASTER_PLOT_ACTUAL_PLATFORM_EROSION_TITLE = "Actual (constrained) shore platform erosion depth";
string const RASTER_PLOT_AVG_SEA_DEPTH_TITLE = "Average sea depth";
string const RASTER_PLOT_AVG_SUSPENDED_SEDIMENT_TITLE = "Average depth of suspended sediment";
string const RASTER_PLOT_AVG_WAVE_HEIGHT_TITLE = "Average wave height";
string const RASTER_PLOT_AVG_WAVE_ANGLE_TITLE = "Average wave orientation";
string const RASTER_PLOT_BASEMENT_ELEVATION_TITLE = "Basement elevation";
string const RASTER_PLOT_BEACH_DEPOSITION_TITLE = "Beach deposition depth";
string const RASTER_PLOT_BEACH_MASK_TITLE = "Beach mask";
string const RASTER_PLOT_BEACH_PROTECTION_TITLE = "Beach protection factor";
string const RASTER_PLOT_CLIFF_COLLAPSE_DEPOSITION_COARSE_TITLE = "Depth of coarse talus from cliff collapse";
string const RASTER_PLOT_CLIFF_COLLAPSE_DEPOSITION_SAND_TITLE = "Depth of sand talus from cliff collapse";
string const RASTER_PLOT_CLIFF_COLLAPSE_EROSION_COARSE_TITLE = "Cliff collapse depth of erosion, coarse sediment";
string const RASTER_PLOT_CLIFF_COLLAPSE_EROSION_FINE_TITLE = "Cliff collapse depth of erosion, fine sediment";
string const RASTER_PLOT_CLIFF_COLLAPSE_EROSION_SAND_TITLE = "Cliff collapse depth of erosion, sand sediment";
string const RASTER_PLOT_CLIFF_COLLAPSE_TIMESTEP_TITLE = "Timestep at which cliff collapse occurred";
string const RASTER_PLOT_CLIFF_NOTCH_ALL_TITLE = "All cliff notch incision";
string const RASTER_PLOT_CLIFF_TOE_TITLE = "Cliff toe cells";
string const RASTER_PLOT_COARSE_CONSOLIDATED_SEDIMENT_TITLE = "Consolidated coarse sediment depth";
string const RASTER_PLOT_COARSE_UNCONSOLIDATED_SEDIMENT_TITLE = "Unconsolidated coarse sediment depth";
string const RASTER_PLOT_COAST_TITLE = "Rasterized coastline";
string const RASTER_PLOT_DEEP_WATER_WAVE_HEIGHT_TITLE = "Deep water wave height";
string const RASTER_PLOT_DEEP_WATER_WAVE_ANGLE_TITLE = "Deep water wave orientation";
string const RASTER_PLOT_DEEP_WATER_WAVE_PERIOD_TITLE = "Deep water wave period";
string const RASTER_PLOT_FINE_CONSOLIDATED_SEDIMENT_TITLE = "Consolidated fine sediment depth";
string const RASTER_PLOT_FINE_UNCONSOLIDATED_SEDIMENT_TITLE = "Unconsolidated fine sediment depth";
string const RASTER_PLOT_INTERVENTION_CLASS_TITLE = "Intervention class";
string const RASTER_PLOT_INTERVENTION_HEIGHT_TITLE = "Intervention height";
string const RASTER_PLOT_INUNDATION_MASK_TITLE = "Inundated area mask";
string const RASTER_PLOT_LANDFORM_TITLE = "Landform class";
string const RASTER_PLOT_NORMAL_PROFILE_TITLE = "Rasterized normal profiles";
string const RASTER_PLOT_POLYGON_GAIN_OR_LOSS_TITLE = "Polygon gain or loss of unconsolidated sediment";
string const RASTER_PLOT_POLYGON_TITLE = "Rasterized polygon boundaries";
string const RASTER_PLOT_POLYGON_UPDRIFT_OR_DOWNDRIFT_TITLE = "Polygon updrift or downdrift movement of unconsolidated sediment";
string const RASTER_PLOT_POTENTIAL_BEACH_EROSION_TITLE = "Potential (unconstrained) beach erosion depth";
string const RASTER_PLOT_POTENTIAL_PLATFORM_EROSION_MASK_TITLE = "Potential (unconstrained) shore platform erosion binary mask";
string const RASTER_PLOT_POTENTIAL_PLATFORM_EROSION_TITLE = "Potential (unconstrained) shore platform erosion depth";
string const RASTER_PLOT_SAND_CONSOLIDATED_SEDIMENT_TITLE = "Consolidated sand sediment depth";
string const RASTER_PLOT_SAND_UNCONSOLIDATED_SEDIMENT_TITLE = "Unconsolidated sand sediment depth";
string const RASTER_PLOT_SEA_DEPTH_TITLE = "Sea depth";
string const RASTER_PLOT_SEDIMENT_INPUT_EVENT_TITLE = "Sediment input event(s) since last GIS save";
string const RASTER_PLOT_SED_TOP_INC_TALUS_ELEV_TITLE = "Elevation of sediment top inc talus";
// string const RASTER_PLOT_SETUP_SURGE_FLOOD_MASK_TITLE = "Mask of setup-surge flood";
// string const RASTER_PLOT_SETUP_SURGE_RUNUP_FLOOD_MASK_TITLE = "Mask of setup-surge-runup flood";
string const RASTER_PLOT_SHADOW_DOWNDRIFT_ZONE_TITLE = "Downdrift of wave shadow zones";
string const RASTER_PLOT_SHADOW_ZONE_TITLE = "Wave shadow zones";
string const RASTER_PLOT_SLICE_TITLE = "Slice though layers at elevation = ";
string const RASTER_PLOT_SLOPE_FOR_CLIFF_TOE_TITLE = "Slope";
string const RASTER_PLOT_CONS_SED_SLOPE_TITLE = "Local slope of consolidated sediment";
string const RASTER_PLOT_SUSPENDED_SEDIMENT_TITLE = "Suspended sediment depth";
string const RASTER_PLOT_TALUS_TITLE = "Talus from cliff collapse";
string const RASTER_PLOT_TOP_ELEV_INC_SEA_TITLE = "Topmost elevation (sediment plus intervention plus se";
string const RASTER_PLOT_TOTAL_ACTUAL_BEACH_EROSION_TITLE = "Total actual (constrained) beach erosion depth";
string const RASTER_PLOT_TOTAL_ACTUAL_PLATFORM_EROSION_TITLE = "Total actual (constrained) shore platform erosion depth";
string const RASTER_PLOT_TOTAL_BEACH_DEPOSITION_TITLE = "Total beach deposition depth";
string const RASTER_PLOT_TOTAL_CLIFF_COLLAPSE_DEPOSITION_COARSE_TITLE = "Total depth of coarse talus from cliff collapse";
string const RASTER_PLOT_TOTAL_CLIFF_COLLAPSE_DEPOSITION_SAND_TITLE = "Total depth of sand talus from cliff collapse";
string const RASTER_PLOT_TOTAL_CLIFF_COLLAPSE_EROSION_COARSE_TITLE = "Total of cliff collapse erosion depth, coarse";
string const RASTER_PLOT_TOTAL_CLIFF_COLLAPSE_EROSION_FINE_TITLE = "Total of cliff collapse erosion depth, fine";
string const RASTER_PLOT_TOTAL_CLIFF_COLLAPSE_EROSION_SAND_TITLE = "Total of cliff collapse erosion depth, sand";
string const RASTER_PLOT_TOTAL_POTENTIAL_BEACH_EROSION_TITLE = "Total potential (unconstrained) beach erosion depth";
string const RASTER_PLOT_TOTAL_POTENTIAL_PLATFORM_EROSION_TITLE = "Total potential (unconstrained) shore platform erosion depth";
string const RASTER_PLOT_WAVE_FLOOD_LINE_TITLE = "Wave flood line";
string const RASTER_PLOT_WAVE_HEIGHT_TITLE = "Wave height";
string const RASTER_PLOT_WAVE_ANGLE_TITLE = "Wave orientation";

// GIS vector output user codes
string const VECTOR_ALL_OUTPUT_CODE = "all";
string const VECTOR_ALL_RIVER_FLOOD_OUTPUT_CODE = "all";
string const VECTOR_AVG_WAVE_ANGLE_AND_HEIGHT_CODE = "avg_wave_angle";
string const VECTOR_AVG_WAVE_ANGLE_AND_HEIGHT_NAME = "avg_wave_angle";
string const VECTOR_BREAKING_WAVE_HEIGHT_CODE = "breaking_wave_height";
string const VECTOR_BREAKING_WAVE_HEIGHT_NAME = "breaking_wave_height";
// string const VECTOR_CLIFF_EDGE_CODE = "cliff_edge";
// string const VECTOR_CLIFF_EDGE_NAME = "cliff_edge";
string const VECTOR_CLIFF_NOTCH_ACTIVE_CODE = "cliff_notch_active";
string const VECTOR_CLIFF_NOTCH_ACTIVE_NAME = "cliff_notch_active";
string const VECTOR_COAST_CODE = "coast";
string const VECTOR_COAST_CURVATURE_CODE = "coast_curvature";
string const VECTOR_COAST_CURVATURE_NAME = "coast_curvature";
string const VECTOR_COAST_NAME = "coast";
string const VECTOR_COAST_SWL_HIGHEST_NAME = "coast_SWL_highest";
string const VECTOR_COAST_SWL_LOWEST_NAME = "coast_SWL_lowest";
string const VECTOR_DEEP_WATER_WAVE_ANGLE_AND_HEIGHT_CODE = "deep_water_wave_angle";
string const VECTOR_DEEP_WATER_WAVE_ANGLE_AND_HEIGHT_NAME = "deep_water_wave_angle";
string const VECTOR_DOWNDRIFT_ZONE_BOUNDARY_CODE = "downdrift_boundary";
string const VECTOR_DOWNDRIFT_ZONE_BOUNDARY_NAME = "downdrift_boundary";
string const VECTOR_FLOOD_LINE_CODE = "flood_line";
string const VECTOR_FLOOD_LINE_NAME = "flood_line";
string const VECTOR_FLOOD_SWL_SETUP_LINE_CODE = "setup";
string const VECTOR_FLOOD_SWL_SETUP_LINE_NAME = "setup";
string const VECTOR_FLOOD_SWL_SETUP_SURGE_LINE_CODE = "setup_surge";
string const VECTOR_FLOOD_SWL_SETUP_SURGE_LINE_NAME = "setup_surge";
string const VECTOR_FLOOD_SWL_SETUP_SURGE_RUNUP_LINE_CODE = "setup_surge_runup";
string const VECTOR_FLOOD_SWL_SETUP_SURGE_RUNUP_LINE_NAME = "setup_surge_runup";
string const VECTOR_INVALID_NORMALS_CODE = "invalid_normals";
string const VECTOR_INVALID_NORMALS_NAME = "invalid_normals";
string const VECTOR_MEAN_WAVE_ENERGY_CODE = "mean_wave_energy";
string const VECTOR_MEAN_WAVE_ENERGY_NAME = "mean_wave_energy";
string const VECTOR_NORMALS_CODE = "normals";
string const VECTOR_NORMALS_NAME = "normals";
string const VECTOR_POLYGON_BOUNDARY_CODE = "polygon";
string const VECTOR_POLYGON_BOUNDARY_NAME = "polygon";
string const VECTOR_POLYGON_NODE_CODE = "polygon_node";
string const VECTOR_POLYGON_NODE_NAME = "polygon_node";
string const VECTOR_RUN_UP_CODE = "run_up";
string const VECTOR_RUN_UP_NAME = "run_up";
string const VECTOR_SHADOW_ZONE_BOUNDARY_CODE = "shadow_boundary";
string const VECTOR_SHADOW_ZONE_BOUNDARY_NAME = "shadow_boundary";
string const VECTOR_STORM_SURGE_CODE = "storm_surge";
string const VECTOR_STORM_SURGE_NAME = "storm_surge";
string const VECTOR_USUAL_OUTPUT_CODE = "usual";
string const VECTOR_WAVE_ANGLE_AND_HEIGHT_CODE = "wave_angle";
string const VECTOR_WAVE_ANGLE_AND_HEIGHT_NAME = "wave_angle";
string const VECTOR_WAVE_ENERGY_SINCE_COLLAPSE_CODE = "wave_energy";
string const VECTOR_WAVE_ENERGY_SINCE_COLLAPSE_NAME = "wave_energy";
string const VECTOR_WAVE_SETUP_CODE = "wave_setup";
string const VECTOR_WAVE_SETUP_NAME = "wave_setup";
string const VECTOR_WAVE_TRANSECT_POINTS_CODE = "wave_transect_points";
string const VECTOR_WAVE_TRANSECT_POINTS_NAME = "wave_transect_points";

// GIS vector output titles
string const VECTOR_PLOT_AVG_WAVE_ANGLE_AND_HEIGHT_TITLE = "Average wave orientation and height";
string const VECTOR_PLOT_BREAKING_WAVE_HEIGHT_TITLE = "Breaking wave height";
// string const VECTOR_PLOT_CLIFF_EDGE_TITLE = "Cliff edge";
string const VECTOR_PLOT_CLIFF_NOTCH_ACTIVE_TITLE = "Active cliff notch incision";
string const VECTOR_PLOT_COAST_CURVATURE_TITLE = "Coastline curvature";
string const VECTOR_PLOT_COAST_SWL_HIGHEST_TITLE = "Coastline at highest SWL";
string const VECTOR_PLOT_COAST_SWL_LOWEST_TITLE = "Coastline at lowest SWL";
string const VECTOR_PLOT_COAST_TITLE = "Coastline";
string const VECTOR_PLOT_DEEP_WATER_WAVE_ANGLE_AND_HEIGHT_TITLE = "Deep water wave orientation and height";
string const VECTOR_PLOT_DOWNDRIFT_ZONE_BOUNDARY_TITLE = "Downdrift zone boundary";
string const VECTOR_PLOT_FLOOD_LINE_TITLE = "Flood ";
string const VECTOR_PLOT_FLOOD_SWL_SETUP_LINE_TITLE = "SWL-Setup line";
string const VECTOR_PLOT_FLOOD_SWL_SETUP_SURGE_LINE_TITLE = "SWL-Setup-Surge line";
string const VECTOR_PLOT_FLOOD_SWL_SETUP_SURGE_RUNUP_LINE_TITLE = "SWL-Setup-Surge-Runup line";
string const VECTOR_PLOT_INVALID_NORMALS_TITLE = "INVALID coastline-normal profiles";
string const VECTOR_PLOT_MEAN_WAVE_ENERGY_TITLE = "Mean wave energy";
string const VECTOR_PLOT_NORMALS_TITLE = "Coastline-normal profiles";
string const VECTOR_PLOT_POLYGON_BOUNDARY_TITLE = "Polygons";
string const VECTOR_PLOT_POLYGON_NODES_TITLE = "Polygon nodes";
string const VECTOR_PLOT_RUN_UP_TITLE = "Run up";
string const VECTOR_PLOT_SHADOW_ZONE_BOUNDARY_TITLE = "Shadow zone boundary";
string const VECTOR_PLOT_STORM_SURGE_TITLE = "Storm surge";
string const VECTOR_PLOT_WAVE_ANGLE_AND_HEIGHT_TITLE = "Wave orientation and height";
string const VECTOR_PLOT_WAVE_ENERGY_SINCE_COLLAPSE_TITLE = "Wave energy since collapse";
string const VECTOR_PLOT_WAVE_SETUP_TITLE = "Wave setup";
string const VECTOR_PLOT_WAVE_TRANSECT_POINTS_TITLE = "Wave transect points (real and synthetic)";

// Time series codes
string const TIME_SERIES_BEACH_CHANGE_NET_CODE = "beach_change_net";
string const TIME_SERIES_BEACH_CHANGE_NET_NAME = "beach_change_net";
string const TIME_SERIES_BEACH_DEPOSITION_CODE = "beach_deposition";
string const TIME_SERIES_BEACH_DEPOSITION_NAME = "beach_deposition";
string const TIME_SERIES_BEACH_EROSION_CODE = "beach_erosion";
string const TIME_SERIES_BEACH_EROSION_NAME = "beach_erosion";
string const TIME_SERIES_CLIFF_COLLAPSE_DEPOSITION_CODE = "cliff_collapse_deposition";
string const TIME_SERIES_CLIFF_COLLAPSE_DEPOSITION_NAME = "cliff_collapse_deposition";
string const TIME_SERIES_CLIFF_COLLAPSE_EROSION_CODE = "cliff_collapse_erosion";
string const TIME_SERIES_CLIFF_COLLAPSE_EROSION_NAME = "cliff_collapse_erosion";
string const TIME_SERIES_CLIFF_NOTCH_ELEV_CODE = "cliff_notch";
string const TIME_SERIES_CLIFF_NOTCH_ELEV_NAME = "cliff_notch";
string const TIME_SERIES_FLOOD_SETUP_SURGE_CODE = "flood_setup_surge";
string const TIME_SERIES_FLOOD_SETUP_SURGE_NAME = "flood_setup_surge";
string const TIME_SERIES_FLOOD_SETUP_SURGE_RUNUP_CODE = "flood_setup_surge_runup";
string const TIME_SERIES_FLOOD_SETUP_SURGE_RUNUP_NAME = "flood_setup_surge_runup";
string const TIME_SERIES_PLATFORM_EROSION_CODE = "platform_erosion";
string const TIME_SERIES_PLATFORM_EROSION_NAME = "platform_erosion";
string const TIME_SERIES_SEA_AREA_CODE = "sea_area";
string const TIME_SERIES_SEA_AREA_NAME = "sea_area";
string const TIME_SERIES_SUSPENDED_SEDIMENT_CODE = "suspended";
string const TIME_SERIES_SUSPENDED_SEDIMENT_NAME = "suspended_sediment";
string const TIME_SERIES_SWL_CODE = "SWL";
string const TIME_SERIES_SWL_NAME = "SWL";

// CShore stuff
string const WAVE_ENERGY_FLUX = "wave_energy_flux";
string const WAVE_HEIGHT_X_FILENAME = "wave_height_x.csv";
string const WAVE_HEIGHT_Y_FILENAME = "wave_height_y.csv";
string const ACTIVE_ZONE_FILENAME = "activezone.csv";

//================================================ Globally-available functions =================================================
template <class T>
T tMax(T a, T b)
{
   return ((a > b) ? a : b);
}

template <class T>
T tMax(T a, T b, T c)
{
   T max = (a < b) ? b : a;
   return ((max < c) ? c : max);
}

template <class T>
T tMin(T a, T b)
{
   return ((a < b) ? a : b);
}

template <class T>
T tMin(T a, T b, T c)
{
   return (a < b ? (a < c ? a : c) : (b < c ? b : c));
}

template <class T>
T tAbs(T a)
{
   // From a posting dated 18 Nov 93 by rmartin@rcmcon.com (Robert Martin), archived in cpp_tips
   return ((a < 0) ? -a : a);
}

// template <class T>
// bool bIsBetween(T a, T b, T c)
// {
//    // Assumes b > c
//    return ((a >= b) && (a <= c));
// }

template <typename T>
string strDblToStr(const T &t)
{
   // From http://stackoverflow.com/questions/2125880/convert-float-to-stdstring-in-c
   ostringstream os;
   os << t;
   return os.str();
}

// ==============================================================================================================================
// For comparison of two floating-point numbers, with a specified accuracy
// ==============================================================================================================================
template <class T>
bool bFPIsEqual(const T d1, const T d2, const T dEpsilon)
{
   // Since the accuracy of floating-point numbers varies with their magnitude, we must compare them by using an accuracy threshold which is relative to the magnitude of the two numbers being compared. This is a blend of an example from Knuth's 'The Art of Computer Programming. Volume 1. Fundamental Algorithms' and a posting dated 18 Nov 93 by rmartin@rcmcon.com (Robert Martin), archived in cpp_tips

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"

   if ((0 == d1) && (tAbs(d2) < dEpsilon))
      return true;
   else if ((0 == d2) && (tAbs(d1) < dEpsilon))
      return true;
   else
      return ((tAbs(d1 - d2) < (dEpsilon * tAbs(d1))) ? true : false);

#pragma GCC diagnostic pop
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
// Definitions are in utilsglobal.cpp
double dRound(double const);
int nRound(double const);
// bool bIsWhole(double const);
bool bIsStringValidDouble(string &);
bool bIsStringValidInt(string &);

struct FillToWidth
{
   FillToWidth(char f, int w) : chFill(f), nWidth(w) {}
   char chFill;
   int nWidth;
};

//! Struct to hold wave data for a single transect/profile
struct TransectWaveData
{
   std::vector<double> VdX;
   std::vector<double> VdY;
   std::vector<double> VdHeightX;
   std::vector<double> VdHeightY;
   std::vector<bool> VbBreaking;
   int nCoastID;
   int nProfileID;
   bool bIsStartOrEndOfCoast;
};

ostream &operator<<(ostream &, const FillToWidth &);

// string strDbl(double const, int const);
string strDblRight(double const, int const, int const, bool const = true);
string strIntRight(int const, int const);
string strCentre(const char *, int const);
string strCentre(const string &, int const);
string strRight(const string &, int const);
string strRight(const char *, int const);
string strLeft(const string &, int const);
string strLeft(const char *, int const);
string strRightPerCent(double const, double const, int const, int const,
                       bool const = true);
#endif

//================================================= debugging stuff =============================================================
// #define CLOCKCHECK          // Uncomment to check CPU clock rollover settings

#endif // CME_H
