/*!
   \file utils.cpp
   \brief Utility routines
   \details TODO 001 A more detailed description of this routine.
   \author David Favis-Mortlock
   \author Andres Payo
   \author Wilf Chun
   \date 2026
   \copyright GNU General Public License
*/

/* ==============================================================================================================================
   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public  License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
==============================================================================================================================*/
#include <assert.h>

#ifdef _WIN32
#include <windows.h>       // Needed for CalcProcessStats()
#include <psapi.h>
#include <io.h>            // For isatty()
#elif defined __GNUG__
#include <sys/resource.h>  // Needed for CalcProcessStats()
#include <unistd.h>        // For isatty()
#include <sys/wait.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include <stdio.h>

#include <cstdlib>

#include <float.h>

#include <cctype>
using std::tolower;

#include <cmath>
using std::floor;

#include <ctime>
using std::clock;
using std::clock_t;
using std::difftime;
using std::localtime;
using std::time;

#include <ios>
using std::fixed;

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
using std::ios;

#include <iomanip>
using std::put_time;
using std::setprecision;
using std::setw;

#include <string>
using std::to_string;

#include <sstream>
using std::stringstream;

#include <algorithm>
using std::transform;

#include <gdal.h>

#include "cme.h"
#include "simulation.h"
#include "coast.h"
#include "2di_point.h"

//===============================================================================================================================
//! Returns the missing value code for double precision
//===============================================================================================================================
double CSimulation::dGetMissingValue(void) const
{
   return m_dMissingValue;
}

//===============================================================================================================================
//! Returns the still water level (SWL)
//===============================================================================================================================
double CSimulation::dGetThisIterSWL(void) const
{
   return m_dThisIterSWL;
}

//===============================================================================================================================
//! Returns the this-iteration total water level TODO 007 Finish surge and runup stuff
//===============================================================================================================================
double CSimulation::dGetThisIterTotWaterLevel(void) const
{
   return m_dThisIterDiffTotWaterLevel;
}

// //===============================================================================================================================
// //! Returns the max elevation of the beach above SWL
// //===============================================================================================================================
// double CSimulation::dGetMaxBeachElevAboveSWL (void) const
// {
// return m_dMaxBeachElevAboveSWL;
// }

//===============================================================================================================================
// Returns the cell side length
//===============================================================================================================================
// double CSimulation::dGetCellSide(void) const
// {
// return m_dCellSide;
// }

//===============================================================================================================================
//! Returns X grid max
//===============================================================================================================================
int CSimulation::nGetGridXMax(void) const
{
   return m_nXGridSize;
}

//===============================================================================================================================
//! Returns Y grid max
//===============================================================================================================================
int CSimulation::nGetGridYMax(void) const
{
   return m_nYGridSize;
}

//===============================================================================================================================
//! Returns D50 for fine sediment
//===============================================================================================================================
double CSimulation::dGetD50Fine(void) const
{
   return m_dD50Fine;
}

//===============================================================================================================================
//! Returns D50 for sand sediment
//===============================================================================================================================
double CSimulation::dGetD50Sand(void) const
{
   return m_dD50Sand;
}

//===============================================================================================================================
//! Returns D50 for coarse sediment
//===============================================================================================================================
double CSimulation::dGetD50Coarse(void) const
{
   return m_dD50Coarse;
}

//===============================================================================================================================
//! Handles command-line parameters
//===============================================================================================================================
int CSimulation::nHandleCommandLineParams(int nArg, char const* pcArgv[])
{
   if ((!isatty(fileno(stdout))) || (!isatty(fileno(stderr))))
      // Running with stdout or stderr not a tty, so either redirected or running as a background job. Ignore all command line parameters
      return RTN_OK;

   // Process the parameters following the name of the executable
   for (int i = 1; i < nArg; i++)
   {
      string strArg = pcArgv[i];
      strArg = strTrim(&strArg);

#ifdef _WIN32
      // Swap any forward slashes to backslashes
      strArg = pstrChangeToBackslash(&strArg);
#endif

      if (strArg.find("--gdal") != string::npos)
      {
         // User wants to know what GDAL raster drivers are available
         cout << GDAL_DRIVERS << endl
              << endl;

         for (int j = 0; j < GDALGetDriverCount(); j++)
         {
            GDALDriverH hDriver = GDALGetDriver(j);

            string strTmp(GDALGetDriverShortName(hDriver));
            strTmp.append("          ");
            strTmp.append(GDALGetDriverLongName(hDriver));

            cout << strTmp << endl;
         }

         return (RTN_HELP_ONLY);
      }

      else
      {
         if (strArg.find("--about") != string::npos)
         {
            // User wants information about CoastalME
            cout << ABOUT << endl;
            cout << THANKS << endl;

            return (RTN_HELP_ONLY);
         }

         else
         {
            if (strArg.find("--yaml") != string::npos)
            {
               // User wants to use YAML format for input datafile
               m_bYamlInputFormat = true;
            }

            else if (strArg.find("--home") != string::npos)
            {
               // Read in user defined runtime directory
               // string strTmp;

               // Find the position of '='
               size_t const pos = strArg.find('=');

               // Was '=' found?
               if (pos != string::npos)
               {
                  // Yes, so get the substring after '=' and assign it to the global variable
                  m_strCMEIni = strArg.substr(pos + 1);
               }

               else
               {
                  // No
                  cout << "No '=' found in the input string" << endl;
               }

               return (RTN_OK);
            }

            // TODO 049 Handle other command line parameters e.g. path to .ini file, path to datafile
            else
            {
               // Display usage information
               cout << USAGE << endl;
               cout << USAGE1 << endl;
               cout << USAGE2 << endl;
               cout << USAGE3 << endl;
               cout << USAGE4 << endl;
               cout << USAGE5 << endl;
               cout << USAGE6 << endl;

               return (RTN_HELP_ONLY);
            }
         }
      }
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Tells the user that we have started the simulation
//===============================================================================================================================
void CSimulation::AnnounceStart(void)
{
   cout << endl << PROGRAM_NAME << " for " << PLATFORM << " " << strGetBuild() << endl;
}

//===============================================================================================================================
//! Starts the clock ticking
//===============================================================================================================================
void CSimulation::StartClock(void)
{
   // First start the 'CPU time' clock ticking
   if (static_cast<clock_t>(-1) == clock())
   {
      // There's a problem with the clock, but continue anyway
      LogStream << NOTE << "CPU time not available" << endl;
      m_dCPUClock = -1;
   }
   else
   {
      // All OK, so get the time in m_dClkLast (this is needed to check for clock rollover on long runs)
      m_dClkLast = static_cast<double>(clock());
      m_dClkLast -= CLOCK_T_MIN; // necessary if clock_t is signed to make m_dClkLast unsigned
   }

   // And now get the actual time we started
   m_tSysStartTime = time(nullptr);
}

//===============================================================================================================================
//! Finds the folder (directory) in which the CoastalME executable is located
//===============================================================================================================================
bool CSimulation::bFindExeDir(char const* pcArg)
{
   string strTmp;
   char szBuf[BUF_SIZE] = "";

#ifdef _WIN32

   if (0 != GetModuleFileName(NULL, szBuf, BUF_SIZE))
      strTmp = szBuf;
   else
      // It failed, so try another approach
      strTmp = pcArg;
#else
   // char* pResult = getcwd(szBuf, BUF_SIZE);          // Used to use this, but what if cwd is not the same as the CoastalME dir?

   if (-1 != readlink("/proc/self/exe", szBuf, BUF_SIZE))
      strTmp = szBuf;
   else
      // It failed, so try another approach
      strTmp = pcArg;

#endif

   // Neither approach has worked, so give up
   if (strTmp.empty())
      return false;

   // It's OK, so trim off the executable's name
   int const nPos = static_cast<int>(strTmp.find_last_of(PATH_SEPARATOR));
   m_strCMEDir = strTmp.substr(0, nPos + 1); // Note that this must be terminated with a backslash

   return true;
}

//===============================================================================================================================
//! Tells the user about the licence, and the simulation's starting conditions
//===============================================================================================================================
void CSimulation::AnnounceLicenceAndStartConditions(void)
{
   cout << COPYRIGHT << endl << endl;
   cout << LINE << endl;
   cout << DISCLAIMER1 << endl;
   cout << DISCLAIMER2 << endl;
   cout << DISCLAIMER3 << endl;
   cout << DISCLAIMER4 << endl;
   cout << DISCLAIMER5 << endl;
   cout << DISCLAIMER6 << endl;
   cout << LINE << endl << endl;

   cout << START_NOTICE << strGetComputerName() << " at " << put_time(localtime(&m_tSysStartTime), "%T on %A %d %B %Y") << endl;

#ifdef _OPENMP
   cout << " - OpenMP is enabled: " << omp_get_max_threads() << " threads available";
   #pragma omp parallel
   {
      #pragma omp single
      cout << ", " << omp_get_num_threads() << " threads used" << endl;
   }
#else
   cout << " - OpenMP is not enabled" << endl;
#endif

   cout << INITIALIZING_NOTICE << endl;

}

//===============================================================================================================================
//! Given a string containing time units, this returns the appropriate multiplier
//===============================================================================================================================
double CSimulation::dGetTimeMultiplier(string const* strIn)
{
   // First decide what the time units are
   int const nTimeUnits = nDoTimeUnits(strIn);

   // Then return the correct multiplier, since m_dTimeStep is in hours
   switch (nTimeUnits)
   {
   case TIME_UNKNOWN:
      return TIME_UNKNOWN;
      break;

   case TIME_HOURS:
      return 1;               // Multiplier for hours
      break;

   case TIME_DAYS:
      return 24;              // Multiplier for days -> hours
      break;

   case TIME_WEEKS:
      return 24 * 7;          // Multiplier for days -> weeks
      break;

   case TIME_MONTHS:
      return 24 * 30.416667;  // Multiplier for months -> hours (assume 30 + 5/12 day months, no leap years)
      break;

   case TIME_YEARS:
      return 24 * 365.25;     // Multiplier for years -> hours
      break;
   }

   return 0;
}

//===============================================================================================================================
//! Given a string containing time units, this sets up the appropriate multiplier and display units for the simulation
//===============================================================================================================================
int CSimulation::nDoSimulationTimeMultiplier(string const* strIn)
{
   // First decide what the time units are
   int const nTimeUnits = nDoTimeUnits(strIn);

   // Next set up the correct multiplier, since m_dTimeStep is in hours
   switch (nTimeUnits)
   {
   case TIME_UNKNOWN:
      return RTN_ERR_TIME_UNITS;
      break;

   case TIME_HOURS:
      m_dDurationUnitsMult = 1;              // Multiplier for hours
      m_strDurationUnits = "hours";
      break;

   case TIME_DAYS:
      m_dDurationUnitsMult = 24;             // Multiplier for days -> hours
      m_strDurationUnits = "days";
      break;

   case TIME_MONTHS:
      m_dDurationUnitsMult = 24 * 30.416667; // Multiplier for months -> hours (assume 30 + 5/12 day months, no leap years)
      m_strDurationUnits = "months";
      break;

   case TIME_YEARS:
      m_dDurationUnitsMult = 24 * 365.25;    // Multiplier for years -> hours
      m_strDurationUnits = "years";
      break;
   }

   return RTN_OK;
}

//===============================================================================================================================
//! This finds time units in a string
//===============================================================================================================================
int CSimulation::nDoTimeUnits(string const* strIn)
{
   if (strIn->find("hour") != string::npos)
      return TIME_HOURS;
   else if (strIn->find("day") != string::npos)
      return TIME_DAYS;
   else if (strIn->find("week") != string::npos)
      return TIME_WEEKS;
   else if (strIn->find("month") != string::npos)
      return TIME_MONTHS;
   else if (strIn->find("year") != string::npos)
      return TIME_YEARS;
   else
      return TIME_UNKNOWN;
}

//===============================================================================================================================
//! Opens the log file
//===============================================================================================================================
bool CSimulation::bOpenLogFile(void)
{
   if (m_nLogFileDetail == 0)
   {
      LogStream.open("/dev/null", ios::out | ios::trunc);
      cout << "Warning: log file is not writting" << endl;
   }
   else
      LogStream.open(m_strLogFile.c_str(), ios::out | ios::trunc);

   if (!LogStream)
   {
      // Error, cannot open log file
      cerr << ERR << "cannot open " << m_strLogFile << " for output" << endl;
      return false;
   }

   return true;
}

//===============================================================================================================================
//! Tells the user that we are now reading the DEM file
//===============================================================================================================================
void CSimulation::AnnounceReadBasementDEM(void) const
{
   // Tell the user what is happening
#ifdef _WIN32
   cout << READING_BASEMENT << pstrChangeToForwardSlash(&m_strInitialBasementDEMFile) << endl;
#else
   cout << READING_BASEMENT << m_strInitialBasementDEMFile << endl;
#endif
}

//===============================================================================================================================
//! Tells the user that we are now allocating memory
//===============================================================================================================================
void CSimulation::AnnounceAllocateMemory(void)
{
   cout << ALLOCATE_MEMORY << endl;
}

//===============================================================================================================================
//! Tells the user that we are now adding layers
//===============================================================================================================================
void CSimulation::AnnounceAddLayers(void)
{
   // Tell the user what is happening
   cout << ADD_LAYERS << endl;
}

//===============================================================================================================================
//! Now reading raster GIS files
//===============================================================================================================================
void CSimulation::AnnounceReadRasterFiles(void)
{
   cout << READING_RASTER_FILES << endl;
}

//===============================================================================================================================
//! Now reading vector GIS files
//===============================================================================================================================
void CSimulation::AnnounceReadVectorFiles(void)
{
   cout << READING_VECTOR_FILES << endl;
}

//===============================================================================================================================
//! Tells the user that we are now reading the Landscape category GIS file
//===============================================================================================================================
void CSimulation::AnnounceReadLGIS(void) const
{
   // Tell the user what is happening
   if (! m_strInitialLandformFile.empty())
#ifdef _WIN32
      cout << READING_LANDFORM_FILE << pstrChangeToForwardSlash(&m_strInitialLandformFile) << endl;
#else
      cout << READING_LANDFORM_FILE << m_strInitialLandformFile << endl;
#endif
}

//===============================================================================================================================
//! Tells the user that we are now reading the Intervention class GIS file
//===============================================================================================================================
void CSimulation::AnnounceReadICGIS(void) const
{
   // Tell the user what is happening
   if (! m_strInterventionClassFile.empty())
#ifdef _WIN32
      cout << READING_INTERVENTION_CLASS_FILE << pstrChangeToForwardSlash(&m_strInterventionClassFile) << endl;
#else
      cout << READING_INTERVENTION_CLASS_FILE << m_strInterventionClassFile << endl;
#endif
}

//===============================================================================================================================
//! Tells the user that we are now reading the Intervention height GIS file
//===============================================================================================================================
void CSimulation::AnnounceReadIHGIS(void) const
{
   // Tell the user what is happening
   if (! m_strInterventionHeightFile.empty())
#ifdef _WIN32
      cout << READING_INTERVENTION_HEIGHT_FILE << pstrChangeToForwardSlash(&m_strInterventionHeightFile) << endl;
#else
      cout << READING_INTERVENTION_HEIGHT_FILE << m_strInterventionHeightFile << endl;
#endif
}

//===============================================================================================================================
//! Tells the user that we are now reading the deep water wave values GIS file
//===============================================================================================================================
void CSimulation::AnnounceReadDeepWaterWaveValuesGIS(void) const
{
   // Tell the user what is happening
   if (! m_strDeepWaterWavesInputFile.empty())
#ifdef _WIN32
      cout << READING_DEEP_WATER_WAVE_FILE << pstrChangeToForwardSlash(&m_strDeepWaterWavesInputFile) << endl;
#else
      cout << READING_DEEP_WATER_WAVE_FILE << m_strDeepWaterWavesInputFile << endl;
#endif
}

//===============================================================================================================================
//! Tells the user that we are now reading the sediment input events GIS file
//===============================================================================================================================
void CSimulation::AnnounceReadSedimentEventInputValuesGIS(void) const
{
   // Tell the user what is happening
   if (! m_strSedimentInputEventFile.empty())
#ifdef _WIN32
      cout << READING_SED_INPUT_EVENT_FILE << pstrChangeToForwardSlash(&m_strSedimentInputEventFile) << endl;
#else
      cout << READING_SED_INPUT_EVENT_FILE << m_strSedimentInputEventFile << endl;
#endif
}

//===============================================================================================================================
//! Tells the user that we are now reading the flood location GIS file
//===============================================================================================================================
void CSimulation::AnnounceReadFloodLocationGIS(void) const
{
   // Tell the user what is happening
   if (! m_strFloodLocationShapefile.empty())
#ifdef _WIN32
      cout << READING_FLOOD_LOCATION << pstrChangeToForwardSlash(&m_strFloodLocationShapefile) << endl;
#else
      cout << READING_FLOOD_LOCATION << m_strFloodLocationShapefile << endl;
#endif
}

//===============================================================================================================================
//! Tells the user that we are now reading the initial suspended sediment depth GIS file
//===============================================================================================================================
void CSimulation::AnnounceReadInitialSuspSedGIS(void) const
{
   // Tell the user what is happening
#ifdef _WIN32
   cout << READING_SUSPENDED_SEDIMENT_FILE << pstrChangeToForwardSlash(&m_strInitialSuspSedimentFile) << endl;
#else
   cout << READING_SUSPENDED_SEDIMENT_FILE << m_strInitialSuspSedimentFile << endl;
#endif
}

//===============================================================================================================================
//! Tells the user that we are now reading the initial fine unconsolidated sediment depth GIS file
//===============================================================================================================================
void CSimulation::AnnounceReadInitialFineUnconsSedGIS(int const nLayer) const
{
   // Tell the user what is happening
#ifdef _WIN32
   cout << READING_UNCONS_FINE_SEDIMENT_FILE << nLayer + 1 << "): " << pstrChangeToForwardSlash(&m_VstrInitialFineUnconsSedimentFile[nLayer]) << endl;
#else
   cout << READING_UNCONS_FINE_SEDIMENT_FILE << nLayer + 1 << "): " << m_VstrInitialFineUnconsSedimentFile[nLayer] << endl;
#endif
}

//===============================================================================================================================
//! Tells the user that we are now reading the initial sand unconsolidated sediment depth GIS file
//===============================================================================================================================
void CSimulation::AnnounceReadInitialSandUnconsSedGIS(int const nLayer) const
{
   // Tell the user what is happening
#ifdef _WIN32
   cout << READING_UNCONS_SAND_SEDIMENT_FILE << nLayer + 1 << "): " << pstrChangeToForwardSlash(&m_VstrInitialSandUnconsSedimentFile[nLayer]) << endl;
#else
   cout << READING_UNCONS_SAND_SEDIMENT_FILE << nLayer + 1 << "): " << m_VstrInitialSandUnconsSedimentFile[nLayer] << endl;
#endif
}

//===============================================================================================================================
//! Tells the user that we are now reading the initial coarse unconsolidated sediment depth GIS file
//===============================================================================================================================
void CSimulation::AnnounceReadInitialCoarseUnconsSedGIS(int const nLayer) const
{
   // Tell the user what is happening
#ifdef _WIN32
   cout << READING_UNCONS_COARSE_SEDIMENT_FILE << nLayer + 1 << "): " << pstrChangeToForwardSlash(&m_VstrInitialCoarseUnconsSedimentFile[nLayer]) << endl;
#else
   cout << READING_UNCONS_COARSE_SEDIMENT_FILE << nLayer + 1 << "): " << m_VstrInitialCoarseUnconsSedimentFile[nLayer] << endl;
#endif
}

//===============================================================================================================================
//! Tells the user that we are now reading the initial fine consolidated sediment depth GIS file
//===============================================================================================================================
void CSimulation::AnnounceReadInitialFineConsSedGIS(int const nLayer) const
{
   // Tell the user what is happening
#ifdef _WIN32
   cout << READING_CONS_FINE_SEDIMENT_FILE << nLayer + 1 << "): " << pstrChangeToForwardSlash(&m_VstrInitialFineConsSedimentFile[nLayer]) << endl;
#else
   cout << READING_CONS_FINE_SEDIMENT_FILE << nLayer + 1 << "): " << m_VstrInitialFineConsSedimentFile[nLayer] << endl;
#endif
}

//===============================================================================================================================
//! Tells the user that we are now reading the initial sand consolidated sediment depth GIS file
//===============================================================================================================================
void CSimulation::AnnounceReadInitialSandConsSedGIS(int const nLayer) const
{
   // Tell the user what is happening
#ifdef _WIN32
   cout << READING_CONS_SAND_SEDIMENT_FILE << nLayer + 1 << "): " << pstrChangeToForwardSlash(&m_VstrInitialSandConsSedimentFile[nLayer]) << endl;
#else
   cout << READING_CONS_SAND_SEDIMENT_FILE << nLayer + 1 << "): " << m_VstrInitialSandConsSedimentFile[nLayer] << endl;
#endif
}

//===============================================================================================================================
//! Tells the user that we are now reading the initial coarse consolidated sediment depth GIS file
//===============================================================================================================================
void CSimulation::AnnounceReadInitialCoarseConsSedGIS(int const nLayer) const
{
   // Tell the user what is happening
#ifdef _WIN32
   cout << READING_CONS_COARSE_SEDIMENT_FILE << nLayer + 1 << "): " << pstrChangeToForwardSlash(&m_VstrInitialCoarseConsSedimentFile[nLayer]) << endl;
#else
   cout << READING_CONS_COARSE_SEDIMENT_FILE << nLayer + 1 << "): " << m_VstrInitialCoarseConsSedimentFile[nLayer] << endl;
#endif
}

//===============================================================================================================================
//! Now reading tide data file
//===============================================================================================================================
void CSimulation::AnnounceReadTideData(void) const
{
#ifdef _WIN32
   cout << READING_TIDE_DATA_FILE << pstrChangeToForwardSlash(&m_strTideDataFile) << endl;
#else
   cout << READING_TIDE_DATA_FILE << m_strTideDataFile << endl;
#endif
}

//===============================================================================================================================
//! Now reading the SCAPE shape function file
//===============================================================================================================================
void CSimulation::AnnounceReadSCAPEShapeFunctionFile(void)
{
   cout << READING_SCAPE_SHAPE_FUNCTION_FILE << endl;
}

//===============================================================================================================================
//! Tells the user that we are now initializing
//===============================================================================================================================
void CSimulation::AnnounceFinalInitialization(void)
{
   // Tell the user what is happening
   cout << INITIALIZING_FINAL << endl;
}

//===============================================================================================================================
//! Tell the user that the simulation is now running
//===============================================================================================================================
void CSimulation::AnnounceIsRunning(void)
{
   cout << RUN_NOTICE << endl;
}

//===============================================================================================================================
//! Return a space-separated string containing the names of the raster GIS output files
//===============================================================================================================================
string CSimulation::strListRasterFiles(void) const
{
   string strTmp;

   if (m_bBasementElevSave)
   {
      strTmp.append(RASTER_BASEMENT_ELEVATION_CODE);
      strTmp.append(", ");
   }

   if (m_bSedIncTalusTopSurfSave)
   {
      strTmp.append(RASTER_SEDIMENT_TOP_CODE);
      strTmp.append(", ");
   }

   if (m_bTopSurfIncSeaSave)
   {
      strTmp.append(RASTER_TOP_ELEVATION_INC_SEA_CODE);
      strTmp.append(", ");
   }

   if (m_bTalusSave)
   {
      strTmp.append(RASTER_TALUS_CODE);
      strTmp.append(", ");
   }

   if (m_bSeaDepthSave)
   {
      strTmp.append(RASTER_SEA_DEPTH_CODE);
      strTmp.append(", ");
   }

   if (m_bAvgSeaDepthSave)
   {
      strTmp.append(RASTER_AVG_SEA_DEPTH_CODE);
      strTmp.append(", ");
   }

   if (m_bSeaMaskSave)
   {
      strTmp.append(RASTER_INUNDATION_MASK_CODE);
      strTmp.append(", ");
   }

   if (m_bWaveHeightSave)
   {
      strTmp.append(RASTER_WAVE_HEIGHT_CODE);
      strTmp.append(", ");
   }

   if (m_bWaveAngleSave)
   {
      strTmp.append(RASTER_WAVE_ANGLE_CODE);
      strTmp.append(", ");
   }

   if (m_bAvgWaveHeightSave)
   {
      strTmp.append(RASTER_AVG_WAVE_HEIGHT_CODE);
      strTmp.append(", ");
   }

   if (m_bBeachProtectionSave)
   {
      strTmp.append(RASTER_BEACH_PROTECTION_CODE);
      strTmp.append(", ");
   }

   if (m_bPotentialPlatformErosionSave)
   {
      strTmp.append(RASTER_POTENTIAL_PLATFORM_EROSION_CODE);
      strTmp.append(", ");
   }

   if (m_bPotentialPlatformErosionMaskSave)
   {
      strTmp.append(RASTER_POTENTIAL_PLATFORM_EROSION_MASK_CODE);
      strTmp.append(", ");
   }

   if (m_bBeachDepositionSave)
   {
      strTmp.append(RASTER_BEACH_DEPOSITION_CODE);
      strTmp.append(", ");
   }

   if (m_bTotalBeachDepositionSave)
   {
      strTmp.append(RASTER_TOTAL_BEACH_DEPOSITION_CODE);
      strTmp.append(", ");
   }

   if (m_bBeachMaskSave)
   {
      strTmp.append(RASTER_BEACH_MASK_CODE);
      strTmp.append(", ");
   }

   if (m_bActualPlatformErosionSave)
   {
      strTmp.append(RASTER_ACTUAL_PLATFORM_EROSION_CODE);
      strTmp.append(", ");
   }

   if (m_bTotalPotentialPlatformErosionSave)
   {
      strTmp.append(RASTER_TOTAL_POTENTIAL_PLATFORM_EROSION_CODE);
      strTmp.append(", ");
   }

   if (m_bTotalActualPlatformErosionSave)
   {
      strTmp.append(RASTER_TOTAL_ACTUAL_PLATFORM_EROSION_CODE);
      strTmp.append(", ");
   }

   if (m_bPotentialBeachErosionSave)
   {
      strTmp.append(RASTER_POTENTIAL_BEACH_EROSION_CODE);
      strTmp.append(", ");
   }

   if (m_bActualBeachErosionSave)
   {
      strTmp.append(RASTER_ACTUAL_BEACH_EROSION_CODE);
      strTmp.append(", ");
   }

   if (m_bTotalPotentialBeachErosionSave)
   {
      strTmp.append(RASTER_TOTAL_POTENTIAL_BEACH_EROSION_CODE);
      strTmp.append(", ");
   }

   if (m_bTotalActualBeachErosionSave)
   {
      strTmp.append(RASTER_TOTAL_ACTUAL_BEACH_EROSION_CODE);
      strTmp.append(", ");
   }

   if (m_bLandformSave)
   {
      strTmp.append(RASTER_LANDFORM_CODE);
      strTmp.append(", ");
   }

   if (m_bConsSedSlopeSave)
   {
      strTmp.append(RASTER_SLOPE_OF_CONSOLIDATED_SEDIMENT_CODE);
      strTmp.append(", ");
   }

   if (m_bInterventionClassSave)
   {
      strTmp.append(RASTER_INTERVENTION_CLASS_CODE);
      strTmp.append(", ");
   }

   if (m_bInterventionHeightSave)
   {
      strTmp.append(RASTER_INTERVENTION_HEIGHT_CODE);
      strTmp.append(", ");
   }

   if (m_bHaveFineSediment && m_bSuspSedSave)
   {
      strTmp.append(RASTER_SUSP_SED_CODE);
      strTmp.append(", ");
   }

   if (m_bHaveFineSediment && m_bAvgSuspSedSave)
   {
      strTmp.append(RASTER_AVG_SUSP_SED_CODE);
      strTmp.append(", ");
   }

   if (m_bHaveFineSediment && m_bFineUnconsSedSave)
   {
      strTmp.append(RASTER_FINE_UNCONS_CODE);
      strTmp.append(", ");
   }

   if (m_bHaveSandSediment && m_bSandUnconsSedSave)
   {
      strTmp.append(RASTER_SAND_UNCONS_CODE);
      strTmp.append(", ");
   }

   if (m_bHaveCoarseSediment && m_bCoarseUnconsSedSave)
   {
      strTmp.append(RASTER_COARSE_UNCONS_CODE);
      strTmp.append(", ");
   }

   if (m_bHaveFineSediment && m_bFineConsSedSave)
   {
      strTmp.append(RASTER_FINE_CONS_CODE);
      strTmp.append(", ");
   }

   if (m_bHaveSandSediment && m_bSandConsSedSave)
   {
      strTmp.append(RASTER_SAND_CONS_CODE);
      strTmp.append(", ");
   }

   if (m_bHaveCoarseSediment && m_bCoarseConsSedSave)
   {
      strTmp.append(RASTER_COARSE_CONS_CODE);
      strTmp.append(", ");
   }

   if (m_bRasterCoastlineSave)
   {
      strTmp.append(RASTER_COAST_CODE);
      strTmp.append(", ");
   }

   if (m_bRasterNormalProfileSave)
   {
      strTmp.append(RASTER_COAST_NORMAL_CODE);
      strTmp.append(", ");
   }

   if (m_bActiveZoneSave)
   {
      strTmp.append(RASTER_ACTIVE_ZONE_CODE);
      strTmp.append(", ");
   }

   if (m_bRasterPolygonSave)
   {
      strTmp.append(RASTER_POLYGON_CODE);
      strTmp.append(", ");
   }

   if (m_bPotentialPlatformErosionMaskSave)
   {
      strTmp.append(RASTER_POTENTIAL_PLATFORM_EROSION_MASK_CODE);
      strTmp.append(", ");
   }

   if (m_bSedimentInputEventSave)
   {
      strTmp.append(RASTER_SEDIMENT_INPUT_EVENT_CODE);
      strTmp.append(", ");
   }

   // Remove the trailing comma and space
   if (strTmp.size() > 2)
      strTmp.resize(strTmp.size() - 2);

   return strTmp;
}

//===============================================================================================================================
//! Return a space-separated string containing the names of the vector GIS output files
//===============================================================================================================================
string CSimulation::strListVectorFiles(void) const
{
   string strTmp;

   if (m_bCoastSave)
   {
      strTmp.append(VECTOR_COAST_CODE);
      strTmp.append(", ");
   }

   if (m_bNormalsSave)
   {
      strTmp.append(VECTOR_NORMALS_CODE);
      strTmp.append(", ");
   }

   if (m_bInvalidNormalsSave)
   {
      strTmp.append(VECTOR_INVALID_NORMALS_CODE);
      strTmp.append(", ");
   }

   if (m_bWaveAngleAndHeightSave)
   {
      strTmp.append(VECTOR_WAVE_ANGLE_AND_HEIGHT_CODE);
      strTmp.append(", ");
   }

   if (m_bAvgWaveAngleAndHeightSave)
   {
      strTmp.append(VECTOR_AVG_WAVE_ANGLE_AND_HEIGHT_CODE);
      strTmp.append(", ");
   }

   if (m_bCoastCurvatureSave)
   {
      strTmp.append(VECTOR_COAST_CURVATURE_CODE);
      strTmp.append(", ");
   }

   if (m_bWaveEnergySinceCollapseSave)
   {
      strTmp.append(VECTOR_WAVE_ENERGY_SINCE_COLLAPSE_CODE);
      strTmp.append(", ");
   }

   if (m_bMeanWaveEnergySave)
   {
      strTmp.append(VECTOR_MEAN_WAVE_ENERGY_CODE);
      strTmp.append(", ");
   }

   if (m_bBreakingWaveHeightSave)
   {
      strTmp.append(VECTOR_BREAKING_WAVE_HEIGHT_CODE);
      strTmp.append(", ");
   }

   if (m_bPolygonNodeSave)
   {
      strTmp.append(VECTOR_POLYGON_NODE_CODE);
      strTmp.append(", ");
   }

   if (m_bPolygonBoundarySave)
   {
      strTmp.append(VECTOR_POLYGON_BOUNDARY_CODE);
      strTmp.append(", ");
   }

   if (m_bCliffNotchSave)
   {
      strTmp.append(VECTOR_CLIFF_NOTCH_ACTIVE_CODE);
      strTmp.append(", ");
   }

   if (m_bShadowBoundarySave)
   {
      strTmp.append(VECTOR_SHADOW_ZONE_BOUNDARY_CODE);
      strTmp.append(", ");
   }

   if (m_bShadowDowndriftBoundarySave)
   {
      strTmp.append(VECTOR_DOWNDRIFT_ZONE_BOUNDARY_CODE);
      strTmp.append(", ");
   }

   if (m_bDeepWaterWaveAngleAndHeightSave)
   {
      strTmp.append(VECTOR_DEEP_WATER_WAVE_ANGLE_AND_HEIGHT_CODE);
      strTmp.append(", ");
   }

   // Remove the trailing comma and space
   if (strTmp.size() > 2)
      strTmp.resize(strTmp.size() - 2);

   return strTmp;
}

//===============================================================================================================================
//! Return a space-separated string containing the names of the time series output files
//===============================================================================================================================
string CSimulation::strListTSFiles(void) const
{
   string strTmp;

   if (m_bSeaAreaTSSave)
   {
      strTmp.append(TIME_SERIES_SEA_AREA_CODE);
      strTmp.append(", ");
   }

   if (m_bSWLTSSave)
   {
      strTmp.append(TIME_SERIES_SWL_CODE);
      strTmp.append(", ");
   }

   if (m_bActualPlatformErosionTSSave)
   {
      strTmp.append(TIME_SERIES_PLATFORM_EROSION_CODE);
      strTmp.append(", ");
   }

   if (m_bCliffCollapseErosionTSSave)
   {
      strTmp.append(TIME_SERIES_CLIFF_COLLAPSE_EROSION_CODE);
      strTmp.append(", ");
   }

   if (m_bCliffCollapseDepositionTSSave)
   {
      strTmp.append(TIME_SERIES_CLIFF_COLLAPSE_DEPOSITION_CODE);
      strTmp.append(", ");
   }

   if (m_bCliffCollapseNetTSSave)
   {
      strTmp.append(TIME_SERIES_CLIFF_COLLAPSE_NET_CODE);
      strTmp.append(", ");
   }

   if (m_bBeachErosionTSSave)
   {
      strTmp.append(TIME_SERIES_BEACH_EROSION_CODE);
      strTmp.append(", ");
   }

   if (m_bBeachDepositionTSSave)
   {
      strTmp.append(TIME_SERIES_BEACH_DEPOSITION_CODE);
      strTmp.append(", ");
   }

   if (m_bBeachSedimentChangeNetTSSave)
   {
      strTmp.append(TIME_SERIES_BEACH_CHANGE_NET_CODE);
      strTmp.append(", ");
   }

   if (m_bSuspSedTSSave)
   {
      strTmp.append(TIME_SERIES_SUSPENDED_SEDIMENT_CODE);
      strTmp.append(", ");
   }

   if (m_bFloodSetupSurgeTSSave)
   {
      strTmp.append(TIME_SERIES_FLOOD_SETUP_SURGE_CODE);
      strTmp.append(", ");
   }

   if (m_bFloodSetupSurgeRunupTSSave)
   {
      strTmp.append(TIME_SERIES_FLOOD_SETUP_SURGE_RUNUP_CODE);
      strTmp.append(", ");
   }

   if (m_bCliffNotchElevTSSave)
   {
      strTmp.append(TIME_SERIES_CLIFF_NOTCH_ELEV_CODE);
      strTmp.append(", ");
   }

   // Remove the trailing comma and space
   if (strTmp.size() > 2)
      strTmp.resize(strTmp.size() - 2);

   return strTmp;
}

//===============================================================================================================================
//! This member function intialises the time series files
//===============================================================================================================================
bool CSimulation::bSetUpTSFiles(void)
{
   string strTSFile;

   if (m_bSeaAreaTSSave)
   {
      // Start with wetted area
      strTSFile = m_strOutPath;
      strTSFile.append(TIME_SERIES_SEA_AREA_NAME);
      strTSFile.append(CSVEXT);

      // Open sea area time-series CSV file
      SeaAreaTSStream.open(strTSFile.c_str(), ios::out | ios::trunc);
      if (! SeaAreaTSStream)
      {
         // Error, cannot open wetted area  time-series file
         cerr << ERR << "cannot open " << strTSFile << " for output" << endl;
         return false;
      }
   }

   if (m_bSWLTSSave)
   {
      // Now SWL
      strTSFile = m_strOutPath;
      strTSFile.append(TIME_SERIES_SWL_NAME);
      strTSFile.append(CSVEXT);

      // Open SWL time-series CSV file
      SWLTSStream.open(strTSFile.c_str(), ios::out | ios::trunc);
      if (! SWLTSStream)
      {
         // Error, cannot open SWL time-series file
         cerr << ERR << "cannot open " << strTSFile << " for output" << endl;
         return false;
      }
   }

   if (m_bActualPlatformErosionTSSave)
   {
      // Erosion (fine, sand, coarse)
      strTSFile = m_strOutPath;
      strTSFile.append(TIME_SERIES_PLATFORM_EROSION_NAME);
      strTSFile.append(CSVEXT);

      // Open erosion time-series CSV file
      PlatformErosionTSStream.open(strTSFile.c_str(), ios::out | ios::trunc);
      if (! PlatformErosionTSStream)
      {
         // Error, cannot open erosion time-series file
         cerr << ERR << "cannot open " << strTSFile << " for output" << endl;
         return false;
      }
   }

   if (m_bCliffCollapseErosionTSSave)
   {
      // Erosion due to cliff collapse
      strTSFile = m_strOutPath;
      strTSFile.append(TIME_SERIES_CLIFF_COLLAPSE_EROSION_NAME);
      strTSFile.append(CSVEXT);

      // Open cliff collapse erosion time-series CSV file
      CliffCollapseErosionTSStream.open(strTSFile.c_str(), ios::out | ios::trunc);
      if (! CliffCollapseErosionTSStream)
      {
         // Error, cannot open cliff collapse erosion time-series file
         cerr << ERR << "cannot open " << strTSFile << " for output" << endl;
         return false;
      }
   }

   if (m_bCliffCollapseDepositionTSSave)
   {
      // Deposition due to cliff collapse
      strTSFile = m_strOutPath;
      strTSFile.append(TIME_SERIES_CLIFF_COLLAPSE_DEPOSITION_NAME);
      strTSFile.append(CSVEXT);

      // Open cliff collapse deposition time-series CSV file
      CliffCollapseDepositionTSStream.open(strTSFile.c_str(), ios::out | ios::trunc);
      if (! CliffCollapseDepositionTSStream)
      {
         // Error, cannot open cliff collapse deposition time-series file
         cerr << ERR << "cannot open " << strTSFile << " for output" << endl;
         return false;
      }
   }

   if (m_bCliffCollapseNetTSSave)
   {
      // Net change in unconsolidated sediment due to cliff collapse
      strTSFile = m_strOutPath;
      strTSFile.append(TIME_SERIES_CLIFF_COLLAPSE_NET_NAME);
      strTSFile.append(CSVEXT);

      // Open net cliff collapse time-series CSV file
      CliffCollapseNetChangeTSStream.open(strTSFile.c_str(), ios::out | ios::trunc);
      if (! CliffCollapseNetChangeTSStream)
      {
         // Error, cannot open net cliff collapse time-series file
         cerr << ERR << "cannot open " << strTSFile << " for output" << endl;
         return false;
      }
   }

   if (m_bBeachErosionTSSave)
   {
      // Beach erosion
      strTSFile = m_strOutPath;
      strTSFile.append(TIME_SERIES_BEACH_EROSION_NAME);
      strTSFile.append(CSVEXT);

      // Open beach erosion time-series CSV file
      BeachErosionTSStream.open(strTSFile.c_str(), ios::out | ios::trunc);
      if (! BeachErosionTSStream)
      {
         // Error, cannot open beach erosion time-series file
         cerr << ERR << "cannot open " << strTSFile << " for output" << endl;
         return false;
      }
   }

   if (m_bBeachDepositionTSSave)
   {
      // Beach deposition
      strTSFile = m_strOutPath;
      strTSFile.append(TIME_SERIES_BEACH_DEPOSITION_NAME);
      strTSFile.append(CSVEXT);

      // Open beach deposition time-series CSV file
      BeachDepositionTSStream.open(strTSFile.c_str(), ios::out | ios::trunc);
      if (! BeachDepositionTSStream)
      {
         // Error, cannot open beach deposition time-series file
         cerr << ERR << "cannot open " << strTSFile << " for output" << endl;
         return false;
      }
   }

   if (m_bBeachSedimentChangeNetTSSave)
   {
      // Beach sediment change
      strTSFile = m_strOutPath;
      strTSFile.append(TIME_SERIES_BEACH_CHANGE_NET_NAME);
      strTSFile.append(CSVEXT);

      // Open net beach sediment change time-series CSV file
      BeachSedimentNetChangeTSStream.open(strTSFile.c_str(), ios::out | ios::trunc);
      if (! BeachSedimentNetChangeTSStream)
      {
         // Error, cannot open beach sediment change time-series file
         cerr << ERR << "cannot open " << strTSFile << " for output" << endl;
         return false;
      }
   }

   if (m_bSuspSedTSSave)
   {
      // Sediment load
      strTSFile = m_strOutPath;
      strTSFile.append(TIME_SERIES_SUSPENDED_SEDIMENT_NAME);
      strTSFile.append(CSVEXT);

      // Open sediment load time-series CSV file
      FineSedSuspensionTSStream.open(strTSFile.c_str(), ios::out | ios::trunc);
      if (! FineSedSuspensionTSStream)
      {
         // Error, cannot open sediment load time-series file
         cerr << ERR << "cannot open " << strTSFile << " for output" << endl;
         return false;
      }
   }

   if (m_bFloodSetupSurgeTSSave)
   {
      // Sediment load
      strTSFile = m_strOutPath;
      strTSFile.append(TIME_SERIES_FLOOD_SETUP_SURGE_CODE);
      strTSFile.append(CSVEXT);

      // Open sediment load time-series CSV file
      FloodSetupSurgeTSStream.open(strTSFile.c_str(), ios::out | ios::trunc);
      if (! FloodSetupSurgeTSStream)
      {
         // Error, cannot open sediment load time-series file
         cerr << ERR << "cannot open " << strTSFile << " for output" << endl;
         return false;
      }
   }

   if (m_bFloodSetupSurgeRunupTSSave)
   {
      // Sediment load
      strTSFile = m_strOutPath;
      strTSFile.append(TIME_SERIES_FLOOD_SETUP_SURGE_RUNUP_CODE);
      strTSFile.append(CSVEXT);

      // Open sediment load time-series CSV file
      FloodSetupSurgeRunupTSStream.open(strTSFile.c_str(), ios::out | ios::trunc);
      if (! FloodSetupSurgeRunupTSStream)
      {
         // Error, cannot open sediment load time-series file
         cerr << ERR << "cannot open " << strTSFile << " for output" << endl;
         return false;
      }
   }

   if (m_bCliffNotchElevTSSave)
   {
      // Elevation of cliff notch
      strTSFile = m_strOutPath;
      strTSFile.append(TIME_SERIES_CLIFF_NOTCH_ELEV_NAME);
      strTSFile.append(CSVEXT);

      // Open cliff notch elevation time-series CSV file
      CliffNotchElevTSStream.open(strTSFile.c_str(), ios::out | ios::trunc);
      if (! CliffNotchElevTSStream)
      {
         // Error, cannot open cliff notch elevation time-series file
         cerr << ERR << "cannot open " << strTSFile << " for output" << endl;
         return false;
      }
   }

   return true;
}

//===============================================================================================================================
//! Checks to see if the simulation has gone on too long, amongst other things
//===============================================================================================================================
bool CSimulation::bTimeToQuit(void)
{
   // Add timestep to the total time simulated so far
   m_dSimElapsed += m_dTimeStep;

   if (m_dSimElapsed >= m_dSimDuration)
   {
      // It is time to quit
      m_dSimElapsed = m_dSimDuration;
      AnnounceProgress();
      return true;
   }

   // Not quitting, so increment the timestep count, and recalc total timesteps
   m_ulIter++;
   m_ulTotTimestep = static_cast<unsigned long>(dRound(m_dSimDuration / m_dTimeStep));

   // Check to see if we have done CLOCK_CHECK_ITERATION timesteps: if so, it is time to reset the CPU time running total in case the clock() function later rolls over
   if (0 == m_ulIter % CLOCK_CHECK_ITERATION)
      DoCPUClockReset();

   // Not yet time to quit
   return false;
}

//===============================================================================================================================
//! Returns a string, hopefully giving the name of the computer on which the simulation is running
//===============================================================================================================================
string CSimulation::strGetComputerName(void)
{
   string strComputerName;

#ifdef _WIN32
   // Being compiled to run under Windows, either by MS VC++, Borland C++, or Cygwin
   strComputerName = getenv("COMPUTERNAME");
#else
   // Being compiled for another platform; assume for Linux-Unix
   char szHostName[BUF_SIZE] = "";
   gethostname(szHostName, BUF_SIZE);

   strComputerName = szHostName;

   if (strComputerName.empty())
      strComputerName = "Unknown Computer";

#endif

   return strComputerName;
}

//===============================================================================================================================
//! Resets the CPU clock timer to prevent it 'rolling over', as can happen during long runs. This is a particularly problem under Unix systems where the value returned by clock() is defined in microseconds (for compatibility with systems that have CPU clocks with much higher resolution) i.e. CLOCKS_PER_SEC is 1000000 rather than the more usual 1000. In this case, the value returned from clock() will wrap around after accumulating only 2147 seconds of CPU time (about 36 minutes).
//===============================================================================================================================
void CSimulation::DoCPUClockReset(void)
{
   if (static_cast<clock_t>(-1) == clock())
   {
      // Error
      LogStream << "CPU time not available" << endl;
      m_dCPUClock = -1;
      return;
   }

   // OK, so carry on
   double dClkThis = static_cast<double>(clock());
   dClkThis -= CLOCK_T_MIN; // necessary when clock_t is signed, to make dClkThis unsigned

   if (dClkThis < m_dClkLast)
   {
      // Clock has 'rolled over'
      m_dCPUClock += (CLOCK_T_RANGE + 1 - m_dClkLast); // this elapsed before rollover
      m_dCPUClock += dClkThis;                         // this elapsed after rollover

#ifdef CLOCKCHECK
      // For debug purposes
      LogStream << "Rolled over: dClkThis=" << dClkThis << " m_dClkLast=" << m_dClkLast << endl
                << "\t"
                << " before rollover=" << (CLOCK_T_RANGE + 1 - m_dClkLast) << endl
                << "\t"
                << " after rollover=" << dClkThis << endl
                << "\t"
                << " ADDED=" << (CLOCK_T_RANGE + 1 - m_dClkLast + dClkThis) << endl;
#endif
   }

   else
   {
      // No rollover
      m_dCPUClock += (dClkThis - m_dClkLast);

#ifdef CLOCKCHECK
      // For debug purposes
      LogStream << "No rollover: dClkThis=" << dClkThis << " m_dClkLast=" << m_dClkLast << " ADDED=" << dClkThis - m_dClkLast << endl;
#endif
   }

   // Reset for next time
   m_dClkLast = dClkThis;
}

//===============================================================================================================================
//! Announce the end of the simulation
//===============================================================================================================================
void CSimulation::AnnounceSimEnd(void)
{
   cout << endl << FINAL_OUTPUT << endl;
}

// //===============================================================================================================================
// //! Calculates and displays time elapsed in terms of CPU time and real time, also calculates time per timestep in terms of both CPU time and real time
// //===============================================================================================================================
// void CSimulation::CalcTime(double const dRunLength)
// {
//    // Reset CPU count for last time
//    DoCPUClockReset();
//
//    if (! bFPIsEqual(m_dCPUClock, -1.0, TOLERANCE))
//    {
//       // Calculate CPU time in secs
//       double const dDuration = m_dCPUClock / CLOCKS_PER_SEC;
//
//       // And write CPU time out to OutStream and LogStream
//       OutStream << "CPU time elapsed: " << strDispTime(dDuration, false, true);
//       LogStream << "CPU time elapsed: " << strDispTime(dDuration, false, true);
//
//       // Calculate CPU time per timestep
//       double const dPerTimestep = dDuration / static_cast<double>(m_ulTotTimestep);
//
//       // And write CPU time per timestep to OutStream and LogStream
//       OutStream << fixed << setprecision(4) << " (" << dPerTimestep << " per timestep)" << endl;
//       LogStream << fixed << setprecision(4) << " (" << dPerTimestep << " per timestep)" << endl;
//
//       // Calculate ratio of CPU time to time simulated
//       OutStream << resetiosflags(ios::floatfield);
//       OutStream << fixed << setprecision(0) << "In terms of CPU time, this is ";
//       LogStream << resetiosflags(ios::floatfield);
//       LogStream << fixed << setprecision(0) << "In terms of CPU time, this is ";
//
//       if (dDuration > dRunLength)
//       {
//          OutStream << dDuration / dRunLength << " x slower than reality" << endl;
//          LogStream << dDuration / dRunLength << " x slower than reality" << endl;
//       }
//
//       else
//       {
//          OutStream << dRunLength / dDuration << " x faster than reality" << endl;
//          LogStream << dRunLength / dDuration << " x faster than reality" << endl;
//       }
//    }
//
//    // Calculate run time
//    double const dDuration = difftime(m_tSysEndTime, m_tSysStartTime);
//
//    // And write run time out to OutStream and LogStream
//    OutStream << "Run time elapsed: " << strDispTime(dDuration, false, false);
//    LogStream << "Run time elapsed: " << strDispTime(dDuration, false, false);
//
//    // Calculate run time per timestep
//    double const dPerTimestep = dDuration / static_cast<double>(m_ulTotTimestep);
//
//    // And write run time per timestep to OutStream and LogStream
//    OutStream << resetiosflags(ios::floatfield);
//    OutStream << " (" << fixed << setprecision(4) << dPerTimestep << " per timestep)" << endl;
//    LogStream << resetiosflags(ios::floatfield);
//    LogStream << " (" << fixed << setprecision(4) << dPerTimestep << " per timestep)" << endl;
//
//    // Calculate ratio of run time to time simulated
//    OutStream << fixed << setprecision(0) << "In terms of run time, this is ";
//    LogStream << fixed << setprecision(0) << "In terms of run time, this is ";
//
//    if (dDuration > dRunLength)
//    {
//       OutStream << dDuration / dRunLength << " x slower than reality" << endl;
//       LogStream << dDuration / dRunLength << " x slower than reality" << endl;
//    }
//
//    else
//    {
//       OutStream << dRunLength / dDuration << " x faster than reality" << endl;
//       LogStream << dRunLength / dDuration << " x faster than reality" << endl;
//    }
// }

//===============================================================================================================================
//! strDispSimTime returns a string formatted as year Julian_day hour, given a parameter in hours
//===============================================================================================================================
string CSimulation::strDispSimTime(const double dTimeIn)
{
   // Make sure no negative times
   double dTmpTime = tMax(dTimeIn, 0.0);

   string strTime;

   // Constants
   double const dHoursInYear = 24 * 365; // it was 365.25
   double const dHoursInDay = 24;

   // Display years
   if (dTmpTime >= dHoursInYear)
   {
      double const dYears = floor(dTmpTime / dHoursInYear);
      dTmpTime -= (dYears * dHoursInYear);

      strTime = to_string(static_cast<int>(dYears));
      strTime.append("y ");
   }
   else
      strTime = "0y ";

   // Display Julian days
   if (dTmpTime >= dHoursInDay)
   {
      double const dJDays = floor(dTmpTime / dHoursInDay);
      dTmpTime -= (dJDays * dHoursInDay);

      stringstream ststrTmp;
      ststrTmp << FillToWidth('0', 3) << static_cast<int>(dJDays);
      strTime.append(ststrTmp.str());
      strTime.append("d ");
   }
   else
      strTime.append("000d ");

   // Display hours
   stringstream ststrTmp;
   ststrTmp << FillToWidth('0', 2) << static_cast<int>(dTmpTime);
   strTime.append(ststrTmp.str());
   strTime.append("h");

   return strTime;
}

//===============================================================================================================================
//! strDispTime returns a string formatted as h:mm:ss, given a parameter in seconds, with rounding and fractions of a second if desired
//===============================================================================================================================
string CSimulation::strDispTime(const double dTimeIn, const bool bRound, const bool bFrac)
{
   // Make sure no negative times
   double dTime = tMax(dTimeIn, 0.0);

   string strTime;

   if (bRound)
      dTime = dRound(dTime);

   unsigned long ulTimeIn = static_cast<unsigned long>(floor(dTime));
   dTime -= static_cast<double>(ulTimeIn);

   // Hours
   if (ulTimeIn >= 3600)
   {
      // Display some hours
      unsigned long const ulHours = ulTimeIn / 3600ul;
      ulTimeIn -= (ulHours * 3600ul);

      strTime = to_string(ulHours);
      strTime.append(":");
   }
   else
      strTime = "0:";

   // Minutes
   if (ulTimeIn >= 60)
   {
      // display some minutes
      unsigned long const ulMins = ulTimeIn / 60ul;
      ulTimeIn -= (ulMins * 60ul);

      stringstream ststrTmp;
      ststrTmp << FillToWidth('0', 2) << ulMins;
      strTime.append(ststrTmp.str());
      strTime.append(":");
   }
   else
      strTime.append("00:");

   // Seconds
   stringstream ststrTmp;
   ststrTmp << FillToWidth('0', 2) << ulTimeIn;
   strTime.append(ststrTmp.str());

   if (bFrac)
   {
      // Fractions of a second
      strTime.append(".");
      ststrTmp.clear();
      ststrTmp.str(string());
      ststrTmp << FillToWidth('0', 2) << static_cast<unsigned long>(dTime * 100);
      strTime.append(ststrTmp.str());
   }

   return strTime;
}

//===============================================================================================================================
//! Returns the date and time on which the program was compiled
//===============================================================================================================================
string CSimulation::strGetBuild(void)
{
   string strBuild("(");
   strBuild.append(__TIME__);
   strBuild.append(" ");
   strBuild.append(__DATE__);
#ifdef _DEBUG
   strBuild.append(" DEBUG");
#endif
   strBuild.append(" build)");

   return strBuild;
}

//===============================================================================================================================
//! Displays information regarding the progress of the simulation
//===============================================================================================================================
void CSimulation::AnnounceProgress(void)
{
   if (isatty(fileno(stdout)))
   {
      // Stdout is connected to a tty, so not running as a background job
      static double sdElapsed = 0;
      static double sdToGo = 0;
      static double sdExpectedRuntime = 0;
      time_t const tNow = time(nullptr);

      // Calculate time elapsed and remaining
      sdElapsed = difftime(tNow, m_tSysStartTime);
      sdExpectedRuntime = (sdElapsed * m_dSimDuration / m_dSimElapsed);
      sdToGo = sdExpectedRuntime  - sdElapsed;

      // Tell the user about progress (note need to make several separate calls to cout here, or MS VC++ compiler appears to get confused)
      cout << SIMULATING << strDispSimTime(m_dSimElapsed);
      cout << fixed << setprecision(3) << setw(9) << 100 * m_dSimElapsed / m_dSimDuration;
      cout << "%   (elapsed " << strDispTime(sdElapsed, false, false) ;

      cout << " remaining " << strDispTime(sdToGo, false, false);
      cout << " total expected " << strDispTime(sdExpectedRuntime, false, false) << ")  ";

      // Add a marker for GIS saves etc.
      cout << setw(6);
      if (m_bSaveGISThisIter)
         cout << "GIS" + to_string(m_nGISSave) + "    ";
      else if (m_bSedimentInputThisIter)
         cout << "SED INPUT ";
      else
         cout << setw(10) << SPACE;

      cout.flush();
   }
}

//===============================================================================================================================
//! This calculates and displays process statistics
//===============================================================================================================================
void CSimulation::CalcProcessStats(void)
{
   string const NA = "Not available";

   OutStream << endl;
   OutStream << "Process statistics" << endl;
   OutStream << "------------------" << endl;

#ifdef _WIN32
   // First, find out which version of Windows we are running under
   OSVERSIONINFOEX osvi;
   BOOL bOsVersionInfoEx;

   ZeroMemory(&osvi, sizeof(OSVERSIONINFOEX)); // fill this much memory with zeros
   osvi.dwOSVersionInfoSize = sizeof(OSVERSIONINFOEX);

   if (!(bOsVersionInfoEx = GetVersionEx((OSVERSIONINFO*)&osvi)))
   {
      // OSVERSIONINFOEX didn't work so try OSVERSIONINFO instead
      osvi.dwOSVersionInfoSize = sizeof(OSVERSIONINFO);

      if (!GetVersionEx((OSVERSIONINFO*)&osvi))
      {
         // That didn't work either, too risky to proceed so give up
         OutStream << NA << endl;
         return;
      }
   }

   // OK, we have Windows version so display it
   OutStream << "Running under                                \t: ";

   switch (osvi.dwPlatformId)
   {
   case VER_PLATFORM_WIN32_NT:
      if (osvi.dwMajorVersion <= 4)
         OutStream << "Windows NT ";

      else if (5 == osvi.dwMajorVersion && 0 == osvi.dwMinorVersion)
         OutStream << "Windows 2000 ";

      else if (5 == osvi.dwMajorVersion && 1 == osvi.dwMinorVersion)
         OutStream << "Windows XP ";

      else if (6 == osvi.dwMajorVersion && 0 == osvi.dwMinorVersion)
         OutStream << "Windows Vista ";

      else if (6 == osvi.dwMajorVersion && 1 == osvi.dwMinorVersion)
         OutStream << "Windows 7 ";

      else if (6 == osvi.dwMajorVersion && 2 == osvi.dwMinorVersion)
         OutStream << "Windows 8 ";

      else if (6 == osvi.dwMajorVersion && 3 == osvi.dwMinorVersion)
         OutStream << "Windows 8.1 ";

      else if (10 == osvi.dwMajorVersion && 0 == osvi.dwMinorVersion)
         OutStream << "Windows 10 ";

      else if (11 == osvi.dwMajorVersion && 0 == osvi.dwMinorVersion)
         OutStream << "Windows 11 ";

      else
         OutStream << "unknown Windows version ";

      // Display version, service pack (if any), and build number
      if (osvi.dwMajorVersion <= 4)
         OutStream << "version " << osvi.dwMajorVersion << "." << osvi.dwMinorVersion << " " << osvi.szCSDVersion << " (Build " << (osvi.dwBuildNumber & 0xFFFF) << ")" << endl;

      else
         OutStream << osvi.szCSDVersion << " (Build " << (osvi.dwBuildNumber & 0xFFFF) << ")" << endl;

      break;

   case VER_PLATFORM_WIN32_WINDOWS:
      if (4 == osvi.dwMajorVersion && 0 == osvi.dwMinorVersion)
      {
         OutStream << "Windows 95";

         if ('C' == osvi.szCSDVersion[1] || 'B' == osvi.szCSDVersion[1])
            OutStream << " OSR2";

         OutStream << endl;
      }

      else if (4 == osvi.dwMajorVersion && 10 == osvi.dwMinorVersion)
      {
         OutStream << "Windows 98";

         if ('A' == osvi.szCSDVersion[1])
            OutStream << "SE";

         OutStream << endl;
      }

      else if (4 == osvi.dwMajorVersion && 90 == osvi.dwMinorVersion)
         OutStream << "Windows Me" << endl;

      else
         OutStream << "unknown 16-bit Windows version " << endl;

      break;

   case VER_PLATFORM_WIN32s:
      OutStream << "Win32s" << endl;
      break;
   }

   // Now get process timimgs: this only works under 32-bit windows
   if (VER_PLATFORM_WIN32_NT == osvi.dwPlatformId)
   {
      FILETIME ftCreate, ftExit, ftKernel, ftUser;

      if (GetProcessTimes(GetCurrentProcess(), &ftCreate, &ftExit, &ftKernel, &ftUser))
      {
         ULARGE_INTEGER ul;
         ul.LowPart = ftUser.dwLowDateTime;
         ul.HighPart = ftUser.dwHighDateTime;
         OutStream << "Time spent executing user code               \t: " << strDispTime(static_cast<double>(ul.QuadPart) * 1e-7, false) << endl;
         ul.LowPart = ftKernel.dwLowDateTime;
         ul.HighPart = ftKernel.dwHighDateTime;
         OutStream << "Time spent executing kernel code             \t: " << strDispTime(static_cast<double>(ul.QuadPart) * 1e-7, false) << endl;
      }
   }

   else
      OutStream << "Process timings                              \t: " << NA << endl;

   // Finally get more process statistics: this needs psapi.dll, so only proceed if it is present on this system
   HINSTANCE hDLL = LoadLibrary("psapi.dll");

   if (hDLL != NULL)
   {
      // The dll has been found
      typedef BOOL(__stdcall * DLLPROC)(HANDLE, PPROCESS_MEMORY_COUNTERS, DWORD);
      DLLPROC ProcAdd;

      // Try to get the address of the function we will call
      ProcAdd = (DLLPROC)GetProcAddress(hDLL, "GetProcessMemoryInfo");

      if (ProcAdd)
      {
         // Address was found
         PROCESS_MEMORY_COUNTERS pmc;

         // Now call the function
         if ((ProcAdd)(GetCurrentProcess(), &pmc, sizeof(pmc)))
         {
            OutStream << "Peak working set size                        \t: " << pmc.PeakWorkingSetSize / (1024.0 * 1024.0) << " Mb" << endl;
            OutStream << "Current working set size                     \t: " << pmc.WorkingSetSize / (1024.0 * 1024.0) << " Mb" << endl;
            OutStream << "Peak paged pool usage                        \t: " << pmc.QuotaPeakPagedPoolUsage / (1024.0 * 1024.0) << " Mb" << endl;
            OutStream << "Current paged pool usage                     \t: " << pmc.QuotaPagedPoolUsage / (1024.0 * 1024.0) << " Mb" << endl;
            OutStream << "Peak non-paged pool usage                    \t: " << pmc.QuotaPeakNonPagedPoolUsage / (1024.0 * 1024.0) << " Mb" << endl;
            OutStream << "Current non-paged pool usage                 \t: " << pmc.QuotaNonPagedPoolUsage / (1024.0 * 1024.0) << " Mb" << endl;
            OutStream << "Peak pagefile usage                          \t: " << pmc.PeakPagefileUsage / (1024.0 * 1024.0) << " Mb" << endl;
            OutStream << "Current pagefile usage                       \t: " << pmc.PagefileUsage / (1024.0 * 1024.0) << " Mb" << endl;
            OutStream << "No. of page faults                           \t: " << pmc.PageFaultCount << endl;
         }
      }

      // Free the memory used by the dll
      FreeLibrary(hDLL);
   }

#elif defined __GNUG__
   rusage ru;

   if (getrusage(RUSAGE_SELF, &ru) >= 0)
   {
      OutStream << "Time spent executing user code               \t: " << strDispTime(static_cast<double>(ru.ru_utime.tv_sec), false, true) << endl;
      // OutStream << "ru_utime.tv_usec                             \t: " << ru.ru_utime.tv_usec << endl;
      OutStream << "Time spent executing kernel code             \t: " << strDispTime(static_cast<double>(ru.ru_stime.tv_sec), false, true) << endl;
      // OutStream << "ru_stime.tv_usec                             \t: " << ru.ru_stime.tv_usec << endl;
      // OutStream << "Maximum resident set size                    \t: " << ru.ru_maxrss/1024.0 << " Mb" << endl;
      // OutStream << "ixrss (???)                                  \t: " << ru.ru_ixrss << endl;
      // OutStream << "Sum of rm_asrss (???)                        \t: " << ru.ru_idrss << endl;
      // OutStream << "isrss (???)                                  \t: " << ru.ru_isrss << endl;
      OutStream << "No. of page faults not requiring physical I/O\t: " << ru.ru_minflt << endl;
      OutStream << "No. of page faults requiring physical I/O    \t: " << ru.ru_majflt << endl;
      // OutStream << "No. of times swapped out of main memory      \t: " << ru.ru_nswap << endl;
      // OutStream << "No. of times performed input (read request)  \t: " << ru.ru_inblock << endl;
      // OutStream << "No. of times performed output (write request)\t: " << ru.ru_oublock << endl;
      // OutStream << "No. of signals received                      \t: " << ru.ru_nsignals << endl;
      OutStream << "No. of voluntary context switches            \t: " << ru.ru_nvcsw << endl;
      OutStream << "No. of involuntary context switches          \t: " << ru.ru_nivcsw << endl;
   }

   else
      OutStream << NA << endl;

#else
   OutStream << NA << endl;
#endif

   OutStream << endl;

#ifdef _OPENMP
#pragma omp parallel
   {
      if (0 == omp_get_thread_num())
      {
         OutStream << "Number of OpenMP threads                     \t: " << omp_get_num_threads() << endl;
         OutStream << "Number of OpenMP processors                  \t: " << omp_get_num_procs() << endl;

         LogStream << "Number of OpenMP threads                     \t: " << omp_get_num_threads() << endl;
         LogStream << "Number of OpenMP processors                  \t: " << omp_get_num_procs() << endl;
      }
   }
#endif

   time_t const tRunTime = m_tSysEndTime - m_tSysStartTime;
   struct tm* ptmRunTime = gmtime(&tRunTime);

   OutStream << "Time required for simulation                 \t: " << put_time(ptmRunTime, "%T") << endl;
   LogStream << "Time required for simulation                 \t: " << put_time(ptmRunTime, "%T") << endl;

   double const dSpeedUp = m_dSimDuration * 3600 / static_cast<double>(tRunTime);
   OutStream << setprecision(0);
   OutStream << "Time simulated / time required for simulation\t: " << dSpeedUp << " x faster than reality" << endl;

   LogStream << setprecision(0);
   LogStream << "Time simulated / time required for simulation\t: " << dSpeedUp << " x faster than reality" << endl;
}

//===============================================================================================================================
//! Returns an error message given an error code
//===============================================================================================================================
string CSimulation::strGetErrorText(int const nErr)
{
   string strErr;

   switch (nErr)
   {
   case RTN_USER_ABORT:
      strErr = "run ended by user";
      break;

   case RTN_ERR_BADPARAM:
      strErr = "in command-line parameter";
      break;

   case RTN_ERR_INI:
      strErr = "reading initialisation file";
      break;

   case RTN_ERR_CMEDIR:
      strErr = "in directory name";
      break;

   case RTN_ERR_RUNDATA:
      strErr = "reading run details file";
      break;

   case RTN_ERR_SCAPE_SHAPE_FUNCTION_FILE:
      strErr = "reading SCAPE shape function file";
      break;

   case RTN_ERR_TIDEDATAFILE:
      strErr = "reading tide data file";
      break;

   case RTN_ERR_LOGFILE:
      strErr = "creating log file";
      break;

   case RTN_ERR_OUTFILE:
      strErr = "creating text output file";
      break;

   case RTN_ERR_TSFILE:
      strErr = "creating time series file";
      break;

   case RTN_ERR_DEMFILE:
      strErr = "reading initial DEM file";
      break;

   case RTN_ERR_RASTER_FILE_READ:
      strErr = "reading raster GIS file";
      break;

   case RTN_ERR_VECTOR_FILE_READ:
      strErr = "reading vector GIS file";
      break;

   case RTN_ERR_MEMALLOC:
      strErr = "allocating memory";
      break;

   case RTN_ERR_RASTER_GIS_OUT_FORMAT:
      strErr = "problem with raster GIS output format";
      break;

   case RTN_ERR_VECTOR_GIS_OUT_FORMAT:
      strErr = "problem with vector GIS output format";
      break;

   case RTN_ERR_TEXT_FILE_WRITE:
      strErr = "writing text output file";
      break;

   case RTN_ERR_RASTER_FILE_WRITE:
      strErr = "writing raster GIS output file";
      break;

   case RTN_ERR_VECTOR_FILE_WRITE:
      strErr = "writing vector GIS output file";
      break;

   case RTN_ERR_TIMESERIES_FILE_WRITE:
      strErr = "writing time series output file";
      break;

   case RTN_ERR_LINETOGRID:
      strErr = "putting linear feature onto raster grid";
      break;

   case RTN_ERR_NOSEACELLS:
      strErr = "no sea cells found";
      break;

   case RTN_ERR_GRID_TO_LINE:
      strErr = "when searching grid for linear feature";
      break;

   case RTN_ERR_NO_COAST:
      strErr = "no coastlines found. Is the SWL correct?";
      break;

   case RTN_ERR_PROFILE_WRITE:
      strErr = "writing coastline-normal profiles";
      break;

   case RTN_ERR_TIME_UNITS:
      strErr = "in time units";
      break;

   case RTN_ERR_NO_SOLUTION_FOR_ENDPOINT:
      strErr = "no solution when finding end point for coastline-normal line";
      break;

   case RTN_ERR_CHECK_INTERSECTION:
      strErr = "while checking for profile intersections";
      break;

   case RTN_ERR_PROFILE_ENDPOINT_IS_INLAND:
      strErr = "end point for coastline-normal line is not in the contiguous sea";
      break;

   case RTN_ERR_CLIFF_NOTCH:
      strErr = "cliff notch is above sediment top elevation";
      break;

   case RTN_ERR_CLIFF_CANNOT_DEPOSIT_ALL:
      strErr = "unable to deposit enough unconsolidated sediment (talus) from cliff collapse";
      break;

   case RTN_ERR_PROFILE_SPACING:
      strErr = "coastline-normal profiles are too closely spaced";
      break;

   case RTN_ERR_NO_PROFILES_1:
      strErr = "no coastline-normal profiles created, check the SWL";
      break;

   case RTN_ERR_NO_PROFILES_2:
      strErr = "no coastline-normal profiles created during rasterization";
      break;

   case RTN_ERR_EDGE_OF_GRID:
      strErr = "hit grid edge when eroding beach";
      break;

   case RTN_ERR_NO_SEAWARD_END_OF_PROFILE_BEACH_EROSION:
      strErr = "could not locate seaward end of profile when creating Dean profile for beach erosion";
      break;

   case RTN_ERR_NO_SEAWARD_END_OF_PROFILE_UPCOAST_BEACH_DEPOSITION:
      strErr = "could not locate seaward end of profile when creating Dean profile for up-coast beach deposition";
      break;

   case RTN_ERR_NO_SEAWARD_END_OF_PROFILE_DOWNCOAST_BEACH_DEPOSITION:
      strErr = "could not locate seaward end of profile when creating Dean profile for down-coast beach deposition";
      break;

   case RTN_ERR_LANDFORM_TO_GRID:
      strErr = "updating grid with landforms";
      break;

   case RTN_ERR_NO_TOP_LAYER_DURING_WAVE_CALC:
      strErr = "no top layer of sediment during wave calculations";
      break;

   case RTN_ERR_NO_TOP_LAYER_DURING_BEACH_CALC:
      strErr = "no top layer of sediment during beach erosion/deposition";
      break;

   case RTN_ERR_NO_TOP_LAYER_DURING_PLATFORM_CALC:
      strErr = "no top layer of sediment during shore platform erosion";
      break;

   case RTN_ERR_NO_TOP_LAYER_DURING_CLIFF_COLLAPSE_CALC:
      strErr = "no top layer of sediment during cliff collapse";
      break;

   case RTN_ERR_NO_ADJACENT_POLYGON:
      strErr = "problem with polygon-to-polygon sediment routing sequence";
      break;

   case RTN_ERR_BAD_MULTILINE:
      strErr = "inconsistent multiline";
      break;

   case RTN_ERR_CANNOT_INSERT_POINT:
      strErr = "cannot insert point into multiline";
      break;

   case RTN_ERR_CANNOT_ASSIGN_COASTAL_LANDFORM:
      strErr = "cannot assign coastal landform";
      break;

   case RTN_ERR_SHADOW_ZONE_FLOOD_FILL_NOGRID:
      strErr = "start point for cell-by-cell fill of wave shadow zone is outside grid";
      break;

   case RTN_ERR_SHADOW_ZONE_FLOOD_START_POINT:
      strErr = "could not find start point for cell-by-cell fill of wave shadow zone";
      break;

   case RTN_ERR_CSHORE_EMPTY_PROFILE:
      strErr = "empty profile during during CShore wave propagation";
      break;

   case RTN_ERR_CSHORE_FILE_INPUT:
      strErr = "creating file for CShore input";
      break;

   case RTN_ERR_READING_CSHORE_FILE_OUTPUT:
      strErr = "reading CShore output file";
      break;

   case RTN_ERR_WAVE_INTERPOLATION_LOOKUP:
      strErr = "during wave interpolation lookup";
      break;

   case RTN_ERR_GRIDCREATE:
      strErr = "while running GDALGridCreate()";
      break;

   case RTN_ERR_COAST_CANT_FIND_EDGE_CELL:
      strErr = "cannot find edge cell while constructing grid-edge profile";
      break;

   case RTN_ERR_CSHORE_ERROR:
      strErr = "CShore did not finish correctly";
      break;

   case RTN_ERR_NO_CELL_UNDER_COASTLINE:
      strErr = "Could not find cell under coastline";
      break;

   case RTN_ERR_OPEN_DEEP_WATER_WAVE_DATA:
      strErr = "opening deep sea wave time series file";
      break;

   case RTN_ERR_READING_DEEP_WATER_WAVE_DATA:
      strErr = "reading deep sea wave time series file";
      break;

   case RTN_ERR_BOUNDING_BOX:
      strErr = "finding edges of the bounding box";
      break;

   case RTN_ERR_READING_SEDIMENT_INPUT_EVENT:
      strErr = "reading sediment input event time series file";
      break;

   case RTN_ERR_SEDIMENT_INPUT_EVENT:
      strErr = "simulating sediment input event";
      break;

   case RTN_ERR_SEDIMENT_INPUT_EVENT_LOCATION:
      strErr = "location of sediment input event is outside grod";
      break;

   case RTN_ERR_WAVESTATION_LOCATION:
      strErr = "location of wavestation is outside grid";
      break;

   case RTN_ERR_CLIFF_NOT_IN_POLYGON:
      strErr = "cliff not in polygon";
      break;

   case RTN_ERR_CELL_MARKED_PROFILE_COAST_BUT_NOT_PROFILE:
      strErr = "Cell marked as profile coast but not as profile";
      break;

   case RTN_ERR_TRACING_FLOOD:
      strErr = "tracing flood line on grid";
      break;

   case RTN_ERR_NO_START_FINISH_POINTS_TRACING_COAST:
      strErr = "tracing coastline on grid, no coast start-finish points found";
      break;

   case RTN_ERR_NO_VALID_COAST:
      strErr = "tracing coastline on grid, no valid coast found";
      break;

   case RTN_ERR_COAST_TRACING_REPEATING:
      strErr = "tracing coastline on grid, coast search just repeats";
      break;

   case RTN_ERR_COAST_TRACING_SAME_START_FINISH:
      strErr = "tracing coastline on grid, start cell is the same as the finish cell";
      break;

   case RTN_ERR_COAST_TRACING_SAME_EDGE_START_FINISH:
      strErr = "tracing coastline on grid, TEST cost starts and ends on same edge";
      break;

   case RTN_ERR_COAST_TRACING_ZERO_LENGTH:
      strErr = "tracing coastline on grid, zero length coast";
      break;

   case RTN_ERR_COAST_TRACING_TOO_SHORT:
      strErr = "tracing coastline on grid, coast less than minimum length";
      break;

   case RTN_ERR_COAST_TRACING_OFF_EDGE:
      strErr = "tracing coastline on grid, hit off-edge cell";
      break;

   case RTN_ERR_COAST_TRACING_TOO_LONG:
      strErr = "tracing coastline on grid, too many times round tracing loop";
      break;

   case RTN_ERR_CELL_NOT_FOUND_IN_HIT_PROFILE_DIFFERENT_COASTS:
      strErr = "intersection cell not found in hit profile";
      break;

   case RTN_ERR_POINT_NOT_FOUND_IN_MULTILINE_DIFFERENT_COASTS:
      strErr = "point not found when truncating multiline for different coasts";
      break;

   case RTN_ERR_CELL_NOT_FOUND_IN_HIT_PROFILE:
      strErr = "cell not found in hit profile";
      break;

   case RTN_ERR_CELL_IN_POLY_BUT_NO_POLY_COAST:
      strErr = "cell marked as in polygon, but does not have polygon's coast";
      break;

   case RTN_ERR_CLIFF_TALUS_TO_UNCONS:
      strErr = "cannot find closest point to coast when moving talus to unconsolidated sediment";
      break;

   case RTN_ERR_POLYGON_FILL_START_POINT:
      strErr = "cannot find a start point for cell-by-cell polygon infilling";
      break;

   case RTN_ERR_UNKNOWN:
      strErr = "unknown error";
      break;

   default:
      // should never get here
      strErr = " error";
   }

   return strErr;
}

//===============================================================================================================================
//! Notifies the user that the simulation has ended, asks for keypress if necessary, and if compiled under GNU can send an email
//===============================================================================================================================
void CSimulation::DoSimulationEnd(int const nRtn)
{
   // If we don't know the time that the run ended (e.g. because it did not finish correctly), then get it now
   if (m_tSysEndTime == 0)
      m_tSysEndTime = time(nullptr);

   switch (nRtn)
   {
   case (RTN_OK):
      // normal ending
      cout << RUN_END_NOTICE << put_time(localtime(&m_tSysEndTime), "%T %A %d %B %Y") << endl;
      break;

   case (RTN_HELP_ONLY):
   case (RTN_CHECK_ONLY):
      return;

   default:
      // Aborting because of some error
      cerr << RUN_END_NOTICE << "iteration " << m_ulIter << ERROR_NOTICE << nRtn << ": \"" << strGetErrorText(nRtn) << "\", " << put_time(localtime(&m_tSysEndTime), "%T %A %d %B %Y") << endl;

      if (m_ulIter > 1)
      {
         // If the run has actually started, then output all GIS files: this is very helpful in tracking down problems
         m_bSaveGISThisIter = true;
         m_nGISSave = 998; // Will get incremented to 999 when we write the files
         bSaveAllRasterGISFiles();
         bSaveAllVectorGISFiles();
      }

      // Write the error message to the logfile and to stdout
      if (LogStream && LogStream.is_open())
      {
         LogStream << ERR << strGetErrorText(nRtn) << " (error code " << nRtn << ") on " << put_time(localtime(&m_tSysEndTime), "%T %A %d %B %Y") << endl;
         LogStream.flush();
      }

      if (OutStream && OutStream.is_open())
      {
         OutStream << ERR << strGetErrorText(nRtn) << " (error code " << nRtn << ") on " << put_time(localtime(&m_tSysEndTime), "%T %A %d %B %Y") << endl;
         OutStream.flush();
      }
   }

#ifdef __GNUG__
   if (isatty(fileno(stdout)))
   {
      // Stdout is connected to a tty, so not running as a background job
      // cout << endl
      // << PRESS_KEY;
      // cout.flush();
      // getchar();
   }
   else
   {
      // Stdout is not connected to a tty, so must be running in the background; if we have something entered for the email address, then send an email
      if (! m_strMailAddress.empty())
      {
         cout << SEND_EMAIL << m_strMailAddress << endl;

         string strCmd("echo \"");

         stringstream ststrTmp;
         ststrTmp << put_time(localtime(&m_tSysEndTime), "%T on %A %d %B %Y") << endl;

         // Send an email using Linux/Unix mail command
         if (RTN_OK == nRtn)
         {
            // Finished normally
            strCmd.append("Simulation ");
            strCmd.append(m_strRunName);
            strCmd.append(", running on ");
            strCmd.append(strGetComputerName());
            strCmd.append(", completed normally at ");
            strCmd.append(ststrTmp.str());
            strCmd.append("\" | mail -s \"");
            strCmd.append(PROGRAM_NAME);
            strCmd.append(": normal completion\" ");
            strCmd.append(m_strMailAddress);
         }
         else
         {
            // Error, so give some information to help debugging
            strCmd.append("Simulation ");
            strCmd.append(m_strRunName);
            strCmd.append(", running on ");
            strCmd.append(strGetComputerName());
            strCmd.append(", aborted with error code ");
            strCmd.append(to_string(nRtn));
            strCmd.append(": ");
            strCmd.append(strGetErrorText(nRtn));
            strCmd.append(" at timestep ");
            strCmd.append(to_string(m_ulIter));
            strCmd.append(" (");
            strCmd.append(strDispSimTime(m_dSimElapsed));
            strCmd.append(").\n\nThis message sent at ");
            strCmd.append(ststrTmp.str());
            strCmd.append("\" | mail -s \"");
            strCmd.append(PROGRAM_NAME);
            strCmd.append(": ERROR\" ");
            strCmd.append(m_strMailAddress);
         }

         int const nRet = system(strCmd.c_str());

         if (WEXITSTATUS(nRet) != 0)
            cerr << ERR << EMAIL_ERROR << endl;
      }
   }
#endif
}

//===============================================================================================================================
//! Changes all forward slashes in the input string to backslashes, leaving the original unchanged
//===============================================================================================================================
string CSimulation::pstrChangeToBackslash(string const* strIn)
{
   string strOut(*strIn);
   strOut.replace(strOut.begin(), strOut.end(), '/', '\\');
   return strOut;
}

//===============================================================================================================================
//! Swaps all backslashes in the input string to forward slashes, leaving the original unchanged
//===============================================================================================================================
string CSimulation::pstrChangeToForwardSlash(string const* strIn)
{
   string strOut(*strIn);
   strOut.replace(strOut.begin(), strOut.end(), '\\', '/');
   return strOut;
}

//===============================================================================================================================
//! Trims whitespace from the left side of a string, does not change the original string
//===============================================================================================================================
string CSimulation::strTrimLeft(string const* strIn)
{
   // Trim leading spaces
   size_t const nStartpos = strIn->find_first_not_of(" \t");

   if (nStartpos == string::npos)
      return *strIn;

   else
      return strIn->substr(nStartpos);
}

//===============================================================================================================================
//! Trims whitespace from the right side of a string, does not change the original string
//===============================================================================================================================
string CSimulation::strTrimRight(string const* strIn)
{
   string strTmp(*strIn);

   // Remove any stray carriage returns (can happen if file was edited in Windows)
   strTmp.erase(remove(strTmp.begin(), strTmp.end(), '\r'), strTmp.end());

   // Trim trailing spaces
   size_t const nEndpos = strTmp.find_last_not_of(" \t");

   if (nEndpos == string::npos)
      return strTmp;

   else
      return strTmp.substr(0, nEndpos + 1);
}

//===============================================================================================================================
//! Trims whitespace from both sides of a string, does not change the original string
//===============================================================================================================================
string CSimulation::strTrim(string const* strIn)
{
   string strTmp = *strIn;

   // Remove any stray carriage returns (can happen if file was edited in Windows)
   strTmp.erase(remove(strTmp.begin(), strTmp.end(), '\r'), strTmp.end());

   // Trim trailing spaces
   size_t nPos = strTmp.find_last_not_of(" \t");

   if (nPos != string::npos)
      strTmp.resize(nPos + 1);

   // Trim leading spaces
   nPos = strTmp.find_first_not_of(" \t");

   if (nPos != string::npos)
      strTmp = strTmp.substr(nPos);

   return strTmp;
}

//===============================================================================================================================
//! Returns the lower case version of an string, leaving the original unchanged
//===============================================================================================================================
string CSimulation::strToLower(string const* strIn)
{
   string strOut = *strIn;
   transform(strIn->begin(), strIn->end(), strOut.begin(), tolower);
   return strOut;
}

//===============================================================================================================================
// Returns the upper case version of an string, leaving the original unchanged
//===============================================================================================================================
// string CSimulation::strToUpper(string const* strIn)
// {
// string strOut = *strIn;
// transform(strIn->begin(), strIn->end(), strOut.begin(), toupper);
// return strOut;
// }

//===============================================================================================================================
//! Returns a string with a substring removed, and with whitespace trimmed
//===============================================================================================================================
string CSimulation::strRemoveSubstr(string* pStrIn, string const* pStrSub)
{
   size_t const nPos = pStrIn->find(*pStrSub);

   if (nPos != string::npos)
   {
      // OK, found the substring
      pStrIn->replace(nPos, pStrSub->size(), "");
      return strTrim(pStrIn);
   }

   else
   {
      // If not found, return the string unchanged
      return *pStrIn;
   }
}

//===============================================================================================================================
//! From http://stackoverflow.com/questions/236129/split-a-string-in-c They implement (approximately) Python's split() function. This first version puts the results into a pre-constructed string vector. It ignores empty items
//===============================================================================================================================
vector<string>* CSimulation::VstrSplit(string const* s, char const delim, vector<string>* elems)
{
   stringstream ss(*s);
   string item;

   while (getline(ss, item, delim))
   {
      if (!item.empty())
         elems->push_back(item);
   }

   return elems;
}

//===============================================================================================================================
//! From http://stackoverflow.com/questions/236129/split-a-string-in-c They implement (approximately) Python's split() function. This second version returns a new string vector (it calls the first version)
//===============================================================================================================================
vector<string> CSimulation::VstrSplit(string const* s, char const delim)
{
   vector<string> elems;
   VstrSplit(s, delim, &elems);
   return elems;
}

// //===============================================================================================================================
// //! Calculates the vector cross product of three points
// //===============================================================================================================================
// double CSimulation::dCrossProduct(double const dX1, double const dY1, double const dX2, double const dY2, double const dX3, double const dY3)
// {
//    // Based on code at http://debian.fmi.uni-sofia.bg/~sergei/cgsr/docs/clockwise.htm
// return (dX2 - dX1) * (dY3 - dY2) - ((dY2 - dY1) * (dX3 - dX2));
// }

// //===============================================================================================================================
// //! Calculates the mean of a pointer to a vector of doubles
// //===============================================================================================================================
// double CSimulation::dGetMean(vector<double> const* pV)
// {
// double dSum = accumulate(pV->begin(), pV->end(), 0.0);
// double dMean = dSum / static_cast<double>(pV->size());
// return dMean;
// }

// //===============================================================================================================================
// //! Calculates the standard deviation of a pointer to a vector of doubles. From http://stackoverflow.com/questions/7616511/calculate-mean-and-standard-deviation-from-a-vector-of-samples-in-c-using-boos
// //===============================================================================================================================
// double CSimulation::dGetStdDev(vector<double> const* pV)
// {
// double dSum = accumulate(pV->begin(), pV->end(), 0.0);
// double dMean = dSum / static_cast<double>(pV->size());
//
// double dSqSum = inner_product(pV->begin(), pV->end(), pV->begin(), 0.0);
// double dStdDev = sqrt(dSqSum / static_cast<double>(pV->size()) - dMean * dMean);
//
// return dStdDev;
// }

//===============================================================================================================================
//! Appends a CGeom2DIPoint to a vector<CGeom2DIPoint>, making sure that the new end point touches the previous end point i.e. that there is no gap between the two points
//===============================================================================================================================
void CSimulation::AppendEnsureNoGap(vector<CGeom2DIPoint>* pVPtiPoints, CGeom2DIPoint const* pPti)
{
   int const nX = pPti->nGetX();
   int const nY = pPti->nGetY();
   int const nXLast = pVPtiPoints->back().nGetX();
   int const nYLast = pVPtiPoints->back().nGetY();
   int const nXDiff = nX - nXLast;
   int const nYDiff = nY - nYLast;
   int const nXDiffA = tAbs(nXDiff);
   int const nYDiffA = tAbs(nYDiff);
   int const nDiff = tMax(nXDiffA, nYDiffA);

   if (nDiff > 1)
   {
      // We have a gap
      double
          dXInc = 0,
          dYInc = 0;

      if (nXDiffA > 1)
         dXInc = static_cast<double>(nXDiff) / nDiff;

      if (nYDiffA > 1)
         dYInc = static_cast<double>(nYDiff) / nDiff;

      for (int n = 1; n < nDiff; n++)
      {
         CGeom2DIPoint const Pti(nXLast + nRound(n * dXInc), nYLast + nRound(n * dYInc));
         pVPtiPoints->push_back(Pti);
      }
   }

   pVPtiPoints->push_back(CGeom2DIPoint(nX, nY));
}

//===============================================================================================================================
//! Calculates a Dean equilibrium profile h(y) = A * y^(2/3) where h(y) is the distance below the highest point in the Dean profile at a distance y from the landward start of the profile
//===============================================================================================================================
void CSimulation::CalcDeanProfile(vector<double>* pdVDeanProfile, double const dInc, double const dDeanTopElev, double const dA, bool const bDeposition, int const nSeawardOffset, double const dStartCellElev)
{
   double dDistFromProfileStart = 0;

   if (bDeposition)
   {
      // This Dean profile is for deposition i.e. seaward displacement of the profile
      pdVDeanProfile->at(0) = dStartCellElev;   // Is coast elevation

      for (int n = 1; n < static_cast<int>(pdVDeanProfile->size()); n++)
      {
         if (n <= nSeawardOffset)
            // As we extend the profile seaward, the elevation of any points coastward of the new coast point of the Dean profile are set to the elevation of the original coast point
            pdVDeanProfile->at(n) = dStartCellElev;
         else
         {
            double const dDistBelowTop = dA * pow(dDistFromProfileStart, DEAN_POWER);
            pdVDeanProfile->at(n) = dDeanTopElev - dDistBelowTop;

            dDistFromProfileStart += dInc;
         }
      }
   }
   else
   {
      // This Dean profile is for erosion i.e. landward displacement of the profile
      for (int n = 0; n < static_cast<int>(pdVDeanProfile->size()); n++)
      {
         double const dDistBelowTop = dA * pow(dDistFromProfileStart, DEAN_POWER);
         pdVDeanProfile->at(n) = dDeanTopElev - dDistBelowTop;

         dDistFromProfileStart += dInc;
      }
   }
}

//===============================================================================================================================
//! Calculate the total elevation difference between every point in two elevation profiles (first profile - second profile)
//===============================================================================================================================
double CSimulation::dSubtractProfiles(vector<double> const* pdVFirstProfile, vector<double> const* pdVSecondProfile, vector<bool> const* pbVIsValid)
{
   double dTotElevDiff = 0;

   // Note that this assumes that all three vectors are of equal length, should really check this
   for (int n = 0; n < static_cast<int>(pdVFirstProfile->size()); n++)
   {
      if (pbVIsValid->at(n))
      {
         double const dProfileDiff = pdVFirstProfile->at(n) - pdVSecondProfile->at(n);

         dTotElevDiff += dProfileDiff;
      }
   }

   //    // DEBUG CODE -----------------------------------------------------
   // LogStream << endl;
   // LogStream << "First profile = ";
   // for (int n = 0; n < static_cast<int>(pdVFirstProfile->size()); n++)
   // {
   // LogStream << pdVFirstProfile->at(n) << " ";
   // }
   // LogStream << endl;
   // LogStream << "Second profile = ";
   // for (int n = 0; n < static_cast<int>(pdVFirstProfile->size()); n++)
   // {
   // LogStream << pdVSecondProfile->at(n) << " ";
   // }
   // LogStream << endl;
   // LogStream << "Difference = ";
   // for (int n = 0; n < static_cast<int>(pdVFirstProfile->size()); n++)
   // {
   // LogStream << pdVFirstProfile->at(n) - pdVSecondProfile->at(n) << " ";
   // }
   // LogStream << endl;
   //    // DEBUG CODE -----------------------------------------------------

   return dTotElevDiff;
}

//===============================================================================================================================
//! Calculate the depth of closure
//===============================================================================================================================
void CSimulation::CalcDepthOfClosure(void)
{
   double dDeepWaterWaveHeight;
   double dDeepWaterPeriod;

   if (m_bSingleDeepWaterWaveValues)
   {
      dDeepWaterWaveHeight = m_dAllCellsDeepWaterWaveHeight;
      dDeepWaterPeriod = m_dAllCellsDeepWaterWavePeriod;
   }

   else
   {
      dDeepWaterWaveHeight = m_dMaxUserInputWaveHeight;
      dDeepWaterPeriod = m_dMaxUserInputWavePeriod;
   }

   // TODO 051 Calculate depth of closure using 'average of the maximum values observed during a typical year'
   // dL = 2.28 * Hsx − (68.5 * Hsx^2 / (g * Tsx^2))
   // where:
   // Hsx is the nearshore storm wave height that is exceeded only 12 hours each year
   // Tsx is the associated wave period
   // from Hallermeier, R.J. (1978). Uses for a calculated limit depth to beach erosion. Proc. 16th Coastal Engineering Conf., ASCE, New York. Pp 1493 - 1512
   //
   // For the time being, and since we assume wave height and period constant just use the actual wave height and period to calculate the depth of closure
   // m_dDepthOfClosure = (2.28 * dDeepWaterWaveHeight) - (68.5 * dDeepWaterWaveHeight * dDeepWaterWaveHeight / (m_dG * dDeepWaterPeriod * dDeepWaterPeriod));

   // An alternative (which produces smaller depth of closure estimates) is Birkemeier (1985) TODO 007 Full reference needed
   // dL = 1.75 * Hsx - (57.9 * Hsx^2/ (g * Tsx^2))
   m_dDepthOfClosure = (1.75 * dDeepWaterWaveHeight) - (57.9 * dDeepWaterWaveHeight * dDeepWaterWaveHeight / (m_dG * dDeepWaterPeriod * dDeepWaterPeriod));
}

// //===============================================================================================================================
// //! Tests a reference to a string to see if it is numeric (modified from https://tfetimes.com/c-determine-if-a-string-is-numeric/)
// //===============================================================================================================================
// bool CSimulation::bIsNumeric(string const*strIn)
// {
// return all_of(strIn->begin(), strIn->end(), isdigit);
// }

//===============================================================================================================================
//! Parses a date string into days, months, and years, and checks each of them
//===============================================================================================================================
bool CSimulation::bParseDate(string const* strDate, int& nDay, int& nMonth, int& nYear)
{
   vector<string> VstrTmp = VstrSplit(strDate, SLASH);

   if (VstrTmp.size() < 3)
   {
      cerr << "date string must include day, month, and year '" << strDate << "'" << endl;
      return false;
   }

   // Sort out day
   if (! bIsStringValidInt(VstrTmp[0]))
   {
      cerr << "invalid integer for day in date '" << strDate << "'" << endl;
      return false;
   }

   nDay = stoi(VstrTmp[0]);

   if ((nDay < 1) || (nDay > 31))
   {
      cerr << "day must be between 1 and 31 in date '" << strDate << "'" << endl;
      return false;
   }

   // Sort out month
   if (! bIsStringValidInt(VstrTmp[1]))
   {
      cerr << "invalid integer for month in date '" << strDate << "'" << endl;
      return false;
   }

   nMonth = stoi(VstrTmp[1]);

   if ((nMonth < 1) || (nMonth > 12))
   {
      cerr << "month must be between 1 and 12 in date '" << strDate << "'" << endl;
      return false;
   }

   // Sort out year
   if (! bIsStringValidInt(VstrTmp[2]))
   {
      cerr << "invalid integer for year in date '" << strDate << "'" << endl;
      return false;
   }

   nYear = stoi(VstrTmp[2]);

   if (nYear < 0)
   {
      cerr << "year must be > 0 in date '" << strDate << "'" << endl;
      return false;
   }

   return true;
}

//===============================================================================================================================
//! Parses a time string into hours, minutes, and seconds, and checks each of them
//===============================================================================================================================
bool CSimulation::bParseTime(string const* strTime, int& nHour, int& nMin, int& nSec)
{
   vector<string> VstrTmp = VstrSplit(strTime, DASH);

   if (VstrTmp.size() < 3)
   {
      cerr << "time string must include hours, minutes, and seconds '" << strTime << "'" << endl;
      return false;
   }

   // Sort out hour
   if (! bIsStringValidInt(VstrTmp[0]))
   {
      cerr << "invalid integer for hours in time '" << strTime << "'" << endl;
      return false;
   }

   nHour = stoi(VstrTmp[0]);

   if ((nHour < 0) || (nHour > 23))
   {
      cerr << "hour must be between 0 and 23 in time '" << strTime << "'" << endl;
      return false;
   }

   // Sort out minutes
   if (! bIsStringValidInt(VstrTmp[1]))
   {
      cerr << "invalid integer for minutes in time '" << strTime << "'" << endl;
      return false;
   }

   nMin = stoi(VstrTmp[1]);

   if ((nMin < 0) || (nMin > 59))
   {
      cerr << "minutes must be betwen 0 and 59 in time '" << strTime << "'" << endl;
      return false;
   }

   // Sort out seconds
   if (! bIsStringValidInt(VstrTmp[2]))
   {
      cerr << "invalid integer for seconds in time '" << strTime << "'" << endl;
      return false;
   }

   nSec = stoi(VstrTmp[2]);

   if ((nSec < 0) || (nSec > 59))
   {
      cerr << "seconds must be between 0 and 59 in time '" << strTime << "'" << endl;
      return false;
   }

   return true;
}

//===============================================================================================================================
//! For sediment input events, parses a string that may be relative (a number of hours or days after the start of the simulation), or absolute (a time/date in the format hh-mm-ss dd/mm/yyyy). Returns the timestep in which the sediment input event occurs
//===============================================================================================================================
unsigned long CSimulation::ulConvertToTimestep(string const* pstrIn) const
{
   unsigned long ulTimeStep = 0;

   // Convert to lower case, remove leading and trailing whitespace
   string strDate = strToLower(pstrIn);
   strDate = strTrim(&strDate);

   if (strDate.find("hour") != string::npos)
   {
      // OK, this is a number of hours (a relative time, from the start of simulation)
      vector<string> VstrTmp = VstrSplit(&strDate, SPACE);

      if ((VstrTmp.size() < 2) || (! bIsStringValidInt(VstrTmp[0])))
      {
         cerr << "Error in number of hours '" + strDate + "' for sediment input event" << endl;
         return SEDIMENT_INPUT_EVENT_ERROR;
      }

      double const dHours = stod(strTrim(&VstrTmp[0]));

      if (dHours > m_dSimDuration)
      {
         cerr << "Sediment input event '" + strDate + "' occurs after end of simulation" << endl;
         return SEDIMENT_INPUT_EVENT_ERROR;
      }

      ulTimeStep = static_cast<unsigned long>(dRound(dHours / m_dTimeStep));
   }

   else if (strDate.find("day") != string::npos)
   {
      // OK, this is a number of days (a relative time, from the start of simulation)
      vector<string> VstrTmp = VstrSplit(&strDate, SPACE);

      if ((VstrTmp.size() < 2) || (! bIsStringValidInt(VstrTmp[0])))
      {
         cerr << "Error in number of days '" + strDate + "' for sediment input event" << endl;
         return SEDIMENT_INPUT_EVENT_ERROR;
      }

      double const dHours = stod(strTrim(&VstrTmp[0])) * 24;

      if (dHours > m_dSimDuration)
      {
         cerr << "Sediment input event '" + strDate + "' occurs after end of simulation" << endl;
         return SEDIMENT_INPUT_EVENT_ERROR;
      }

      ulTimeStep = static_cast<unsigned long>(dRound(dHours / m_dTimeStep));
   }

   else
   {
      // This is an absolute time/date in the format hh-mm-ss dd/mm/yyyy
      vector<string> VstrTmp = VstrSplit(&strDate, SPACE);

      if (VstrTmp.size() < 2)
      {
         cerr << "Error in time/date '" + strDate + "' of sediment input event" << endl;
         return SEDIMENT_INPUT_EVENT_ERROR;
      }

      int nHour = 0;
      int nMin = 0;
      int nSec = 0;

      // OK, first sort out the time
      if (! bParseTime(&VstrTmp[0], nHour, nMin, nSec))
      {
         cerr << "Error in time '" + VstrTmp[0] + "' of sediment input event" << endl;
         return SEDIMENT_INPUT_EVENT_ERROR;
      }

      int nDay = 0;
      int nMonth = 0;
      int nYear = 0;

      // Now sort out the time
      if (! bParseDate(&VstrTmp[1], nDay, nMonth, nYear))
      {
         cerr << "Error in date '" + VstrTmp[1] + "' of sediment input event" << endl;
         return SEDIMENT_INPUT_EVENT_ERROR;
      }

      // This is modified from https://stackoverflow.com/questions/14218894/number-of-days-between-two-dates-c
      struct tm tmSimStart = {};
      tmSimStart.tm_sec = m_nSimStartSec;
      tmSimStart.tm_min = m_nSimStartMin;
      tmSimStart.tm_hour = m_nSimStartHour;
      tmSimStart.tm_mday = m_nSimStartDay;
      tmSimStart.tm_mon = m_nSimStartMonth - 1;
      tmSimStart.tm_year = m_nSimStartYear - 1900;

      struct tm tmSimEvent = {};
      tmSimEvent.tm_sec = nSec;
      tmSimEvent.tm_min = nMin;
      tmSimEvent.tm_hour = nHour;
      tmSimEvent.tm_mday = nDay;
      tmSimEvent.tm_mon = nMonth - 1;
      tmSimEvent.tm_year = nYear - 1900;

      time_t const tStart = mktime(&tmSimStart);
      time_t const tEvent = mktime(&tmSimEvent);

      if (tStart == (time_t)(-1))
      {
         cerr << "Error in simulation start time/date" << endl;
         return SEDIMENT_INPUT_EVENT_ERROR;
      }

      if (tEvent == (time_t)(-1))
      {
         cerr << "Error in time/date '" + strDate + "' of sediment input event" << endl;
         return SEDIMENT_INPUT_EVENT_ERROR;
      }

      double const dHours = difftime(tEvent, tStart) / (60 * 60);

      if (dHours < 0)
      {
         cerr << "Sediment input event '" + strDate + "' occurs before start of simulation" << endl;
         return SEDIMENT_INPUT_EVENT_ERROR;
      }

      if (dHours > m_dSimDuration)
      {
         cerr << "Sediment input event '" + strDate + "' occurs after end of simulation" << endl;
         return SEDIMENT_INPUT_EVENT_ERROR;
      }

      ulTimeStep = static_cast<unsigned long>(dHours / m_dTimeStep);
   }

   return ulTimeStep;
}

//===============================================================================================================================
//! Returns true if the cell is an intervention
//===============================================================================================================================
bool CSimulation::bIsInterventionCell(int const nX, int const nY) const
{
   int const nCat = m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory();
   if ((nCat == LF_INTERVENTION_STRUCT) || (nCat == LF_INTERVENTION_NON_STRUCT))
      return true;

   return false;
}

//===============================================================================================================================
//! Do end-of-run memory clearance
//===============================================================================================================================
void CSimulation::DoEndOfRunDeletes(void)
{
   // Clear all vector coastlines, profiles, and polygons
   m_VCoast.clear();

   // m_VFloodWaveSetup.clear();
   m_VFloodWaveSetupSurge.clear();
   m_VFloodWaveSetupSurgeRunup.clear();
}

//===============================================================================================================================
//! Calculate Mean High Water (MHW) elevation for a given duration (in days). This is a tidal datum determined from the arithmetic mean of the high water heights observed each tidal day
//===============================================================================================================================
void CSimulation::CalcMHWElevation(int const nTideDataCount)
{
   // Calculate the number of tide values to read
   int const nTidevaluesPerDay = tMax(nRound(24 / m_dTimeStep), 1);
   int const nTideValuesToRead = nTidevaluesPerDay * NUM_DAYS_FOR_MEAN_HIGH_WATER_CALC;

   int const nNumTideValues = static_cast<int>(m_VdTideData.size());

   // Now read the tide data (note that this assumes that the first line of the tide data is the first tide reading of the day)
   int nThisCount = nTideDataCount;
   double dTotMaxTide = 0;
   for (int n = 0; n < nTideValuesToRead; n++)
   {
      // Read in this days's tide data
      double dDayMaxTide = -DBL_MAX;
      for (int m = 0; m < nTidevaluesPerDay; m++)
      {
         // If necessary, wrap the tide data, i.e. start again with the first line of the tide data if we do not have enough
         if (nThisCount > nNumTideValues - 1)
            nThisCount = 0;

         double const dThisTideData = m_VdTideData[nThisCount];

         if (dThisTideData > dDayMaxTide)
            dDayMaxTide = dThisTideData;

         nThisCount++;
      }

      // We have the max tide for the day, so increment the total of highest daily tides
      dTotMaxTide += dDayMaxTide;
   }

   // Now calculate the average max tide for the next NUM_DAYS_FOR_MEAN_HIGH_WATER_CALC days
   double const dMaxTideAvg = dTotMaxTide / nTideValuesToRead;

   // Finally, calculate MHW for this iteration (includes long-term SWL change)
   m_dThisIterMHWElev = m_dThisIterMeanSWL + dMaxTideAvg;

   // And set the apex elevation of any new cliff notches (i.e. cliff notches which will be created during this timestep) to be at or slightly above MHW level
   m_dThisIterNewNotchApexElev = m_dThisIterMHWElev + m_dNotchApexAboveMHW;

   // LogStream << m_ulIter << ": this-iteration MHW elevation = " << m_dThisIterMHWElev << " elevation of apex of new cliff notches = " << m_dThisIterNewNotchApexElev << endl;
}
