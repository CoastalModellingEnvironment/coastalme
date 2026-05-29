/*!
   \file simulation.cpp
   \brief The start-of-simulation routine
   \details TODO 001 A more detailed description of this routine.
   \author David Favis-Mortlock
   \author Andres Payo
   \author Wilf Chun
   \date 2026
   \copyright GNU General Public License
*/

/* ==============================================================================================================================
   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
==============================================================================================================================*/
#include <assert.h>

#include <cstdio>
#include <unistd.h>
#include <stdio.h>

#include <climits>

#include <ios>
using std::fixed;

#include <iostream>
using std::cerr;
using std::cin;
using std::endl;
using std::ios;

#include <iomanip>
using std::setprecision;

#include <cfloat>

#include <string>
using std::to_string;

#include <filesystem>      // C++17 and later, needed for missing output directory creation
using std::filesystem::is_directory;
using std::filesystem::exists;
using std::filesystem::create_directories;

#include <gdal.h>

#include "cme.h"
#include "simulation.h"
#include "raster_grid.h"
#include "coast.h"

//===============================================================================================================================
//! The CSimulation constructor
//===============================================================================================================================
CSimulation::CSimulation(void)
{
   // Initialization
   m_bHaveFineSediment = false;
   m_bHaveSandSediment = false;
   m_bHaveCoarseSediment = false;
   m_bBasementElevSave = false;
   m_bSedIncTalusTopSurfSave = false;
   m_bTopSurfIncSeaSave = false;
   m_bTalusSave = false;
   m_bSliceSave = false;
   m_bSeaDepthSave = false;
   m_bAvgSeaDepthSave = false;
   m_bWaveHeightSave = false;
   m_bAvgWaveHeightSave = false;
   m_bWaveAngleSave = false;
   m_bAvgWaveAngleSave = false;
   m_bWaveAngleAndHeightSave = false;
   m_bAvgWaveAngleAndHeightSave = false;
   m_bDeepWaterWaveAngleAndHeightSave = false;
   m_bBeachProtectionSave = false;
   m_bWaveEnergySinceCollapseSave = false;
   m_bMeanWaveEnergySave = false;
   m_bBreakingWaveHeightSave = false;
   m_bPotentialPlatformErosionSave = false;
   m_bActualPlatformErosionSave = false;
   m_bTotalPotentialPlatformErosionSave = false;
   m_bTotalActualPlatformErosionSave = false;
   m_bPotentialBeachErosionSave = false;
   m_bActualBeachErosionSave = false;
   m_bTotalPotentialBeachErosionSave = false;
   m_bTotalActualBeachErosionSave = false;
   m_bBeachDepositionSave = false;
   m_bTotalBeachDepositionSave = false;
   m_bLandformSave = false;
   m_bConsSedSlopeSave = false;
   m_bSlopeSaveForCliffToe = false;
   m_bInterventionClassSave = false;
   m_bInterventionHeightSave = false;
   m_bSuspSedSave = false;
   m_bAvgSuspSedSave = false;
   m_bFineUnconsSedSave = false;
   m_bSandUnconsSedSave = false;
   m_bCoarseUnconsSedSave = false;
   m_bFineConsSedSave = false;
   m_bSandConsSedSave = false;
   m_bCoarseConsSedSave = false;
   m_bRasterCoastlineSave = false;
   m_bRasterNormalProfileSave = false;
   m_bActiveZoneSave = false;
   m_bCliffCollapseSave = false;
   m_bTotCliffCollapseSave = false;
   m_bCliffCollapseDepositionSave = false;
   m_bTotCliffCollapseDepositionSave = false;
   m_bCliffNotchAllSave = false;
   m_bCliffCollapseTimestepSave = false;
   m_bRasterPolygonSave = false;
   m_bPotentialPlatformErosionMaskSave = false;
   m_bSeaMaskSave = false;
   m_bBeachMaskSave = false;
   m_bShadowZoneCodesSave = false;
   m_bSaveRegular = false;
   m_bCoastSave = false;
   m_bCliffEdgeSave = false;
   m_bNormalsSave = false;
   m_bInvalidNormalsSave = false;
   m_bCoastCurvatureSave = false;
   m_bPolygonNodeSave = false;
   m_bPolygonBoundarySave = false;
   m_bCliffNotchSave = false;
   m_bWaveTransectPointsSave = false;
   m_bShadowBoundarySave = false;
   m_bShadowDowndriftBoundarySave = false;
   m_bDeepWaterWaveAngleSave = false;
   m_bDeepWaterWaveHeightSave = false;
   m_bDeepWaterWavePeriodSave = false;
   m_bPolygonUnconsSedUpOrDownDriftSave = false;
   m_bPolygonUnconsSedGainOrLossSave = false;
   m_bCliffToeSave = false;
   m_bSeaAreaTSSave = false;
   m_bSWLTSSave = false;
   m_bActualPlatformErosionTSSave = false;
   m_bSuspSedTSSave = false;
   m_bFloodSetupSurgeTSSave = false;
   m_bFloodSetupSurgeRunupTSSave = false;
   m_bCliffCollapseTalusDepositionTSSave = false;
   m_bCliffCollapseErosionTSSave = false;
   m_bBeachErosionTSSave = false;
   m_bBeachDepositionTSSave = false;
   m_bBeachSedimentChangeNetTSSave = false;
   m_bCliffNotchElevTSSave = false;
   m_bSaveGISThisIter = false;
   m_bOutputConsolidatedProfileData = false;
   m_bOutputParallelProfileData = false;
   m_bOutputErosionPotentialData = false;
   m_bOmitSearchNorthEdge = false;
   m_bOmitSearchSouthEdge = false;
   m_bOmitSearchWestEdge = false;
   m_bOmitSearchEastEdge = false;
   m_bDoShorePlatformErosion = false;
   m_bDoCliffCollapse = false;
   m_bDoBeachSedimentTransport = false;
   m_bGDALCanWriteFloat = false;
   m_bGDALCanWriteInt32 = false;
   m_bScaleRasterOutput = false;
   m_bWorldFile = false;
   m_bSingleDeepWaterWaveValues = false;
   m_bHaveWaveStationData = false;
   m_bSedimentInput = false;
   m_bSedimentInputAtPoint = false;
   m_bSedimentInputAtCoast = false;
   m_bSedimentInputAlongLine = false;
   m_bSedimentInputThisIter = false;
   m_bSedimentInputEventSave = false;
   m_bWaveSetupSave = false;
   m_bStormSurgeSave = false;
   m_bRiverineFlooding = false;
   m_bRunUpSave = false;
   m_bSetupSurgeFloodMaskSave = false;
   m_bSetupSurgeRunupFloodMaskSave = false;
   m_bRasterWaveFloodLineSave = false;
   m_bVectorWaveFloodLineSave = false;
   m_bFloodLocationSave = false;
   m_bFloodSWLSetupLineSave = false;
   m_bFloodSWLSetupSurgeLine = false;
   m_bFloodSWLSetupSurgeRunupLineSave = false;
   m_bGISSaveDigitsSequential = false;
   m_bHaveConsolidatedSediment = false;
   m_bGDALOptimisations = false;
   m_bCliffToeLocate = false;
   m_bHighestSWLSoFar = false;
   m_bLowestSWLSoFar = false;

   m_bGDALCanCreate = true;
   m_bCSVPerTimestepResults = true; // Default to CSV output format
   m_bYamlInputFormat = false;      // Default to .dat format

   m_papszGDALRasterOptions = NULL;
   m_papszGDALVectorOptions = NULL;

   m_nLayers = 0;
   m_nCoastSmooth = 0;
   m_nCoastSmoothingWindowSize = 0;
   m_nSavGolCoastPoly = 0;
   m_nCliffEdgeSmooth = 0;
   m_nCliffEdgeSmoothWindow = 0;
   m_nSavGolCliffEdgePoly = 0;
   m_nProfileSmoothWindow = 0;
   m_nCoastProfileSpacing = 0;
   m_nCoastProfileInterventionSpacing = 0;
   m_nCoastCurvatureInterval = 0;
   m_nGISMaxSaveDigits = 0;
   m_nGISSave = 0;
   m_nUSave = 0;
   m_nThisSave = 0;
   m_nXGridSize = 0;
   m_nYGridSize = 0;
   m_nCoastMax = 0;
   m_nCoastMin = 0;
   m_nNumThisIterCliffCollapse = 0;
   m_nNumTotCliffCollapse = 0;
   m_nUnconsSedimentHandlingAtGridEdges = 0;
   m_nBeachErosionDepositionEquation = 0;
   m_nWavePropagationModel = 0;
   m_nSimStartSec = 0;
   m_nSimStartMin = 0;
   m_nSimStartHour = 0;
   m_nSimStartDay = 0;
   m_nSimStartMonth = 0;
   m_nSimStartYear = 0;
   m_nDeepWaterWaveDataNumTimeSteps = 0;
   m_nLogFileDetail = 0;
   m_nRunUpEquation = 0;
   // m_nLevel = 0;
   m_nExtra = 0;

   // TODO 011 May wish to make this a user-supplied value
   m_nGISMissingValue = INT_NODATA;
   m_nMissingValue = INT_NODATA;

   m_nXMinBoundingBox = INT_MAX;
   m_nXMaxBoundingBox = INT_MIN;
   m_nYMinBoundingBox = INT_MAX;
   m_nYMaxBoundingBox = INT_MIN;

   // cppcheck-suppress useInitializationList
   m_GDALWriteIntDataType = GDT_Unknown;
   // cppcheck-suppress useInitializationList
   m_GDALWriteFloatDataType = GDT_Unknown;

   m_lGDALMaxCanWrite = 0;
   m_lGDALMinCanWrite = 0;

   m_ulIter = 0;
   m_ulTotTimestep = 0;
   m_ulThisIterNumPotentialBeachErosionCells = 0;
   m_ulThisIterNumActualBeachErosionCells = 0;
   m_ulThisIterNumBeachDepositionCells = 0;
   m_ulTotPotentialPlatformErosionOnProfiles = 0;
   m_ulTotPotentialPlatformErosionBetweenProfiles = 0;
   m_ulMissingValueBasementCells = 0;
   m_ulNumCells = 0;
   m_ulThisIterNumSeaCells = 0;
   m_ulThisIterNumCoastCells = 0;
   m_ulThisIterNumPotentialPlatformErosionCells = 0;
   m_ulThisIterNumActualPlatformErosionCells = 0;

   m_ulMissingValue = UNSIGNED_LONG_NODATA;

   for (int i = 0; i < NUMBER_OF_RNGS; i++)
      m_ulRandSeed[i] = 0;

   for (int i = 0; i < SAVEMAX; i++)
      m_dUSaveTime[i] = 0;

   m_dDurationUnitsMult = 0;
   m_dNorthWestXExtCRS = 0;
   m_dNorthWestYExtCRS = 0;
   m_dSouthEastXExtCRS = 0;
   m_dSouthEastYExtCRS = 0;
   m_dExtCRSGridArea = 0;
   m_dCellSide = 0;
   m_dCellDiagonal = 0;
   m_dInvCellSide = 0;
   m_dInvCellDiagonal = 0;
   m_dCellArea = 0;
   m_dSimDuration = 0;
   m_dTimeStep = 0;
   m_dSimElapsed = 0;
   m_dRegularSaveTime = 0;
   m_dRegularSaveInterval = 0;
   m_dClkLast = 0;
   m_dCPUClock = 0;
   m_dSeaWaterDensity = 0;
   m_dThisIterSWL = 0;
   m_dThisIterMeanSWL = 0;
   m_dInitialMeanSWL = 0;
   m_dFinalMeanSWL = 0;
   m_dDeltaSWLPerTimestep = 0;
   m_dBreakingWaveHeight = 0;
   m_dC_0 = 0;
   m_dL_0 = 0;
   m_dWaveDepthRatioForWaveCalcs = 0;
   m_dAllCellsDeepWaterWaveHeight = 0;
   m_dAllCellsDeepWaterWaveAngle = 0;
   m_dAllCellsDeepWaterWavePeriod = 0;
   m_dMaxUserInputWaveHeight = 0;
   m_dMaxUserInputWavePeriod = 0;
   m_dR = 0;
   m_dD50Fine = 0;
   m_dD50Sand = 0;
   m_dD50Coarse = 0;
   m_dBeachSedimentDensity = 0;
   m_dBeachSedimentPorosity = 0;
   m_dFineErodibility = 0;
   m_dSandErodibility = 0;
   m_dCoarseErodibility = 0;
   m_dFineErodibilityNormalized = 0;
   m_dSandErodibilityNormalized = 0;
   m_dCoarseErodibilityNormalized = 0;
   m_dKLS = 0;
   m_dKamphuis = 0;
   m_dG = 0;
   m_dInmersedToBulkVolumetric = 0;
   m_dDepthOfClosure = 0;
   m_dCoastProfileSpacing = 0;
   m_dCoastProfileInterventionSpacing = 0;
   m_dCoastNormalLength = 0;
   m_dSyntheticTransectSpacing = 0;
   m_dThisIterTotSeaDepth = 0;
   m_dThisIterPotentialSedLostBeachErosion = 0;
   m_dThisIterLeftGridUnconsFine = 0;        // TODO067
   m_dThisIterLeftGridUnconsSand = 0;
   m_dThisIterLeftGridUnconsCoarse = 0;
   m_dThisIterPotentialPlatformErosion = 0;
   m_dThisIterActualPlatformErosionFineCons = 0;
   m_dThisIterActualPlatformErosionSandCons = 0;
   m_dThisIterActualPlatformErosionCoarseCons = 0;
   m_dThisIterPotentialBeachErosion = 0;
   m_dThisIterBeachErosionFine = 0;
   m_dThisIterBeachErosionSand = 0;
   m_dThisIterBeachErosionCoarse = 0;
   m_dThisIterBeachDepositionSand = 0;
   m_dThisIterBeachDepositionCoarse = 0;
   m_dThisIterFineSedimentToSuspension = 0;
   m_dDepositionSandDiff = 0;
   m_dDepositionCoarseDiff = 0;
   m_dDepthOverDBMax = 0;
   m_dTotPotentialPlatformErosionOnProfiles = 0;
   m_dTotPotentialPlatformErosionBetweenProfiles = 0;
   m_dProfileMaxSlope = 0;
   m_dMaxBeachElevAboveSWL = 0;
   m_dCliffErosionResistance = 0;
   m_dNotchIncisionAtCollapse = 0;
   m_dThisIterNewNotchApexElev = 0;
   m_dNotchApexAboveMHW = 0;
   m_dCliffDepositionA = 0;
   m_dCliffDepositionPlanviewWidth = 0;
   m_dCliffTalusMinDepositionLength = 0;
   m_dMinCliffTalusHeightFrac = 0;
   m_dThisIterCliffCollapseErosionFineUncons = 0;
   m_dThisIterCliffCollapseErosionSandUncons = 0;
   m_dThisIterCliffCollapseErosionCoarseUncons = 0;
   m_dThisIterCliffCollapseErosionFineCons = 0;
   m_dThisIterCliffCollapseErosionSandCons = 0;
   m_dThisIterCliffCollapseErosionCoarseCons = 0;
   m_dThisIterFineCliffTalusDeposition = 0;
   m_dThisIterSandCliffTalusDeposition = 0;
   m_dThisIterCoarseCliffTalusDeposition = 0;
   m_dCoastNormalRandSpacingFactor = 0;
   m_dDeanProfileStartAboveSWL = 0;
   m_dAccumulatedSeaLevelChange = 0;
   m_dBreakingWaveHeightDepthRatio = 0;
   m_dWaveDataWrapHours = 0;
   m_dThisIterTopElevMax = 0;
   m_dThisIterTopElevMin = 0;
   m_dThisiterUnconsFineInput = 0;
   m_dThisiterUnconsSandInput = 0;
   m_dThisiterUnconsCoarseInput = 0;
   m_dStartIterSuspFineAllCells = 0;
   m_dStartIterSuspFineInPolygons = 0;
   m_dStartIterUnconsFineAllCells = 0;
   m_dStartIterUnconsSandAllCells = 0;
   m_dStartIterUnconsCoarseAllCells = 0;
   m_dStartIterConsFineAllCells = 0;
   m_dStartIterConsSandAllCells = 0;
   m_dStartIterConsCoarseAllCells = 0;
   m_dThisIterDiffTotWaterLevel = 0;         // Used in surge stuff TODO 007 Finish surge and runup stuff
   m_dThisIterDiffWaveSetupWaterLevel = 0;
   m_dThisIterDiffWaveSetupSurgeWaterLevel = 0;
   m_dThisIterDiffWaveSetupSurgeRunupWaterLevel = 0;
   m_dTotalFineUnconsInPolygons = 0;
   m_dTotalSandUnconsInPolygons = 0;
   m_dTotalCoarseUnconsInPolygons = 0;
   m_dUnconsSandNotDepositedLastIter = 0;
   m_dUnconsCoarseNotDepositedLastIter = 0;
   m_dTotalFineConsInPolygons = 0;
   m_dTotalSandConsInPolygons = 0;
   m_dTotalCoarseConsInPolygons = 0;
   m_dSlopeThresholdForCliffToe = 0;
   m_dThisIterMHWElev = 0;

   m_dMinSWLSoFar = DBL_MAX;
   m_dMaxSWLSoFar = DBL_MIN;

   for (int i = 0; i < 6; i++)
      m_dGeoTransform[i] = 0;

   // TODO 011 May wish to make this a user-supplied value
   m_dGISMissingValue = DBL_NODATA;
   m_dMissingValue = DBL_NODATA;

   m_ldGTotPotentialPlatformErosion = 0;
   m_ldGTotFineActualPlatformErosion = 0;
   m_ldGTotSandActualPlatformErosion = 0;
   m_ldGTotCoarseActualPlatformErosion = 0;
   m_ldGTotPotentialSedLostBeachErosion = 0;
   m_ldGTotActualFineLostBeachErosion = 0;
   m_ldGTotActualSandLostBeachErosion = 0;
   m_ldGTotActualCoarseLostBeachErosion = 0;
   m_ldGTotSandSedLostCliffCollapse = 0;
   m_ldGTotCoarseSedLostCliffCollapse = 0;
   m_ldGTotCliffCollapseFine = 0;
   m_ldGTotCliffCollapseSand = 0;
   m_ldGTotCliffCollapseCoarse = 0;
   m_ldGTotCliffTalusFineDeposition = 0;
   m_ldGTotCliffTalusSandDeposition = 0;
   m_ldGTotCliffTalusCoarseDeposition = 0;
   m_ldGTotCliffCollapseFineErodedDuringDeposition = 0;
   m_ldGTotCliffCollapseSandErodedDuringDeposition = 0;
   m_ldGTotCliffCollapseCoarseErodedDuringDeposition = 0;
   m_ldGTotPotentialBeachErosion = 0;
   m_ldGTotActualFineBeachErosion = 0;
   m_ldGTotActualSandBeachErosion = 0;
   m_ldGTotActualCoarseBeachErosion = 0;
   m_ldGTotSandBeachDeposition = 0;
   m_ldGTotCoarseBeachDeposition = 0;
   m_ldGTotSuspendedSediment = 0;
   m_ldGTotSandDepositionDiff = 0;
   m_ldGTotCoarseDepositionDiff = 0;
   m_ldGTotFineSedimentInput = 0;
   m_ldGTotSandSedimentInput = 0;
   m_ldGTotCoarseSedimentInput = 0;

   m_tSysStartTime = 0;
   m_tSysEndTime = 0;

   m_pRasterGrid = NULL;
}

//===============================================================================================================================
//! The CSimulation destructor
//===============================================================================================================================
CSimulation::~CSimulation(void)
{
   // Close output files if open
   if (LogStream && LogStream.is_open())
   {
      LogStream.flush();
      LogStream.close();
   }

   if (OutStream && OutStream.is_open())
   {
      OutStream.flush();
      OutStream.close();
   }

   if (SeaAreaTSStream && SeaAreaTSStream.is_open())
   {
      SeaAreaTSStream.flush();
      SeaAreaTSStream.close();
   }

   if (SWLTSStream && SWLTSStream.is_open())
   {
      SWLTSStream.flush();
      SWLTSStream.close();
   }

   if (PlatformErosionTSStream && PlatformErosionTSStream.is_open())
   {
      PlatformErosionTSStream.flush();
      PlatformErosionTSStream.close();
   }

   if (CliffCollapseErosionTSStream && CliffCollapseErosionTSStream.is_open())
   {
      CliffCollapseErosionTSStream.flush();
      CliffCollapseErosionTSStream.close();
   }

   if (CliffCollapseDepositionTSStream && CliffCollapseDepositionTSStream.is_open())
   {
      CliffCollapseDepositionTSStream.flush();
      CliffCollapseDepositionTSStream.close();
   }

   if (CliffCollapseNetChangeTSStream && CliffCollapseNetChangeTSStream.is_open())
   {
      CliffCollapseNetChangeTSStream.flush();
      CliffCollapseNetChangeTSStream.close();
   }

   if (FineSedSuspensionTSStream && FineSedSuspensionTSStream.is_open())
   {
      FineSedSuspensionTSStream.flush();
      FineSedSuspensionTSStream.close();
   }

   if (FloodSetupSurgeTSStream && FloodSetupSurgeTSStream.is_open())
   {
      FloodSetupSurgeTSStream.flush();
      FloodSetupSurgeTSStream.close();
   }

   if (FloodSetupSurgeRunupTSStream && FloodSetupSurgeRunupTSStream.is_open())
   {
      FloodSetupSurgeRunupTSStream.flush();
      FloodSetupSurgeRunupTSStream.close();
   }

   if (CliffNotchElevTSStream && CliffNotchElevTSStream.is_open())
   {
      CliffNotchElevTSStream.flush();
      CliffNotchElevTSStream.close();
   }

   if (m_pRasterGrid)
      delete m_pRasterGrid;
}

//===============================================================================================================================
//! The nDoSimulation member function of CSimulation sets up and runs the simulation
//===============================================================================================================================
int CSimulation::nDoSimulation(int nArg, char const* pcArgv[])
{
   // ================================================== initialisation section ================================================
   // Hello, World!
   AnnounceStart();

   // Start the clock ticking
   StartClock();

   // Deal with command-line parameters
   int nRet = nHandleCommandLineParams(nArg, pcArgv);

   if (nRet != RTN_OK)
      return (nRet);

   // Find out the folder in which the CoastalME executable sits, in order to open the .ini file (they are assumed to be in the same folder)
   if (! bFindExeDir(pcArgv[0]))
      return (RTN_ERR_CMEDIR);

   // OK, we are off, tell the user about the licence and the start time
   AnnounceLicenceAndStartConditions();

   // Read the .ini file and get the name of the run-data file, and path for output etc.
   if (! bReadIniFile())
      return (RTN_ERR_INI);

   // Check if output dir exists
   if ((! is_directory(m_strOutPath.c_str())) || (! exists(m_strOutPath.c_str())))
   {
      // Output dir does not exist
      bool bCreateDir = false;

      if ((isatty(fileno(stdout))) && (isatty(fileno(stderr))))
      {
         // Running with stdout and stderr as a tty, so ask the user if they wish to create it
         char ch;
         cerr << endl << "Output folder '" << m_strOutPath << "' does not exist. Create it? (Y/N) ";
         cerr.flush();
         cin.get(ch);

         if ((ch == 'y') || (ch == 'Y'))
            bCreateDir = true;
      }
      else
      {
         // Running with stdout or stderr not a tty, so create output dir rather than abort
         bCreateDir = true;
      }

      if (bCreateDir)
      {
         // Yes, so create the directory
         create_directories(m_strOutPath.c_str());
         cerr << m_strOutPath << " created" << endl << endl;
      }
      else
         // Nope, just end the run
         return RTN_USER_ABORT;
   }

   // We have the name of the run-data input file, so read it
   if (! bReadRunDataFile())
      return RTN_ERR_RUNDATA;

   // Check raster GIS output format
   if (! bCheckRasterGISOutputFormat())
      return (RTN_ERR_RASTER_GIS_OUT_FORMAT);

   // Check vector GIS output format
   if (! bCheckVectorGISOutputFormat())
      return (RTN_ERR_VECTOR_GIS_OUT_FORMAT);

   // Open log file
   if (! bOpenLogFile())
      return (RTN_ERR_LOGFILE);

   // Set up the time series output files
   if (! bSetUpTSFiles())
      return (RTN_ERR_TSFILE);

   // Initialize the random number generators
   for (int n = 0; n < NUMBER_OF_RNGS; n++)
      m_Rand[n].seed(m_ulRandSeed[n]);

   // If we are doing Savitzky-Golay smoothing of the vector coastline(s), calculate the filter coefficients
   if (m_nCoastSmooth == SMOOTH_SAVITZKY_GOLAY)
      CalcSavitzkyGolayCoeffs();

   // Create the raster grid object
   m_pRasterGrid = new CGeomRasterGrid(this);

   // Read in the basement layer (must have this file), create the raster grid, then read in the basement DEM data to the array
   AnnounceReadBasementDEM();
   nRet = nReadRasterBasementDEM();

   if (nRet != RTN_OK)
      return nRet;

   // Are we simulating cliff collapse?
   if (m_bDoCliffCollapse)
   {
      // We are: now that we have a value for m_dCellSide, we can check some more input parameters. Talus must be more than one cell wide, and since the number of cells must be odd, three cells is the minimum width
      int const nTmp = nConvertMetresToNumCells(m_dCliffDepositionPlanviewWidth);
      if (nTmp < 3)
      {
         string const strErr = ERR + "cliff deposition must have a planview width of at least three cells. The current setting of " + to_string(m_dCliffDepositionPlanviewWidth) + " m gives a planview width of " + to_string(nTmp) + " cells. Please edit " + m_strDataPathName;
         cerr << strErr << endl;
         LogStream << strErr << endl;
         OutStream << strErr << endl;
         return RTN_ERR_RUNDATA;
      }
   }

   // Do some more initialisation
   // cppcheck-suppress truncLongCastAssignment
   m_ulNumCells = m_nXGridSize * m_nYGridSize;

   // Mark edge cells, as defined by the basement layer
   nRet = nMarkBoundingBoxEdgeCells();
   if (nRet != RTN_OK)
      return nRet;

   //    // DEBUG CODE =================================================================================================================
   // for (int n = 0; n < m_VPtiAllEdgeCell.size(); n++)
   // {
   // LogStream << "[" << m_VPtiAllEdgeCell[n].nGetX() << "][" << m_VPtiAllEdgeCell[n].nGetY() << "] = {" << dGridCentroidXToExtCRSX(m_VPtiAllEdgeCell[n].nGetX()) << ", " << dGridCentroidYToExtCRSY(m_VPtiAllEdgeCell[n].nGetY()) << "} " << m_VPtiAllEdgeCellEdge[n] << endl;
   // }
   //    // DEBUG CODE =================================================================================================================

   // If we are using the default cell spacing, then now that we know the size of the raster cells, we can set the size of profile spacing in m
   if (bFPIsEqual(m_dCoastProfileSpacing, 0.0, TOLERANCE))
      m_dCoastProfileSpacing = DEFAULT_PROFILE_SPACING * m_dCellSide;
   else
   {
      // The user specified a profile spacing, is this too small?
      m_nCoastProfileSpacing = nRound(m_dCoastProfileSpacing / m_dCellSide);

      if (m_nCoastProfileSpacing < DEFAULT_PROFILE_SPACING)
      {
         cerr << ERR << "profile spacing was specified as " << m_dCoastProfileSpacing << " m, which is " << m_nCoastProfileSpacing << " cells. Polygon creation works poorly if profile spacing is less than " << DEFAULT_PROFILE_SPACING << " cells, i.e. " << DEFAULT_PROFILE_SPACING * m_dCellSide << " m" << endl;

         LogStream << ERR << "profile spacing was specified as " << m_dCoastProfileSpacing << " m, which is " << m_nCoastProfileSpacing << " cells. Polygon creation works poorly if profile spacing is less than " << DEFAULT_PROFILE_SPACING << " cells, i.e. " << DEFAULT_PROFILE_SPACING * m_dCellSide << " m" << endl;

         return RTN_ERR_PROFILE_SPACING;
      }
   }

   // Set the profile spacing on interventions (both structural and non-structural)
   m_dCoastProfileInterventionSpacing = dRound(m_dCoastProfileSpacing * INTERVENTION_PROFILE_SPACING_FACTOR);
   m_nCoastProfileInterventionSpacing = nRound(m_dCoastProfileInterventionSpacing / m_dCellSide);

   // We have at least one filename for the first layer, so add the correct number of layers. Note the the number of layers does not change during the simulation: however layers can decrease in thickness until they have zero thickness
   AnnounceAddLayers();

   for (int nX = 0; nX < m_nXGridSize; nX++)
      for (int nY = 0; nY < m_nYGridSize; nY++)
         m_pRasterGrid->m_Cell[nX][nY].AppendLayers(m_nLayers);

   // Tell the user what is happening then read in the layer files
   AnnounceReadRasterFiles();

   for (int nLayer = 0; nLayer < m_nLayers; nLayer++)
   {
      if (! m_VstrInitialFineUnconsSedimentFile[nLayer].empty())
      {
         // Read in the initial fine unconsolidated sediment depth file(s)
         AnnounceReadInitialFineUnconsSedGIS(nLayer);
         nRet = nReadRasterGISFile(FINE_UNCONS_RASTER, nLayer);
         if (nRet != RTN_OK)
            return (nRet);
      }

      if (! m_VstrInitialSandUnconsSedimentFile[nLayer].empty())
      {
         // Read in the initial sand unconsolidated sediment depth file
         AnnounceReadInitialSandUnconsSedGIS(nLayer);
         nRet = nReadRasterGISFile(SAND_UNCONS_RASTER, nLayer);
         if (nRet != RTN_OK)
            return (nRet);
      }

      if (! m_VstrInitialCoarseUnconsSedimentFile[nLayer].empty())
      {
         // Read in the initial coarse unconsolidated sediment depth file
         AnnounceReadInitialCoarseUnconsSedGIS(nLayer);
         nRet = nReadRasterGISFile(COARSE_UNCONS_RASTER, nLayer);
         if (nRet != RTN_OK)
            return (nRet);
      }

      if (! m_VstrInitialFineConsSedimentFile[nLayer].empty())
      {
         // Read in the initial fine consolidated sediment depth file
         AnnounceReadInitialFineConsSedGIS(nLayer);
         nRet = nReadRasterGISFile(FINE_CONS_RASTER, nLayer);
         if (nRet != RTN_OK)
            return (nRet);
      }

      if (! m_VstrInitialSandConsSedimentFile[nLayer].empty())
      {
         // Read in the initial sand consolidated sediment depth file
         AnnounceReadInitialSandConsSedGIS(nLayer);
         nRet = nReadRasterGISFile(SAND_CONS_RASTER, nLayer);
         if (nRet != RTN_OK)
            return (nRet);
      }

      if (! m_VstrInitialCoarseConsSedimentFile[nLayer].empty())
      {
         // Read in the initial coarse consolidated sediment depth file
         AnnounceReadInitialCoarseConsSedGIS(nLayer);
         nRet = nReadRasterGISFile(COARSE_CONS_RASTER, nLayer);
         if (nRet != RTN_OK)
            return (nRet);
      }
   }

   if (! m_strInitialSuspSedimentFile.empty())
   {
      // Read in the initial suspended sediment depth file
      AnnounceReadInitialSuspSedGIS();
      nRet = nReadRasterGISFile(SUSP_SED_RASTER, 0);
      if (nRet != RTN_OK)
         return (nRet);
   }

   // Maybe read in the landform class data, otherwise calculate this during the first timestep using identification rules
   if (! m_strInitialLandformFile.empty())
   {
      AnnounceReadLGIS();
      nRet = nReadRasterGISFile(LANDFORM_RASTER, 0);
      if (nRet != RTN_OK)
         return (nRet);
   }

   // Maybe read in intervention data
   if (! m_strInterventionClassFile.empty())
   {
      AnnounceReadICGIS();
      nRet = nReadRasterGISFile(INTERVENTION_CLASS_RASTER, 0);
      if (nRet != RTN_OK)
         return (nRet);

      AnnounceReadIHGIS();
      nRet = nReadRasterGISFile(INTERVENTION_HEIGHT_RASTER, 0);
      if (nRet != RTN_OK)
         return (nRet);
   }

   // Maybe read in the tide data
   if (! m_strTideDataFile.empty())
   {
      AnnounceReadTideData();
      nRet = nReadTideDataFile();
      if (nRet != RTN_OK)
         return (nRet);
   }

   // Read in the erosion potential shape function data
   AnnounceReadSCAPEShapeFunctionFile();
   nRet = nReadShapeFunctionFile();
   if (nRet != RTN_OK)
      return (nRet);

   // Do we want to output the erosion potential look-up values, for checking purposes?
   if (m_bOutputErosionPotentialData)
      WriteLookUpData();

   // OK, now read in the vector files (if any)
   if (m_bHaveWaveStationData || m_bSedimentInput)
      AnnounceReadVectorFiles();

   // Maybe read in deep water wave station data
   if (m_bHaveWaveStationData)
   {
      // We are reading deep water wave height, orientation and period from a file of vector points and file time series
      AnnounceReadDeepWaterWaveValuesGIS();

      // Read in vector points
      nRet = nReadVectorGISFile(DEEP_WATER_WAVE_STATIONS_VEC);
      if (nRet != RTN_OK)
         return (nRet);

      int const nWaveStations = static_cast<int>(m_VnDeepWaterWaveStationID.size());

      if (nWaveStations == 1)
         m_bSingleDeepWaterWaveValues = true;

      // Read in time series values, and initialise the vector which stores each timestep's deep water wave height, orientation and period
      nRet = nReadWaveStationInputFile(nWaveStations);
      if (nRet != RTN_OK)
         return (nRet);
   }

   // Maybe read in sediment input event data
   if (m_bSedimentInput)
   {
      // We are reading sediment input event data
      AnnounceReadSedimentEventInputValuesGIS();

      // Read in vector points for sediment input events
      nRet = nReadVectorGISFile(SEDIMENT_INPUT_EVENT_LOCATION_VEC);
      if (nRet != RTN_OK)
         return (nRet);

      // Read in the time series values for sediment input events
      nRet = nReadSedimentInputEventFile();
      if (nRet != RTN_OK)
         return (nRet);
   }

   // Maybe read in flood input location
   if (m_bFloodLocationSave)
   {
      // We are reading sediment input event data
      AnnounceReadFloodLocationGIS();

      // Read in vector points for sediment input events
      nRet = nReadVectorGISFile(FLOOD_LOCATION_VEC);
      if (nRet != RTN_OK)
         return (nRet);
   }

   // Open the main output
   OutStream.open(m_strOutFile.c_str(), ios::out | ios::trunc);

   if (!OutStream)
   {
      // Error, cannot open Out file
      cerr << ERR << "cannot open " << m_strOutFile << " for output" << endl;
      return (RTN_ERR_OUTFILE);
   }

   // Write beginning-of-run information to Out and Log files
   WriteStartRunDetails();

   // Final stage of initialization
   AnnounceFinalInitialization();

   // Misc initialisation calcs
   m_nCoastMax = COAST_LENGTH_MAX_CONST * tMax(m_nXGridSize, m_nYGridSize);                     // Arbitrary but probably OK
   m_nCoastMin = COAST_LENGTH_MIN_CONST * nRound(m_dCoastProfileSpacing);                        // Arbitrary but probably OK
   m_nCoastCurvatureInterval = tMax(nRound(m_dCoastProfileSpacing / (m_dCellSide * 2)), 2);      // Arbitrary but probably OK

   // For beach erosion/deposition, conversion from immersed weight to bulk volumetric (sand and voids) transport rate (Leo Van Rijn) TODO 007 need full reference
   m_dInmersedToBulkVolumetric = 1 / ((m_dBeachSedimentDensity - m_dSeaWaterDensity) * (1 - m_dBeachSedimentPorosity) * m_dG);

   m_bConsSedChangedThisIter.resize(m_nLayers, false);
   m_bUnconsChangedThisIter.resize(m_nLayers, false);

   // Normalize sediment erodibility values, so that none are > 1
   double const dTmp = m_dFineErodibility + m_dSandErodibility + m_dCoarseErodibility;
   m_dFineErodibilityNormalized = m_dFineErodibility / dTmp;
   m_dSandErodibilityNormalized = m_dSandErodibility / dTmp;
   m_dCoarseErodibilityNormalized = m_dCoarseErodibility / dTmp;

   // Intialise SWL
   m_dThisIterSWL = m_dInitialMeanSWL;

   // If SWL changes during the simulation, calculate the per-timestep increment (could be -ve)
   if (! bFPIsEqual(m_dFinalMeanSWL, m_dInitialMeanSWL, TOLERANCE))
   {
      m_dDeltaSWLPerTimestep = (m_dTimeStep * (m_dFinalMeanSWL - m_dInitialMeanSWL)) / m_dSimDuration;
      m_dAccumulatedSeaLevelChange -= m_dDeltaSWLPerTimestep;
   }

   // Calculate the sea-to-land processing direction for each side of the grid
   CalcGridEdgeSeaToLandDirection();

   // ===================================================== The main loop ======================================================
   // Tell the user what is happening
   AnnounceIsRunning();

   while (true)
   {
      // Check that we haven't gone on too long: if not then update timestep number etc.
      if (bTimeToQuit())
         break;

      // Tell the user how the simulation is progressing
      AnnounceProgress();

      if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
         LogStream << "TIMESTEP " << m_ulIter << " " << string(154, '=') << endl;

      LogStream << fixed << setprecision(3);

      // Note: m_prSlumpDirtyCells is NOT cleared here - it accumulates across timesteps. It will be cleared after GIS output is written (when saving at intervals)

      // Check to see if there is a new intervention in place: if so, update it on the RasterGrid array
      nRet = nUpdateIntervention();
      if (nRet != RTN_OK)
         return nRet;

      // Calculate changes due to external forcing (change in still water level, tide level and deep water waves height, orientation and period)
      nRet = nCalcExternalForcing();
      if (nRet != RTN_OK)
         return nRet;

      // Do per-timestep initialisation: set up the grid cells ready for this timestep, also initialise per-timestep totals. Note that in the first timestep, all cells (including hinterland cells) are given the deep water wave values
      nRet = nInitGridAndCalcStillWaterLevel();
      if (nRet != RTN_OK)
         return nRet;

      // Next find out which cells are inundated and locate the coastline(s). This also gives to all sea cells, wave values which are the same as the deep water values. For shallow water sea cells, these wave values will be changed later, in nDoAllPropagateWaves()
      nRet = nLocateSeaAndCoasts();
      if (nRet != RTN_OK)
         return nRet;

      // Tell the user how the simulation is progressing
      AnnounceProgress();

      // Locate estuaries TODO someday...

      // if (m_bHaveConsolidatedSediment && m_bDoCliffCollapse && m_bCliffToeLocate)
      // {
      //    // Locate and trace cliff toe
      //    nRet = nLocateCliffToe();
      //    if (nRet != RTN_OK)
      //       return nRet;
      // }

      // For all cells, use classification rules to assign sea and hinterland landform categories
      nRet = nAssignLandformsForAllCells();
      if (nRet != RTN_OK)
         return nRet;

      // For every coastline, use classification rules to assign landform categories
      nRet = nAssignLandformsForAllCoasts();
      if (nRet != RTN_OK)
         return nRet;

      // For every coastline, create all coastline-normal profiles
      nRet = nCreateAllProfiles();
      if (nRet != RTN_OK)
         return nRet;

      // For every coastline, check all coastline-normal profiles for intersection, then modify intersecting profiles so that the sections of each profile seaward of the point of intersection are 'shared' i.e. are multi-lines. This creates the boundaries of the triangular polygons
      CheckAllProfilesForIntersection();

      // For every coastline, do further checks on profiles
      nRet = nFurtherCheckAndMarkAllProfiles();
      if (nRet != RTN_OK)
         return nRet;

      if (m_VCoast.size() > 1)
      {
         // We have multiple coastlines
         nRet = nDoMultipleCoastlines();
         if (nRet != RTN_OK)
         return nRet;
      }

      // Mark valid profiles on the raster grid
      nRet = nMarkProfilesOnGrid();
      if (nRet != RTN_OK)
         return nRet;

      // Tell the user how the simulation is progressing
      AnnounceProgress();

      // Create the coast polygons
      nRet = nCreateAllPolygons();
      if (nRet != RTN_OK)
         return nRet;

      // Mark cells of the raster grid that are within each polygon, and do some polygon initialisation
      nRet = nMarkPolygonCells();
      if (nRet != RTN_OK)
         return nRet;

      // // DEBUG CODE ================
      // m_nGISSave++;
      // if (! bWriteVectorGISFile(VECTOR_PLOT_COAST, &VECTOR_PLOT_COAST_TITLE))
      //    return false;
      // if (! bWriteVectorGISFile(VECTOR_PLOT_NORMALS, &VECTOR_PLOT_NORMALS_TITLE))
      //    return false;
      // if (! bWriteVectorGISFile(VECTOR_PLOT_INVALID_NORMALS, &VECTOR_PLOT_INVALID_NORMALS_TITLE))
      //    return false;
      // if (! bWriteRasterGISFile(RASTER_PLOT_NORMAL_PROFILE, &RASTER_PLOT_NORMAL_PROFILE_TITLE))
      //    return false;
      // if (! bWriteRasterGISFile(RASTER_PLOT_COAST, &RASTER_PLOT_COAST_TITLE))
      //    return false;
      // if (! bWriteRasterGISFile(RASTER_PLOT_POLYGON, &RASTER_PLOT_POLYGON_TITLE))
      //    return false;
      // if (! bWriteVectorGISFile(VECTOR_PLOT_POLYGON_BOUNDARY, &VECTOR_PLOT_POLYGON_BOUNDARY_TITLE))
      //    return false;
      // // DEBUG CODE ================

      // Calculate the length of the shared normal between each polygon and the adjacent polygon(s)
      nRet = nDoPolygonSharedBoundaries();
      if (nRet != RTN_OK)
         return nRet;

      // Tell the user how the simulation is progressing
      AnnounceProgress();

      // // DEBUG CODE =========================================================================================================
      // nNODATA = 0;
      // nPoly0 = 0;
      // nPoly24 = 0;
      // for (int nX = 0; nX < m_nXGridSize; nX++)
      // {
      // for (int nY = 0; nY < m_nYGridSize; nY++)
      // {
      // int nTmp = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonID();
      // if (nTmp == INT_NODATA)
      // nNODATA++;
      //
      // if (nTmp == 0)
      // nPoly0++;
      //
      // if (nTmp == 24)
      // nPoly24++;
      // }
      // }
      // LogStream << "After marking polygon cells, N cells with NODATA polygon ID = " << nNODATA << endl;
      // // LogStream << "After marking polygon cells, N cells with zero polygon ID = " << nPoly0 << endl;
      // LogStream << "After marking polygon cells, N cells with 24 polygon ID = " << nPoly24 << endl;
      // // DEBUG CODE =========================================================================================================
      // PropagateWind();

      // Give every coast point a value for deep water wave height and direction
      nRet = nSetAllCoastpointDeepWaterWaveValues();
      if (nRet != RTN_OK)
         return nRet;

      // // DEBUG CODE ===============
      // for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
      // {
      // LogStream << "====================" << endl;
      //
      // for (int nProfile = 0; nProfile < m_VCoast[nCoast].nGetNumProfiles(); nProfile++)
      // {
      // CGeomProfile* pProfile = m_VCoast[nCoast].pGetProfile(nProfile);
      // int nCell = pProfile->nGetNumCellsInProfile();
      // LogStream << "Profile " << pProfile->nGetProfileID() << " nGetNumCellsInProfile() = " << nCell << endl;
      // }
      //
      // LogStream << endl;
      //
      // for (int nProfile = 0; nProfile < m_VCoast[nCoast].nGetNumProfiles(); nProfile++)
      // {
      // CGeomProfile* pProfile = m_VCoast[nCoast].pGetProfileWithDownCoastSeq(nProfile);
      // int nCell = pProfile->nGetNumCellsInProfile();
      // LogStream << "Profile " << pProfile->nGetProfileID() << " nGetNumCellsInProfile() = " << nCell << endl;
      // }
      //
      // LogStream << "====================" << endl;
      // }
      // // DEBUG CODE =====================

      // Change the wave properties in all shallow water sea cells: propagate waves and define the active zone, also locate wave shadow zones
      nRet = nDoAllPropagateWaves();
      if (nRet != RTN_OK)
         return nRet;

      // Output polygon share table and pre-existing sediment table to log file
      if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
      {
         WritePolygonInfoTable();
         WritePolygonPreExistingSedimentTable();
      }

      // Tell the user how the simulation is progressing
      AnnounceProgress();

      //       // DEBUG CODE ===========================================================================================================
      // string strOutFile = m_strOutPath;
      // strOutFile += "sea_wave_height_CHECKPOINT_";
      // strOutFile += to_string(m_ulIter);
      // strOutFile += ".tif";
      //
      // GDALDriver* pDriver = GetGDALDriverManager()->GetDriverByName("gtiff");
      // GDALDataset* pDataSet = pDriver->Create(strOutFile.c_str(), m_nXGridSize, m_nYGridSize, 1, GDT_Float64, m_papszGDALRasterOptions);
      // pDataSet->SetProjection(m_strGDALBasementDEMProjection.c_str());
      // pDataSet->SetGeoTransform(m_dGeoTransform);
      //
      // int nn = 0;
      // double* pdRaster = new double[m_nXGridSize * m_nYGridSize];
      // for (int nY = 0; nY < m_nYGridSize; nY++)
      // {
      // for (int nX = 0; nX < m_nXGridSize; nX++)
      // {
      // pdRaster[nn++] = m_pRasterGrid->m_Cell[nX][nY].dGetWaveHeight();
      // }
      // }
      //
      // GDALRasterBand* pBand = pDataSet->GetRasterBand(1);
      // pBand->SetNoDataValue(m_dMissingValue);
      // int nRet = pBand->RasterIO(GF_Write, 0, 0, m_nXGridSize, m_nYGridSize, pdRaster, m_nXGridSize, m_nYGridSize, GDT_Float64, 0, 0, NULL);
      //
      // if (nRet == CE_Failure)
      // return RTN_ERR_GRIDCREATE;
      //
      // GDALClose(pDataSet);
      // delete[] pdRaster;
      //       // DEBUG CODE ===========================================================================================================
      //
      //       // DEBUG CODE ===========================================================================================================
      // strOutFile = m_strOutPath;
      // strOutFile += "sea_wave_angle_CHECKPOINT_";
      // strOutFile += to_string(m_ulIter);
      // strOutFile += ".tif";
      //
      // pDriver = GetGDALDriverManager()->GetDriverByName("gtiff");
      // pDataSet = pDriver->Create(strOutFile.c_str(), m_nXGridSize, m_nYGridSize, 1, GDT_Float64, m_papszGDALRasterOptions);
      // pDataSet->SetProjection(m_strGDALBasementDEMProjection.c_str());
      // pDataSet->SetGeoTransform(m_dGeoTransform);
      //
      // nn = 0;
      // pdRaster = new double[m_nXGridSize * m_nYGridSize];
      // for (int nY = 0; nY < m_nYGridSize; nY++)
      // {
      // for (int nX = 0; nX < m_nXGridSize; nX++)
      // {
      // pdRaster[nn++] = m_pRasterGrid->m_Cell[nX][nY].dGetWaveAngle();
      // }
      // }
      //
      // pBand = pDataSet->GetRasterBand(1);
      // pBand->SetNoDataValue(m_dMissingValue);
      // nRet = pBand->RasterIO(GF_Write, 0, 0, m_nXGridSize, m_nYGridSize, pdRaster, m_nXGridSize, m_nYGridSize, GDT_Float64, 0, 0, NULL);
      //
      // if (nRet == CE_Failure)
      // return RTN_ERR_GRIDCREATE;
      //
      // GDALClose(pDataSet);
      // delete[] pdRaster;
      //       // DEBUG CODE ===========================================================================================================

      // Save the not-deposited values, to be shown in the logfile after we've finished beach sediment movement
      m_dUnconsSandNotDepositedLastIter = m_dDepositionSandDiff;
      m_dUnconsCoarseNotDepositedLastIter = m_dDepositionCoarseDiff;

      if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
         if (m_dDepositionSandDiff > MASS_BALANCE_TOLERANCE)
         {
            LogStream << m_ulIter << ": AT ITERATION START m_dDepositionSandDiff = " << m_dDepositionSandDiff * m_dCellArea << " m_dUnconsSandNotDepositedLastIter = " << m_dUnconsSandNotDepositedLastIter << endl;
            LogStream << m_ulIter << ": AT ITERATION START m_dDepositionCoarseDiff = " << m_dDepositionCoarseDiff * m_dCellArea << " m_dUnconsCoarseNotDepositedLastIter = " << m_dUnconsCoarseNotDepositedLastIter << endl;
         }

      if (m_bDoShorePlatformErosion)
      {
         // Calculate elevation change on the consolidated sediment which comprises the coastal platform
         nRet = nDoAllShorePlatFormErosion();
         if (nRet != RTN_OK)
            return nRet;
      }

      // Output shore platform erosion table to log file
      if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
         WritePolygonShorePlatformErosion();

      // Are we considering cliff collapse?
      if (m_bHaveConsolidatedSediment && m_bDoCliffCollapse)
      {
         // Distribute wave energy to coast landforms (currently, this is only relevant to cliffs), maybe incise cliff notches, maybe do cliff collapses
         nRet = nDoAllWaveEnergyToCoastLandforms();
         if (nRet != RTN_OK)
            return nRet;

         // Output cliff collapse table to log file
         if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
            WritePolygonCliffCollapseErosion();

         // Move some cliff collapse talus to unconsolidated sediment
         nRet = nMoveCliffTalusToUnconsolidatedOrSuspension();
         if (nRet != RTN_OK)
            return nRet;
      }

      // Tell the user how the simulation is progressing
      AnnounceProgress();

      if (m_bDoBeachSedimentTransport)
      {
         // Next simulate beach erosion and deposition i.e. simulate alongshore transport of unconsolidated sediment (longshore drift) between polygons. First calculate potential sediment movement between polygons
         DoAllPotentialBeachErosion();

         // Do within-polygon redistribution of unconsolidated sediment, constraining potential sediment movement to give actual (i.e. supply-limited) sediment movement to/from each polygon in three size classes
         nRet = nDoAllActualBeachErosionAndDeposition();
         if (nRet != RTN_OK)
            return nRet;
      }

      // If we have sediment input events, then check to see whether this is time for an event to occur. If it is, then do it
      if (m_bSedimentInput)
      {
         nRet = nCheckForSedimentInputEvent();
         if (nRet != RTN_OK)
            return nRet;

         // If we have had at least one sediment input event this iteration, then output the sediment event per polygon table to the log file
         if (m_bSedimentInputThisIter && (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL))
            WritePolygonSedimentInputEventTable();
      }

      // Do slumping of unconsolidated sediment on dune areas, for cells that changed this timestep TODO make this an option
      nRet = nDoSedimentSlumping();
      if (nRet != RTN_OK)
         return nRet;

      // Tell the user how the simulation is progressing
      AnnounceProgress();

      //       // DEBUG CODE ===========================================================================================================
      // string strOutFile = m_strOutPath;
      // strOutFile += "sea_wave_height_CHECKPOINT_";
      // strOutFile += to_string(m_ulIter);
      // strOutFile += ".tif";
      //
      // GDALDriver* pDriver = GetGDALDriverManager()->GetDriverByName("gtiff");
      // GDALDataset* pDataSet = pDriver->Create(strOutFile.c_str(), m_nXGridSize, m_nYGridSize, 1, GDT_Float64, m_papszGDALRasterOptions);
      // pDataSet->SetProjection(m_strGDALBasementDEMProjection.c_str());
      // pDataSet->SetGeoTransform(m_dGeoTransform);
      //
      // int nn = 0;
      // double* pdRaster = new double[m_nXGridSize * m_nYGridSize];
      // for (int nY = 0; nY < m_nYGridSize; nY++)
      // {
      // for (int nX = 0; nX < m_nXGridSize; nX++)
      // {
      // pdRaster[nn++] = m_pRasterGrid->m_Cell[nX][nY].dGetWaveHeight();
      // }
      // }
      //
      // GDALRasterBand* pBand = pDataSet->GetRasterBand(1);
      // pBand->SetNoDataValue(m_dMissingValue);
      // int nRet = pBand->RasterIO(GF_Write, 0, 0, m_nXGridSize, m_nYGridSize, pdRaster, m_nXGridSize, m_nYGridSize, GDT_Float64, 0, 0, NULL);
      //
      // if (nRet == CE_Failure)
      // return RTN_ERR_GRIDCREATE;
      //
      // GDALClose(pDataSet);
      // delete[] pdRaster;
      //       // DEBUG CODE ===========================================================================================================
      //
      //       // DEBUG CODE ===========================================================================================================
      // strOutFile = m_strOutPath;
      // strOutFile += "sea_wave_angle_CHECKPOINT_";
      // strOutFile += to_string(m_ulIter);
      // strOutFile += ".tif";
      //
      // pDriver = GetGDALDriverManager()->GetDriverByName("gtiff");
      // pDataSet = pDriver->Create(strOutFile.c_str(), m_nXGridSize, m_nYGridSize, 1, GDT_Float64, m_papszGDALRasterOptions);
      // pDataSet->SetProjection(m_strGDALBasementDEMProjection.c_str());
      // pDataSet->SetGeoTransform(m_dGeoTransform);
      //
      // nn = 0;
      // pdRaster = new double[m_nXGridSize * m_nYGridSize];
      // for (int nY = 0; nY < m_nYGridSize; nY++)
      // {
      // for (int nX = 0; nX < m_nXGridSize; nX++)
      // {
      // pdRaster[nn++] = m_pRasterGrid->m_Cell[nX][nY].dGetWaveAngle();
      // }
      // }
      //
      // pBand = pDataSet->GetRasterBand(1);
      // pBand->SetNoDataValue(m_dMissingValue);
      // nRet = pBand->RasterIO(GF_Write, 0, 0, m_nXGridSize, m_nYGridSize, pdRaster, m_nXGridSize, m_nYGridSize, GDT_Float64, 0, 0, NULL);
      //
      // if (nRet == CE_Failure)
      // return RTN_ERR_GRIDCREATE;
      //
      // GDALClose(pDataSet);
      // delete[] pdRaster;
      //       // DEBUG CODE ===========================================================================================================

      // Do some end-of-timestep updates to the raster grid, also update per-timestep and running totals
      nRet = nEndOfTimestepUpdateGrid();
      if (nRet != RTN_OK)
         return nRet;

      // // Make water level inundation on grid
      // if (m_bFloodSWLSetupSurgeLine || m_bSetupSurgeFloodMaskSave)
      // {
      //    m_nLevel = 0;
      //
      //    nRet = nLocateFloodAndCoasts();
      //    if (nRet != RTN_OK)
      //       return nRet;
      // }
      //
      // if (m_bFloodSWLSetupSurgeRunupLineSave || m_bSetupSurgeRunupFloodMaskSave)
      // {
      //    // TODO 007 Finish surge and runup stuff
      //    m_nLevel = 1;
      //
      //    nRet = nLocateFloodAndCoasts();
      //    if (nRet != RTN_OK)
      //       return nRet;
      // }

      // Now save results, first the raster and vector GIS files if required
      m_bSaveGISThisIter = false;

      if ((m_bSaveRegular && (m_dSimElapsed >= m_dRegularSaveTime) && (m_dSimElapsed < m_dSimDuration)) || (! m_bSaveRegular && (m_dSimElapsed >= m_dUSaveTime[m_nThisSave])))
      {
         m_bSaveGISThisIter = true;

         // Save the values from the RasterGrid array into raster GIS files
         if (! bSaveAllRasterGISFiles())
            return (RTN_ERR_RASTER_FILE_WRITE);

         // Tell the user how the simulation is progressing
         AnnounceProgress();

         // Save the vector GIS files
         if (! bSaveAllVectorGISFiles())
            return (RTN_ERR_VECTOR_FILE_WRITE);

         // Clear dirty cells now that GIS output has been written. This allows dirty cells to accumulate across timesteps when saving at intervals
         m_prSlumpDirtyCells.clear();

         // Tell the user how the simulation is progressing
         AnnounceProgress();
      }

      // Output per-timestep results to the .out file
      if (! bWritePerTimestepResults())
         return (RTN_ERR_TEXT_FILE_WRITE);

      // Now output time series CSV stuff
      if (! bWriteTSFiles())
         return (RTN_ERR_TIMESERIES_FILE_WRITE);

      // Tell the user how the simulation is progressing
      AnnounceProgress();

      // Update grand totals
      DoEndOfTimestepTotals();

   } // ================================================ End of main loop ======================================================

   // =================================================== post-loop tidying =====================================================
   // Tell the user what is happening
   AnnounceSimEnd();

   // Write end-of-run information to Out, Log and time-series files
   nRet = nWriteEndRunDetails();
   if (nRet != RTN_OK)
      return (nRet);

   // Do end-of-run memory clearance
   DoEndOfRunDeletes();
   return RTN_OK;
}
