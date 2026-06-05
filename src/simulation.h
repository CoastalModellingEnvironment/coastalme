/*!
   \class CSimulation
   \brief This class runs CoastalME simulations
   \details TODO 001 This is a more detailed description of the CSimulation class
   \author David Favis-Mortlock
   \author Andres Payo
   \author Wilf Chun
   \date 2026
   \copyright GNU General Public License
   \file simulation.h
   \brief Contains CSimulation definitions
*/

#ifndef SIMULATION_H
#define SIMULATION_H
/* ===============================================================================================================================
   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
===============================================================================================================================*/
#include <vector>
using std::vector;

#include <ctime>
using std::localtime;
using std::time;
using std::time_t;

#include <fstream>
using std::ofstream;

#include <string>
using std::string;

#include <set>
using std::set;

#include <utility>
using std::pair;

#include <stack>
using std::stack;

#include <random>
using std::default_random_engine;
using std::normal_distribution;

#include <gdal.h>
using ::GDALDataType;

#include "yaml_parser.h"
#include "configuration.h"
#include "2d_point.h"
#include "2di_point.h"
#include "line.h"
#include "cme.h"
#include "line.h"
#include "configuration.h"

class CGeomRasterGrid; // Forward declarations
class CRWCoast;
class CGeomProfile;
class CGeomCoastPolygon;
class CRWCliff;
class CRWSedInputEvent;
class CRWCellLandform;
class CConfiguration;

class CSimulation
{
 private:
   //! Does this simulation consider fine-sized sediment?
   bool m_bHaveFineSediment;

   //! Does this simulation consider sand-sized sediment?
   bool m_bHaveSandSediment;

   //! Does this simulation consider coarse-sized sediment?
   bool m_bHaveCoarseSediment;

   //! Save basement raster DEMs?
   bool m_bBasementElevSave;

   //! Save sediment (inc talus) top surface raster DEMs?
   bool m_bSedIncTalusTopSurfSave;

   //! Save top surface (sediment, talus, and sea) raster DEMs?
   bool m_bTopSurfIncSeaSave;

   //! Save talus depth?
   bool m_bTalusSave;

   //! Save slope raster grids (used for cliff toe location)?
   bool m_bSlopeSaveForCliffToe;

   //! Save cliff toe raster grids?
   bool m_bCliffToeSave;

   //! Save slices?
   bool m_bSliceSave;

   //! Save sea depth raster GIS files?
   bool m_bSeaDepthSave;

   //! Save average sea depth raster GIS files?
   bool m_bAvgSeaDepthSave;

   //! Save wave height raster GIS files?
   bool m_bWaveHeightSave;

   //! Save wave height raster GIS files?
   bool m_bAvgWaveHeightSave;

   //! Save wave angle raster GIS files?
   bool m_bWaveAngleSave;

   //! Save average wave angle raster GIS files?
   bool m_bAvgWaveAngleSave;

   //! Save wave angle and wave height raster GIS files?
   bool m_bWaveAngleAndHeightSave;

   //! Save average wave angle and average wave height raster GIS files?
   bool m_bAvgWaveAngleAndHeightSave;

   //! Save deep water wave angle and wave height raster GIS files?
   bool m_bDeepWaterWaveAngleAndHeightSave;

   //! Save wave energy since cliff collapse raster GIS files?
   bool m_bWaveEnergySinceCollapseSave;

   //! Save mean wave energy raster GIS files?
   bool m_bMeanWaveEnergySave;

   //! Save breaking wave height raster GIS files?
   bool m_bBreakingWaveHeightSave;

   //! Save beach protection raster GIS files>
   bool m_bBeachProtectionSave;

   //! Save potential shore platform erosion raster GIS files?
   bool m_bPotentialPlatformErosionSave;

   //! Save actual (supply-limited) shore platform erosion raster GIS files?
   bool m_bActualPlatformErosionSave;

   //! Save total potential shore platform erosion raster GIS files?
   bool m_bTotalPotentialPlatformErosionSave;

   //! Save total actual (supply-limited) shore platform erosion raster GIS files?
   bool m_bTotalActualPlatformErosionSave;

   //! Save potential beach (unconsolidated sediment) erosion raster GIS files?
   bool m_bPotentialBeachErosionSave;

   //! Save actual (supply-limited) beach (unconsolidated sediment) erosion raster GIS files?
   bool m_bActualBeachErosionSave;

   //! Save total potential beach (unconsolidated sediment) erosion raster GIS files?
   bool m_bTotalPotentialBeachErosionSave;

   //! Save total actual (supply-limited) beach (unconsolidated sediment) erosion raster GIS files?
   bool m_bTotalActualBeachErosionSave;

   //! Save beach (unconsolidated sediment) deposition raster GIS files?
   bool m_bBeachDepositionSave;

   //! Save total beach (unconsolidated sediment) deposition raster GIS files?
   bool m_bTotalBeachDepositionSave;

   //! Save coast landform raster GIS files?
   bool m_bLandformSave;

   //! Save slope of consolidated sediment raster GIS files?
   bool m_bConsSedSlopeSave;

   //! Save intervention class raster GIS files?
   bool m_bInterventionClassSave;

   //! Save intervention height raster GIS files?
   bool m_bInterventionHeightSave;

   //! Save suspended sediment raster GIS files?
   bool m_bSuspSedSave;

   //! Save average suspended sediment raster GIS files?
   bool m_bAvgSuspSedSave;

   //! Save fine unconsolidated sediment raster GIS files?
   bool m_bFineUnconsSedSave;

   //! Save sand unconsolidated sediment raster GIS files?
   bool m_bSandUnconsSedSave;

   //! Save coarse unconsolidated sediment raster GIS files?
   bool m_bCoarseUnconsSedSave;

   //! Save fine consolidated sediment raster GIS files?
   bool m_bFineConsSedSave;

   //! Save sand consolidated sediment raster GIS files?
   bool m_bSandConsSedSave;

   //! Save coarse consolidated sediment raster GIS files?
   bool m_bCoarseConsSedSave;

   //! Save rasterized coastline GIS files?
   bool m_bRasterCoastlineSave;

   //! Save rasterized coastline-normal profiles GIS files?
   bool m_bRasterNormalProfileSave;

   //! Save active zone raster GIS files?
   bool m_bActiveZoneSave;

   //! Save cliff collapse raster GIS files?
   bool m_bCliffCollapseSave;

   //! Save total cliff collapse raster GIS files?
   bool m_bTotCliffCollapseSave;

   //! Save cliff collapse deposition raster GIS files?
   bool m_bCliffCollapseDepositionSave;

   //! Save total cliff collapse deposition raster GIS files?
   bool m_bTotCliffCollapseDepositionSave;

   //! Save raster polygon raster GIS files?
   bool m_bRasterPolygonSave;

   //! Save potential platform erosion mask raster GIS files?
   bool m_bPotentialPlatformErosionMaskSave;

   //! Save sea mask raster GIS files?
   bool m_bSeaMaskSave;

   //! Save beach mask raster GIS files?
   bool m_bBeachMaskSave;

   //! Save wave shadow zones raster GIS files?
   bool m_bShadowZoneCodesSave;

   //! Save deep water wave angle raster GIS files?
   bool m_bDeepWaterWaveAngleSave;

   //! Save deep water wave height raster GIS files?
   bool m_bDeepWaterWaveHeightSave;

   //! Save deep water wave period raster GIS files?
   bool m_bDeepWaterWavePeriodSave;

   //! Save polygon unconsolidated sediment up- or down-drift raster GIS files?
   bool m_bPolygonUnconsSedUpOrDownDriftSave;

   //! Save polygon unconsolidated sediment gain or loss raster GIS files?
   bool m_bPolygonUnconsSedGainOrLossSave;

   //! Save GIS files at regular intervals?
   bool m_bSaveRegular;

   //! Save coastline as vector GIS file?
   bool m_bCoastSave;

   //! Save cliff edge vector GIS files?
   bool m_bCliffEdgeSave;

   //! Save coastline-normal vector GIS files?
   bool m_bNormalsSave;

   //! Save invalid coastline-normal vector GIS files?
   bool m_bInvalidNormalsSave;

   //! Save coastline-curvature vector GIS files?
   bool m_bCoastCurvatureSave;

   //! Save polygon node vector GIS files?
   bool m_bPolygonNodeSave;

   //! Save polygon boundary vector GIS files?
   bool m_bPolygonBoundarySave;

   //! Save cliff notch incision depth vector GIS files?
   bool m_bCliffNotchSave;

   //! Save wave transect points vector GIS files?
   bool m_bWaveTransectPointsSave;

   //! Save wave shadow boundary vector GIS files?
   bool m_bShadowBoundarySave;

   //! Save wave shadow downdrift boundary vector GIS files?
   bool m_bShadowDowndriftBoundarySave;

   //! Save the sea area time series file?
   bool m_bSeaAreaTSSave;

   //! Save the SWL (still water level) time series file?
   bool m_bSWLTSSave;

   //! Save the actual (supply-limited) shore platform erosion time series file?
   bool m_bActualPlatformErosionTSSave;

   //! Save the cliff collapse erosion time series file?
   bool m_bCliffCollapseErosionTSSave;

   //! Save the cliff collapse deposition time series file?
   bool m_bCliffCollapseTalusDepositionTSSave;

   //! Save the beach (unconsolidated sediment) erosion time series file?
   bool m_bBeachErosionTSSave;

   //! Save the beach (unconsolidated sediment) deposition time series file?
   bool m_bBeachDepositionTSSave;

   //! Save the beach (unconsolidated sediment) net change time series file?
   bool m_bBeachSedimentChangeNetTSSave;

   //! Save the suspended sediment time series file?
   bool m_bSuspSedTSSave;

   //! Save the flood setup surge time series file? TODO 007 Finish surge and runup stuff
   bool m_bFloodSetupSurgeTSSave;

   //! Save the flood setup surge runup time series file? TODO 007 Finish surge and runup stuff
   bool m_bFloodSetupSurgeRunupTSSave;

   //! Save the cliff notch elevation time series file?
   bool m_bCliffNotchElevTSSave;

   //! Save GIS files this iteration?
   bool m_bSaveGISThisIter;

   //! Output profile data?
   bool m_bOutputConsolidatedProfileData;

   //! Output parallel profile data?
   bool m_bOutputParallelProfileData;

   //! Output erosion potential data?
   bool m_bOutputErosionPotentialData;

   //! Output per-timestep results in CSV format instead of fixed-width?
   bool m_bCSVPerTimestepResults;

   //! Use YAML format for input datafile instead of .dat format?
   bool m_bYamlInputFormat;

   //! Omit the north edge of the grid from coast-end searches?
   bool m_bOmitSearchNorthEdge;

   //! Omit the south edge of the grid from coast-end searches?
   bool m_bOmitSearchSouthEdge;

   //! Omit the west edge of the grid from coast-end searches?
   bool m_bOmitSearchWestEdge;

   //! Omit the east edge of the grid from coast-end searches?
   bool m_bOmitSearchEastEdge;

   //! Simulate shore platform erosion?
   bool m_bDoShorePlatformErosion;

   //! Simulate cliff collapse?
   bool m_bDoCliffCollapse;

   //! Simulate unconsolidated sediment (beach) transport?
   bool m_bDoBeachSedimentTransport;

   //! Is the selected GDAL output file format capable of writing files?
   bool m_bGDALCanCreate;

   //! Is the selected GDAL output file format capable of writing floating-point values to files?
   bool m_bGDALCanWriteFloat;

   //! Is the selected GDAL output file format capable of writing 32-bit integers to files?
   bool m_bGDALCanWriteInt32;

   //! Scale raster output?
   bool m_bScaleRasterOutput;

   //! Write a GIS World file?
   bool m_bWorldFile;

   //! Do we have just a point source for (i.e. only a single measurement of) deep water wave values
   bool m_bSingleDeepWaterWaveValues;

   //! Do we have wave station data?
   bool m_bHaveWaveStationData;

   //! Do we have sediment input events?
   bool m_bSedimentInput;

   //! Do we have sediment input at a point?
   bool m_bSedimentInputAtPoint;

   //! Do we have sediment input at the coast?
   bool m_bSedimentInputAtCoast;

   //! Do we have sediment input along a line?
   bool m_bSedimentInputAlongLine;

   //! Save sediment input data?
   bool m_bSedimentInputEventSave;

   //! Do we have a sediment input event this iteration?
   bool m_bSedimentInputThisIter;

   //! Are we doing riverine flooding?
   bool m_bRiverineFlooding;

   //! Are we saving the wave setup? TODO 007 Finish surge and runup stuff
   bool m_bWaveSetupSave;

   //! Are we saving the storm surge? TODO 007 Finish surge and runup stuff
   bool m_bStormSurgeSave;

   //! Are we saving runup? TODO 007 Finish surge and runup stuff
   bool m_bRunUpSave;

   //! Are we saving the setup surge flood mask? TODO 007 Finish surge and runup stuff
   bool m_bSetupSurgeFloodMaskSave;

   //! Are we saving the setup surge runup flood mask? TODO 007 Finish surge and runup stuff
   bool m_bSetupSurgeRunupFloodMaskSave;

   //! Are we saving the raster wave flood line? TODO 007 Finish surge and runup stuff
   bool m_bRasterWaveFloodLineSave;

   //! Are we saving the vector wave flood line? TODO 007 Finish surge and runup stuff
   bool m_bVectorWaveFloodLineSave;

   //! Are we saving the flood location? TODO 007 Finish surge and runup stuff
   bool m_bFloodLocationSave;

   //! Are we saving the flood still water level setup line? TODO 007 Finish surge and runup stuff
   bool m_bFloodSWLSetupLineSave;

   //! Are we saving the flood still water level setup surge line? TODO 007 Finish surge and runup stuff
   bool m_bFloodSWLSetupSurgeLine;

   //! Are we saving all cliff notches?
   bool m_bCliffNotchAllSave;

   //! Are we saving the timestep at which each cliff occurred?
   bool m_bCliffCollapseTimestepSave;

   //! Are we saving the flood still water level setup surge runup line? TODO 007 Finish surge and runup stuff
   bool m_bFloodSWLSetupSurgeRunupLineSave;

   //! Are the GIS save digits (which are part of each GIS file name) sequential, or are they the iteration number?
   bool m_bGISSaveDigitsSequential;

   //! Does this simulation consider consolidated sediment, or is it an unconsolidated sediment only simulation?>
   bool m_bHaveConsolidatedSediment;

   //! GDAL optimisations enabled?
   bool m_bGDALOptimisations;

   //! Cliff to location?
   bool m_bCliffToeLocate;

   //! Do we have the highest SWL so far?
   bool m_bHighestSWLSoFar;

   //! Do we have the lowest SWL so far?
   bool m_bLowestSWLSoFar;

   //! In order to go from low elevation (sea) to high elevation (land), process the north edge cells forward (i.e. ascending indices)?
   bool m_bSearchNorthEdgeForward = false;

   //! In order to go from low elevation (sea) to high elevation (land), process the south edge cells forward (i.e. ascending indices)?
   bool m_bSearchSouthEdgeForward = false;

   //! In order to go from low elevation (sea) to high elevation (land), process the east edge cells forward (i.e. ascending indices)?
   bool m_bSearchEastEdgeForward = false;

   //! In order to go from low elevation (sea) to high elevation (land), process the west edge cells forward (i.e. ascending indices)?
   bool m_bSearchWestEdgeForward = false;

   //! Are we simulating slumping?
   bool m_bSlumping = false;

   //! Options for GDAL when handling raster files
   char **m_papszGDALRasterOptions;

   //! Options for GDAL when handling vector files
   char **m_papszGDALVectorOptions;

   //! The size of the grid in the x direction
   int m_nXGridSize;

   //! The size of the grid in the y direction
   int m_nYGridSize;

   //! The number of sediment layers
   int m_nLayers;

   //! Which method to use for coast smoothing
   int m_nCoastSmooth;

   //! The size of the window used for coast smoothing. Must be an odd number
   int m_nCoastSmoothingWindowSize;

   //! The order of the coastline profile smoothing polynomial if Savitzky-Golay smoothing is used (usually 2 or 4, max is 6)
   int m_nSavGolCoastPoly;

   //! Which method to use for cliff edge smoothing
   int m_nCliffEdgeSmooth;

   //! The size of the window used for cliff edge smoothing. Must be an odd number
   int m_nCliffEdgeSmoothWindow;

   //! The order of the cliff edge smoothing polynomial if Savitzky-Golay smoothing is used (usually 2 or 4, max is 6)
   int m_nSavGolCliffEdgePoly;

   //! The size of the window used for running-mean coast-normal profile smoothing (must be odd)
   int m_nProfileSmoothWindow;

   //! Average spacing between coastline normals, measured in cells
   int m_nCoastProfileSpacing;

   //! Average spacing between coastline normals on interventions, measured in cells
   int m_nCoastProfileInterventionSpacing;

   //! Coast curvature interval is a length, measured in coastline points
   int m_nCoastCurvatureInterval;

   //! The maximum number of digits in GIS filenames. These can be sequential, or the iteration number
   int m_nGISMaxSaveDigits;

   //! The save number for GIS files (can be sequential, or the iteration number)
   int m_nGISSave;

   //! If user-defined GIS save intervals, the number of these
   int m_nUSave;

   //! Used in calculations of GIS save intervals
   int m_nThisSave;

   //! Maximum valid coast length when searching for coasts
   int m_nCoastMax;

   //! Minimum valid coast length when searching for coasts
   int m_nCoastMin;

   //! The number of cells with cliff collapse this iteration
   int m_nNumThisIterCliffCollapse;

   //! The total number of cells with cliff collapse since the start of the simulation
   int m_nNumTotCliffCollapse;

   //! How sediment which moves off an edge of the grid is handled. Possible values are GRID_EDGE_CLOSED, GRID_EDGE_OPEN, GRID_EDGE_RECIRCULATE
   int m_nUnconsSedimentHandlingAtGridEdges;

   //! Which beach erosion-deposition equation is used. Possible values are UNCONS_SEDIMENT_EQUATION_CERC and UNCONS_SEDIMENT_EQUATION_KAMPHUIS
   int m_nBeachErosionDepositionEquation;

   //! The value for integer missing values, as read from GIS input files
   int m_nGISMissingValue;

   //! Used by CoastalME for integer missing values
   int m_nMissingValue;

   //! The minimum x value of the bounding box
   int m_nXMinBoundingBox;

   //! The maximum x value of the bounding box
   int m_nXMaxBoundingBox;

   //! The minimum y value of the bounding box
   int m_nYMinBoundingBox;

   //! The maximum y value of the bounding box
   int m_nYMaxBoundingBox;

   //! The wave propagation model used. Possible values are WAVE_MODEL_CSHORE and WAVE_MODEL_COVE
   int m_nWavePropagationModel;

   //! Start time of the simulation (seconds)
   int m_nSimStartSec;

   //! Start time of the simulation (minutes)
   int m_nSimStartMin;

   //! Start time of the simulation (hours)
   int m_nSimStartHour;

   //! Start date of the simulation (day)
   int m_nSimStartDay;

   //! Start date of the simulation (month)
   int m_nSimStartMonth;

   //! Start date of the simulation (year)
   int m_nSimStartYear;

   //! The duration of data for deep water waves, expressed as a number of timesteps
   int m_nDeepWaterWaveDataNumTimeSteps;

   //! The level of detail in the log file output. Can be LOG_FILE_LOW_DETAIL, LOG_FILE_MIDDLE_DETAIL, LOG_FILE_HIGH_DETAIL, or LOG_FILE_ALL
   int m_nLogFileDetail;

   //! The run-up equation used TODO 007 Finish surge and runup stuff
   int m_nRunUpEquation;

   // //! TODO 007 Used in WAVESETUP + SURGE + RUNUP Finish surge and runup stuff
   // int m_nLevel;

   //! Used in debugging, if we need to output multiple sets of vector files within the same timestep
   int m_nExtra;

   //! The data type used by GDAL for integer operations, can be GDT_Byte, GDT_Int16, GDT_UInt16, GDT_Int32, or GDT_UInt32
   GDALDataType m_GDALWriteIntDataType;

   //! Thw data type used by GDAL for floating point operations, can be GDT_Byte, GDT_Int16, GDT_UInt16, GDT_Int32, GDT_UInt32, or GDT_Float32
   GDALDataType m_GDALWriteFloatDataType;

   //! The maximum integer value which GDAL can write, can be UINT8_MAX, INT16_MAX, UINT16_MAX, INT32_MAX, or UINT32_MAX,
   long m_lGDALMaxCanWrite;

   //! The minimum integer value which GDAL can write, can be zero, INT16_MIN, INT32_MIN
   long m_lGDALMinCanWrite;

   //! The number of the current iteration (time step)
   unsigned long m_ulIter;

   //! The target number of iterations
   unsigned long m_ulTotTimestep;

   //! A seed for each of the random number generators
   unsigned long m_ulRandSeed[NUMBER_OF_RNGS];

   //! The number of cells in the grid
   unsigned long m_ulNumCells;

   //! The number of grid cells which are marked as sea, for this iteration
   unsigned long m_ulThisIterNumSeaCells;

   //! The number of grid cells which are marked as coast, for this iteration
   unsigned long m_ulThisIterNumCoastCells;

   //! The number of grid cells on which potential platform erosion occurs, for this iteration
   unsigned long m_ulThisIterNumPotentialPlatformErosionCells;

   //! The number of grid cells on which actual platform erosion occurs, for this iteration
   unsigned long m_ulThisIterNumActualPlatformErosionCells;

   //! The number of grid cells on which potential beach (unconsolidated sediment) erosion occurs, for this iteration
   unsigned long m_ulThisIterNumPotentialBeachErosionCells;

   //! The number of grid cells on which actual beach (unconsolidated sediment) erosion occurs, for this iteration
   unsigned long m_ulThisIterNumActualBeachErosionCells;

   //! The number of grid cells on which beach (unconsolidated sediment) deposition occurs, for this iteration
   unsigned long m_ulThisIterNumBeachDepositionCells;

   //! The number of cells on which on-profile average potential shore platform erosion occurs
   unsigned long m_ulTotPotentialPlatformErosionOnProfiles;

   //! The number of cells on which between-profile average potential shore platform erosion occurs
   unsigned long m_ulTotPotentialPlatformErosionBetweenProfiles;

   //! The number of basement cells marked with as missing value
   unsigned long m_ulMissingValueBasementCells;

   //! Used by CoastalME for unsigned long integer missing values
   unsigned long m_ulMissingValue;

   //! Multiplier for duration units, to convert to hours
   double m_dDurationUnitsMult;

   //! The north-west x coordinate, in the external coordinate reference system (CRS)
   double m_dNorthWestXExtCRS;

   //! The north-west y coordinate, in the external coordinate reference system (CRS)
   double m_dNorthWestYExtCRS;

   //! The south-east x coordinate, in the external coordinate reference system (CRS)
   double m_dSouthEastXExtCRS;

   //! The south-east y coordinate, in the external coordinate reference system (CRS)
   double m_dSouthEastYExtCRS;

   //! The area of the grid (in external CRS units)
   double m_dExtCRSGridArea;

   //! Length of a cell side (in external CRS units)
   double m_dCellSide;

   //! Area of a cell (in external CRS units)
   double m_dCellArea;

   //! Length of a cell's diagonal (in external CRS units)
   double m_dCellDiagonal;

   //! Inverse of m_dCellSide
   double m_dInvCellSide;

   //! Inverse of m_dCellDiagonal
   double m_dInvCellDiagonal;

   //! Duration of simulation, in hours
   double m_dSimDuration;

   //! The length of an iteration (a timestep) in hours
   double m_dTimeStep;

   //! Time simulated so far, in hours
   double m_dSimElapsed;

   //! The time of the next save, in hours from the start of the simulation, if we are saving regularly
   double m_dRegularSaveTime;

   //! The interval between regular saves, in hours
   double m_dRegularSaveInterval;

   //! Save time, in hours from the start of the simukation, if we are not saving regularly
   double m_dUSaveTime[SAVEMAX];

   //! Last value returned by clock()
   double m_dClkLast;

   //! Total elapsed CPU time
   double m_dCPUClock;

   //! GDAL geotransformation info (see http://www.gdal.org/classGDALDataset.html)
   double m_dGeoTransform[6];

   //! Density of sea water in kg/m**3
   double m_dSeaWaterDensity;

   //! The start-of-simulation still water level (m)
   double m_dInitialMeanSWL;

   //! The end-of-simulation still water (m), is same as m_dInitialMeanSWL unless SWL changes
   double m_dFinalMeanSWL;

   //! If long-term SWL changes, the increment per timestep
   double m_dDeltaSWLPerTimestep;

   //! The still water level for this timestep (this includes tidal changes and any long-term SWL change)
   double m_dThisIterSWL;

   //! The mean still water level for this timestep (does not include tidal changes, but includes any long-term SWL change)
   double m_dThisIterMeanSWL;

   //! If long-term SWL changes, the total change (m) so far since the start of simulation
   double m_dAccumulatedSeaLevelChange;

   //! Minimum still water level (m)
   double m_dMinSWLSoFar;

   //! Maximum still water level (m)
   double m_dMaxSWLSoFar;

   //! TODO 007 Finish surge and runup stuff
   double m_dThisIterDiffTotWaterLevel;

   //! TODO 007 Finish surge and runup stuff
   double m_dThisIterDiffWaveSetupWaterLevel;

   //! TODO 007 Finish surge and runup stuff
   double m_dThisIterDiffWaveSetupSurgeWaterLevel;

   //! TODO 007 Finish surge and runup stuff
   double m_dThisIterDiffWaveSetupSurgeRunupWaterLevel;

   //! The height of breaking waves (m)
   double m_dBreakingWaveHeight;

   //! Deep water wave speed (m/s)
   double m_dC_0;

   //! Deep water wave length (m)
   double m_dL_0;

   //! Start depth for wave calculations
   double m_dWaveDepthRatioForWaveCalcs;

   //! Breaking wave height-to-depth ratio
   double m_dBreakingWaveHeightDepthRatio;

   //! Deep water wave height (m) for all sea cells
   double m_dAllCellsDeepWaterWaveHeight;

   //! Deep water wave angle for all sea cells
   double m_dAllCellsDeepWaterWaveAngle;

   //! Deep water wave period for all sea cells
   double m_dAllCellsDeepWaterWavePeriod;

   //! Maximum deep water wave height
   double m_dMaxUserInputWaveHeight;

   //! Used to constrain depth of closure
   double m_dMaxUserInputWavePeriod;

   //! Coast platform resistance to erosion R, see Walkden & Hall, 2011
   double m_dR;

   //! The D50 for fine sediment
   double m_dD50Fine;

   //! The D50 for sand sediment
   double m_dD50Sand;

   //! The D50 for coarse sediment
   double m_dD50Coarse;

   //! The density of unconsolidated beach sediment (kg/m**3)
   double m_dBeachSedimentDensity;

   //! The porosity of unconsolidated beach sediment (0-1)
   double m_dBeachSedimentPorosity;

   //! The relative erodibility (0-1) of fine unconsolidated beach sediment
   double m_dFineErodibility;

   //! The relative erodibility (0-1) of sand unconsolidated beach sediment
   double m_dSandErodibility;

   //! The relative erodibility (0-1) of coarse unconsolidated beach sediment
   double m_dCoarseErodibility;

   //! Relative erodibility of fine unconsolidated beach sediment, normalized
   double m_dFineErodibilityNormalized;

   //! Relative erodibility of sand unconsolidated beach sediment, normalized
   double m_dSandErodibilityNormalized;

   //! Relative erodibility of coarse unconsolidated beach sediment, normalized
   double m_dCoarseErodibilityNormalized;

   //! Transport parameter KLS in the CERC equation
   double m_dKLS;

   //! Transport parameter for the Kamphuis equation
   double m_dKamphuis;

   //! Gravitational acceleration (m**2/sec)
   double m_dG;

   //! For beach erosion/deposition, conversion from immersed weight to bulk volumetric (sand and voids) transport rate (Leo Van Rijn) TODO 007 Need date of reference
   double m_dInmersedToBulkVolumetric;

   //! Depth of closure (in m) TODO 007 can be calculated using Hallermeier, R.J. (1978) or Birkemeier (1985) TODO 045 This needs to be a user decision
   double m_dDepthOfClosure;

   //! Average spacing of the coastline-normal profiles, in m
   double m_dCoastProfileSpacing;

   //! Average spacing of the coastline-normal profiles on interventions, in m
   double m_dCoastProfileInterventionSpacing;

   //! Length of the coastline-normal profiles, in m
   double m_dCoastNormalLength;

   //! Approximate minimum spacing (m) between wave transects (real and synthetic) for wave interpolation densification
   double m_dSyntheticTransectSpacing;

   //! Total sea depth (m) for this iteration
   double m_dThisIterTotSeaDepth;

   //! Total potential platform erosion (all size classes of consolidated sediment) for this iteration (depth in m)
   double m_dThisIterPotentialPlatformErosion;

   //! Total actual platform erosion (fine consolidated sediment) for this iteration (depth in m)
   double m_dThisIterActualPlatformErosionFineCons;

   //! Total actual platform erosion (sand consolidated sediment) for this iteration (depth in m)
   double m_dThisIterActualPlatformErosionSandCons;

   //! Total actual platform erosion (coarse consolidated sediment) for this iteration (depth in m)
   double m_dThisIterActualPlatformErosionCoarseCons;

   //! Total potential beach erosion (all size classes of unconsolidated sediment) for this iteration (depth in m)
   double m_dThisIterPotentialBeachErosion;

   //! Total actual beach erosion (fine unconsolidated sediment) for this iteration (depth in m)
   double m_dThisIterBeachErosionFine;

   //! Total actual beach erosion (sand unconsolidated sediment) for this iteration (depth in m)
   double m_dThisIterBeachErosionSand;

   //! Total actual  beach erosion (coarse unconsolidated sediment) for this iteration (depth in m)
   double m_dThisIterBeachErosionCoarse;

   //! Total beach deposition (sand unconsolidated sediment) for this iteration (depth in m)
   double m_dThisIterBeachDepositionSand;

   //! Total beach deposition (coarse unconsolidated sediment) for this iteration (depth in m)
   double m_dThisIterBeachDepositionCoarse;

   //! Total fine unconsolidated sediment in suspension for this iteration (depth in m)
   double m_dThisIterFineSedimentToSuspension;

   //! Total unconsolidated sediment from beach erosion (all size classes) lost from the grid this iteration (depth in m)
   double m_dThisIterPotentialSedLostBeachErosion;

   //! Total fine unconsolidated sediment lost from the grid this iteration (depth in m)
   double m_dThisIterLeftGridUnconsFine;

   //! Total sand unconsolidated sediment lost from the grid this iteration (depth in m)
   double m_dThisIterLeftGridUnconsSand;

   //! Total coarse unconsolidated sediment lost from the grid this iteration (depth in m)
   double m_dThisIterLeftGridUnconsCoarse;

   //! Total depth of fine talus for this iteration
   double m_dThisIterTotFineTalusDepth;

   //! Total depth of sand talus for this iteration
   double m_dThisIterTotSandTalusDepth;

   //! Total depth of coarse talus for this iteration
   double m_dThisIterTotCoarseTalusDepth;

   //! Error term: if we are unable to deposit enough unconslidated sand on polygon(s), this is held over to be deposited the next iteration
   double m_dDepositionSandDiff;

   //! Error term: if we are unable to deposit enough unconslidated coarse on polygon(s), this is held over to be deposited the next iteration
   double m_dDepositionCoarseDiff;

   //! Maximum value of deoth over DB, is used in erosion potential look-up function
   double m_dDepthOverDBMax;

   //! Total potential platform erosion on profiles
   double m_dTotPotentialPlatformErosionOnProfiles;

   //! Total potential platform erosion between profiles
   double m_dTotPotentialPlatformErosionBetweenProfiles;

   //! Maximum slope on coastline-normal profiles
   double m_dProfileMaxSlope;

   //! Maximum elevation of beach above SWL (m)
   double m_dMaxBeachElevAboveSWL;

   //! Resistance of cliff to notch erosion
   double m_dCliffErosionResistance;

   //! Notch overhang (i.e. length of horizontal incision at the apex elevation) to initiate collapse (m)
   double m_dNotchIncisionAtCollapse;

   //! Elevation (m) of the apex of any cliff notches created during this iteration
   double m_dThisIterNewNotchApexElev;

   //! Distance of notch base below SWL (m)
   double m_dNotchApexAboveMHW;

   //! Scale parameter A for cliff deposition (m^(1/3)), may be zero for auto-calculation
   double m_dCliffDepositionA;

   //! Planview width of cliff collapse talus (m)
   double m_dCliffDepositionPlanviewWidth;

   //! Planview length of cliff deposition talus (m)
   double m_dCliffTalusMinDepositionLength;

   //! Minimum height of the landward end of cliff collapse talus, as a fraction of cliff elevation
   double m_dMinCliffTalusHeightFrac;

   //! This-iteration total of fine unconsolidated sediment produced by cliff collapse (m^3)
   double m_dThisIterCliffCollapseErosionFineUncons;

   //! This-iteration total of sand unconsolidated sediment produced by cliff collapse (m^3)
   double m_dThisIterCliffCollapseErosionSandUncons;

   //! This-iteration total of coarse unconsolidated sediment produced by cliff collapse (m^3)
   double m_dThisIterCliffCollapseErosionCoarseUncons;

   //! This-iteration total of fine consolidated sediment produced by cliff collapse (m^3)
   double m_dThisIterCliffCollapseErosionFineCons;

   //! This-iteration total of sand consolidated sediment produced by cliff collapse (m^3)
   double m_dThisIterCliffCollapseErosionSandCons;

   //! This-iteration total of coarse consolidated sediment produced by cliff collapse (m^3)
   double m_dThisIterCliffCollapseErosionCoarseCons;

   //! This-iteration total of fine unconsolidated sediment deposited due to cliff collapse (m^3)
   double m_dThisIterFineCliffTalusDeposition;

   //! This-iteration total of sand unconsolidated sediment deposited due to cliff collapse (m^3)
   double m_dThisIterSandCliffTalusDeposition;

   //! This-iteration total of coarse unconsolidated sediment deposited due to cliff collapse (m^3)
   double m_dThisIterCoarseCliffTalusDeposition;

   //! Random factor for spacing of along-coast normals
   double m_dCoastNormalRandSpacingFactor;

   //! Berm height i.e. height above SWL of start of depositional Dean profile
   double m_dDeanProfileStartAboveSWL;

   //! The value for floating-point missing values, as read from GIS input files
   double m_dGISMissingValue;

   //! Used by CoastalME for floating-point missing values
   double m_dMissingValue;

   //! Number of hours after which deep water wave data wraps
   double m_dWaveDataWrapHours;

   //! This-iteration highest elevation of DEM
   double m_dThisIterTopElevMax;

   //! This-iteration lowest elevation of DEM
   double m_dThisIterTopElevMin;

   //! Depth (m) of fine unconsolidated sediment added, at this iteration
   double m_dThisiterUnconsFineInput;

   //! Depth (m) of sand unconsolidated sediment added, at this iteration
   double m_dThisiterUnconsSandInput;

   //! Depth (m) of coarse unconsolidated sediment added, at this iteration
   double m_dThisiterUnconsCoarseInput;

   //! Depth (m) of fine suspended sediment at the start of the simulation, all cells (both inside and outside polygons)
   double m_dStartIterSuspFineAllCells;

   //! Depth (m) of fine suspended sediment at the start of the simulation (only cells in polygons)
   double m_dStartIterSuspFineInPolygons;

   //! Depth (m) of fine unconsolidated sediment at the start of the simulation, all cells (both inside and outside polygons)
   double m_dStartIterUnconsFineAllCells;

   //! Depth (m) of sand unconsolidated sediment at the start of the simulation, all cells (both inside and outside polygons)
   double m_dStartIterUnconsSandAllCells;

   //! Depth (m) of coarse unconsolidated sediment at the start of the simulation, all cells (both inside and outside polygons)
   double m_dStartIterUnconsCoarseAllCells;

   //! Depth (m) of fine consolidated sediment at the start of the simulation, all cells (both inside and outside polygons)
   double m_dStartIterConsFineAllCells;

   //! Depth (m) of sand consolidated sediment at the start of the simulation, all cells (both inside and outside polygons)
   double m_dStartIterConsSandAllCells;

   //! Depth (m) of coarse consolidated sediment at the start of the simulation, all cells (both inside and outside polygons)
   double m_dStartIterConsCoarseAllCells;

   //! Total fine unconsolidated sediment in all polygons, before polygon-to-polygon movement (only cells in polygons)
   double m_dTotalFineUnconsInPolygons;

   //! Total sand unconsolidated sediment in all polygons, before polygon-to-polygon movement (only cells in polygons)
   double m_dTotalSandUnconsInPolygons;

   //! Total coarse unconsolidated sediment in all polygons, before polygon-to-polygon movement (only cells in polygons)
   double m_dTotalCoarseUnconsInPolygons;

   //! Total fine consolidated sediment in all polygons, before polygon-to-polygon movement (only cells in polygons)
   double m_dTotalFineConsInPolygons;

   //! Total sand consolidated sediment in all polygons, before polygon-to-polygon movement (only cells in polygons)
   double m_dTotalSandConsInPolygons;

   //! Total coarse consolidated sediment in all polygons, before polygon-to-polygon movement (only cells in polygons)
   double m_dTotalCoarseConsInPolygons;

   //! Depth of unconsolidated sand sediment that could not be deposited during the last iteration, carried forward to this iteration
   double m_dUnconsSandNotDepositedLastIter;

   //! Depth of unconsolidated coarse sediment that could not be deposited during the last iteration, carried forward to this iteration
   double m_dUnconsCoarseNotDepositedLastIter;

   //! This iteration's Mean High Water (MHW) elevation, calculated using a moving time window. If we have no tide data, this is set to mean SWL
   double m_dThisIterMHWElev;

   //! Slope limit for cliff toe detection
   double m_dSlopeThresholdForCliffToe;

   // These grand totals are all long doubles. The aim is to minimize rounding errors when many very small numbers are added to a single much larger number, see e.g. http://www.ddj.com/cpp/184403224
   //! All-simulation total of potential platform erosion (m), all size classes
   long double m_ldGTotPotentialPlatformErosion;

   //! All-simulation total of fine sediment actual platform erosion (m)
   long double m_ldGTotFineActualPlatformErosion;

   //! All-simulation total of sand sediment actual platform erosion (m)
   long double m_ldGTotSandActualPlatformErosion;

   //! All-simulation total of coarse sediment actual platform erosion (m)
   long double m_ldGTotCoarseActualPlatformErosion;

   //! All-simulation total of potential sediment lost via beach (unconsolidated) sediment movement (m), all size classes
   long double m_ldGTotPotentialSedLostBeachErosion;

   //! All-simulation total of fine sediment lost via beach (unconsolidated) sediment movement (m)
   long double m_ldGTotActualFineLostBeachErosion;

   //! All-simulation total of sand sediment lost via beach (unconsolidated) sediment movement (m)
   long double m_ldGTotActualSandLostBeachErosion;

   //! All-simulation total of coarse sediment lost via beach (unconsolidated) sediment movement (m)
   long double m_ldGTotActualCoarseLostBeachErosion;

   //! All-simulation total of sand sediment lost via cliff collapse (m)
   long double m_ldGTotSandSedLostCliffCollapse;

   //! All-simulation total of coarse sediment lost via cliff collapse (m)
   long double m_ldGTotCoarseSedLostCliffCollapse;

   //! All-simulation total of fine sediment from cliff collapse (m)
   long double m_ldGTotCliffCollapseFine;

   //! All-simulation total of sand sediment from cliff collapse (m)
   long double m_ldGTotCliffCollapseSand;

   //! All-simulation total of coarse sediment from cliff collapse (m)
   long double m_ldGTotCliffCollapseCoarse;

   //! All-simulation total of fine sediment moved to suspension, due to cliff collapse (m)
   long double m_ldGTotCliffTalusFineDeposition;

   //! All-simulation total of sand sediment deposited as talus following cliff collapse (m)
   long double m_ldGTotCliffTalusSandDeposition;

   //! All-simulation total of coarse sediment deposited as talus following cliff collapse (m)
   long double m_ldGTotCliffTalusCoarseDeposition;

   //! All-simulation total of fine sediment eroded during talus deposition following cliff collapse (m)
   long double m_ldGTotCliffCollapseFineErodedDuringDeposition;

   //! All-simulation total of sand sediment eroded during talus deposition following cliff collapse (m)
   long double m_ldGTotCliffCollapseSandErodedDuringDeposition;

   //! All-simulation total of coarse sediment eroded during talus deposition following cliff collapse (m)
   long double m_ldGTotCliffCollapseCoarseErodedDuringDeposition;

   //! All-simulation total of potential beach erosion (m), all size classes
   long double m_ldGTotPotentialBeachErosion;

   //! All-simulation total of fine sediment eroded during beach (unconsolidated sediment) movement (m)
   long double m_ldGTotActualFineBeachErosion;

   //! All-simulation total of sand sediment eroded during beach (unconsolidated sediment) movement (m)
   long double m_ldGTotActualSandBeachErosion;

   //! All-simulation total of coarse sediment eroded during beach (unconsolidated sediment) movement (m)
   long double m_ldGTotActualCoarseBeachErosion;

   //! All-simulation total of sand sediment deposited during beach (unconsolidated sediment) movement (m)
   long double m_ldGTotSandBeachDeposition;

   //! All-simulation total of coarse sediment deposited during beach (unconsolidated sediment) movement (m)
   long double m_ldGTotCoarseBeachDeposition;

   //! All-simulation total of suspended sediment (m)
   long double m_ldGTotSuspendedSediment;

   //! All-simulation total of shortfall in unconsolidated sand sediment deposition (m, not currently used)
   long double m_ldGTotSandDepositionDiff;

   //! All-simulation total of shortfall in unconsolidated coarse sediment deposition (m, not currently used)
   long double m_ldGTotCoarseDepositionDiff;

   //! All-simulation total of fine sediment input (m)
   long double m_ldGTotFineSedimentInput;

   //! All-simulation total of sand sediment input (m)
   long double m_ldGTotSandSedimentInput;

   //! All-simulation total of coarse sediment input (m)
   long double m_ldGTotCoarseSedimentInput;

   //! The CME folder
   string m_strCMEDir;

   //! Folder for the CME .ini file
   string m_strCMEIni;

   //! An email addresx to which to send end-of-simulation messages
   string m_strMailAddress;

   //! Folder in which the CME data file is found
   string m_strDataPathName;

   //! Raster GIS output format
   string m_strRasterGISOutFormat;

   //! Vector GIS output format
   string m_strVectorGISOutFormat;

   //! Name of initial basement DEM file
   string m_strInitialBasementDEMFile;

   //! Name of initial landform file
   string m_strInitialLandformFile;

   //! Name of intervention class file
   string m_strInterventionClassFile;

   //! Name of intervention height file
   string m_strInterventionHeightFile;

   //! Name of initial suspended sediment file
   string m_strInitialSuspSedimentFile;

   //! Name of SCAPE shape function file
   string m_strSCAPEShapeFunctionFile;

   //! Name of tide data file
   string m_strTideDataFile;

   //! Name of output log file
   string m_strLogFile;

   //! Path for all output files
   string m_strOutPath;

   //! Name of main output file
   string m_strOutFile;

   //! GDAL code for the basement DEM raster file type
   string m_strGDALBasementDEMDriverCode;

   //! GDAL description of the basement DEM raster file type
   string m_strGDALBasementDEMDriverDesc;

   //! GDAL projection string for the basement DEM raster file
   string m_strGDALBasementDEMProjection;

   //! GDAL data type for the basement DEM raster file
   string m_strGDALBasementDEMDataType;

   //! GDAL code for the for the initial landform class raster file
   string m_strGDALLDriverCode;

   //! GDAL description of the initial landform class raster file
   string m_strGDALLDriverDesc;

   //! GDAL projection string for the initial landform class raster file
   string m_strGDALLProjection;

   //! GDAL data type for the initial landform class raster file
   string m_strGDALLDataType;

   //! GDAL code for the initial intervention class raster file
   string m_strGDALICDriverCode;

   //! GDAL description of the initial intervention class raster file
   string m_strGDALICDriverDesc;

   //! GDAL projection string for the initial intervention class raster file
   string m_strGDALICProjection;

   //! GDAL data type of the initial intervention class raster file
   string m_strGDALICDataType;

   //! GDAL code for the initial intervention height raster file
   string m_strGDALIHDriverCode;

   //! GDAL description for the initial intervention height raster file
   string m_strGDALIHDriverDesc;

   //! GDAL projection string for the initial intervention height raster file
   string m_strGDALIHProjection;

   //! GDAL data type for the initial intervention height raster file
   string m_strGDALIHDataType;

   //! GDAL code for the initial water depth raster file
   string m_strGDALIWDriverCode;

   //! GDAL description for the initial water depth raster file
   string m_strGDALIWDriverDesc;

   //! GDAL projection string for the initial water depth raster file
   string m_strGDALIWProjection;

   //! GDAL data type for the initial water depth raster file
   string m_strGDALIWDataType;

   //! GDAL code for the initial suspended sediment raster file
   string m_strGDALISSDriverCode;

   //! GDAL description for the initial suspended sediment raster file
   string m_strGDALISSDriverDesc;

   //! GDAL projection string for the initial suspended sediment raster file
   string m_strGDALISSProjection;

   //! GDAL data type for the initial suspended sediment raster file
   string m_strGDALISSDataType;

   //! GDAL code for the deep water wave stations vector file
   string m_strOGRDWWVDriverCode;

   // GDAL description for the deep water wave stations vector file
   string m_strOGRDWWVDriverDesc;

   //! GDAL geometry for the deep water wave stations vector file
   string m_strOGRDWWVGeometry;

   //! GDAL data type for the deep water wave stations vector file
   string m_strOGRDWWVDataType;

   //! GDAL code for the sediment input event locations vector file
   string m_strOGRSedInputDriverCode;

   // GDAL description for the sediment input event locations vector file
   string m_strOGRSedInputDriverDesc;

   //! GDAL geometry for the sediment input event locations vector file
   string m_strOGRSedInputGeometry;

   //! GDAL data type for the sediment input event locations vector file
   string m_strOGRSedInputDataType;

   //! GDAL code for the flood input locations point or vector file
   string m_strOGRFloodDriverCode;

   // GDAL description for the flood input locations point or vector file
   string m_strOGRFloodDriverDesc;

   //! GDAL geometry for the flood input locations point or vector file
   string m_strOGRFloodGeometry;

   //! GDAL data type for the flood input locations point or vector file
   string m_strOGRFloodDataType;

   //! GDAL raster output driver long name
   string m_strGDALRasterOutputDriverLongname;

   //! GDAL raster output driver file extension
   string m_strGDALRasterOutputDriverExtension;

   //! GDAL-OGR vector output drive file extension
   string m_strOGRVectorOutputExtension;

   //! The name of this simulation
   string m_strRunName;

   //! The duration units for this simulation
   string m_strDurationUnits;

   //! The name of the deep water wave stations shape file
   string m_strDeepWaterWaveStationsShapefile;

   //! The name of the deep water wave stations time series file
   string m_strDeepWaterWavesInputFile;

   //! The name of the sediment input events shape file
   string m_strSedimentInputEventShapefile;

   //! The name of the sediment input events time series file
   string m_strSedimentInputEventFile;

   //! The name of the flood loction events shape file
   string m_strFloodLocationShapefile;

   //! System start-simulation time
   time_t m_tSysStartTime;

   //! System finish-simulation time
   time_t m_tSysEndTime;

   //! The main output file stream
   ofstream OutStream;

   //! Sea area time series file output stream
   ofstream SeaAreaTSStream;

   //! SWL time series file output stream
   ofstream SWLTSStream;

   //! Shore platform erosion time series file output stream
   ofstream PlatformErosionTSStream;

   //! Cliff collapse erosion time series file output stream
   ofstream CliffCollapseErosionTSStream;

   //! Cliff collapse deposition time series file output stream
   ofstream CliffCollapseDepositionTSStream;

   //! Cliff collapse net change (erosion - deposition) time series file output stream
   ofstream CliffCollapseNetChangeTSStream;

   //! Beach sediment erosion time series file output stream
   ofstream BeachErosionTSStream;

   //! Beach sediment deposition time series file output stream
   ofstream BeachDepositionTSStream;

   //! Beach sediment net change (erosion - deposition) time series file output stream
   ofstream BeachSedimentNetChangeTSStream;

   //! Fine sediment in suspension time series file output stream
   ofstream FineSedSuspensionTSStream;

   //! Flood setup surge time series file output stream
   ofstream FloodSetupSurgeTSStream;

   //! Flood setup surge runup time series file output stream
   ofstream FloodSetupSurgeRunupTSStream;

   //! Cliff notch elevation time series file output stream
   ofstream CliffNotchElevTSStream;

   //! One element per layer: has the consolidated sediment of this layer been changed during this iteration?
   vector<bool> m_bConsSedChangedThisIter;

   //! One element per layer: has the consolidated sediment of this layer been changed during this iteration?
   vector<bool> m_bUnconsChangedThisIter;

   //! The numbers of the profiles which are to be saved
   vector<int> m_VnProfileToSave;

   //! ID for deep water wave station, this corresponds with the ID in the wave time series file
   vector<int> m_VnDeepWaterWaveStationID;

   //! ID for sediment input location, this corresponds with the ID in the sediment input time series file
   vector<int> m_VnSedimentInputLocationID;

   //! ID for flood location
   vector<int> m_VnFloodLocationID;

   //! Savitzky-Golay shift index for the coastline vector(s)
   vector<int> m_VnSavGolIndexCoast;

   //! Timesteps at which to save profiles
   vector<unsigned long> m_VulProfileTimestep;

   //! Calculate deep water wave values at these timesteps
   vector<unsigned long> m_VlDeepWaterWaveValuesAtTimestep;

   //! Elevations for raster slice output
   vector<double> m_VdSliceElev;

   //! For erosion potential lookup
   vector<double> m_VdErosionPotential;

   //! For erosion potential lookup
   vector<double> m_VdDepthOverDB;

   //! Savitzky-Golay filter coefficients for the coastline vector(s)
   vector<double> m_VdSavGolFCRWCoast;

   //! Savitzky-Golay filter coefficients for the profile vectors
   vector<double> m_VdSavGolFCGeomProfile;

   //! Tide data: one record per timestep, is the change (m) from still water level for that timestep
   vector<double> m_VdTideData;

   //! X coordinate (grid CRS) for deep water wave station
   vector<double> m_VdDeepWaterWaveStationX;

   //! Y coordinate (grid CRS) for deep water wave station
   vector<double> m_VdDeepWaterWaveStationY;

   //! This-iteration wave height at deep water wave station
   vector<double> m_VdThisIterDeepWaterWaveStationHeight;

   //! This-iteration wave orientation at deep water wave station
   vector<double> m_VdThisIterDeepWaterWaveStationAngle;

   //! This-iteration wave period at deep water wave station
   vector<double> m_VdThisIterDeepWaterWaveStationPeriod;

   //! Time series of wave heights at deep water wave station
   vector<double> m_VdTSDeepWaterWaveStationHeight;

   //! Time series of wave orientation at deep water wave station
   vector<double> m_VdTSDeepWaterWaveStationAngle;

   //! Time series of wave period at deep water wave station
   vector<double> m_VdTSDeepWaterWaveStationPeriod;

   //! Storage for wave transect points (real and synthetic) for debug output
   vector<TransectWaveData> m_VAllTransectsWithSynthetic;

   //! X coordinate (grid CRS) for sediment input event
   vector<double> m_VdSedimentInputLocationX;

   //! X coordinate (grid CRS) for sediment input event
   vector<double> m_VdSedimentInputLocationY;

   //! X coordinate (grid CRS) for total water level flooding
   vector<double> m_VdFloodLocationX;

   //! X coordinate (grid CRS) for total water level flooding
   vector<double> m_VdFloodLocationY;

   //! The name of the initial fine-sized unconsolidated sediment GIS file
   vector<string> m_VstrInitialFineUnconsSedimentFile;

   //! The name of the initial sand-sized unconsolidated sediment GIS file
   vector<string> m_VstrInitialSandUnconsSedimentFile;

   //! The name of the initial coarse-sized unconsolidated sediment GIS file
   vector<string> m_VstrInitialCoarseUnconsSedimentFile;

   //! The name of the initial fine-sized consolidated sediment GIS file
   vector<string> m_VstrInitialFineConsSedimentFile;

   //! The name of the initial sand-sized consolidated sediment GIS file
   vector<string> m_VstrInitialSandConsSedimentFile;

   //! The name of the initial coarse-sized consolidated sediment GIS file
   vector<string> m_VstrInitialCoarseConsSedimentFile;

   //! GDAL driver code for the initial unconsolidated fine sediment GIS data
   vector<string> m_VstrGDALIUFDriverCode;

   //! GDAL driver description for the initial unconsolidated fine sediment GIS data
   vector<string> m_VstrGDALIUFDriverDesc;

   //! GDAL projection  for the initial unconsolidated fine sediment GIS data
   vector<string> m_VstrGDALIUFProjection;

   //! GDAL data type for the initial unconsolidated fine sediment GIS data
   vector<string> m_VstrGDALIUFDataType;

   //! GDAL driver code for the initial unconsolidated sand sediment GIS data
   vector<string> m_VstrGDALIUSDriverCode;

   //! GDAL driver description for the initial unconsolidated sand sediment GIS data
   vector<string> m_VstrGDALIUSDriverDesc;

   //! GDAL projection for the initial unconsolidated sand sediment GIS data
   vector<string> m_VstrGDALIUSProjection;

   //! GDAL data type for the initial unconsolidated sand sediment GIS data
   vector<string> m_VstrGDALIUSDataType;

   //! GDAL driver code for the initial unconsolidated coarse sediment GIS data
   vector<string> m_VstrGDALIUCDriverCode;

   //! GDAL driver description for the initial unconsolidated coarse sediment GIS data
   vector<string> m_VstrGDALIUCDriverDesc;

   //! GDAL projection for the initial unconsolidated coarse sediment GIS data
   vector<string> m_VstrGDALIUCProjection;

   //! GDAL data type for the initial unconsolidated coarse sediment GIS data
   vector<string> m_VstrGDALIUCDataType;

   //! GDAL driver code for the initial consolidated fine sediment GIS data
   vector<string> m_VstrGDALICFDriverCode;

   //! GDAL driver description for the initial consolidated fine sediment GIS data
   vector<string> m_VstrGDALICFDriverDesc;

   //! GDAL projection for the initial consolidated fine sediment GIS data
   vector<string> m_VstrGDALICFProjection;

   //! GDAL data type for the initial consolidated fine sediment GIS data
   vector<string> m_VstrGDALICFDataType;

   //! GDAL driver code for the initial consolidated sand sediment GIS data
   vector<string> m_VstrGDALICSDriverCode;

   //! GDAL driver description for the initial consolidated sand sediment GIS data
   vector<string> m_VstrGDALICSDriverDesc;

   //! GDAL dprojection for the initial consolidated sand sediment GIS data
   vector<string> m_VstrGDALICSProjection;

   //! GDAL data type for the initial consolidated sand sediment GIS data
   vector<string> m_VstrGDALICSDataType;

   //! GDAL driver code for the initial consolidated coarse sediment GIS data
   vector<string> m_VstrGDALICCDriverCode;

   //! GDAL driver decription for the initial consolidated coarse sediment GIS data
   vector<string> m_VstrGDALICCDriverDesc;

   //! GDAL projection for the initial consolidated coarse sediment GIS data
   vector<string> m_VstrGDALICCProjection;

   //! GDAL data type for the initial consolidated coarse sediment GIS data
   vector<string> m_VstrGDALICCDataType;

   //! Pointer to the raster grid object
   CGeomRasterGrid* m_pRasterGrid;

   //! The coastline objects
   vector<CRWCoast> m_VCoast;

   //! Coastline (external CRS) at the highest SWL so far during this simulation
   vector<CGeomLine> m_VHighestSWLCoastLine;

   //! Coastline (external CRS) at the lowest SWL so far during this simulation
   vector<CGeomLine> m_VLowestSWLCoastLine;

   //! The traced cliff toe lines (in external CRS)
   vector<CGeomLine> m_VCliffToe;

   //! TODO 007 Finish surge and runup stuff
   vector<CRWCoast> m_VFloodWaveSetupSurge;

   //! TODO 007 Finish surge and runup stuff
   vector<CRWCoast> m_VFloodWaveSetupSurgeRunup;

   //! Corners of the bounding box
   vector<CGeom2DIPoint> m_VPtiBoundingBoxCorner;

   //! North edge cells
   vector<CGeom2DIPoint> m_VPtiNorthEdgeCell;

   //! South edge cells
   vector<CGeom2DIPoint> m_VPtiSouthEdgeCell;

   //! West edge cells
   vector<CGeom2DIPoint> m_VPtiWestEdgeCell;

   //! East edge cells
   vector<CGeom2DIPoint> m_VPtiEastEdgeCell;

   //! All edge cells
   vector<CGeom2DIPoint> m_VPtiAllEdgeCell;

   //! The location to compute the total water level for flooding
   vector<int> m_VCellFloodLocation;

   //! Sediment input events
   vector<CRWSedInputEvent*> m_pVSedInputEvent;

   //! Cells that had sediment changes this timestep (for sediment slump processing)
   set<pair<int, int>> m_prSlumpDirtyCells;

   //! The c++11 random number generators
   default_random_engine m_Rand[NUMBER_OF_RNGS];

   //! c++11 unit normal distribution (mean = 0, stdev = 1)
   normal_distribution<double> m_dGetFromUnitNormalDist{0.0, 1.0};

   //! Vector-of-vectors to hold this-iteration per-polygon fine talus amounts
   vector<vector<double>> m_VVdFineTalus;

   //! Vector-of-vectors to hold this-iteration per-polygon sand talus amounts
   vector<vector<double>> m_VVdSandTalus;

   //! Vector-of-vectors to hold this-iteration per-polygon coarse talus amounts
   vector<vector<double>> m_VVdCoarseTalus;

   //! Vector-of-vectors to hold this-iteration per-polygon fine talus amounts added
   vector<vector<double>> m_VVdFineTalusAdded;

   //! Vector-of-vectors to hold this-iteration per-polygon sand talus amounts added
   vector<vector<double>> m_VVdSandTalusAdded;

   //! Vector-of-vectors to hold this-iteration per-polygon coarse talus amounts added
   vector<vector<double>> m_VVdCoarseTalusAdded;

   //! Vector-of-vectors to hold this-iteration per-polygon fine talus amounts removed
   vector<vector<double>> m_VVdFineTalusRemoved;

   //! Vector-of-vectors to hold this-iteration per-polygon sand talus amounts removed
   vector<vector<double>> m_VVdSandTalusRemoved;

   //! Vector-of-vectors to hold this-iteration per-polygon coarse talus amounts removed
   vector<vector<double>> m_VVdCoarseTalusRemoved;

private:
   // Input and output routines
   int nHandleCommandLineParams(int, char const*[]);
   bool bReadIniFile(void);
   bool bReadIniYamlFile(void);
   bool bReadRunDataFile(void);
   bool bReadYamlFile(void);
   bool bFileIsYamlFormat(string const& strFileName);
   bool bConfigureFromDatFile(CConfiguration& config);
   bool bConfigureFromYamlFile(CConfiguration& config);
   void ApplyConfiguration(CConfiguration const& config);
   string const* processFilePath(string*);
   bool bOpenLogFile(void);
   bool bSetUpTSFiles(void);
   void WriteStartRunDetails(void);
   bool bWritePerTimestepResults(void);
   bool bWritePerTimestepResultsFixedWidth(void);
   bool bWritePerTimestepResultsCSV(void);
   bool bWriteTSFiles(void);
   int nWriteEndRunDetails(void);
   int nReadShapeFunctionFile(void);
   int nReadWaveStationInputFile(int const);
   int nReadSedimentInputEventFile(void);
   int nReadTideDataFile(void);
   int nSaveProfile(int const, CGeomProfile const*, int const, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<CGeom2DIPoint>*const, vector<double> const*) const;
   bool bWriteProfileData(int const, CGeomProfile const*, int const, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<CGeom2DIPoint>*const, vector<double> const*) const;
   int nSaveParProfile(int const, CGeomProfile const*, int const, int const, int const, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<CGeom2DIPoint>* const, vector<double> const*) const;
   bool bWriteParProfileData(int const, int const, int const, int const, int const, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<CGeom2DIPoint>* const, vector<double> const*) const;
   void WriteLookUpData(void);

   // GIS input and output stuff
   void InitializeGDALPerformance(void);
   int nReadRasterBasementDEM(void);
   int nReadRasterGISFile(int const, int const);
   int nReadVectorGISFile(int const);
   bool bWriteRasterGISFile(int const, string const*, int const = 0, double const = 0, string = "");
   bool bWriteVectorGISFile(int const, string const*, string = "");
   void GetRasterOutputMinMax(int const, double&, double&, int const, double const);
   void SetRasterFileCreationDefaults(void);
   int nInterpolateWavesToPolygonCells(vector<TransectWaveData> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*);

   // Initialization
   bool bCreateErosionPotentialLookUp(vector<double>*, vector<double>*, vector<double>*);
   void CalcMHWAndNewCliffNotchElevation(int const);

   // Top-level simulation routines
   static int nUpdateIntervention(void);
   int nCheckForSedimentInputEvent(void);
   int nCalcExternalForcing(void);
   int nInitGridAndCalcStillWaterLevel(void);
   int nLocateSeaAndCoasts(void);
   // int nLocateFloodAndCoasts(void);
   int nAssignLandformsForAllCoasts(void);
   int nAssignLandformsForAllCells(void);
   int nDoAllPropagateWaves(void);
   void GenerateSyntheticTransects(vector<TransectWaveData> const*, vector<TransectWaveData>*);
   int nDoAllShorePlatFormErosion(void);
   int nDoAllWaveEnergyToCoastLandforms(void);
   int nDoCliffCollapse(int const, int const, int const, double const);
   void DoCliffCollapseTalusDeposition(int const, int const, double const, double const, double const, int const);
   int nMoveCliffTalusToUnconsolidatedOrSuspension(void);
   double dCalcSlopeForUnconsSlumping(int const, int const, int const, int const) const;
   int nEndOfTimestepUpdateGrid(void);

   // Lower-level simulation routines
   void FindAllSeaCells(void);
   int FindAllInundatedCells(void);
   void CellByCellFillSea(int const, int const);
   bool nAdjacentCellIsSea(int const, int const);
   void FloodFillLand(int const, int const);
   int nTraceCoastLine(int const, vector<CGeom2DIPoint> const*, vector<bool>*, int const, int const, int const);
   int nTraceAllCoasts(void);
   // int nTraceFloodCoastLine(unsigned int const, int const, int const, vector<bool>*, vector<CGeom2DIPoint> const*);
   // int nTraceAllFloodCoasts(void);
   void DoCoastCurvature(int const, int const);
   int nFurtherCheckAndMarkAllProfiles(void);
   int nCreateAllProfiles(void);
   void LocateAndCreateProfiles(int const, int&, vector<bool>*, vector<pair<int, double>> const*);
   int nCreateProfile(int const, int const, int const, int const, bool const, CGeom2DIPoint const*);
   int nLocateAndCreateGridEdgeProfile(bool const, int const, int&);
   int nMarkProfilesOnGrid(void);
   void CheckAllProfilesForIntersection(void);
   static bool bCheckForIntersection(CGeomProfile* const, CGeomProfile* const, int&, int&, double&, double&, double&, double&);
   void MergeProfilesAtFinalLineSegments(int const, CGeomProfile*, CGeomProfile*, int const, int const, double const, double const, double const, double const);
   void TruncateOneProfileRetainOtherProfile(int const, CGeomProfile*, CGeomProfile*, double, double, int, int);
   int nInsertPointIntoProfilesIfNeededThenUpdate(int const, CGeomProfile*, double const, double const, int const, CGeomProfile*, int const, bool const);
   void TruncateProfileAndAppendNew(int const, CGeomProfile*, int const, vector<CGeom2DPoint> const*, vector<vector<pair<int, int>>> const*);
   void CreateRasterizedProfile(int const, CGeomProfile*, vector<CGeom2DIPoint>*, vector<bool>*, bool&, /*bool const&,*/ bool&, bool&, bool&/*, bool&*/);
   static void CalcDeanProfile(vector<double>*, double const, double const, double const, bool const, int const, double const);
   static double dSubtractProfiles(vector<double> const*, vector<double> const*, vector<bool> const*);
   int nCalcPotentialPlatformErosionOnProfile(int const, CGeomProfile*);
   int nCalcPotentialPlatformErosionBetweenProfiles(int const, CGeomProfile*, int);
   void ConstructParallelProfile(int const, int const, int const, int const, int const, vector<CGeom2DIPoint>* const, vector<CGeom2DIPoint>*, vector<CGeom2DPoint> *);
   double dCalcBeachProtectionFactor(int const, int const, double const);
   void FillInBeachProtectionHolesAndRemoveLegacyCliffs(void);
   void FillPotentialPlatformErosionHoles(void);
   void DoActualPlatformErosionOnCell(int const, int const);
   double dLookUpErosionPotential(double const);
   static CGeom2DPoint PtChooseEndPoint(int const, CGeom2DPoint const*, CGeom2DPoint const*, double const, double const, double const, double const);
   int nGetCoastNormalEndPoint(int const, int const, int const, CGeom2DPoint const*, double const, CGeom2DPoint *, CGeom2DIPoint *, bool const);
   int nLandformToGrid(int const, int const);
   int nCalcWavePropertiesOnProfile(int const, int const, CGeomProfile*, vector<double>*, vector<double>*, vector<double>*, vector<double>*, vector<bool>*);
   int nGetThisProfileElevationsForCShore(int const, CGeomProfile*, int const, vector<double>*, vector<double>*, vector<double>*);
   int nCreateCShoreInfile(int const, int const, int const, int const, int const, int const, int const, int const, int const, int const, int const, int const, int const, double const, double const, double const, double const, double const, double const, double const, double const, vector<double> const*, vector<double> const*, vector<double> const*);
   int nReadCShoreOutput(int const, string const*, int const, int const, vector<double> const*, vector<double>*);
   /*   static */ void InterpolateCShoreOutput(vector<double> const*, int const, int const, vector<double>*, vector<double>*, vector<double>*, vector<double>*, vector<double>*, vector<double>*, vector<double>*, vector<double>*, vector<double>*);
   static double dCalcWaveAngleToCoastNormal(double const, double const, int const);
   void CalcCoastTangents(int const);
   void InterpolateWavePropertiesBetweenProfiles(int const, int const);
   void InterpolateWaveHeightToCoastPoints(int const);
   // void InterpolateWavePropertiesToCells(int const, int const, int const);
   void ModifyBreakingWavePropertiesWithinShadowZoneToCoastline(int const, int const);
   static double dCalcCurvature(int const, CGeom2DPoint const*, CGeom2DPoint const*, CGeom2DPoint const*);
   void CalcD50AndFillWaveCalcHoles(void);
   int nDoAllShadowZones(void);
   static bool bOnOrOffShoreAndUpOrDownCoast(double const, double const, int const, bool&);
   static CGeom2DIPoint PtiFollowWaveAngle(CGeom2DIPoint const*, double const, double&);
   // int nFindAllShadowZones(void);
   int nCellByCellFillShadowZone(int const, int const, CGeom2DIPoint const*, CGeom2DIPoint const*, CGeom2DIPoint const*);
   void DoShadowZoneAndDownDriftZone(int const, int const, int const, int const);
   void ProcessDownDriftCell(int const, int const, int const, double const, int const);
   void ProcessShadowZoneCell(int const, int const, int const, CGeom2DIPoint const*, int const, int const, int const);
   int nCreateAllPolygons(void);
   void RasterizePolygonJoiningLine(int const, CGeom2DIPoint const*, CGeom2DIPoint const*, int const);
   int nMarkPolygonCells(void);
   int nDoPolygonSharedBoundaries(void);
   void DoAllPotentialBeachErosion(void);
   int nDoAllActualBeachErosionAndDeposition(void);
   int nDoParallelProfileUnconsErosion(CGeomCoastPolygon *, int const, int const, int const, int const, int const, int const, vector<CGeom2DIPoint> const*, vector<double> const*, double&, double&, double&);
   void ErodeCellBeachSedimentSupplyLimited(int const, int const, int const, int const, double const, double&);
   int nDoUnconsErosionOnPolygon(int const, CGeomCoastPolygon *, int const, double const, double&);
   int nDoUnconsDepositionOnPolygon(int const, CGeomCoastPolygon *, int const, double, double&);
   void CalcDepthOfClosure(void);
   int nInterpolateAllDeepWaterWaveValues(void);
   int nSetAllCoastpointDeepWaterWaveValues(void);
   int nDoSedimentInputEvent(int const);
   void AllPolygonsUpdateStoredUncons(int const);
   bool bIsInterventionCell(int const, int const) const;
   bool bSurroundedByDriftCells(int const, int const);
   bool bElevAboveDeanElev(int const, int const, double const, CRWCellLandform const*);
   // void CreatePolygonIndexIDSeq(int const);
   int nDoMultipleCoastlines(void);
   int nTruncateProfilesDifferentCoasts(int const, int const, int const, int const, int const, int const);
   int nTruncateProfileHitDifferentCoast(int const, int const, int const, int const);
   int nTruncateProfileMultiLineDifferentCoasts(CGeomProfile*, double const, double const);
   bool bIncreaseCliffNotchIncision(int const, int const, int const, CRWCliff*, double const);
   bool bCreateNotchInland(int const, int const, /*int const, int const,*/ double const, double const);
   double dCalculateSlumpInstability(int const, int const) const;
   set<pair<int, int>> prDoSlumpRedistributeSediment(int const, int const);
   void SlumpMarkCellDirty(int const, int const);
   int nDoSedimentSlumping(void);
   bool bIdentifyPossibleCoastStart(int const, int const, int const, int const, vector<CGeom2DIPoint>*);

   // GIS utility routines
   int nMarkBoundingBoxEdgeCells(void);
   void CalcGridEdgeSeaToLandDirection(void);
   bool bCheckRasterGISOutputFormat(void);
   bool bCheckVectorGISOutputFormat(void);
   bool bSaveAllRasterGISFiles(void);
   bool bSaveAllVectorGISFiles(void);
   bool bIsWithinValidGrid(int const, int const) const;
   bool bIsWithinValidGrid(CGeom2DIPoint const*) const;
   double dGridCentroidXToExtCRSX(int const) const;
   double dGridCentroidYToExtCRSY(int const) const;
   double dGridXToExtCRSX(double const) const;
   double dGridYToExtCRSY(double const) const;
   // double dExtCRSXToGridCentroidX(double const) const;
   // double dExtCRSYToGridCentroidY(double const) const;
   CGeom2DPoint PtGridCentroidToExt(CGeom2DIPoint const*) const;
   double dExtCRSXToGridX(double const) const;
   double dExtCRSYToGridY(double const) const;
   static double dGetDistanceBetween(CGeom2DPoint const*, CGeom2DPoint const*);
   static double dGetDistanceBetween(CGeom2DIPoint const*, CGeom2DIPoint const*);
   static double dGetDistanceBetween(double const, double const, double const, double const);
   static double dTriangleAreax2(CGeom2DPoint const*, CGeom2DPoint const*, CGeom2DPoint const*);
   void KeepWithinValidGrid(int&, int&) const;
   void KeepWithinValidGrid(int, int, int&, int&) const;
   void KeepWithinValidGrid(CGeom2DIPoint const*, CGeom2DIPoint *) const;
   void KeepWithinValidGrid(double&, double&, int&, int&) const;
   static double dKeepWithin360(double const);
   // vector<CGeom2DPoint> VGetPerpendicular(CGeom2DPoint const*, CGeom2DPoint
   // const*, double const, int const);
   // static CGeom2DPoint PtGetPerpendicular(CGeom2DPoint const*, CGeom2DPoint const*, double const, int const);
   static CGeom2DIPoint PtiGetPerpendicular(CGeom2DIPoint const*, CGeom2DIPoint const*, double const, int const);
   static CGeom2DIPoint PtiGetPerpendicular(int const, int const, int const, int const, double const, int const);
   static CGeom2DPoint PtAverage(CGeom2DPoint const*, CGeom2DPoint const*);
   static CGeom2DPoint PtAverage(vector<CGeom2DPoint>*);
   static CGeom2DIPoint PtiAverage(CGeom2DIPoint const*, CGeom2DIPoint const*);
   static CGeom2DIPoint PtiAverage(vector<CGeom2DIPoint>*);
   static CGeom2DIPoint PtiWeightedAverage(CGeom2DIPoint const*, CGeom2DIPoint const*, double const);
   static CGeom2DIPoint PtiPolygonCentroid(vector<CGeom2DIPoint>*);
   static double dAngleSubtended(CGeom2DIPoint const*, CGeom2DIPoint const*, CGeom2DIPoint const*);
   static int nGetOppositeDirection(int const);
   // static void GetSlopeAndInterceptFromPoints(CGeom2DIPoint const*,
   // CGeom2DIPoint const*, double&, double&);
   CGeom2DIPoint PtiFindClosestCoastPoint(int const, int const, int&);
   int nFindClosestCoastPoint(int const, int const, int&);
   int nConvertMetresToNumCells(double const) const;
   bool bIsAdjacentEdgeCell(CGeom2DIPoint const*, CGeom2DIPoint const*);
   void FindClosestPointOnStraightLine(double const, double const, double const, double const, double const, double const, double&, double&);

   // Interpolation routines
   double dGetInterpolatedValue(vector<double> const*, vector<double> const*, double, bool);
   double dGetInterpolatedValue(vector<int> const*, vector<double> const*, int, bool);
   int nFindIndex(vector<double> const*, double const);
   vector<double> VdInterpolateCShoreProfileOutput(vector<double> const*, vector<double> const*, vector<double> const*);

   // Utility routines
   static void AnnounceStart(void);
   void AnnounceLicenceAndStartConditions(void);
   void AnnounceReadBasementDEM(void) const;
   static void AnnounceAddLayers(void);
   static void AnnounceReadRasterFiles(void);
   static void AnnounceReadVectorFiles(void);
   void AnnounceReadLGIS(void) const;
   void AnnounceReadICGIS(void) const;
   void AnnounceReadIHGIS(void) const;
   static void AnnounceFinalInitialization(void);
   void AnnounceReadInitialSuspSedGIS(void) const;
   void AnnounceReadInitialFineUnconsSedGIS(int const) const;
   void AnnounceReadInitialSandUnconsSedGIS(int const) const;
   void AnnounceReadInitialCoarseUnconsSedGIS(int const) const;
   void AnnounceReadInitialFineConsSedGIS(int const) const;
   void AnnounceReadInitialSandConsSedGIS(int const) const;
   void AnnounceReadInitialCoarseConsSedGIS(int const) const;
   void AnnounceReadDeepWaterWaveValuesGIS(void) const;
   void AnnounceReadSedimentEventInputValuesGIS(void) const;
   void AnnounceReadFloodLocationGIS(void) const;
   void AnnounceReadTideData(void) const;
   static void AnnounceReadSCAPEShapeFunctionFile(void);
   static void AnnounceAllocateMemory(void);
   static void AnnounceIsRunning(void);
   static void AnnounceSimEnd(void);
   void StartClock(void);
   bool bFindExeDir(char const*);
   bool bTimeToQuit(void);
   static int nDoTimeUnits(string const*);
   int nDoSimulationTimeMultiplier(string const*);
   static double dGetTimeMultiplier(string const*);
   static bool bParseDate(string const*, int&, int&, int&);
   static bool bParseTime(string const*, int&, int&, int&);
   void DoEndOfTimestepTotals(void);
   static string strGetBuild(void);
   static string strGetComputerName(void);
   void DoCPUClockReset(void);
   // void CalcTime(double const);
   static string strDispTime(double const, bool const, bool const);
   static string strDispSimTime(double const);
   void AnnounceProgress(void);
   static string strGetErrorText(int const);
   string strListRasterFiles(void) const;
   string strListVectorFiles(void) const;
   string strListTSFiles(void) const;
   void CalcProcessStats(void);
   void CalcSavitzkyGolayCoeffs(void);
   CGeomLine LSmoothCoastSavitzkyGolay(CGeomLine *, int const, int const) const;
   CGeomLine LSmoothCoastRunningMean(CGeomLine *) const;
   vector<double> dVSmoothProfileSlope(vector<double>*) const;
   // vector<double> dVCalCGeomProfileSlope(vector<CGeom2DPoint>*, vector<double>*);         // TODO 007 Why was this removed? vector<double>
   // dVSmoothProfileSavitzkyGolay(vector<double>*, vector<double>*);         //
   // TODO 007 was this removed? vector<double>
   // dVSmoothProfileRunningMean(vector<double>*);                            //
   // TODO 007 was this removed?
   static void CalcSavitzkyGolay(double[], int const, int const, int const, int const, int const);
   static string pstrChangeToBackslash(string const*);
   static string pstrChangeToForwardSlash(string const*);
   static string strTrim(string const*);
   static string strTrimLeft(string const*);
   static string strTrimRight(string const*);
   // static string strToUpper(string const*);
   static string strRemoveSubstr(string*, string const*);
   static vector<string>* VstrSplit(string const*, char const, vector<string>*);
   static vector<string> VstrSplit(string const*, char const);
   // static double dCrossProduct(double const, double const, double const,
   // double const, double const, double const); static double
   // dGetMean(vector<double> const*); static double dGetStdDev(vector<double>
   // const*);
   static void AppendEnsureNoGap(vector<CGeom2DIPoint>*, CGeom2DIPoint const*);
   // static bool bIsNumeric(string const*);
   unsigned long ulConvertToTimestep(string const*) const;
   void WritePolygonInfoTable(void);
   void WritePolygonPreExistingSedimentTable(void);
   void WritePolygonSedimentInputEventTable(void);
   void WritePolygonShorePlatformErosion(void);
   void WritePolygonCliffCollapseErosion(void);
   void WritePolygonSedimentBeforeMovement(void);
   void WritePolygonPotentialErosion(void);
   // void WritePolygonUnconsErosion(int const);
   void WritePolygonUnsortedSequence(vector<vector<vector<int>>> const&);
   void WritePolygonSortedSequence(vector<vector<vector<int>>> const&);
   void WritePolygonActualMovement(vector<vector<vector<int>>> const&);
   void DoEndOfRunDeletes(void);
   void DoRand1Shuffle(vector<int>*);

// #ifdef _DEBUG
//    void DEBUG_PrintProfileDetails(CGeomProfile*, CGeomProfile*);
//    void DEBUG_PrintPolygonDetails(int const, CGeomCoastPolygon*);
// #endif

 protected:
 public:
   ofstream LogStream;

   CSimulation(void);
   ~CSimulation(void);

   //! Returns the NODATA value
   double dGetMissingValue(void) const;

   //! Returns this timestep's SWL
   double dGetThisIterSWL(void) const;

   //! Returns this timestep's total water level TODO 007 Finish surge and runup stuff
   double dGetThisIterTotWaterLevel(void) const;

   // //! Returns the vertical tolerance for beach cells to be included in smoothing
   // double dGetMaxBeachElevAboveSWL(void) const;

   //! Returns the cell size
   // double dGetCellSide(void) const;

   //! Returns the size of the grid in the X direction
   int nGetGridXMax(void) const;

   //! Returns the size of the grid in the Y direction
   int nGetGridYMax(void) const;

   //! Returns the global d50 value for fine sediment
   double dGetD50Fine(void) const;

   //! Returns the global d50 value for sand sediment
   double dGetD50Sand(void) const;

   //! Returns the global d50 value for coarse sediment
   double dGetD50Coarse(void) const;

   //! Runs the simulation
   int nDoSimulation(int, char const*[]);

   //! Carries out end-of-simulation tidying (error messages etc.)
   void DoSimulationEnd(int const);

   //! Returns a pointer to the lower case version of a string
   static string strToLower(string const*);

   //! Transforms a pointer to a CGeom2DPoint in the external CRS to the equivalent CGeom2DIPoint in the raster grid CRS (both values rounded). The result may be outside the grid
   CGeom2DIPoint PtiExtCRSToGridRound(CGeom2DPoint const*) const;

   //! Constrains the supplied point (in the grid CRS) to be a valid cell within the raster grid
   void KeepWithinValidGrid(CGeom2DIPoint*) const;

   //! Determines whether a point (grid CRS) is within the polygon (but is not on a polygon edge)
   bool bIsWithinPolygon(CGeom2DIPoint const*, vector<CGeom2DIPoint> const*);
};
#endif // SIMULATION_H
