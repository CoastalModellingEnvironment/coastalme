/*!
   \file configuration.h
   \brief Unified configuration class for CoastalME simulation parameters
   \details Provides a single interface for accessing simulation parameters regardless of input format (.dat or YAML)
   \author Wilf Chun
   \author David Favis-Mortlock
   \author Andres Payo
   \date 2026
   \copyright GNU General Public License
*/

#ifndef CONFIGURATION_H
#define CONFIGURATION_H
/* ==============================================================================================================================
   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public  License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
==============================================================================================================================*/

#include <algorithm>

#include <cctype>

#include <string>
using std::string;

#include <vector>
using std::vector;

#include "simulation.h"

//! Unified configuration class for CoastalME simulation parameters
class CConfiguration
{
private:
   // Run information
   string m_strRunName;
   int m_nLogFileDetail;
   bool m_bCSVPerTimestepResults;

   // Simulation timing
   string m_strStartDateTime;
   string m_strDuration;
   string m_strTimestep;
   vector<string> m_VstrSaveTimes;
   int m_nRandomSeed;
   bool m_bUseSystemTimeForSeed;

   // GIS output
   int m_nMaxSaveDigits;
   string m_strSaveDigitsMode;
   vector<string> m_VstrRasterFiles;
   string m_strRasterFormat;
   bool m_bWorldFile;
   bool m_bScaleValues;
   vector<double> m_VdSliceElevations;
   vector<string> m_VstrVectorFiles;
   string m_strVectorFormat;
   vector<string> m_VstrTimeSeriesFiles;

   // Grid and coastline
   int m_nCoastlineSmoothing;
   int m_nCoastlineSmoothingWindow;
   int m_nPolynomialOrder;
   string m_strOmitGridEdges;
   int m_nProfileSmoothingWindow;
   double m_dMaxLocalSlope;
   double m_dMaxBeachElevation;

   // Layers and files
   int m_nNumLayers;
   string m_strBasementDEMFile;
   vector<string> m_VstrUnconsFineFiles;
   vector<string> m_VstrUnconsSandFiles;
   vector<string> m_VStrUnconsCoarseFiles;
   vector<string> m_VstrConsFineFiles;
   vector<string> m_VstrConsSandFiles;
   vector<string> m_VstrConsCoarseFiles;
   string m_strSuspendedSedFile;
   string m_strLandformFile;
   string m_strInterventionClassFile;
   string m_strInterventionHeightFile;

   // Hydrology
   int m_nWavePropagationModel;
   double m_dSeawaterDensity;
   double m_dInitialWaterLevel;
   double m_dFinalWaterLevel;
   bool m_bbHasFinalWaterLevel;

   // Waves
   string m_strWaveInputMode;
   string m_strWaveHeightTimeSeries;
   string m_strWaveStationDataFile;
   double m_dDeepWaterWaveHeight;
   double m_dDeepWaterWaveOrientation;
   double m_dWavePeriod;

   // tides
   string m_strTideDataFile;
   double m_dBreakingWaveRatio;

   // Sediment, erosion, and deposition
   bool m_bCoastPlatformErosion;
   double m_dPlatformErosionResistance;
   bool m_bBeachSedimentTransport;
   int m_nBeachTransportAtEdges;
   int m_nBeachErosionEquation;
   double m_dFineMedianSize;
   double m_dSandMedianSize;
   double m_dCoarseMedianSize;
   double m_dSedimentDensity;
   double m_dBeachSedimentPorosity;
   double m_dFineErosivity;
   double m_dSandErosivity;
   double m_dCoarseErosivity;
   double m_dTransportKLS;
   double m_dKamphuis;
   double m_dBermHeight;

   // Cliff collapse
   bool m_bCliffCollapse;
   double m_dCliffErosionResistance;
   double m_dNotchOverhang;
   double m_dNotchBase;
   double m_dTalusWidth;
   double m_dMinTalusLength;
   double m_dMinTalusHeight;
   double m_dCliffCollapseTalusErodibility;
   double m_dCliffCollapseFineTalusRemovalRate;
   double m_dCliffCollapseSandTalusRemovalRate;
   double m_dCliffCollapseCoarseTalusRemovalRate;

   // Riverine flooding
   bool m_bFloodInput;
   string m_strFloodCoastline;
   int m_nRunUpEquation;
   string m_strFloodLocations;
   string m_strFloodInputLocation;
   vector<string> m_VstrFloodFiles;

   // Sediment input
   bool m_bSedimentInput;
   string m_strSedimentInputLocation;
   string m_strSedimentInputType;
   string m_strSedimentInputDetails;

   // Slumping
   bool m_bSlumping;

   // Barrier
   bool m_bWaveUprush;

   // Physics and geometry
   double m_dGravitationalAcceleration;
   double m_dNormalSpacing;
   double m_dRandomFactor;
   double m_dNormalLength;
   double m_dStartDepthRatio;
   double m_dTemporaryProfileSpacing;

   // Profile and output options
   bool m_bSaveProfileData;
   vector<int> m_VnProfileNumbers;
   vector<unsigned long> m_VulProfileTimesteps;
   bool m_bSaveParallelProfiles;
   bool m_bOutputErosionPotential;
   int m_nCurvatureWindow;
   int m_nSmoothProfileMethod;

public:
   CConfiguration(void);
   ~CConfiguration(void);

   // Setters for all parameters
   void SetRunName(string const* pStr)
   {
      m_strRunName = *pStr;
   }

   void SetLogFileDetail(int n)
   {
      m_nLogFileDetail = n;
   }

   void SetCSVPerTimestepResults(bool b)
   {
      m_bCSVPerTimestepResults = b;
   }

   void SetStartDateTime(string const* pStr)
   {
      m_strStartDateTime = *pStr;
   }

   void SetDuration(string const* pStr)
   {
      m_strDuration = *pStr;
   }

   void SetTimestep(string const* pStr)
   {
      m_strTimestep = *pStr;
   }

   void SetSaveTimes(vector<string> const* pVStr)
   {
      m_VstrSaveTimes = *pVStr;
   }

   void SetRandomSeed(int n)
   {
      m_nRandomSeed = n;
      m_bUseSystemTimeForSeed = false;
   }
   void UseSystemTimeForSeed(void)
   {
      m_bUseSystemTimeForSeed = false;
   }

   void SetMaxSaveDigits(int n)
   {
      m_nMaxSaveDigits = n;
   }

   void SetSaveDigitsMode(string const* pStr)
   {
      m_strSaveDigitsMode = *pStr;
   }

   void SetRasterFiles(vector<string> const* pVStr)
   {
      m_VstrRasterFiles = *pVStr;
   }

   void SetRasterFormat(string const* pStr)
   {
      m_strRasterFormat = *pStr;
   }

   void SetWorldFile(bool b)
   {
      m_bWorldFile = b;
   }

   void SetScaleValues(bool b)
   {
      m_bScaleValues = b;
   }

   void SetSliceElevations(vector<double> const* pVStr)
   {
      m_VdSliceElevations = *pVStr;
   }

   void SetVectorFiles(vector<string> const* pVStr)
   {
      m_VstrVectorFiles = *pVStr;
   }

   void SetVectorFormat(string const* pStr)
   {
      m_strVectorFormat = *pStr;
   }

   void SetTimeSeriesFiles(vector<string> const &vec)
   {
      m_VstrTimeSeriesFiles = vec;
   }

   void SetCoastlineSmoothing(int n)
   {
      m_nCoastlineSmoothing = n;
   }

   void SetCoastlineSmoothingWindow(int n)
   {
      m_nCoastlineSmoothingWindow = n;
   }

   void SetPolynomialOrder(int n)
   {
      m_nPolynomialOrder = n;
   }

   void SetOmitGridEdges(string const* pStr)
   {
      m_strOmitGridEdges = *pStr;
   }

   void SetProfileSmoothingWindow(int n)
   {
      m_nProfileSmoothingWindow = n;
   }

   void SetMaxLocalSlope(double d)
   {
      m_dMaxLocalSlope = d;
   }
   void SetMaxBeachElevation(double d)
   {
      m_dMaxBeachElevation = d;
   }

   void SetNumLayers(int n)
   {
      m_nNumLayers = n;
   }

   void SetBasementDEMFile(string const* pStr)
   {
      m_strBasementDEMFile = *pStr;
   }

   void SetUnconsFineFiles(vector<string> const &vec)
   {
      m_VstrUnconsFineFiles = vec;
   }

   void SetUnconsSandFiles(vector<string> const &vec)
   {
      m_VstrUnconsSandFiles = vec;
   }

   void SetUnconsCoarseFiles(vector<string> const &vec)
   {
      m_VStrUnconsCoarseFiles = vec;
   }

   void SetConsFineFiles(vector<string> const &vec)
   {
      m_VstrConsFineFiles = vec;
   }

   void SetConsSandFiles(vector<string> const &vec)
   {
      m_VstrConsSandFiles = vec;
   }

   void SetConsCoarseFiles(vector<string> const &vec)
   {
      m_VstrConsCoarseFiles = vec;
   }

   void SetSuspendedSedFile(string const* pStr)
   {
      m_strSuspendedSedFile = *pStr;
   }

   void SetLandformFile(string const* pStr)
   {
      m_strLandformFile = *pStr;
   }

   void SetInterventionClassFile(string const* pStr)
   {
      m_strInterventionClassFile = *pStr;
   }

   void SetInterventionHeightFile(string const* pStr)
   {
      m_strInterventionHeightFile = *pStr;
   }

   void SetWavePropagationModel(int n)
   {
      m_nWavePropagationModel = n;
   }

   void SetSeawaterDensity(double d)
   {
      m_dSeawaterDensity = d;
   }

   void SetInitialWaterLevel(double d)
   {
      m_dInitialWaterLevel = d;
   }

   void SetFinalWaterLevel(double d)
   {
      m_dFinalWaterLevel = d;
      m_bbHasFinalWaterLevel = true;
   }

   // Wave height Data
   void SetWaveInputMode(string const* pStr)
   {
      m_strWaveInputMode = *pStr;
   }

   void SetWaveHeightTimeSeries(string const* pStr)
   {
      m_strWaveHeightTimeSeries = *pStr;
   }

   void SetWaveStationDataFile(string const* pStr)
   {
      m_strWaveStationDataFile = *pStr;
   }

   void SetDeepWaterWaveHeight(double d)
   {
      m_dDeepWaterWaveHeight = d;
   }

   void SetDeepWaterWaveOrientation(double d)
   {
      m_dDeepWaterWaveOrientation = d;
   }

   void SetWavePeriod(double d)
   {
      m_dWavePeriod = d;
   }

   void SetTideDataFile(string const* pStr)
   {
      m_strTideDataFile = *pStr;
   }

   void SetBreakingWaveRatio(double d)
   {
      m_dBreakingWaveRatio = d;
   }

   // Additional setters for comprehensive YAML support
   void SetCoastPlatformErosion(bool b)
   {
      m_bCoastPlatformErosion = b;
   }

   void SetPlatformErosionResistance(double d)
   {
      m_dPlatformErosionResistance = d;
   }

   void SetBeachSedimentTransport(bool b)
   {
      m_bBeachSedimentTransport = b;
   }

   void SetBeachTransportAtEdges(int n)
   {
      m_nBeachTransportAtEdges = n;
   }

   void SetBeachErosionEquation(int n)
   {
      m_nBeachErosionEquation = n;
   }

   void SetFineMedianSize(double d)
   {
      m_dFineMedianSize = d;
   }

   void SetSandMedianSize(double d)
   {
      m_dSandMedianSize = d;
   }

   void SetCoarseMedianSize(double d)
   {
      m_dCoarseMedianSize = d;
   }

   void SetSedimentDensity(double d)
   {
      m_dSedimentDensity = d;
   }

   void SetBeachSedimentPorosity(double d)
   {
      m_dBeachSedimentPorosity = d;
   }

   void SetFineErosivity(double d)
   {
      m_dFineErosivity = d;
   }

   void SetSandErosivity(double d)
   {
      m_dSandErosivity = d;
   }

   void SetCoarseErosivity(double d)
   {
      m_dCoarseErosivity = d;
   }
   void SetTransportKLS(double d)
   {
      m_dTransportKLS = d;
   }

   void SetKamphuis(double d)
   {
      m_dKamphuis = d;
   }

   void SetBermHeight(double d)
   {
      m_dBermHeight = d;
   }

   void SetCliffCollapse(bool b)
   {
      m_bCliffCollapse = b;
   }

   void SetCliffErosionResistance(double d)
   {
      m_dCliffErosionResistance = d;
   }

   void SetNotchOverhang(double d)
   {
      m_dNotchOverhang = d;
   }

   void SetNotchBase(double d)
   {
      m_dNotchBase = d;
   }

   void SetTalusWidth(double d)
   {
      m_dTalusWidth = d;
   }

   void SetMinTalusLength(double d)
   {
      m_dMinTalusLength = d;
   }

   void SetMinTalusHeight(double d)
   {
      m_dMinTalusHeight = d;
   }

   void SetFloodInput(bool b)
   {
      m_bFloodInput = b;
   }

   void SetFloodFiles(vector<string> const*v)
   {
      m_VstrFloodFiles = *v;
   }

   void SetFloodCoastline(string const* pStr)
   {
      m_strFloodCoastline = *pStr;
   }

   void SetRunUpEquation(string const* pStr)
   {
      if ((*pStr == "") or (*pStr == " "))
      {
         m_nRunUpEquation = 0;
      }
      else
      {
         m_nRunUpEquation = stoi(*pStr);
      }
   }

   void SetFloodLocations(string const* pStr)
   {
      m_strFloodLocations = *pStr;
   }

   void SetFloodInputLocation(string const* pStr)
   {
      m_strFloodInputLocation = *pStr;
   }

   void SetSedimentInput(bool b)
   {
      m_bSedimentInput = b;
   }

   void SetSedimentInputLocation(string const* pStr)
   {
      m_strSedimentInputLocation = *pStr;
   }

   void SetSedimentInputType(string const* pStr)
   {
      m_strSedimentInputType = *pStr;
   }

   void SetSedimentInputDetails(string const* pStr)
   {
      m_strSedimentInputDetails = *pStr;
   }

   void SetSlumping(bool b)
   {
      m_bSlumping = b;
   }

   void SetBarrier(bool b)
   {
      m_bWaveUprush = b;
   }

   void SetGravitationalAcceleration(double d)
   {
      m_dGravitationalAcceleration = d;
   }

   void SetNormalSpacing(double d)
   {
      m_dNormalSpacing = d;
   }

   void SetRandomFactor(double d)
   {
      m_dRandomFactor = d;
   }

   void SetNormalLength(double d)
   {
      m_dNormalLength = d;
   }

   void SetStartDepthRatio(double d)
   {
      m_dStartDepthRatio = d;
   }

   void SetSyntheticTransectSpacing(double d)
   {
      m_dTemporaryProfileSpacing = d;
   }

   void SetSaveProfileData(bool b)
   {
      m_bSaveProfileData = b;
   }

   void SetProfileNumbers(vector<int> const &vec)
   {
      m_VnProfileNumbers = vec;
   }

   void SetProfileTimesteps(vector<unsigned long> const &vec)
   {
      m_VulProfileTimesteps = vec;
   }

   void SetSaveParallelProfiles(bool b)
   {
      m_bSaveParallelProfiles = b;
   }

   void SetOutputErosionPotential(bool b)
   {
      m_bOutputErosionPotential = b;
   }

   void SetCurvatureWindow(int n)
   {
      m_nCurvatureWindow = n;
   }

   void SetProfileSmoothingMethod(int i)
   {
      m_nSmoothProfileMethod = i;
   }

   void SetCliffCollapseTalusErodibility(double d)
   {
      m_dCliffCollapseTalusErodibility = d;
   }

   void SetCliffCollapseFineTalusRemovalRate(double d)
   {
      m_dCliffCollapseFineTalusRemovalRate = d;
   }

   void SetCliffCollapseSandTalusRemovalRate(double d)
   {
      m_dCliffCollapseSandTalusRemovalRate = d;
   }

   void SetCliffCollapseCoarseTalusRemovalRate(double d)
   {
      m_dCliffCollapseCoarseTalusRemovalRate = d;
   }

   // Getters for all parameters
   string const* pstrGetRunName(void) const
   {
      return &m_strRunName;
   }

   int nGetLogFileDetail(void) const
   {
      return m_nLogFileDetail;
   }

   bool bGetCSVPerTimestepResults(void) const
   {
      return m_bCSVPerTimestepResults;
   }

   string const* pstrGetStartDateTime(void) const
   {
      return &m_strStartDateTime;
   }

   string const* pstrGetDuration(void) const
   {
      return &m_strDuration;
   }

   string const* pstrGetTimestep(void) const
   {
      return &m_strTimestep;
   }

   vector<string> const* pVstrGetSaveTimes(void) const
   {
      return &m_VstrSaveTimes;
   }

   int nGetRandomSeed(void) const
   {
      return m_nRandomSeed;
   }

   bool bUseSystemTimeForRandomSeed(void) const
   {
      return m_bUseSystemTimeForSeed;
   }

   int nGetMaxSaveDigits(void) const
   {
      return m_nMaxSaveDigits;
   }

   string const* pstrGetSaveDigitsMode(void) const
   {
      return &m_strSaveDigitsMode;
   }

   void GetRasterFiles(vector<string>*) const;

   string const* pstrGetRasterFormat(void) const
   {
      return &m_strRasterFormat;
   }

   bool bGetWorldFile(void) const
   {
      return m_bWorldFile;
   }

   bool bGetScaleValues(void) const
   {
      return m_bScaleValues;
   }

   vector<double> const* pVdGetSliceElevations(void) const
   {
      return &m_VdSliceElevations;
   }

   void GetVectorFiles(vector<string>*) const;

   string const* pstrGetVectorFormat(void) const
   {
      return &m_strVectorFormat;
   }

   void GetTimeSeriesFiles(vector<string>*) const;

   int nGetCoastlineSmoothing(void) const
   {
      return m_nCoastlineSmoothing;
   }

   int nGetCoastlineSmoothingWindow(void) const
   {
      return m_nCoastlineSmoothingWindow;
   }

   int nGetPolynomialOrder(void) const
   {
      return m_nPolynomialOrder;
   }

   string const strGetOmitGridEdges(void) const;

   int nGetProfileSmoothingWindow(void) const
   {
      return m_nProfileSmoothingWindow;
   }

   double nGetMaxLocalSlope(void) const
   {
      return m_dMaxLocalSlope;
   }

   double dGetMaxBeachElevation(void) const
   {
      return m_dMaxBeachElevation;
   }

   int nGetNumLayers(void) const
   {
      return m_nNumLayers;
   }

   string const* pstrGetBasementDEMFile(void) const
   {
      return &m_strBasementDEMFile;
   }

   vector<string> const* pVstrGetUnconsFineFiles(void) const
   {
      return &m_VstrUnconsFineFiles;
   }

   vector<string> const* pVstrGetUnconsSandFiles(void) const
   {
      return &m_VstrUnconsSandFiles;
   }

   vector<string> const* pVstrGetUnconsCoarseFiles(void) const
   {
      return &m_VStrUnconsCoarseFiles;
   }

   vector<string> const* pVstrGetConsFineFiles(void) const
   {
      return &m_VstrConsFineFiles;
   }

   vector<string> const* pVstrGetConsSandFiles(void) const
   {
      return &m_VstrConsSandFiles;
   }

   vector<string> const* pVstrGetConsCoarseFiles(void) const
   {
      return &m_VstrConsCoarseFiles;
   }

   string const* pstrGetSuspendedSedFile(void) const
   {
      return &m_strSuspendedSedFile;
   }

   string const* pstrGetLandformFile(void) const
   {
      return &m_strLandformFile;
   }

   string const* pstrGetInterventionClassFile(void) const
   {
      return &m_strInterventionClassFile;
   }

   string const* pstrGetInterventionHeightFile(void) const
   {
      return &m_strInterventionHeightFile;
   }

   int nGetWavePropagationModel(void) const
   {
      return m_nWavePropagationModel;
   }

   double dGetSeawaterDensity(void) const
   {
      return m_dSeawaterDensity;
   }

   double dGetInitialWaterLevel(void) const
   {
      return m_dInitialWaterLevel;
   }

   double dGetFinalWaterLevel(void) const
   {
      return m_dFinalWaterLevel;
   }

   bool bHasFinalWaterLevel(void) const
   {
      return m_bbHasFinalWaterLevel;
   }

   string const* pstrGetWaveInputMode(void) const
   {
      return &m_strWaveInputMode;
   }

   // Wave data configuration getters (Cases 37-40)
   string const* pstrGetWaveHeightTimeSeries(void) const
   {
      return &m_strWaveHeightTimeSeries;
   }

   string const* pstrGetWaveStationDataFile(void) const
   {
      return &m_strWaveStationDataFile;
   }

   double dGetDeepWaterWaveHeight(void) const
   {
      return m_dDeepWaterWaveHeight;
   }

   double dGetDeepWaterWaveOrientation(void) const
   {
      return m_dDeepWaterWaveOrientation;
   }

   double dGetWavePeriod(void) const
   {
      return m_dWavePeriod;
   }

   string const* pstrGetTideDataFile(void) const
   {
      return &m_strTideDataFile;
   }

   double dGetBreakingWaveRatio(void) const
   {
      return m_dBreakingWaveRatio;
   }

   // Sediment and Erosion parameters
   bool bGetCoastPlatformErosion(void) const
   {
      return m_bCoastPlatformErosion;
   }

   double dGetPlatformErosionResistance(void) const
   {
      return m_dPlatformErosionResistance;
   }

   bool bGetBeachSedimentTransport(void) const
   {
      return m_bBeachSedimentTransport;
   }

   int nGetBeachTransportAtEdges(void) const
   {
      return m_nBeachTransportAtEdges;
   }

   int nGetBeachErosionEquation(void) const
   {
      return m_nBeachErosionEquation;
   }

   double dGetFineMedianSize(void) const
   {
      return m_dFineMedianSize;
   }

   double dGetSandMedianSize(void) const
   {
      return m_dSandMedianSize;
   }

   double dGetCoarseMedianSize(void) const
   {
      return m_dCoarseMedianSize;
   }

   double dGetSedimentDensity(void) const
   {
      return m_dSedimentDensity;
   }

   double dGetBeachSedimentPorosity(void) const
   {
      return m_dBeachSedimentPorosity;
   }

   double dGetFineErosivity(void) const
   {
      return m_dFineErosivity;
   }

   double dGetSandErosivity(void) const
   {
      return m_dSandErosivity;
   }

   double dGetCoarseErosivity(void) const
   {
      return m_dCoarseErosivity;
   }

   double dGetTransportKLS(void) const
   {
      return m_dTransportKLS;
   }

   double dGetKamphuis(void) const
   {
      return m_dKamphuis;
   }

   double dGetBermHeight(void) const
   {
      return m_dBermHeight;
   }

   // Cliff parameters
   bool bGetCliffCollapse(void) const
   {
      return m_bCliffCollapse;
   }

   double dGetCliffErosionResistance(void) const
   {
      return m_dCliffErosionResistance;
   }

   double dGetNotchOverhang(void) const
   {
      return m_dNotchOverhang;
   }

   double dGetNotchBase(void) const
   {
      return m_dNotchBase;
   }

   double dGetTalusWidth(void) const
   {
      return m_dTalusWidth;
   }

   double dGetMinTalusLength(void) const
   {
      return m_dMinTalusLength;
   }

   double dGetMinTalusHeight(void) const
   {
      return m_dMinTalusHeight;
   }

   // Flood parameters
   bool bGetFloodInput(void) const
   {
      return m_bFloodInput;
   }

   void GetFloodFiles(vector<string>*) const;

   string const* pstrGetFloodCoastline(void) const
   {
      return &m_strFloodCoastline;
   }

   int nGetRunUpEquation(void) const
   {
      return m_nRunUpEquation;
   }

   string const* pstrGetFloodLocations(void) const
   {
      return &m_strFloodLocations;
   }

   string const* pstrGetFloodInputLocation(void) const
   {
      return &m_strFloodInputLocation;
   }

   // Sediment Input parameters
   bool bGetSedimentInput(void) const
   {
      return m_bSedimentInput;
   }

   string const* pstrGetSedimentInputLocation(void) const
   {
      return &m_strSedimentInputLocation;
   }

   string const* pstrGetSedimentInputType(void) const
   {
      return &m_strSedimentInputType;
   }

   string const* pstrGetSedimentInputDetails(void) const
   {
      return &m_strSedimentInputDetails;
   }

   // Slumping
   bool bGetSlumping(void) const
   {
      return m_bSlumping;
   }

   // Barriers
   bool bGetBarrier(void) const
   {
      return m_bWaveUprush;
   }

   // Physics and geometry parameters
   double dGetGravitationalAcceleration(void) const
   {
      return m_dGravitationalAcceleration;
   }

   double dGetNormalSpacing(void) const
   {
      return m_dNormalSpacing;
   }

   double dGetRandomFactor(void) const
   {
      return m_dRandomFactor;
   }

   double dGetNormalLength(void) const
   {
      return m_dNormalLength;
   }

   double dGetStartDepthRatio(void) const
   {
      return m_dStartDepthRatio;
   }

   double dGetSyntheticTransectSpacing(void) const
   {
      return m_dTemporaryProfileSpacing;
   }

   // Profile and Output Options
   bool bGetSaveProfileData(void) const
   {
      return m_bSaveProfileData;
   }

   vector<int> const* pVnGetProfileNumbers(void) const
   {
      return &m_VnProfileNumbers;
   }

   vector<unsigned long> const* pVulGetProfileTimesteps(void) const
   {
      return &m_VulProfileTimesteps;
   }

   bool bGetSaveParallelProfiles(void) const
   {
      return m_bSaveParallelProfiles;
   }

   bool bGetOutputErosionPotential(void) const
   {
      return m_bOutputErosionPotential;
   }

   int nGetCurvatureWindow(void) const
   {
      return m_nCurvatureWindow;
   }

   int nGetProfileSmoothingMethod(void) const
   {
      return m_nSmoothProfileMethod;
   }

   double dGetCliffCollapseTalusErodibility(void) const
   {
      return m_dCliffCollapseTalusErodibility;
   }

   double dGetCliffCollapseFineTalusRemovalRate(void) const
   {
      return m_dCliffCollapseFineTalusRemovalRate;
   }

   double dGetCliffCollapseSandTalusRemovalRate(void) const
   {
      return m_dCliffCollapseSandTalusRemovalRate;
   }

   double dGetCliffCollapseCoarseTalusRemovalRate(void) const
   {
      return m_dCliffCollapseCoarseTalusRemovalRate;
   }

   // Initialize with default values
   void InitializeDefaults(void);
};
#endif      //CONFIGURATION_H
