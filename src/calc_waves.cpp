/*!
   \file calc_waves.cpp
   \brief Simulates wave propagation using CShore or the COVE approach
   \details TODO 001 A more detailed description of these routines.
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

#include <cmath>
using std::sin;
using std::cos;
using std::pow;
using std::fmod;

#include <string>
using std::to_string;

#include <iostream>
using std::endl;
using std::ios;

#include <fstream>
using std::ifstream;

#include <algorithm>
using std::reverse_copy;

#include <utility>
using std::move;

#include "cme.h"
#include "coast.h"
#include "simulation.h"
#include "cshore/cshore.h"
#include "2d_point.h"

//===============================================================================================================================
//! Give every coast point a value for deep water wave height and direction TODO 005 This may not be realistic, maybe better to use end-of-profile value instead (how?)
//===============================================================================================================================
int CSimulation::nSetAllCoastpointDeepWaterWaveValues(void)
{
   LogStream << endl << m_ulIter << ": Calculating waves" << endl;

   // For each coastline, put a value for deep water wave height and direction at each coastline point
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      bool bCoastEnd = false;
      int nDistFromPrevProfile = 0;
      int nDistToNextProfile = 0;

      double dPrevProfileDeepWaterWaveHeight = 0;
      double dPrevProfileDeepWaterWaveAngle = 0;
      double dPrevProfileDeepWaterWavePeriod = 0;
      double dNextProfileDeepWaterWaveHeight = 0;
      double dNextProfileDeepWaterWaveAngle = 0;
      double dNextProfileDeepWaterWavePeriod = 0;
      double dDist = 0;

      for (int nPoint = 0; nPoint < m_VCoast[nCoast].nGetCoastlineSize(); nPoint++)
      {
         // We are going down-coast (i.e. along the coast in the direction of increasing coastline point numbers)
         if (m_VCoast[nCoast].bIsProfileAtCoastPoint(nPoint))
         {
            // OK, a coastline-normal profile begins at this coastline point, so set the deep water wave values at this coastline point to be the values at the seaward end of the coastline normal
            CGeomProfile const* pProfile = m_VCoast[nCoast].pGetProfileAtCoastPoint(nPoint);
            // int nProfile = pProfile->nGetProfileID();

            double const dThisDeepWaterWaveHeight = pProfile->dGetProfileDeepWaterWaveHeight();
            double const dThisDeepWaterWaveAngle = pProfile->dGetProfileDeepWaterWaveAngle();
            double const dThisDeepWaterWavePeriod = pProfile->dGetProfileDeepWaterWavePeriod();

            m_VCoast[nCoast].SetCoastDeepWaterWaveHeight(nPoint, dThisDeepWaterWaveHeight);
            m_VCoast[nCoast].SetCoastDeepWaterWaveAngle(nPoint, dThisDeepWaterWaveAngle);
            m_VCoast[nCoast].SetCoastDeepWaterWavePeriod(nPoint, dThisDeepWaterWavePeriod);

            // Reset for next time
            nDistFromPrevProfile = 0;
            dPrevProfileDeepWaterWaveHeight = dThisDeepWaterWaveHeight;
            dPrevProfileDeepWaterWaveAngle = dThisDeepWaterWaveAngle;
            dPrevProfileDeepWaterWavePeriod = dThisDeepWaterWavePeriod;

            // Find the next profile
            CGeomProfile const* pNextProfile = pProfile->pGetDownCoastAdjacentProfile();

            if (pNextProfile == NULL)
            {
               // We are at the end of the coast
               bCoastEnd = true;

               break;
            }

            // And the distance (in along-coast points) to the next profile
            nDistToNextProfile = pNextProfile->nGetCoastPoint() - nPoint;
            dDist = nDistToNextProfile;

            // And the next profile's deep water wave values
            dNextProfileDeepWaterWaveHeight = pNextProfile->dGetProfileDeepWaterWaveHeight();
            dNextProfileDeepWaterWaveAngle = pNextProfile->dGetProfileDeepWaterWaveAngle();
            dNextProfileDeepWaterWavePeriod = pNextProfile->dGetProfileDeepWaterWavePeriod();

            // LogStream << m_ulIter << ": coast point = " << nPoint << " is start of profile " << pProfile->nGetProfileID() << ", next profile is " << pNextProfile->nGetProfileID() << ", which starts at coast point " << pNextProfile->nGetCoastPoint() << ", dThisDeepWaterWaveHeight = " << dThisDeepWaterWaveHeight << ", dThisDeepWaterWaveAngle = " << dThisDeepWaterWaveAngle << " nDistToNextProfile = " << nDistToNextProfile << endl;
         }
         else
         {
            // This coast point is not the start of a coastline normal, so set the deep water wave values at this coastline point to be a weighted average of those from the up-coast and down-coast profiles
            nDistFromPrevProfile++;
            nDistToNextProfile--;

            double const dPrevWeight = (dDist - nDistFromPrevProfile) / dDist;
            double const dNextWeight = (dDist - nDistToNextProfile) / dDist;
            double const dThisDeepWaterWaveHeight = (dPrevWeight * dPrevProfileDeepWaterWaveHeight) + (dNextWeight * dNextProfileDeepWaterWaveHeight);
            double const dThisDeepWaterWaveAngle = dKeepWithin360((dPrevWeight * dPrevProfileDeepWaterWaveAngle) + (dNextWeight * dNextProfileDeepWaterWaveAngle));
            double const dThisDeepWaterWavePeriod = (dPrevWeight * dPrevProfileDeepWaterWavePeriod) + (dNextWeight * dNextProfileDeepWaterWavePeriod);

            m_VCoast[nCoast].SetCoastDeepWaterWaveHeight(nPoint, dThisDeepWaterWaveHeight);
            m_VCoast[nCoast].SetCoastDeepWaterWaveAngle(nPoint, dThisDeepWaterWaveAngle);
            m_VCoast[nCoast].SetCoastDeepWaterWavePeriod(nPoint, dThisDeepWaterWavePeriod);

            // LogStream << m_ulIter << ": coast point = " << nPoint << " dThisDeepWaterWaveHeight = " << dThisDeepWaterWaveHeight << " dThisDeepWaterWaveAngle = " << dThisDeepWaterWaveAngle << " nDistFromPrevProfile = " << nDistFromPrevProfile << " nDistToNextProfile = " << nDistToNextProfile << endl;
         }
      }

      if (bCoastEnd)
         break;
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Creates temporary profiles between existing coast-normal profiles, to increase the density densify of data points for wave interpolation. The number of synthetic profiles created between each pair depends on both the distance between the, and the desired spacing. Uses OpenMP parallelization to concurrently process multiple temporary profiles. OPTIMIZED: Pre-computes external CRS positions to avoid redundant coordinate conversions
//===============================================================================================================================
void CSimulation::CreateTemporaryProfiles(vector<ProfileWaveData> const* pVRealProfiles, vector<ProfileWaveData>* pVAllProfiles)
{
   // Start with all the real profiles
   *pVAllProfiles = *pVRealProfiles;

   // If no spacing specified or spacing is too large, just return the real ones
   if (m_dTemporaryProfileSpacing <= 0)
      return;

   int const nNumRealProfiles = static_cast<int>(pVRealProfiles->size());
   if (nNumRealProfiles < 2)
      return;

   // OPTIMIZATION: Pre-compute external CRS positions for all profileData start points. This avoids redundant coordinate conversions (each profileData appears in 2 pairs)
   vector<double> VdProfileExtX(nNumRealProfiles);
   vector<double> VdProfileExtY(nNumRealProfiles);

   for (int i = 0; i < nNumRealProfiles; i++)
   {
      ProfileWaveData const& profileData = (*pVRealProfiles)[i];
      if (!profileData.VdX.empty())
      {
         VdProfileExtX[i] = dGridCentroidXToExtCRSX(static_cast<int>(profileData.VdX[0]));
         VdProfileExtY[i] = dGridCentroidYToExtCRSY(static_cast<int>(profileData.VdY[0]));
      }
      else
      {
         VdProfileExtX[i] = 0.0;
         VdProfileExtY[i] = 0.0;
      }
   }

   // First pass: calculate distances and determine how many temporary profiles are needed for each pair
   vector<int> VnNumTempProfilesPerPair(nNumRealProfiles - 1, 0);
   int nTotalTempProfiles = 0;

   for (int nPair = 0; nPair < nNumRealProfiles - 1; nPair++)
   {
      ProfileWaveData const& profileData1 = (*pVRealProfiles)[nPair];
      ProfileWaveData const& profileData2 = (*pVRealProfiles)[nPair + 1];

      // Only interpolate between profiles on the same coast
      if (profileData1.nCoastID != profileData2.nCoastID)
         continue;

      // Skip if either profile has no points
      if (profileData1.VdX.empty() || profileData2.VdX.empty())
         continue;

      // Use pre-computed external CRS positions
      double const dX1 = VdProfileExtX[nPair];
      double const dY1 = VdProfileExtY[nPair];
      double const dX2 = VdProfileExtX[nPair + 1];
      double const dY2 = VdProfileExtY[nPair + 1];

      // Calculate distance using pre-computed positions
      double const dDX = dX2 - dX1;
      double const dDY = dY2 - dY1;
      double const dDistance = sqrt(dDX * dDX + dDY * dDY);

      // Calculate number of temporary profiles needed: distance / spacing, rounded to nearest integer. If distance <= spacing, we don't add any temporary profiles
      int nNumTemporarys = 0;
      if (dDistance > m_dTemporaryProfileSpacing)
      {
         nNumTemporarys = static_cast<int>(round(dDistance / m_dTemporaryProfileSpacing));
         // Subtract 1 because we already have the two endpoints (real profiles)
         nNumTemporarys = std::max(0, nNumTemporarys - 1);
      }

      VnNumTempProfilesPerPair[nPair] = nNumTemporarys;
      nTotalTempProfiles += nNumTemporarys;
   }

   if (nTotalTempProfiles <= 0)
   {
      LogStream << m_ulIter << ": No temporary profiles needed (all profile spacings <= " << m_dTemporaryProfileSpacing << " m)" << endl;
      return;
   }

   // Pre-allocate space for temporary profiles
   vector<ProfileWaveData> VTemporaryProfilesWaveData(nTotalTempProfiles);

   // Second pass: generate the temporary profiles. Use OpenMP to parallelize, each thread handles one pair
   // int nCurrentIndex = 0;

   #pragma omp parallel for schedule(dynamic)
   for (int nPair = 0; nPair < nNumRealProfiles - 1; nPair++)
   {
      int const nNumTemporarys = VnNumTempProfilesPerPair[nPair];
      if (nNumTemporarys == 0)
         continue;

      ProfileWaveData const& profileData1 = (*pVRealProfiles)[nPair];
      ProfileWaveData const& profileData2 = (*pVRealProfiles)[nPair + 1];

      // Calculate the starting index for this pair's temporary profiles
      int nPairStartIndex = 0;
      for (int i = 0; i < nPair; i++)
         nPairStartIndex += VnNumTempProfilesPerPair[i];

      // Create temporary profiles for this pair
      for (int nSynth = 1; nSynth <= nNumTemporarys; nSynth++)
      {
         int const nSynthIndex = nPairStartIndex + (nSynth - 1);
         ProfileWaveData& tempProfile = VTemporaryProfilesWaveData[nSynthIndex];

         // Calculate interpolation weight (0 < alpha < 1)
         double const dAlpha = static_cast<double>(nSynth) / (nNumTemporarys + 1);
         double const dOneMinusAlpha = 1.0 - dAlpha;

         // Set metadata for temporary profile
         tempProfile.nCoastID = profileData1.nCoastID;
         tempProfile.nProfileID = -1; // Mark as synthetic
         tempProfile.bIsStartOrEndOfCoast = false;

         // Interpolate points - use the minimum length to avoid extrapolation
         size_t const nMinLength = std::min(profileData1.VdX.size(), profileData2.VdX.size());

         // Reserve space for efficiency
         tempProfile.VdX.reserve(nMinLength);
         tempProfile.VdY.reserve(nMinLength);
         tempProfile.VdHeightX.reserve(nMinLength);
         tempProfile.VdHeightY.reserve(nMinLength);
         tempProfile.VbBreaking.reserve(nMinLength);

         // Interpolate each point along the profileData
         for (size_t i = 0; i < nMinLength; i++)
         {
            // Linear interpolation of position
            tempProfile.VdX.push_back(dOneMinusAlpha * profileData1.VdX[i] + dAlpha * profileData2.VdX[i]);
            tempProfile.VdY.push_back(dOneMinusAlpha * profileData1.VdY[i] + dAlpha * profileData2.VdY[i]);

            // Linear interpolation of wave height components
            tempProfile.VdHeightX.push_back(dOneMinusAlpha * profileData1.VdHeightX[i] + dAlpha * profileData2.VdHeightX[i]);
            tempProfile.VdHeightY.push_back(dOneMinusAlpha * profileData1.VdHeightY[i] + dAlpha * profileData2.VdHeightY[i]);

            // Breaking status: true if either parent profileData has breaking at this point
            tempProfile.VbBreaking.push_back(profileData1.VbBreaking[i] || profileData2.VbBreaking[i]);
         }
      }
   }

   // Append all temporary profiles to the output vector
   pVAllProfiles->insert(pVAllProfiles->end(), VTemporaryProfilesWaveData.begin(), VTemporaryProfilesWaveData.end());

   LogStream << m_ulIter << ": Generated " << nTotalTempProfiles << " temporary profiles between " << nNumRealProfiles << " real profiles (target spacing: " << m_dTemporaryProfileSpacing << " m)" << endl;
}

//===============================================================================================================================
//! Simulates wave propagation along all coastline-normal profiles, on all coasts
//===============================================================================================================================
int CSimulation::nDoAllPropagateWaves(void)
{
   // Set up vector to hold wave data for each coast-normal profile
   vector<ProfileWaveData> VAllProfilesWaveData;

   // Calculate wave properties for every coast
   // bool bSomeNonStartOrEndOfCoastProfiles = false;

   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      int const nCoastSize = m_VCoast[nCoast].nGetCoastlineSize();
      int const nNumProfiles = m_VCoast[nCoast].nGetNumProfiles();

      static bool bDownCoast = true;

      // Calculate wave properties at every point along each valid coast-normal profile, and for the cells under the profiles. Do this alternately in up-coast and down-coast sequence
      for (int nn = 0; nn < nNumProfiles; nn++)
      {
         ProfileWaveData profileData;
         CGeomProfile* pProfile;

         if (bDownCoast)
            pProfile = m_VCoast[nCoast].pGetProfileWithDownCoastSeq(nn);
         else
            pProfile = m_VCoast[nCoast].pGetProfileWithUpCoastSeq(nn);

         int const nRet = nCalcWavePropertiesOnProfile(nCoast, nCoastSize, pProfile, &profileData.VdX, &profileData.VdY, &profileData.VdHeightX, &profileData.VdHeightY, &profileData.VbBreaking);
         if (nRet != RTN_OK)
         {
            if (nRet == RTN_ERR_CSHORE_ERROR)
            {
               // Abandon calculations on this profile, and flag the profile
               pProfile->SetCShoreProblem();

               // Move on to next profile
            }
            else
            {
               // A serious CShore error, so abort the run
               return nRet;
            }
         }
         else
         {
            // CShore worked fine on this profile. Is the profile a start of coast or end of coast profile?
            if (! pProfile->bIsStartOrEndOfCoast())
            {
               // It is neither a start of coast or an end of coast profile, so set switch
               // bSomeNonStartOrEndOfCoastProfiles = true;
            }

            // Are the waves on this profile off-shore? If so, do nothing more with this profile, the wave values for cells have already been given the off-shore value
            if (profileData.VbBreaking.empty())
            {
               continue;
            }

            // Store metadata about this profileData
            profileData.nCoastID = nCoast;
            profileData.nProfileID = pProfile->nGetProfileID();
            profileData.bIsStartOrEndOfCoast = pProfile->bIsStartOrEndOfCoast();

            // Add this profileData to the collection
            VAllProfilesWaveData.push_back(move(profileData));
         }
      }

      bDownCoast = ! bDownCoast;
   }

   // OK, do we have some profiles other than start of coast or end of coast profiles in the all-profile vectors? We need to check this, because GDALGridCreate() in nInterpolateWavePropertiesToWithinPolygonCells() does not work if we give it only a start-of-coast or an end-of-coast profile to work with TODO 006 Is this still true?
   // if (! bSomeNonStartOrEndOfCoastProfiles)
   // {
   //    LogStream << m_ulIter << ":\t waves are on-shore only, for start and/or end of coast profiles" << endl;
   //
   //    return RTN_OK;
   // }

   // We need to also send the deepwater points from the edge of the grid to nInterpolateWavePropertiesToWithinPolygonCells(), this is necessary to prevent GDALGridCreate() leaving holes in the interpolated grid when the polygons are far from regular
   // Store deep water grid edge points separately (not as a profile)
   vector<double> VdDeepWaterX;
   vector<double> VdDeepWaterY;
   vector<double> VdDeepWaterHeightX;
   vector<double> VdDeepWaterHeightY;

   double dDeepWaterWaveX;
   double dDeepWaterWaveY;

   if (m_bSingleDeepWaterWaveValues)
   {
      // Just using the same value of deep water wave hehght and angle for all cells
      dDeepWaterWaveX = m_dAllCellsDeepWaterWaveHeight * sin(m_dAllCellsDeepWaterWaveAngle * PI / 180),
      dDeepWaterWaveY = m_dAllCellsDeepWaterWaveHeight * cos(m_dAllCellsDeepWaterWaveAngle * PI / 180);
   }

   for (int nX = 0; nX < m_nXGridSize; nX++)
   {
      if (m_pRasterGrid->m_Cell[nX][0].bIsInContiguousSea())
      {
         int const nPolyID = m_pRasterGrid->m_Cell[nX][0].nGetPolygonID();

         if (nPolyID == INT_NODATA)
         {
            // Not in a polygon
            VdDeepWaterX.push_back(nX);
            VdDeepWaterY.push_back(0);

            if (! m_bSingleDeepWaterWaveValues)
            {
               // Not using the same value of deep water height and angle for all cells, so get this cell's deep water height and angle values
               dDeepWaterWaveX = m_pRasterGrid->m_Cell[nX][0].dGetCellDeepWaterWaveHeight() * sin(m_pRasterGrid->m_Cell[nX][0].dGetCellDeepWaterWaveAngle() * PI / 180);
               dDeepWaterWaveY = m_pRasterGrid->m_Cell[nX][0].dGetCellDeepWaterWaveHeight() * cos(m_pRasterGrid->m_Cell[nX][0].dGetCellDeepWaterWaveAngle() * PI / 180);
            }

            VdDeepWaterHeightX.push_back(dDeepWaterWaveX);
            VdDeepWaterHeightY.push_back(dDeepWaterWaveY);
         }
      }

      if (m_pRasterGrid->m_Cell[nX][m_nYGridSize - 1].bIsInContiguousSea())
      {
         int const nPolyID = m_pRasterGrid->m_Cell[nX][m_nYGridSize - 1].nGetPolygonID();

         if (nPolyID == INT_NODATA)
         {
            // Not in a polygon
            VdDeepWaterX.push_back(nX);
            VdDeepWaterY.push_back(m_nYGridSize - 1);

            if (! m_bSingleDeepWaterWaveValues)
            {
               // Not using the same value of deep water height and angle for all cells, so get this cell's deep water height and angle values
               dDeepWaterWaveX = m_pRasterGrid->m_Cell[nX][m_nYGridSize - 1].dGetCellDeepWaterWaveHeight() * sin(m_pRasterGrid->m_Cell[nX][m_nYGridSize - 1].dGetCellDeepWaterWaveAngle() * PI / 180);
               dDeepWaterWaveY = m_pRasterGrid->m_Cell[nX][m_nYGridSize - 1].dGetCellDeepWaterWaveHeight() * cos(m_pRasterGrid->m_Cell[nX][m_nYGridSize - 1].dGetCellDeepWaterWaveAngle() * PI / 180);
            }

            VdDeepWaterHeightX.push_back(dDeepWaterWaveX);
            VdDeepWaterHeightY.push_back(dDeepWaterWaveY);
         }
      }
   }

   for (int nY = 0; nY < m_nYGridSize; nY++)
   {
      if (m_pRasterGrid->m_Cell[0][nY].bIsInContiguousSea())
      {
         int const nPolyID = m_pRasterGrid->m_Cell[0][nY].nGetPolygonID();

         if (nPolyID == INT_NODATA)
         {
            // Not in a polygon
            VdDeepWaterX.push_back(0);
            VdDeepWaterY.push_back(nY);

            if (! m_bSingleDeepWaterWaveValues)
            {
               // Not using the same value of deep water height and angle for all cells, so get this cell's deep water height and angle values
               dDeepWaterWaveX = m_pRasterGrid->m_Cell[0][nY].dGetCellDeepWaterWaveHeight() * sin(m_pRasterGrid->m_Cell[0][nY].dGetCellDeepWaterWaveAngle() * PI / 180);
               dDeepWaterWaveY = m_pRasterGrid->m_Cell[0][nY].dGetCellDeepWaterWaveHeight() * cos(m_pRasterGrid->m_Cell[0][nY].dGetCellDeepWaterWaveAngle() * PI / 180);
            }

            VdDeepWaterHeightX.push_back(dDeepWaterWaveX);
            VdDeepWaterHeightY.push_back(dDeepWaterWaveY);
         }
      }

      if (m_pRasterGrid->m_Cell[m_nXGridSize - 1][nY].bIsInContiguousSea())
      {
         int const nPolyID = m_pRasterGrid->m_Cell[m_nXGridSize - 1][nY].nGetPolygonID();

         if (nPolyID == INT_NODATA)
         {
            // Not in a polygon
            VdDeepWaterX.push_back(m_nXGridSize - 1);
            VdDeepWaterY.push_back(nY);

            if (! m_bSingleDeepWaterWaveValues)
            {
               // Not using the same value of deep water height and angle for all cells, so get this cell's deep water height and angle values
               dDeepWaterWaveX = m_pRasterGrid->m_Cell[m_nXGridSize - 1][nY].dGetCellDeepWaterWaveHeight() * sin(m_pRasterGrid->m_Cell[m_nXGridSize - 1][nY].dGetCellDeepWaterWaveAngle() * PI / 180);
               dDeepWaterWaveY = m_pRasterGrid->m_Cell[m_nXGridSize - 1][nY].dGetCellDeepWaterWaveHeight() * cos(m_pRasterGrid->m_Cell[m_nXGridSize - 1][nY].dGetCellDeepWaterWaveAngle() * PI / 180);
            }

            VdDeepWaterHeightX.push_back(dDeepWaterWaveX);
            VdDeepWaterHeightY.push_back(dDeepWaterWaveY);
         }
      }
   }

   // Are the waves off-shore for every profile? If so, do nothing more. So check if any profile has breaking waves
   bool bHasBreakingWaves = false;
   for (const auto &profileData : VAllProfilesWaveData)
   {
      if (! profileData.VbBreaking.empty())
      {
         bHasBreakingWaves = true;
         break;
      }
   }

   if (! bHasBreakingWaves)
   {
      LogStream << m_ulIter << ": waves off-shore for all profiles" << endl;
      return RTN_OK;
   }

   // Generate temporary profiles to create more points for better interpolation
   vector<ProfileWaveData> VAllProfilesIncTemporary;
   CreateTemporaryProfiles(&VAllProfilesWaveData, &VAllProfilesIncTemporary);

   // Store profileDatas for potential debug output
   m_VAllProfilesIncTemporary = VAllProfilesIncTemporary;

   // Some waves are on-shore, so interpolate the wave attributes from all profile points to all within-polygon sea cells, also update the active zone status for each cell
   int nRet = nInterpolateWavesToPolygonCells(&VAllProfilesIncTemporary, &VdDeepWaterX, &VdDeepWaterY, &VdDeepWaterHeightX, &VdDeepWaterHeightY);

   if (nRet != RTN_OK)
      return nRet;

   // // DEBUG CODE ===========================================================================================================
   // string strOutFile = m_strOutPath;
   // strOutFile += "sea_wave_height_CHECKPOINT_1_";
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
   // nRet = pBand->RasterIO(GF_Write, 0, 0, m_nXGridSize, m_nYGridSize, pdRaster, m_nXGridSize, m_nYGridSize, GDT_Float64, 0, 0, NULL);
   //
   // if (nRet == CE_Failure)
   // return RTN_ERR_GRIDCREATE;
   //
   // GDALClose(pDataSet);
   // delete[] pdRaster;
   // // DEBUG CODE ===========================================================================================================

   // // DEBUG CODE ===========================================================================================================
   // strOutFile = m_strOutPath;
   // strOutFile += "sea_wave_angle_CHECKPOINT_1_";
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
   // // DEBUG CODE ===========================================================================================================

   // Find any shadow zones and then modify waves in and adjacent to them
   nRet = nDoAllShadowZones();
   if (nRet != RTN_OK)
      return nRet;

   // Calculate the D50 for each polygon. Also fill in any artefactual 'holes' in active zone and wave property patterns
   CalcD50AndFillWaveCalcHoles();

   // // DEBUG CODE ===========================================================================================================
   // strOutFile = m_strOutPath;
   // strOutFile += "sea_wave_height_CHECKPOINT_2_";
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
   // pdRaster[nn++] = m_pRasterGrid->m_Cell[nX][nY].dGetWaveHeight();
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
   // // DEBUG CODE ===========================================================================================================

   // // DEBUG CODE ===========================================================================================================
   // strOutFile = m_strOutPath;
   // strOutFile += "sea_wave_angle_CHECKPOINT_2_";
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
   // // DEBUG CODE ===========================================================================================================

   // Modify the wave breaking properties (wave height, wave dir, breaking depth, breaking distance) for coastline points within the shadow zone
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      int const nNumProfiles = m_VCoast[nCoast].nGetNumProfiles();

      for (int nProfile = 0; nProfile < nNumProfiles; nProfile++)
         ModifyBreakingWavePropertiesWithinShadowZoneToCoastline(nCoast, nProfile);
   }

   // // DEBUG CODE ===========================================================================================================
   // strOutFile = m_strOutPath;
   // strOutFile += "sea_wave_height_CHECKPOINT_3_";
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
   // pdRaster[nn++] = m_pRasterGrid->m_Cell[nX][nY].dGetWaveHeight();
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
   // // DEBUG CODE ===========================================================================================================

   // // DEBUG CODE ===========================================================================================================
   // strOutFile = m_strOutPath;
   // strOutFile += "sea_wave_angle_CHECKPOINT_3_";
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
   // // DEBUG CODE ===========================================================================================================

   // Interpolate these wave properties for all remaining coastline points. Do this in along-coastline sequence
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      int const nCoastSize = m_VCoast[nCoast].nGetCoastlineSize();
      int const nNumProfiles = m_VCoast[nCoast].nGetNumProfiles();

      // Interpolate these wave properties for all remaining coastline points. Do this in along-coastline sequence, but do not do this for the end-of-coastline profile (which is the final one)
      for (int n = 0; n < nNumProfiles - 1; n++)
         InterpolateWavePropertiesBetweenProfiles(nCoast, n);

      // And do the same for the coastline cells
      InterpolateWaveHeightToCoastPoints(nCoast);

      // Calculate wave energy at breaking for every point on the coastline using equation 4 from Walkden & Hall, 2005
      for (int nCoastPoint = 0; nCoastPoint < nCoastSize; nCoastPoint++)
      {
         double const dBreakingWaveHeight = m_VCoast[nCoast].dGetBreakingWaveHeight(nCoastPoint);
         double const dCoastPointWavePeriod = m_VCoast[nCoast].dGetCoastDeepWaterWavePeriod(nCoastPoint);
         double dWaveEnergy = 0;

         if (bFPIsEqual(dBreakingWaveHeight, DBL_NODATA, TOLERANCE))
         {
            m_VCoast[nCoast].SetBreakingWaveHeight(nCoastPoint, 0);
         }
         else if (dBreakingWaveHeight <= 0)
         {
            // Safety check
            m_VCoast[nCoast].SetBreakingWaveHeight(nCoastPoint, 0);
         }
         else
         {
            double const dErosiveWaveForce = pow(dBreakingWaveHeight, WALKDEN_HALL_PARAM_1) * pow(dCoastPointWavePeriod, WALKDEN_HALL_PARAM_2);

            // Calculate total wave energy at this coast point during this timestep
            dWaveEnergy = dErosiveWaveForce * m_dTimeStep * 3600;
            m_VCoast[nCoast].SetWaveEnergyAtBreaking(nCoastPoint, dWaveEnergy);
         }

         if (m_nLogFileDetail >= LOG_FILE_ALL)
         {
            int const nX = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint)->nGetX();
            int const nY = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint)->nGetY();

            LogStream << m_ulIter << ": calc wave energy at breaking, nCoastPoint = " << nCoastPoint << " [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "} breaking wave height = " << m_VCoast[nCoast].dGetBreakingWaveHeight(nCoastPoint) << " dCoastPointWavePeriod = " << dCoastPointWavePeriod << " wave energy = " << m_VCoast[nCoast].dGetWaveEnergyAtBreaking(nCoastPoint) << endl;
         }
      }
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Calculates the angle between the wave direction and a normal to the coastline tangent. If wave direction has a component which is down-coast (i.e. in the direction with increasing coast point numbers), then the angle returned is +ve. If wave direction has a component which is up-coast (i.e. in the direction with decreasing coast point numbers), then the angle returned is -ve. If waves are in an off-shore direction, DBL_NODATA is returned
//===============================================================================================================================
double CSimulation::dCalcWaveAngleToCoastNormal(double const dCoastAngle, double const dWaveAngle, int const nSeaHand)
{
   double dWaveToNormalAngle = 0;

   if (nSeaHand == LEFT_HANDED)
      // Left-handed coast
      dWaveToNormalAngle = fmod((dWaveAngle - dCoastAngle + 360), 360) - 90;

   else
      // Right-handed coast
      dWaveToNormalAngle = fmod((dWaveAngle - dCoastAngle + 360), 360) - 270;

   if ((dWaveToNormalAngle >= 90) || (dWaveToNormalAngle <= -90))
      dWaveToNormalAngle = DBL_NODATA;

   return dWaveToNormalAngle;
}

//===============================================================================================================================
//! Calculates wave properties along a coastline-normal profile using either the COVE linear wave theory approach or the external CShore model
//===============================================================================================================================
int CSimulation::nCalcWavePropertiesOnProfile(int const nCoast, int const nCoastSize, CGeomProfile* pProfile, vector<double>* pVdX, vector<double>* pVdY, vector<double>* pVdHeightX, vector<double>* pVdHeightY, vector<bool>* pVbBreaking)
{
   // Only work on this profile if it is either a start-of-coast or an end-of-coast profile, or it is problem-free
   if (! pProfile->bProfileOKOrStartEndProfile())
   {
      if (m_nLogFileDetail >= LOG_FILE_ALL)
         LogStream << m_ulIter << ":\t coast " << nCoast << ", profile " << pProfile->nGetProfileID() << " is not OK, will not calc wave properties on this profile" << endl;

      return RTN_OK;
   }

   // Calculate some wave properties based on the wave period following Airy wave theory
   double const dDeepWaterWavePeriod = pProfile->dGetProfileDeepWaterWavePeriod();

   m_dC_0 = (m_dG * dDeepWaterWavePeriod) / (2 * PI); // Deep water (offshore) wave celerity (m/s)
   m_dL_0 = m_dC_0 * dDeepWaterWavePeriod;            // Deep water (offshore) wave length (m)

   int const nSeaHand = m_VCoast[nCoast].nGetSeaHandedness();
   int const nCoastPoint = pProfile->nGetCoastPoint();

   // Get the flux orientation (the orientation of a line which is tangential to the coast) at adjacent coastline points. Note special treatment for the coastline end points
   double const dFluxOrientationThis = m_VCoast[nCoast].dGetFluxOrientation(nCoastPoint);
   double dFluxOrientationPrev = 0;
   double dFluxOrientationNext = 0;

   if (nCoastPoint == 0)
   {
      dFluxOrientationPrev = dFluxOrientationThis;
      dFluxOrientationNext = m_VCoast[nCoast].dGetFluxOrientation(1);
   }
   else if (nCoastPoint == nCoastSize - 1)
   {
      dFluxOrientationPrev = m_VCoast[nCoast].dGetFluxOrientation(nCoastPoint - 2);
      dFluxOrientationNext = dFluxOrientationThis;
   }
   else
   {
      dFluxOrientationPrev = m_VCoast[nCoast].dGetFluxOrientation(nCoastPoint - 1);
      dFluxOrientationNext = m_VCoast[nCoast].dGetFluxOrientation(nCoastPoint + 1);
   }

   // Get the deep water wave orientation for this profile
   double const dDeepWaterWaveAngle = pProfile->dGetProfileDeepWaterWaveAngle();

   // Calculate the angle between the deep water wave direction and a normal to the coast tangent
   double dWaveToNormalAngle = dCalcWaveAngleToCoastNormal(dFluxOrientationThis, dDeepWaterWaveAngle, nSeaHand);

   // Are the waves off-shore?
   if (bFPIsEqual(dWaveToNormalAngle, DBL_NODATA, TOLERANCE))
   {
      // They are, so do nothing (each cell under the profile has already been initialised with deep water wave height and wave direction)

      // LogStream << m_ulIter << ": profile " << nProfile << " has sea to " << (m_VCoast[nCoast].nGetSeaHandedness() == RIGHT_HANDED ? "right" : "left") << " dWaveToNormalAngle = " << dWaveToNormalAngle << " which is off-shore" << endl;

      return RTN_OK;
   }

   // Calculate the angle between the deep water wave direction and a normal to the coast tangent for the previous coast point
   double dWaveToNormalAnglePrev;

   if (nCoastPoint > 0)
   {
      // Get the deep water wave orientation for the up-coast point
      double const dPrevDeepWaterWaveAngle = m_VCoast[nCoast].dGetCoastDeepWaterWaveAngle(nCoastPoint - 1);

      dWaveToNormalAnglePrev = dCalcWaveAngleToCoastNormal(dFluxOrientationPrev, dPrevDeepWaterWaveAngle, nSeaHand);
   }
   else
   {
      dWaveToNormalAnglePrev = dWaveToNormalAngle;
   }

   // Calculate the angle between the deep water wave direction and a normal to the coast tangent for the next coast point
   double dWaveToNormalAngleNext;

   if (nCoastPoint < nCoastSize - 1)
   {
      // Get the deep water wave orientation for the down-coast point
      double const dNextDeepWaterWaveAngle = m_VCoast[nCoast].dGetCoastDeepWaterWaveAngle(nCoastPoint + 1);

      dWaveToNormalAngleNext = dCalcWaveAngleToCoastNormal(dFluxOrientationNext, dNextDeepWaterWaveAngle, nSeaHand);
   }
   else
   {
      dWaveToNormalAngleNext = dWaveToNormalAngle;
   }

   // Following Ashton and Murray (2006), if we have high-angle waves then use the flux orientation of the previous (up-coast) profile, if transitioning from diffusive to antidiffusive use flux maximizing angle (45 degrees)
   if ((dWaveToNormalAngle > 0) && (! bFPIsEqual(dWaveToNormalAnglePrev, DBL_NODATA, TOLERANCE)) && (dWaveToNormalAnglePrev > 0))
   {
      if (dWaveToNormalAngle > 45)
      {
         if (dWaveToNormalAnglePrev < 45)
         {
            dWaveToNormalAngle = 45;
            // LogStream << "\tA1" << endl;
         }
         else
         {
            dWaveToNormalAngle = dWaveToNormalAnglePrev;
            // LogStream << "\tA2" << endl;
         }
      }
   }
   else if ((dWaveToNormalAngle < 0) && (! bFPIsEqual(dWaveToNormalAngleNext, DBL_NODATA, TOLERANCE)) && (dWaveToNormalAngleNext < 0))
   {
      if (dWaveToNormalAngle < -45)
      {
         if (dWaveToNormalAngleNext > -45)
         {
            dWaveToNormalAngle = -45;
            // LogStream << "\tB1" << endl;
         }
         else
         {
            dWaveToNormalAngle = dWaveToNormalAngleNext;
            // LogStream << "\tB2" << endl;
         }
      }
   }
   else if ((dWaveToNormalAngle > 45) && (! bFPIsEqual(dWaveToNormalAnglePrev, DBL_NODATA, TOLERANCE)) && (dWaveToNormalAnglePrev > 0))
   {
      // The wave direction here has an up-coast (decreasing indices) component: so for high-angle waves use the orientation from the up-coast (previous) profile
      // LogStream << "\tCCC" << endl;

      dWaveToNormalAngle = dFluxOrientationPrev;
   }
   else if ((dWaveToNormalAngle < -45) && (! bFPIsEqual(dWaveToNormalAngleNext, DBL_NODATA, TOLERANCE)) && (dWaveToNormalAngleNext < 0))
   {
      // The wave direction here has a down-coast (increasing indices) component: so for high-angle waves use the orientation from the down-coast (next) profile
      // LogStream << "\tDDD" << endl;

      dWaveToNormalAngle = dFluxOrientationNext;
   }

   // Initialize the wave properties at breaking for this profile
   bool bBreaking = false;

   int const nProfileSize = pProfile->nGetNumCellsInProfile();
   int nProfileBreakingDist = 0;

   double dProfileBreakingWaveHeight = DBL_NODATA;
   double dProfileBreakingWaveAngle = 0;
   double dProfileBreakingDepth = 0;
   double dProfileWaveHeight = DBL_NODATA;
   double const dProfileDeepWaterWaveHeight = pProfile->dGetProfileDeepWaterWaveHeight();
   double dProfileWaveAngle = DBL_NODATA;
   double const dProfileDeepWaterWaveAngle = pProfile->dGetProfileDeepWaterWaveAngle();

   vector<bool> VbWaveIsBreaking(nProfileSize, 0);

   vector<double> VdWaveHeight(nProfileSize, 0);
   vector<double> VdWaveSetupSurge(nProfileSize, 0);
   // vector<double> VdStormSurge(nProfileSize, 0);
   vector<double> const VdWaveSetupRunUp(nProfileSize, 0);
   vector<double> VdWaveDirection(nProfileSize, 0);

   if (m_nWavePropagationModel == WAVE_MODEL_CSHORE)
   {
      // We are using CShore to propagate the waves
      double const dCShoreTimeStep = 3600;      // In seconds, not important because we are not using CShore to erode the profile, just to get the hydrodynamics
      double const dSurgeLevel = CSHORE_SURGE_LEVEL;

      // Set up vectors for the coastline-normal profile elevations. The length of this vector line is given by the number of cells 'under' the profile. Thus each point on the vector relates to a single cell in the grid. This assumes that all points on the profile vector are equally spaced (not quite true, depends on the orientation of the line segments which comprise the profile)
      vector<double> VdProfileZ;                // Initial (pre-erosion) elevation of both consolidated and unconsolidated sediment for cells 'under' the profile, in CShore units
      vector<double> VdProfileDistXY;           // Along-profile distance measured from the seaward limit, in CShore units
      vector<double> VdProfileFrictionFactor;   // Along-profile friction factor from seaward limit

      // The elevation of each of these profile points is the elevation of the centroid of the cell that is 'under' the point. However we cannot always be confident that this is the 'true' elevation of the point on the vector since (unless the profile runs planview N-S or W-E) the vector does not always run exactly through the centroid of the cell
      int nRet = nGetThisProfileElevationsForCShore(nCoast, pProfile, nProfileSize, &VdProfileDistXY, &VdProfileZ, &VdProfileFrictionFactor);
      if (nRet != RTN_OK)
      {
         // Could not create the profile elevation vectors
         LogStream << m_ulIter << ":\t coast " << nCoast << " could not create CShore profile elevation vectors for profile " << pProfile->nGetProfileID() << endl;

         return nRet;
      }

      assert(static_cast<int>(VdProfileDistXY.size()) == nProfileSize);
      assert(static_cast<int>(VdProfileZ.size()) == nProfileSize);
      assert(static_cast<int>(VdProfileFrictionFactor.size()) == nProfileSize);

      if (VdProfileDistXY.empty())
      {
         // The profile elevation vector was created, but was not populated
         LogStream << m_ulIter << ":\t coast " << nCoast << " could not populate CShore profile elevation vector for profile " << pProfile->nGetProfileID() << endl;

         return RTN_ERR_CSHORE_EMPTY_PROFILE;
      }

      // Constrain the wave to normal angle to be between -80 and 80 degrees, this is a requirement of CShore
      dWaveToNormalAngle = tMax(dWaveToNormalAngle, -80.0);
      dWaveToNormalAngle = tMin(dWaveToNormalAngle, 80.0);

      int const nProfileDistXYSize = static_cast<int>(VdProfileDistXY.size());
      vector<double> VdFreeSurfaceStd(nProfileDistXYSize, 0);        // This is converted to Hrms by Hrms = sqr(8)*FreeSurfaceStd
      vector<double> VdSinWaveAngleRadians(nProfileDistXYSize, 0);   // This is converted to deg by asin(VdSinWaveAngleRadians)*(180/pi)
      vector<double> VdFractionBreakingWaves(nProfileDistXYSize, 0); // Is 0 if no wave breaking, and 1 if all waves breaking

      // Now define the other values that CShore requires
      int const nILine = 1;  // This is the number of cross-shore lines i.e. the number of CoastalME profiles. Only one at a time, at present
      int const nIProfl = 0; // 0 for fixed bottom profile, 1 for profile evolution computation
      int const nIPerm = 0;  // 0 for impermeable bottom, 1 for permeable bottom of stone structure
      int const nIOver = 0;  // 0 for no wave overtopping and overflow on crest, 1 for wave overtopping and overflow
      int const nIWCInt = 0; // 0 for no wave/current interaction, 1 for wave/current interaction in frequency dispersion, momentum and wave action equations
      int const nIRoll = 0;  // 0 for no roller effects, 1 for roller effects in governing equations
      int const nIWind = 0;  // 0 for no wind effects, 1 for wind shear stresses on momentum equations
      int const nITide = 0;  // 0 for no tidal effect on currents, 1 for longshore and cross-shore tidal currents
      int const nILab = 0;   // 0 for field data set, 1 for laboratory data set
      int const nNWave = 1;  // Number of waves at x = 0 starting from time = 0
      int const nNSurge = 1; // Number of water levels at x = 0 from time = 0
      double const dDX = VdProfileDistXY.back() / (static_cast<double>(VdProfileDistXY.size() - 1));
      double const dWaveInitTime = 0;  // CShore wave start time
      double const dSurgeInitTime = 0; // CShore surge start time

#if defined CSHORE_FILE_INOUT || CSHORE_BOTH
      // Move to the CShore folder
      nRet = chdir(CSHORE_DIR.c_str());

      if (nRet != RTN_OK)
         return nRet;

#endif

#if defined CSHORE_FILE_INOUT
      // We are communicating with CShore using ASCII files, so create an input file for this profile which will be read by CShore
      nRet = nCreateCShoreInfile(nCoast, nProfile, nILine, nIProfl, nIPerm, nIOver, nIWCInt, nIRoll, nIWind, nITide, nILab, nNWave, nNSurge, dDX, dCShoreTimeStep, dWaveInitTime, dDeepWaterWavePeriod, dProfileDeepWaterWaveHeight, dWaveToNormalAngle, dSurgeInitTime, dSurgeLevel, &VdProfileDistXY, &VdProfileZ, &VdProfileFrictionFactor);

      if (nRet != RTN_OK)
         return nRet;

      // Set the error flag: this will be changed to 0 within CShore if CShore returns correctly
      nRet = -1;

      // Run CShore for this profile
      CShore(&nRet);
      // CShore(&nRet, &m_ulIter, &nProfile, &nProfile);

      // Check return code for error
      if (nRet != 0)
      {
         string strErr = to_string(m_Uliter) + ":\t ";

         switch (nRet)
         {
            case -1:
               strErr = "CShore WARNING 1 negative depth at the first node, set to zero";
               break;

            case 2:
               strErr = "CShore WARNING 2 negative value at end of landward marching computation, set to 1.0D-9";
               break;

            case 3:
               strErr = "CShore WARNING 3 large energy gradients at the first node: small waves with short period at sea boundary";
               break;

            case 4:
               strErr = "CShore WARNING 4 zero energy at the first node";
               break;

            case 5:
               strErr = "CShore WARNING 5 at end of landward marching computation, insufficient water depth";
               break;

            case 7:
               strErr = "CShore WARNING 7 did not reach convergence";
               break;

            default:
               strErr = "CShore UNKOWN WARNING";
         }

         strErr += " on coast " + to_string(nCoast) + " profile " + to_string(pProfile->nGetProfileID()) + " profile length " + to_string(nOutSize);

         // OK, give up for this profile
         LogStream << strErr << endl;

         // return RTN_ERR_CSHORE_ERROR;
      }

      // Fetch the CShore results by reading files written by CShore
      string strOSETUP = "OSETUP";
      string strOYVELO = "OYVELO";
      string strOPARAM = "OPARAM";

      nRet = nReadCShoreOutput(nProfile, &strOSETUP, 4, 4, &VdProfileDistXY, &VdFreeSurfaceStd);

      if (nRet != RTN_OK)
         return nRet;

      nRet = nReadCShoreOutput(nProfile, &strOYVELO, 4, 2, &VdProfileDistXY, &VdSinWaveAngleRadians);

      if (nRet != RTN_OK)
         return nRet;

      nRet = nReadCShoreOutput(nProfile, &strOPARAM, 4, 3, &VdProfileDistXY, &VdFractionBreakingWaves);

      if (nRet != RTN_OK)
         return nRet;

      // Read surge outputs
      // VdTSurg = {dSurgeInitTime, dCShoreTimeStep},                           // Ditto
      // VdSWLin = {dSurgeLevel, dSurgeLevel},                                  // Ditto

      // Clean up the CShore outputs
#ifdef _WIN32
      nRet = system("./clean.bat");
#else

      if (SAVE_CSHORE_OUTPUT)
      {
         string strCommand = "./save_CShore_output.sh ";
         strCommand += to_string(m_ulIter);
         strCommand += " ";
         strCommand += to_string(nCoast);
         strCommand += " ";
         strCommand += to_string(nProfile);

         nRet = system(strCommand.c_str());

         if (nRet != RTN_OK)
            return nRet;
      }

      nRet = system("./clean.sh");
#endif

      if (nRet != RTN_OK)
         return nRet;

      // And return to the CoastalME folder
      nRet = chdir(m_strCMEDir.c_str());

      if (nRet != RTN_OK)
         return nRet;

#endif

#if defined CSHORE_ARG_INOUT || CSHORE_BOTH
      // We are communicating with CShore by passing arguments
      int nOutSize = 0; // CShore will return the size of the output vectors

      // Set the error flag: this will be changed within CShore if there is a problem
      nRet = 0;

      vector<double> VdInitTime = {dWaveInitTime, dCShoreTimeStep};                          // Size is nNwave+1, value 1 is for the start of the CShore run, value 2 for end of CShore run
      vector<double> VdTPIn = {dDeepWaterWavePeriod, dDeepWaterWavePeriod};                  // Ditto
      vector<double> VdHrmsIn = {dProfileDeepWaterWaveHeight, dProfileDeepWaterWaveHeight};  // Ditto
      vector<double> VdWangIn = {dWaveToNormalAngle, dWaveToNormalAngle};                    // Ditto
      vector<double> VdTSurg = {dSurgeInitTime, dCShoreTimeStep};                            // Ditto
      vector<double> VdSWLin = {dSurgeLevel, dSurgeLevel};                                   // Ditto
      vector<double> VdFPInp = VdProfileFrictionFactor;                                      // Set the value for wave friction at every point of the normal profile
      vector<double> VdXYDistFromCShoreOut(CSHORE_ARRAY_OUT_SIZE, 0);                        // Output from CShore
      vector<double> VdFreeSurfaceStdOut(CSHORE_ARRAY_OUT_SIZE, 0);                          // Ditto
      vector<double> VdWaveSetupSurgeOut(CSHORE_ARRAY_OUT_SIZE, 0);                          // Ditto
      // vector<double> VdStormSurgeOut(CSHORE_ARRAY_OUT_SIZE, 0);                           // Ditto
      vector<double> const VdWaveSetupRunUpOut(CSHORE_ARRAY_OUT_SIZE, 0);                    // Ditto
      vector<double> VdSinWaveAngleRadiansOut(CSHORE_ARRAY_OUT_SIZE, 0);                     // Ditto
      vector<double> VdFractionBreakingWavesOut(CSHORE_ARRAY_OUT_SIZE, 0);                   // Ditto

      // Call CShore using the argument-passing wrapper
      // long lIter = static_cast<long>(m_ulIter);    // Bodge to get round compiler 'invalid conversion' error

      CShoreWrapper(&nILine,                          /* In_ILINE */
                    &nIProfl,                         /* In_IPROFL */
                    &nIPerm,                          /* In_IPERM */
                    &nIOver,                          /* In_IOVER */
                    &nIWCInt,                         /* In_IWCINT */
                    &nIRoll,                          /* In_IROLL */
                    &nIWind,                          /* In_IWIND */
                    &nITide,                          /* In_ITIDE */
                    &nILab,                           /* In_ILAB */
                    &nNWave,                          /* In_NWAVE */
                    &nNSurge,                         /* In_NSURG */
                    &dDX,                             /* In_DX */
                    &m_dBreakingWaveHeightDepthRatio, /* In_GAMMA */
                    &VdInitTime[0],                   /* In_TWAVE */
                    &VdTPIn[0],                       /* In_TPIN */
                    &VdHrmsIn[0],                     /* In_HRMSIN */
                    &VdWangIn[0],                     /* In_WANGIN */
                    &VdTSurg[0],                      /* In_TSURG */
                    &VdSWLin[0],                      /* In_SWLIN */
                    &nProfileDistXYSize,              /* In_NBINP */
                    &VdProfileDistXY[0],              /* In_XBINP */
                    &VdProfileZ[0],                   /* In_ZBINP */
                    &VdFPInp[0],                      /* In_FBINP */
                    &nRet,                            /* Out_IError */
                    &nOutSize,                        /* Out_nOutSize */
                    &VdXYDistFromCShoreOut[0],        /* Out_XYDist */
                    &VdFreeSurfaceStdOut[0],          /* Out_FreeSurfaceStd */
                    &VdWaveSetupSurgeOut[0],          /* Out_WaveSetupSurge */
                    &VdSinWaveAngleRadiansOut[0],     /* Out_SinWaveAngleRadians */
                    &VdFractionBreakingWavesOut[0]);  /* Out_FractionBreakingWaves */

      // OK, now check for warnings and errors
      if (nOutSize < 2)
      {
         // CShore sometimes returns only one row of results, which contains data only for the seaward point of the profile. This happens when all other (more coastward) points give an invalid result during CShore's calculations. This is a problem. We don't want to abandon the simulation just because of this, so instead we just put some dummy data into the second row, and carry on with these two rows. The profile will get ignored later, since it is too small to be useful
         if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
            LogStream << m_ulIter << ":\t" << WARN << "for coast " << nCoast << " profile " << pProfile->nGetProfileID() << ", only " << nOutSize << " CShore output rows, abandoning this profile" << endl;

         // // Set dummy data in the second row
         // VdXYDistFromCShoreOut[1] = 1e-5;    // Dummy data, must not be the same as VdXYDistFromCShoreOut[0] tho', or get crash in linear interpolation routine
         // VdFreeSurfaceStdOut[1] = VdFreeSurfaceStdOut[0];
         // VdSinWaveAngleRadiansOut[1] = VdSinWaveAngleRadiansOut[0];
         // VdFractionBreakingWavesOut[1] = VdFractionBreakingWavesOut[0];
         // VdWaveSetupSurgeOut[1] = VdWaveSetupSurgeOut[0];
         // // VdStormSurgeOut[1] = VdStormSurgeOut[0];
         // // VdWaveSetupRunUpOut[1] = VdWaveSetupRunUpOut[0];
         //
         // // And increase the expected number of rows
         // nOutSize = 2;

         return RTN_ERR_CSHORE_ERROR;
      }

      if (nRet != 0)
      {
         string strErr = to_string(m_ulIter) + ":\t ";

         switch (nRet)
         {
            case -1:
               strErr += "CShore WARNING 1 negative depth at the first node, set to zero";
               break;

            case 2:
               strErr += "CShore WARNING 2 negative value at end of landward marching computation, set to 1.0D-9";
               break;

            case 3:
               strErr += "CShore WARNING 3 large energy gradients at the first node: small waves with short period at sea boundary";
               break;

            case 4:
               strErr += "CShore WARNING 4 zero energy at the first node";
               break;

            case 5:
               strErr += "CShore WARNING 5 at end of landward marching computation, insufficient water depth";
               break;

            case 7:
               strErr += "CShore WARNING 7 did not reach convergence";
               break;
         }

         strErr += " on coast " + to_string(nCoast) + " profile " + to_string(pProfile->nGetProfileID()) + " profile length " + to_string(nOutSize);

         // OK, give up for this profile
         LogStream << strErr << endl;

         // TEST
         // return RTN_ERR_CSHORE_ERROR;
      }

      // LogStream << m_ulIter << ": interpolating profile " << nProfile << endl;

      // All OK, so interpolate the CShore output, and convert from the CShore convention (cross-shore distance has its origin at the seaward end) to the CoastalME convention (origin at the shoreline)
      InterpolateCShoreOutput(&VdProfileDistXY, nOutSize, nProfileSize, &VdXYDistFromCShoreOut, &VdFreeSurfaceStdOut, &VdWaveSetupSurgeOut, &VdSinWaveAngleRadiansOut, &VdFractionBreakingWavesOut, &VdFreeSurfaceStd, &VdWaveSetupSurge, &VdSinWaveAngleRadians, &VdFractionBreakingWaves);
#endif

#if defined CSHORE_BOTH
#if !defined _WIN32

      if (SAVE_CSHORE_OUTPUT)
      {
         string strCommand = "./save_CShore_output.sh ";
         strCommand += to_string(m_ulIter);
         strCommand += " ";
         strCommand += to_string(nCoast);
         strCommand += " ";
         strCommand += to_string(nProfile);

         nRet = system(strCommand.c_str());

         if (nRet != RTN_OK)
            return nRet;
      }

#endif

      // Return to the CoastalME folder
      nRet = chdir(m_strCMEDir.c_str());

      if (nRet != RTN_OK)
         return nRet;

#endif

      // Do more safety checks, then convert the CShore output to wave height and wave direction, and update wave profile attributes
      for (int nProfilePoint = (nProfileSize - 1); nProfilePoint >= 0; nProfilePoint--)
      {
         int const nX = pProfile->pPtiGetCellInProfile(nProfilePoint)->nGetX();
         int const nY = pProfile->pPtiGetCellInProfile(nProfilePoint)->nGetY();

         // Safety check
         if (nProfilePoint > static_cast<int>(VdFreeSurfaceStd.size()) - 1)
            continue;

         // Safety check: deal with NaN values
         if (isnan(VdFreeSurfaceStd[nProfilePoint]))
            VdFreeSurfaceStd[nProfilePoint] = 0;

         VdWaveHeight[nProfilePoint] = sqrt(8) * VdFreeSurfaceStd[nProfilePoint];

         // Another safety check: deal with NaN values
         if (isnan(VdSinWaveAngleRadians[nProfilePoint]))
         {
            VdSinWaveAngleRadians[nProfilePoint] = 0;
            VdWaveHeight[nProfilePoint] = 0;
         }

         // More safety checks: constrain to the interval -1 to +1 to keep asin() happy
         if (VdSinWaveAngleRadians[nProfilePoint] < -1)
            VdSinWaveAngleRadians[nProfilePoint] = -1;

         if (VdSinWaveAngleRadians[nProfilePoint] > 1)
            VdSinWaveAngleRadians[nProfilePoint] = 1;

         double const dAlpha = asin(VdSinWaveAngleRadians[nProfilePoint]) * (180 / PI);

         if (nSeaHand == LEFT_HANDED)
            VdWaveDirection[nProfilePoint] = dKeepWithin360(dAlpha + 90 + dFluxOrientationThis);
         else
            VdWaveDirection[nProfilePoint] = dKeepWithin360(dAlpha + 270 + dFluxOrientationThis);

         // Yet another safety check: deal with NaN values
         if (isnan(VdFractionBreakingWaves[nProfilePoint]))
         {
            VdFractionBreakingWaves[nProfilePoint] = 0;
            VdWaveHeight[nProfilePoint] = 0;
         }

         // if ((VdFractionBreakingWaves[nProfilePoint] >= 0.10) && (! bBreaking)) // Sometimes is possible that waves break again
         if ((VdFractionBreakingWaves[nProfilePoint] >= 0.10) && (m_dDepthOfClosure >= m_pRasterGrid->m_Cell[nX][nY].dGetSeaDepth()) && (! bBreaking))
         {
            bBreaking = true;
            assert(VdWaveHeight[nProfilePoint] >= 0);
            dProfileBreakingWaveHeight = VdWaveHeight[nProfilePoint];
            dProfileBreakingWaveAngle = VdWaveDirection[nProfilePoint];
            dProfileBreakingDepth = m_pRasterGrid->m_Cell[nX][nY].dGetSeaDepth(); // Water depth for the cell 'under' this point in the profile
            nProfileBreakingDist = nProfilePoint + 1;                             // At the nearest point nProfilePoint = 0, so, plus one

            // LogStream << m_ulIter << ":\t coast " << nCoast << " CShore breaking at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "} nProfile = " << nProfile << ", nProfilePoint = " << nProfilePoint << ", dBreakingWaveHeight = " << dBreakingWaveHeight << ", dBreakingWaveAngle = " << dBreakingWaveAngle << ", dProfileBreakingDepth = " << dProfileBreakingDepth << ", nProfileBreakingDist = " << nProfileBreakingDist << endl;
         }

         VbWaveIsBreaking[nProfilePoint] = bBreaking;
      }

      // Safety check
      if (dProfileBreakingWaveHeight >= dProfileDeepWaterWaveHeight)
         dProfileBreakingWaveHeight = dProfileDeepWaterWaveHeight;
   }
   else if (m_nWavePropagationModel == WAVE_MODEL_COVE)
   {
      // We are using COVE's linear wave theory to propagate the waves
      double const dDepthLookupMax = m_dWaveDepthRatioForWaveCalcs * dProfileDeepWaterWaveHeight;

      // Go landwards along the profile, calculating wave height and wave angle for every inundated point on the profile (don't do point zero, this is on the coastline) until the waves start to break  after breaking wave height is assumed to decrease linearly to zero at the shoreline and wave angle is equalt to wave angle at breaking
      for (int nProfilePoint = (nProfileSize - 1); nProfilePoint >= 0; nProfilePoint--)
      {
         int const nX = pProfile->pPtiGetCellInProfile(nProfilePoint)->nGetX();
         int const nY = pProfile->pPtiGetCellInProfile(nProfilePoint)->nGetY();

         // Safety check
         if (! m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea())
            continue;

         double const dSeaDepth = m_pRasterGrid->m_Cell[nX][nY].dGetSeaDepth(); // Water depth for the cell 'under' this point in the profile

         if (dSeaDepth > dDepthLookupMax)
         {
            // Sea depth is too large relative to wave height to feel the bottom, so do nothing since each cell under the profile has already been initialised with deep water wave height and wave direction
            dProfileWaveHeight = dProfileDeepWaterWaveHeight;
            dProfileWaveAngle = dProfileDeepWaterWaveAngle;
         }
         else
         {
            if (! bBreaking)
            {
               // Start calculating wave properties using linear wave theory
               double const dL = m_dL_0 * sqrt(tanh((2 * PI * dSeaDepth) / m_dL_0));                          // Wavelength (m) in intermediate-shallow waters
               double const dC = m_dC_0 * tanh((2 * PI * dSeaDepth) / dL);                                    // Wave speed (m/s) set by dSeaDepth, dL and m_dC_0
               double const dk = 2 * PI / dL;                                                                 // Wave number (1/m)
               double const dn = ((2 * dSeaDepth * dk) / (sinh(2 * dSeaDepth * dk)) + 1) / 2;                 // Shoaling factor
               double const dKs = sqrt(m_dC_0 / (dn * dC * 2));                                               // Shoaling coefficient
               double const dAlpha = (180 / PI) * asin((dC / m_dC_0) * sin((PI / 180) * dWaveToNormalAngle)); // Calculate angle between wave direction and the normal to the coast tangent
               double const dKr = sqrt(cos((PI / 180) * dWaveToNormalAngle) / cos((PI / 180) * dAlpha));      // Refraction coefficient
               dProfileWaveHeight = dProfileDeepWaterWaveHeight * dKs * dKr;                                  // Calculate wave height, based on the previous (more seaward) wave height

               if (nSeaHand == LEFT_HANDED)
                  dProfileWaveAngle = dKeepWithin360(dAlpha + 90 + dFluxOrientationThis);
               else
                  dProfileWaveAngle = dKeepWithin360(dAlpha + 270 + dFluxOrientationThis);

               // Test to see if the wave breaks at this depth
               if (dProfileWaveHeight < (dSeaDepth * m_dBreakingWaveHeightDepthRatio))
               {
                  dProfileBreakingWaveHeight = dProfileWaveHeight;
                  dProfileBreakingWaveAngle = dProfileWaveAngle;
               }
            }
            else
            {
               // It does
               bBreaking = true;

               // Wave has already broken
               dProfileWaveAngle = dProfileBreakingWaveAngle; // Wave orientation remains equal to wave orientation at breaking

               dProfileWaveHeight = dProfileBreakingWaveHeight * (nProfilePoint / nProfileBreakingDist); // Wave height decreases linearly to zero at shoreline
               // dProfileWaveHeight = dSeaDepth * m_dBreakingWaveHeightDepthRatio; // Wave height is limited by depth
               dProfileBreakingDepth = dSeaDepth;
               nProfileBreakingDist = nProfilePoint;
            }
         }

         // Save current wave attributes
         VdWaveDirection[nProfilePoint] = dProfileWaveAngle;
         VdWaveHeight[nProfilePoint] = dProfileWaveHeight;
         VbWaveIsBreaking[nProfilePoint] = bBreaking;
      }
   }

   // Go landwards along the profile, fetching the calculated wave height and wave angle for every inundated point on this profile
   for (int nProfilePoint = (nProfileSize - 1); nProfilePoint >= 0; nProfilePoint--)
   {
      int nX = pProfile->pPtiGetCellInProfile(nProfilePoint)->nGetX();
      int nY = pProfile->pPtiGetCellInProfile(nProfilePoint)->nGetY();

      // Safety check
      if (! m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea())
         continue;

      // Get the wave attributes calculated for this profile: wave height, wave angle, and whether is in the active zone
      double const dWaveHeight = VdWaveHeight[nProfilePoint];
      double const dWaveAngle = VdWaveDirection[nProfilePoint];

      bBreaking = VbWaveIsBreaking[nProfilePoint];

      // And store the wave properties for this point in the all-profiles vectors
      pVdX->push_back(nX);
      pVdY->push_back(nY);
      pVdHeightX->push_back(dWaveHeight * sin(dWaveAngle * PI / 180));
      pVdHeightY->push_back(dWaveHeight * cos(dWaveAngle * PI / 180));
      pVbBreaking->push_back(bBreaking);
   }

   // Obtain the profile nodes near the coast
   int const nX = pProfile->pPtiGetCellInProfile(nProfileSize - 2)->nGetX();
   int const nY = pProfile->pPtiGetCellInProfile(nProfileSize - 2)->nGetY();
   int const nX1 = pProfile->pPtiGetCellInProfile(nProfileSize - 1)->nGetX();
   int const nY1 = pProfile->pPtiGetCellInProfile(nProfileSize - 1)->nGetY();

   // Calculate the horizontal distance between the profile points
   double dXDist = tAbs(dGridCentroidXToExtCRSX(nX1) - dGridCentroidXToExtCRSX(nX));
   double dYDist = tAbs(dGridCentroidYToExtCRSY(nY1) - dGridCentroidYToExtCRSY(nY));

   // Safety check
   if (bFPIsEqual(dXDist, 0.0, TOLERANCE))
      dXDist = 1e-3;

   // Safety check
   if (bFPIsEqual(dYDist, 0.0, TOLERANCE))
      dYDist = 1e-3;

   double const dDiffProfileDistXY = hypot(dXDist, dYDist);

   // Compute the beach slope
   double const dtanBeta = tan(tAbs(m_pRasterGrid->m_Cell[nX1][nY1].dGetSeaDepth() - m_pRasterGrid->m_Cell[nX][nY].dGetSeaDepth()) / dDiffProfileDistXY);

   // Compute the wave run-up using NIELSEN & HANSLOW (1991) & DHI (2004)
   int nValidPointsWaveHeight = 0;
   int nValidPointsWaveSetup = 0;

   for (int nPoint = 0; nPoint < static_cast<int>(VdWaveHeight.size()); nPoint++)
   {
      if (VdWaveHeight[nPoint] > 1e-4)
         nValidPointsWaveHeight += 1;
      else
         break;
   }

   nValidPointsWaveHeight -= 1;

   for (int nPoint = 0; nPoint < static_cast<int>(VdWaveSetupSurge.size()); nPoint++)
   {
      if (tAbs(VdWaveSetupSurge[nPoint]) < 1) // limiting the absolute value of setup + surge if cshore run fails
         nValidPointsWaveSetup += 1;
      else
         break;
   }

   nValidPointsWaveSetup -= 1;

   double dWaveHeight = 0;

   // Safety checks
   if ((nValidPointsWaveHeight >= 0) && (! bFPIsEqual(VdWaveHeight[nValidPointsWaveHeight], DBL_NODATA, TOLERANCE)))
      dWaveHeight = VdWaveHeight[nValidPointsWaveHeight];

   // TODO 060 Remove the following 'magic numbers'
   double dRunUp = 0;

   if (m_nRunUpEquation == RUNUP_EQUATION_NIELSEN_HANSLOW)
   {
      // Compute the run-up using Nielsen, P. & Hanslow, D. J. 1991. Wave Runup Distributions on Natural Beaches. Journal of Coastal Research, 7, 1139-1152.  *** & DHI (2004) ???
      dRunUp = 0.36 * pow(9.81, 0.5) * dtanBeta * pow(dWaveHeight, 0.5) * dDeepWaterWavePeriod;
   }
   else if (m_nRunUpEquation == RUNUP_EQUATION_MASE)
   {
      // Compute the run-up using Mase, H. 1989. Random Wave Runup Height on Gentle Slope. Journal of Waterway, Port, Coastal, and Ocean Engineering, 115, 649-661.
      double const dS0 = 2 * PI * dWaveHeight / (9.81 * dDeepWaterWavePeriod * dDeepWaterWavePeriod);
      dRunUp = 1.86 * dWaveHeight * pow(pow(dtanBeta / dS0, 0.5), 0.71);
   }
   else if (m_nRunUpEquation == RUNUP_EQUATION_STOCKDON)
   {
      // Compute the run-up using Stockdon, H. F., Holman, R. A., Howd, P. A. & Sallenger JR, A. H. 2006. Empirical parameterization of setup, swash, and runup. Coastal Engineering, 53, 573-588. This version coded by Tom Ashby
      double const dS0 = 2 * PI * dWaveHeight / (9.81 * dDeepWaterWavePeriod * dDeepWaterWavePeriod);
      dRunUp = 1.1 * ((0.35 * dtanBeta * pow((1.0 / dS0) * dWaveHeight * dWaveHeight, 0.5)) + (0.5 * pow((1.0 / dS0) * dWaveHeight * dWaveHeight * (0.563 * dtanBeta * dtanBeta + 0.004), 0.5)));

      // double const dS0 = 2 * PI * dWaveHeight / (9.81 * dDeepWaterWavePeriod * dDeepWaterWavePeriod);
      // dRunUp = 1.1 * ((0.35 * dWaveHeight * (pow((1 / dS0) * dWaveHeight * dWaveHeight, 0.5))) + (((((1 / dS0) * dWaveHeight * dWaveHeight) * (0.563 * dWaveHeight * dWaveHeight + 0.0004)), 0.5)) / 2);
      // double const dH0OverL0 = (1 / dS0) * dWaveHeight;
      // double const dTmp1 = 0.35 * dWaveHeight * pow(dH0OverL0, 0.5);
      // double const dTmp2 = pow(dH0OverL0 * ((0.563 * dWaveHeight * dWaveHeight) + 0.0004), 0.5);
      // dRunUp = 1.1 * (dTmp1 + (dTmp2 / 2));
   }

   if ((tAbs(dRunUp) < 1e-4) || (isnan(dRunUp)))
      dRunUp = 0;

   double dWaveSetupSurge = 0;

   // Safety checks
   if ((nValidPointsWaveSetup >= 0) && (! bFPIsEqual(VdWaveSetupSurge[nValidPointsWaveSetup], DBL_NODATA, TOLERANCE)))
      dWaveSetupSurge = VdWaveSetupSurge[nValidPointsWaveSetup];

   if ((isnan(dWaveSetupSurge) || (tAbs(dWaveSetupSurge) < 1e-4)))
      dWaveSetupSurge = 0;

   // Update wave attributes along the coastline object. Wave height at the coast is always calculated (i.e. whether or not waves are breaking)
   // cout << "Wave Height at the coast is " << VdWaveHeight[nProfileSize - 1] << endl;
   m_VCoast[nCoast].SetCoastWaveHeight(nCoastPoint, dWaveHeight);
   m_VCoast[nCoast].SetWaveSetupSurge(nCoastPoint, dWaveSetupSurge);
   m_VCoast[nCoast].SetRunUp(nCoastPoint, dRunUp);

   if (nProfileBreakingDist > 0)
   {
      // This coast point is in the active zone, so set breaking wave height, breaking wave angle, and depth of breaking for the coast point
      m_VCoast[nCoast].SetBreakingWaveHeight(nCoastPoint, dProfileBreakingWaveHeight);
      m_VCoast[nCoast].SetBreakingWaveAngle(nCoastPoint, dProfileBreakingWaveAngle);
      m_VCoast[nCoast].SetDepthOfBreaking(nCoastPoint, dProfileBreakingDepth);
      m_VCoast[nCoast].SetBreakingDistance(nCoastPoint, nProfileBreakingDist);

      // LogStream << m_ulIter << ": nProfile = " << nProfile << ", nCoastPoint = " << nCoastPoint << " in active zone, dBreakingWaveHeight = " << dBreakingWaveHeight << endl;
   }
   else
   {
      // This coast point is not in the active zone, so breaking wave height, breaking wave angle, and depth of breaking are all meaningless
      m_VCoast[nCoast].SetBreakingWaveHeight(nCoastPoint, DBL_NODATA);
      m_VCoast[nCoast].SetBreakingWaveAngle(nCoastPoint, DBL_NODATA);
      m_VCoast[nCoast].SetDepthOfBreaking(nCoastPoint, DBL_NODATA);
      m_VCoast[nCoast].SetBreakingDistance(nCoastPoint, INT_NODATA);

      // LogStream << m_ulIter << ": nProfile = " << nProfile << ", nCoastPoint = " << nCoastPoint << " NOT in active zone" << endl;
   }

   return RTN_OK;
}

#if defined CSHORE_FILE_INOUT
//===============================================================================================================================
//! Create and write to the CShore input file
//===============================================================================================================================
int CSimulation::nCreateCShoreInfile(int const nCoast, int const nProfile, int const nILine, int const nIProfl, int const nIPerm, int const nIOver, int const nIWcint, int const nIRoll, int const nIWind, int const nITide, int const nILab, int const nWave, int const nSurge, double const dX, double const dTimestep, double const dWaveInitTime, double const dWavePeriod, double const dHrms, double const dWaveAngle, double const dSurgeInitTime, double const dSurgeLevel, vector<double> const *pVdXdist, vector<double> const *pVdBottomElevation, vector<double> const *pVdWaveFriction)
{
   // Create the CShore input file
   ofstream CShoreOutStream;
   CShoreOutStream.open(CSHORE_INFILE.c_str(), ios::out | ios::app);

   if (CShoreOutStream.fail())
   {
      // Error, cannot open file for writing
      LogStream << m_ulIter << ": " << ERR << "cannot write to CShore input file '" << CSHORE_INFILE << "'" << endl;
      return RTN_ERR_CSHORE_FILE_INPUT;
   }

   // And write to the file
   CShoreOutStream << 3 << endl; // Number of comment lines
   CShoreOutStream << "------------------------------------------------------------" << endl;
   CShoreOutStream << "CShore input file created by CoastalME for iteration " << m_ulIter << ", coast " << nCoast << ", profile " << nProfile << endl;
   CShoreOutStream << "------------------------------------------------------------" << endl;
   CShoreOutStream << nILine << "                                         -> ILINE" << endl;
   CShoreOutStream << nIProfl << "                                         -> IPROFL" << endl;
   CShoreOutStream << nIPerm << "                                         -> IPERM" << endl;
   CShoreOutStream << nIOver << "                                         -> IOVER" << endl;
   CShoreOutStream << nIWcint << "                                         -> IWCINT" << endl;
   CShoreOutStream << nIRoll << "                                         -> IROLL" << endl;
   CShoreOutStream << nIWind << "                                         -> IWIND" << endl;
   CShoreOutStream << nITide << "                                         -> ITIDE" << endl;
   CShoreOutStream << fixed;
   CShoreOutStream << setw(11) << setprecision(4) << dX << "                               -> DX" << endl;
   CShoreOutStream << setw(11) << m_dBreakingWaveHeightDepthRatio << "                               -> GAMMA" << endl;
   CShoreOutStream << setw(11) << nILab << "                               -> ILAB" << endl;
   CShoreOutStream << setw(11) << nWave << "                               -> NWAVE" << endl;
   CShoreOutStream << setw(11) << nSurge << "                               -> NSURGE" << endl;

   // Line 18 of infile
   CShoreOutStream << setw(11) << setprecision(2) << dWaveInitTime; // TWAVE(1) in CShore
   CShoreOutStream << setw(11) << setprecision(4) << dWavePeriod;   // TPIN(1) in CShore
   CShoreOutStream << setw(11) << dHrms;                            // HRMS(1) in CShore
   CShoreOutStream << setw(11) << dWaveAngle << endl;               // WANGIN(1) in CShore

   // Line 19 of infile
   CShoreOutStream << setw(11) << setprecision(2) << dTimestep;   // TWAVE(2) in CShore
   CShoreOutStream << setw(11) << setprecision(4) << dWavePeriod; // TPIN(2) in CShore
   CShoreOutStream << setw(11) << dHrms;                          // HRMS(2) in CShore
   CShoreOutStream << setw(11) << dWaveAngle << endl;             // WANGIN(2) in CShore

   // Line 20 of infile
   CShoreOutStream << setw(11) << setprecision(2) << dSurgeInitTime;      // TSURG(1) in CShore
   CShoreOutStream << setw(11) << setprecision(4) << dSurgeLevel << endl; // SWLIN(1) in CShore

   // Line 21 of infile
   CShoreOutStream << setw(11) << setprecision(2) << dTimestep;           // TSURG(2) in CShore
   CShoreOutStream << setw(11) << setprecision(4) << dSurgeLevel << endl; // SWLIN(2) in CShore

   // Line 22 of infile
   CShoreOutStream << setw(8) << pVdXdist->size() << "                                  -> NBINP" << endl;

   CShoreOutStream << fixed << setprecision(4);

   for (unsigned int i = 0; i < pVdXdist->size(); i++)
      // These are BINP(J,1), ZBINP(J,1), FBINP(J-1,1) in CShore
      CShoreOutStream << setw(11) << pVdXdist->at(i) << setw(11) << pVdBottomElevation->at(i) << setw(11) << pVdWaveFriction->at(i) << endl;

   CShoreOutStream << endl;

   // File written, so close it
   CShoreOutStream.close();

   return RTN_OK;
}
#endif

//===============================================================================================================================
//! Get profile horizontal distance and bottom elevation vectors in CShore units
//===============================================================================================================================
int CSimulation::nGetThisProfileElevationsForCShore(int const nCoast, CGeomProfile *pProfile, int const nProfSize, vector<double> *VdDistXY, vector<double> *VdVZ, vector<double> *VdFricF)
{
   bool bIsBehindIntervention = false;

   int nX1 = 0;
   int nY1 = 0;

   double dXDist;
   double dYDist;
   double dProfileDistXY = 0;
   double dProfileFricFact;
   double dPrevDist = -1;

   for (int i = nProfSize - 1; i >= 0; i--)
   {
      int const nX = pProfile->pPtiVGetCellsInProfile()->at(i).nGetX();
      int const nY = pProfile->pPtiVGetCellsInProfile()->at(i).nGetY();

      // Calculate the horizontal distance relative to the most seaward point
      if (i == nProfSize - 1)
         dProfileDistXY = 0;
      else
      {
         dXDist = dGridCentroidXToExtCRSX(nX1) - dGridCentroidXToExtCRSX(nX),
         dYDist = dGridCentroidYToExtCRSY(nY1) - dGridCentroidYToExtCRSY(nY),
         dProfileDistXY = dProfileDistXY + hypot(dXDist, dYDist);
      }

      // Before we store the X-Y distance, must check that it is not the same as the previously-stored distance (if it is, we get zero-divide errors in CShore). If they are the same, add on a small distance
      if (bFPIsEqual(dProfileDistXY, dPrevDist, TOLERANCE))
         dProfileDistXY += 0.1; // TODO 084 Improve this

      // Update the cell indexes, the initial cell is now the previous one
      nX1 = nX;
      nY1 = nY;

      // Get the number of the highest layer with non-zero thickness
      int const nTopLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopNonZeroLayerAboveBasement();

      // Safety check
      if (nTopLayer == NO_NONZERO_THICKNESS_LAYERS)
      {
         // Down to basement
         LogStream << WARN << "down to basement during wave calcs at [" << nX << "][" << nY << "] for coast " << pProfile->nGetCoastID() << " profile " << pProfile->nGetProfileID() << endl;

         // return RTN_ERR_BASEMENT_DURING_WAVE_CALC;
         return RTN_OK;
      }

      // Get the elevation for both consolidated and unconsolidated sediment (including any talus) on this cell
      double const dTopElev = m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevIncTalus() + m_pRasterGrid->m_Cell[nX][nY].dGetInterventionHeight();
      double const VdProfileZ = dTopElev - m_dThisIterSWL;

      // Check that landward elevation is greater than SWL
      if (i == 0)
      {
         if (VdProfileZ < 0)
         {
            VdVZ->push_back(0.1); // TODO 053 Set it to a small +ve elevation (compared with SWL). However there must be a better way of doing this

            // Could not create the profile elevation vectors
            LogStream << m_ulIter << ":\t" << WARN << "coast " << nCoast << " profile " << pProfile->nGetProfileID() << " elevation at the landward end is " << dTopElev << " m. This is lower than this-iteration SWL (" << m_dThisIterSWL << " m). For CShore, changing the landward elevation for profile " << pProfile->nGetProfileID() << " to " << m_dThisIterSWL + 0.1 << "m" << endl;
         }

         else
         {
            VdVZ->push_back(VdProfileZ);
         }
      }
      else
      {
         VdVZ->push_back(VdProfileZ);
      }

      // Now store the X-Y plane distance from the start of the profile
      VdDistXY->push_back(dProfileDistXY);

      // Get the landform type at each point along the profile
      double const dInterventionHeight = m_pRasterGrid->m_Cell[nX][nY].dGetInterventionHeight();

      // Modify default friction factor if a structural intervention is found, otherwise use the default
      if (dInterventionHeight > 0 || bIsBehindIntervention)
      {
         // Use an arbitrarily high value if structure is present
         dProfileFricFact = 100 * CSHORE_FRICTION_FACTOR;
         bIsBehindIntervention = true;
      }
      else
         dProfileFricFact = CSHORE_FRICTION_FACTOR;

      // Store the friction factor
      VdFricF->push_back(dProfileFricFact);

      // For next time round
      dPrevDist = dProfileDistXY;
   }

   return RTN_OK;
}

#if defined CSHORE_FILE_INOUT
//===============================================================================================================================
//! Reads a CShore output file and creates a vector holding interpolated values
//===============================================================================================================================
int CSimulation::nReadCShoreOutput(int const nProfile, string const *strCShoreFilename, int const nExpectedColumns, int const nCShorecolumn, vector<double> const *pVdProfileDistXYCME, vector<double> *pVdInterpolatedValues)
{
   // Read in the first column (contains XY distance relative to seaward limit) and CShore column from the CShore output file
   ifstream InStream;
   InStream.open(strCShoreFilename->c_str(), ios::in);

   // Did it open OK?
   if (!InStream.is_open())
   {
      // Error: cannot open CShore file for input
      LogStream << m_ulIter << ": " << ERR << "for profile " << nProfile << ", cannot open " << *strCShoreFilename << " for input" << endl;

      return RTN_ERR_READING_CSHORE_FILE_OUTPUT;
   }

   // Opened OK, so set up the vectors to hold the CShore output data
   vector<double> VdXYDistCShore;
   vector<double> VdValuesCShore;

   // And read in the data
   int n = -1;
   int nExpectedRows = 0;
   string strLineIn;

   while (getline(InStream, strLineIn))
   {
      n++;

      if (n == 0)
      {
         // Read in the header line
         vector<string> VstrItems = VstrSplit(&strLineIn, SPACE);

         if (! bIsStringValidInt(VstrItems[1]))
         {
            string strErr = ERR + "invalid integer for number of expected rows '" + VstrItems[1] + "' in " + *strCShoreFilename + "\n";
            cerr << strErr;
            LogStream << strErr;

            return RTN_ERR_READING_CSHORE_FILE_OUTPUT;
         }

         // Get the number of expected rows
         nExpectedRows = stoi(VstrItems[1]);
      }

      else
      {
         // Read in a data line
         vector<string> VstrItems = VstrSplit(&strLineIn, SPACE);

         int nCols = static_cast<int>(VstrItems.size());

         if (nCols != nExpectedColumns)
         {
            // Error: did not read the expected number of CShore output columns
            LogStream << m_ulIter << ": " << ERR << "for profile " << nProfile << ", expected " << nExpectedColumns << " CShore output columns but read " << nCols << " columns from header section of file " << *strCShoreFilename << endl;

            return RTN_ERR_READING_CSHORE_FILE_OUTPUT;
         }

         // Number of columns is OK
         VdXYDistCShore.push_back(strtod(VstrItems[0].c_str(), NULL));
         VdValuesCShore.push_back(strtod(VstrItems[nCShorecolumn - 1].c_str(), NULL));
      }
   }

   // Check that we have read nExpectedRows from the file
   int nReadRows = static_cast<int>(VdXYDistCShore.size());

   if (nReadRows != nExpectedRows)
   {
      // Error: did not get nExpectedRows CShore output rows
      LogStream << m_ulIter << ": " << ERR << "for profile " << nProfile << ", expected " << nExpectedRows << " CShore output rows, but read " << nReadRows << " rows from file " << *strCShoreFilename << endl;

      return RTN_ERR_READING_CSHORE_FILE_OUTPUT;
   }

   if (nReadRows < 2)
   {
      // CShore sometimes returns only one row, which contains data for the seaward point of the profile. This happens when all other (more coastward) points give an invalid result during CShore's calculations. This is a problem. We don't want to abandon the simulation just because of this, so instead we just duplicate the row, so that the profile will later get marked as invalid
      if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
         LogStream << m_ulIter << ":\t" << WARN << "for profile " << nProfile << ", only " << nReadRows << " CShore output rows in file " << *strCShoreFilename << endl;

      // Duplicate the data
      VdXYDistCShore.push_back(VdXYDistCShore[0]);
      VdValuesCShore.push_back(VdValuesCShore[0]);

      // And increase the expected number of rows
      nReadRows++;
   }

   // The output is OK, so change the origin of the across-shore distance from the CShore convention to the one used here (i.e. with the origin at the shoreline)
   vector<double> VdXYDistCShoreTmp(nReadRows, 0);

   for (int i = 0; i < nReadRows; i++)
      VdXYDistCShoreTmp[i] = VdXYDistCShore[nReadRows - 1] - VdXYDistCShore[i];

   // Reverse the XY-distance and value vectors (i.e. first point is at the shoreline and must be in strictly ascending order)
   reverse(VdXYDistCShoreTmp.begin(), VdXYDistCShoreTmp.end());

   // Similarly, reverse the CShore output
   reverse(VdValuesCShore.begin(), VdValuesCShore.end());

   // Using a simple linear interpolation approach
   vector<double> VdDistXYCopy(pVdProfileDistXYCME->begin(), pVdProfileDistXYCME->end());

   assertVdXYDistCShoreTmp.size() == VdValuesCShore.size());
   *pVdInterpolatedValues = VdInterpolateCShoreProfileOutput(&VdXYDistCShoreTmp, &VdValuesCShore, &VdDistXYCopy);

   return RTN_OK;
}
#endif

#if defined CSHORE_ARG_INOUT || CSHORE_BOTH
//===============================================================================================================================
//! Interpolates CShore output, and converts from the CShore convention (cross-shore distance has its origin at the seaward end) to the CoastalME convention (origin at the shoreline)
//===============================================================================================================================
void CSimulation::InterpolateCShoreOutput(vector<double> const *pVdProfileDistXYCME, int const nOutSize, int const nProfileSize, vector<double> *pVdXYDistFromCShore, vector<double> *pVdFreeSurfaceStdCShore, vector<double> *pVdWaveSetupSurgeCShore, vector<double> *pVdSinWaveAngleRadiansCShore, vector<double> *pVdFractionBreakingWavesCShore, vector<double> *pVdFreeSurfaceStdCME, vector<double> *pVdWaveSetupSurgeCME, vector<double> *pVdSinWaveAngleRadiansCME, vector<double> *pVdFractionBreakingWavesCME)
{
   // // DEBUG CODE ========================================================================================================
   // LogStream << m_ulIter << ": nOutSize = " << nOutSize << " nProfileSize = " << nProfileSize << " pVdProfileDistXYCME->size() = " << pVdProfileDistXYCME->size() << " pVdXYDistFromCShore->size() = " << pVdXYDistFromCShore->size() << endl;
   // // DEBUG CODE ========================================================================================================

   // Sometimes the profile length returned from CShore is shorter than the CoastalME profile length
   bool bTooShort = false;
   int nTooShort = 0;

   if (nOutSize < nProfileSize)
   {
      bTooShort = true;
      nTooShort = nProfileSize - nOutSize - 1;

      // LogStream << m_ulIter << ": CShore PROFILE IS TOO SHORT" << endl;
   }

   // // DEBUG CODE ========================================================================================================
   // for (int n = 0; n < static_cast<int>(pVdProfileDistXYCME->size()); n++)
   // LogStream << "pVdProfileDistXYCME[" << n << "] = " << pVdProfileDistXYCME->at(n) << endl;
   //
   // LogStream << endl;
   //
   // LogStream << "ORIGINAL" << endl;
   // for (int n = 0; n < nOutSize; n++)
   // LogStream << "pVdXYDistFromCShore[" << n << "] = " << pVdXYDistFromCShore->at(n) << " pVdFreeSurfaceStdCShore[" << n << "] = " << pVdFreeSurfaceStdCShore->at(n) << " pVdWaveSetupSurgeCShore[" << n << "] = " << pVdWaveSetupSurgeCShore->at(n) << " pVdSinWaveAngleRadiansCShore[" << n << "] = " << pVdSinWaveAngleRadiansCShore->at(n) << " pVdFractionBreakingWavesCShore[" << n << "] = " << pVdFractionBreakingWavesCShore->at(n) << endl;
   //
   // LogStream << endl;
   // // DEBUG CODE ========================================================================================================

   if (bTooShort)
   {
      // Add extrapolated value(s) to the end of the valid part of each profile vector so that we have nProfileSize valid values
      double dLastDiff = pVdXYDistFromCShore->at(nOutSize - 1) - pVdXYDistFromCShore->at(nOutSize - 2);

      for (int n = 0; n < nTooShort; n++)
         pVdXYDistFromCShore->at(nOutSize + n) = pVdXYDistFromCShore->at(nOutSize - 1 + n) + dLastDiff;

      for (int n = 0; n < nTooShort; n++)
         pVdFreeSurfaceStdCShore->at(nOutSize + n) = pVdFreeSurfaceStdCShore->at(nOutSize - 1);

      for (int n = 0; n < nTooShort; n++)
         pVdWaveSetupSurgeCShore->at(nOutSize + n) = pVdWaveSetupSurgeCShore->at(nOutSize - 1);

      // TODO 007 Do same for pVdStormSurgeCShore and pVdWaveSetupRunUpCShore ?

      for (int n = 0; n < nTooShort; n++)
         pVdSinWaveAngleRadiansCShore->at(nOutSize + n) = pVdSinWaveAngleRadiansCShore->at(nOutSize - 1);

      dLastDiff = pVdFractionBreakingWavesCShore->at(nOutSize - 1) - pVdFractionBreakingWavesCShore->at(nOutSize - 2);

      for (int n = 0; n < nTooShort; n++)
         pVdFractionBreakingWavesCShore->at(nOutSize + n) = tMin(pVdFractionBreakingWavesCShore->at(nOutSize - 1 + n) + dLastDiff, 1.0);

      // // DEBUG CODE ========================================================================================================
      // LogStream << "EXTENDED" << endl;
      // for (int n = 0; n < nProfileSize; n++)
      // LogStream << "pVdXYDistFromCShore[" << n << "] = " << pVdXYDistFromCShore->at(n) << " pVdFreeSurfaceStdCShore[" << n << "] = " << pVdFreeSurfaceStdCShore->at(n) << " pVdWaveSetupSurgeCShore[" << n << "] = " << pVdWaveSetupSurgeCShore->at(n) << " pVdSinWaveAngleRadiansCShore[" << n << "] = " << pVdSinWaveAngleRadiansCShore->at(n) << " pVdFractionBreakingWavesCShore[" << n << "] = " << pVdFractionBreakingWavesCShore->at(n) << " " << (n > (nOutSize-1) ? "EXTRAPOLATED" : "") << endl;
      //
      // LogStream << endl;
      // // DEBUG CODE ========================================================================================================
   }

   vector<double> VdXYDistCShoreTmp(nProfileSize);
   copy(pVdXYDistFromCShore->begin(), pVdXYDistFromCShore->begin() + nProfileSize, begin(VdXYDistCShoreTmp));

   // The CShore cross-shore distance has its origin at the seaward end, but the CoastalME convention has its origin at the shoreline. So we mst reverse the other vectors to conform with the CoastalME convention
   vector<double> VdFreeSurfaceStdCShoreTmp(nProfileSize);
   reverse_copy(pVdFreeSurfaceStdCShore->begin(), pVdFreeSurfaceStdCShore->begin() + nProfileSize, begin(VdFreeSurfaceStdCShoreTmp));

   vector<double> VdWaveSetupSurgeCShoreTmp(nProfileSize);
   reverse_copy(pVdWaveSetupSurgeCShore->begin(), pVdWaveSetupSurgeCShore->begin() + nProfileSize, begin(VdWaveSetupSurgeCShoreTmp));

   vector<double> VdSinWaveAngleRadiansCShoreTmp(nProfileSize);
   reverse_copy(pVdSinWaveAngleRadiansCShore->begin(), pVdSinWaveAngleRadiansCShore->begin() + nProfileSize, begin(VdSinWaveAngleRadiansCShoreTmp));

   vector<double> VdFractionBreakingWavesCShoreTmp(nProfileSize);
   reverse_copy(pVdFractionBreakingWavesCShore->begin(), pVdFractionBreakingWavesCShore->begin() + nProfileSize, begin(VdFractionBreakingWavesCShoreTmp));

   // // DEBUG CODE ========================================================================================================
   // LogStream << "REVERSED" << endl;
   // for (int n = 0; n < nProfileSize; n++)
   // LogStream << "VdXYDistCShoreTmp[" << n << "] = " << VdXYDistCShoreTmp.at(n) << " VdFreeSurfaceStdCShoreTmp[" << n << "] = " << VdFreeSurfaceStdCShoreTmp.at(n) << " VdWaveSetupSurgeCShoreTmp[" << n << "] = " << VdWaveSetupSurgeCShoreTmp.at(n) << " VdSinWaveAngleRadiansCShoreTmp[" << n << "] = " << VdSinWaveAngleRadiansCShoreTmp.at(n) << " VdFractionBreakingWavesCShoreTmp[" << n << "] = " << VdFractionBreakingWavesCShoreTmp.at(n) << endl;
   //
   // LogStream << endl;
   // // DEBUG CODE ========================================================================================================

   // Finally we linearly interpolate the CShore output vectors to each point on the CoastalME profile. Input parameters x_old, y_old, x_new and the routine returns y_new
   *pVdFreeSurfaceStdCME = VdInterpolateCShoreProfileOutput(&VdXYDistCShoreTmp, pVdFreeSurfaceStdCShore, pVdProfileDistXYCME); // was &VdDistXYCopy);
   *pVdWaveSetupSurgeCME = VdInterpolateCShoreProfileOutput(&VdXYDistCShoreTmp, pVdWaveSetupSurgeCShore, pVdProfileDistXYCME);
   // *pVdStormSurgeCME = VdInterpolateCShoreProfileOutput(&VdXYDistCShoreTmp, pVdStormSurgeCShore, pVdProfileDistXYCME);
   // *pVdWaveSetupRunUpCME = VdInterpolateCShoreProfileOutput(&VdXYDistCShoreTmp, pVdWaveSetupRunUpCShore, pVdProfileDistXYCME);
   *pVdSinWaveAngleRadiansCME = VdInterpolateCShoreProfileOutput(&VdXYDistCShoreTmp, pVdSinWaveAngleRadiansCShore, pVdProfileDistXYCME);
   *pVdFractionBreakingWavesCME = VdInterpolateCShoreProfileOutput(&VdXYDistCShoreTmp, pVdFractionBreakingWavesCShore, pVdProfileDistXYCME);

   // LogStream << "INTERPOLATED" << endl;
   // for (int n = 0; n < nProfileSize; n++)
   // LogStream << "pVdProfileDistXYCME[" << n << "] = " << pVdProfileDistXYCME->at(n) << " pVdFreeSurfaceStdCME[" << n << "] = " << pVdFreeSurfaceStdCME->at(n) << " pVdWaveSetupSurgeCME[" << n << "] = " << pVdWaveSetupSurgeCME->at(n) << " pVdSinWaveAngleRadiansCME[" << n << "] = " << pVdSinWaveAngleRadiansCME->at(n) << " pVdFractionBreakingWavesCME[" << n << "] = " << pVdFractionBreakingWavesCME->at(n) << endl;
   //
   // LogStream << "================================================ " << endl;
}
#endif

//===============================================================================================================================
//! Modifies the wave breaking properties at coastline points of profiles within the shadow zone
//===============================================================================================================================
void CSimulation::ModifyBreakingWavePropertiesWithinShadowZoneToCoastline(int const nCoast, int const nProfile)
{
   CGeomProfile *pProfile = m_VCoast[nCoast].pGetProfile(nProfile);

   // Only work on this profile if it is either a start-of-coast or an end-of-coast profile, or it is problem-free
   if (! pProfile->bProfileOKOrStartEndProfile())
   {
      if (m_nLogFileDetail >= LOG_FILE_ALL)
         LogStream << m_ulIter << ":\t coast " << nCoast << ", profile " << pProfile->nGetProfileID() << " is not OK, will not modify shadow zone breaking wave properties on this profile" << endl;

      return;
   }

   bool bModfiedWaveHeightisBreaking = false;
   bool bProfileIsinShadowZone = false;
   int const nThisCoastPoint = pProfile->nGetCoastPoint();
   int const nProfileSize = pProfile->nGetNumCellsInProfile();
   int nThisBreakingDist = m_VCoast[nCoast].nGetBreakingDistance(nThisCoastPoint);                 // This could be DBL_NODATA
   double dThisBreakingWaveHeight = m_VCoast[nCoast].dGetBreakingWaveHeight(nThisCoastPoint);      // This could be DBL_NODATA
   double dThisBreakingWaveAngle = m_VCoast[nCoast].dGetBreakingWaveAngle(nThisCoastPoint);        // This could be DBL_NODATA
   double dThisBreakingDepth = m_VCoast[nCoast].dGetDepthOfBreaking(nThisCoastPoint);              // This could be DBL_NODATA

   // Traverse the profile landwards, checking if any profile cell is within the shadow zone
   for (int nProfilePoint = (nProfileSize - 1); nProfilePoint >= 0; nProfilePoint--)
   {
      int const nX = pProfile->pPtiGetCellInProfile(nProfilePoint)->nGetX();
      int const nY = pProfile->pPtiGetCellInProfile(nProfilePoint)->nGetY();

      // If there is any cell profile  within the shadow zone and waves are breaking then modify wave breaking properties otherwise continue
      if (m_pRasterGrid->m_Cell[nX][nY].bIsinAnyShadowZone())
      {
         bProfileIsinShadowZone = true;

         // Check if the new wave height is breaking
         double const dWaveHeight = m_pRasterGrid->m_Cell[nX][nY].dGetWaveHeight();

         // Check that wave height at the given point is lower than maximum real wave height. If breaking wave height is expected that no good wave height are obtained, so, do not take it
         if (dWaveHeight > (m_dDepthOfClosure * m_dBreakingWaveHeightDepthRatio) && (! bModfiedWaveHeightisBreaking) && (! bFPIsEqual(dThisBreakingWaveHeight, DBL_NODATA, TOLERANCE)))
         {
            // It is breaking
            bModfiedWaveHeightisBreaking = true;

            dThisBreakingWaveHeight = m_dDepthOfClosure * m_dBreakingWaveHeightDepthRatio;
            dThisBreakingWaveAngle = m_pRasterGrid->m_Cell[nX][nY].dGetWaveAngle();
            dThisBreakingDepth = m_dDepthOfClosure;
            nThisBreakingDist = nProfilePoint;
         }
      }
   }

   // Update breaking wave properties along coastal line object (Wave height, dir, distance). TODO 010 Update the active zone cells
   if (bProfileIsinShadowZone && bModfiedWaveHeightisBreaking) // Modified wave height is still breaking
   {
      // This coast point is in the active zone, so set breaking wave height, breaking wave angle, and depth of breaking for the coast point TODO 007 Where does the 0.78 come from? TODO 011 Should it be an input variable or a named constant?
      if (dThisBreakingWaveHeight > dThisBreakingDepth * 0.78)
      {
         dThisBreakingWaveHeight = dThisBreakingDepth * 0.78; // Likely CShore output wave height is not adequately reproduced due to input profile and wave properties. TODO 007 Finish surge and runup stuff. Does something need to be changed then?
      }

      assert(dThisBreakingWaveHeight >= 0);
      m_VCoast[nCoast].SetBreakingWaveHeight(nThisCoastPoint, dThisBreakingWaveHeight);
      m_VCoast[nCoast].SetBreakingWaveAngle(nThisCoastPoint, dThisBreakingWaveAngle);
      m_VCoast[nCoast].SetDepthOfBreaking(nThisCoastPoint, dThisBreakingDepth);
      m_VCoast[nCoast].SetBreakingDistance(nThisCoastPoint, nThisBreakingDist);

      // LogStream << m_ulIter << ": nProfile = " << nProfile << ", nCoastPoint = " << nCoastPoint << " in active zone, dBreakingWaveHeight = " << dBreakingWaveHeight << endl;
   }

   else if (bProfileIsinShadowZone && (! bModfiedWaveHeightisBreaking))
   {
      // This coast point is no longer in the active zone
      m_VCoast[nCoast].SetBreakingWaveHeight(nThisCoastPoint, DBL_NODATA);
      m_VCoast[nCoast].SetBreakingWaveAngle(nThisCoastPoint, DBL_NODATA);
      m_VCoast[nCoast].SetDepthOfBreaking(nThisCoastPoint, DBL_NODATA);
      m_VCoast[nCoast].SetBreakingDistance(nThisCoastPoint, INT_NODATA);

      // LogStream << m_ulIter << ": nProfile = " << nProfile << ", nCoastPoint = " << nCoastPoint << " NOT in active zone" << endl;
   }

   return;
}

//===============================================================================================================================
//! Interpolates wave properties from profiles to the coastline points between two profiles. Do this by weighting the wave properties so that they change smoothly between the two profiles
//===============================================================================================================================
void CSimulation::InterpolateWavePropertiesBetweenProfiles(int const nCoast, int const nCount)
{
   CGeomProfile const* pProfile = m_VCoast[nCoast].pGetProfileWithDownCoastSeq(nCount);

   // Only do this for profiles without problems, including the start-of-coast profile (but not the end-of-coast profile)
   // if (! pProfile->bOKIncStartOfCoast())
   // return;

   int const nThisCoastPoint = pProfile->nGetCoastPoint();

   // For the breaking wave stuff, to go into the in-between coastline points
   int const nThisBreakingDist = m_VCoast[nCoast].nGetBreakingDistance(nThisCoastPoint);
   double const dThisBreakingWaveHeight = m_VCoast[nCoast].dGetBreakingWaveHeight(nThisCoastPoint);      // This could be DBL_NODATA
   double const dThisBreakingWaveAngle = m_VCoast[nCoast].dGetBreakingWaveAngle(nThisCoastPoint);        // This could be DBL_NODATA
   double const dThisBreakingDepth = m_VCoast[nCoast].dGetDepthOfBreaking(nThisCoastPoint);              // This could be DBL_NODATA
   double const dThisWaveSetupSurge = m_VCoast[nCoast].dGetWaveSetupSurge(nThisCoastPoint);              // This could be DBL_NODATA
   double const dThisRunUp = m_VCoast[nCoast].dGetRunUp(nThisCoastPoint);

   // Get the next profile along the coast, in the down-coast direction. If this next profile has a problem, go to the one after that, etc
   CGeomProfile const* pTmpProfile = pProfile;
   CGeomProfile const* pNextProfile;

   while (true)
   {
      pNextProfile = pTmpProfile->pGetDownCoastAdjacentProfile();

      if (pNextProfile->bProfileOKOrStartEndProfile())
      {
         // The next profile is either OK, or is a start-of-coast or an end-of-coast profile
         break;
      }

      // The next profile is not OK, so prepare to the one after that
      pTmpProfile = pNextProfile;

      if (pTmpProfile->bIsEndOfCoast())
      {
         // Uh-oh, we've reached the down-coast end of the coast without finding an OK down-coast profile. So give up
         return;
      }
   }

   // The next profile is OK
   int const nNextCoastPoint = pNextProfile->nGetCoastPoint();
   int const nDistBetween = nNextCoastPoint - nThisCoastPoint;

   // Safety check
   if (nDistBetween <= 0)
      // Nothing to do
      return;

   int const nNextBreakingDist = m_VCoast[nCoast].nGetBreakingDistance(nNextCoastPoint);                 // This could be DBL_NODATA
   double const dNextBreakingWaveHeight = m_VCoast[nCoast].dGetBreakingWaveHeight(nNextCoastPoint);      // This could be DBL_NODATA
   double const dNextBreakingWaveAngle = m_VCoast[nCoast].dGetBreakingWaveAngle(nNextCoastPoint);        // This could be DBL_NODATA
   double const dNextBreakingDepth = m_VCoast[nCoast].dGetDepthOfBreaking(nNextCoastPoint);              // This could be DBL_NODATA
   double const dNextWaveSetupSurge = m_VCoast[nCoast].dGetWaveSetupSurge(nNextCoastPoint);
   double const dNextRunUp = m_VCoast[nCoast].dGetRunUp(nNextCoastPoint);

   // OK, fill coast point between profiles for setupsurge and runup
   for (int n = nThisCoastPoint; n <= nNextCoastPoint; n++)
   {
      // At this point on the coast, are waves on- or off-shore, and up- or down-coast? First get the flux orientation (a tangent to the coastline)
      double const dFluxOrientation = m_VCoast[nCoast].dGetFluxOrientation(n);

      // If this coast point is in the active zone, use the breaking wave orientation, otherwise use the deep water wave orientation
      double dWaveAngle;
      if (bFPIsEqual(m_VCoast[nCoast].dGetDepthOfBreaking(n), DBL_NODATA, TOLERANCE))
         // Not in active zone
         dWaveAngle = m_VCoast[nCoast].dGetCoastDeepWaterWaveAngle(n);
      else
         // In active zone
         dWaveAngle = m_VCoast[nCoast].dGetBreakingWaveAngle(n);
      bool bDownCoast = false;
      int const nSeaHand = m_VCoast[nCoast].nGetSeaHandedness();
      bool const bOnShore = bOnOrOffShoreAndUpOrDownCoast(dFluxOrientation, dWaveAngle, nSeaHand, bDownCoast);
      m_VCoast[nCoast].SetWavesOnShore(n, bOnShore);
      m_VCoast[nCoast].SetWavesDownCoast(n, bDownCoast);

      // Fill first wave setup and surge
      int const nDist = n - nThisCoastPoint;
      double const dThisWeight = (nDistBetween - nDist) / static_cast<double>(nDistBetween);
      double const dNextWeight = 1 - dThisWeight;
      double dWaveSetupSurge = 0;
      double dRunUp = 0;

      dWaveSetupSurge = (dThisWeight * dThisWaveSetupSurge) + (dNextWeight * dNextWaveSetupSurge);
      m_VCoast[nCoast].SetWaveSetupSurge(n, dWaveSetupSurge);

      dRunUp = (dThisWeight * dThisRunUp) + (dNextWeight * dNextRunUp);
      m_VCoast[nCoast].SetRunUp(n, dRunUp);
   }

   // int const nNextProfile = pNextProfile->nGetProfileID();

   // If both this profile and the next profile are not in the active zone, then do no more
   if ((bFPIsEqual(dThisBreakingWaveHeight, DBL_NODATA, TOLERANCE)) && (bFPIsEqual(dNextBreakingWaveHeight, DBL_NODATA, TOLERANCE)))
   {
      // if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
      // LogStream << m_ulIter << ": both profile " << pProfile->nGetProfileID() << " at coast point " << nThisCoastPoint << ", and profile " << nNextProfile << " at coast point " << nNextCoastPoint << ", are not in the active zone" << endl;

      // Set the breaking wave height, breaking wave angle, and depth of breaking to DBL_NODATA
      for (int n = nThisCoastPoint; n < nNextCoastPoint; n++)
      {
         m_VCoast[nCoast].SetBreakingWaveHeight(nThisCoastPoint, DBL_NODATA);
         m_VCoast[nCoast].SetBreakingWaveAngle(nThisCoastPoint, DBL_NODATA);
         m_VCoast[nCoast].SetDepthOfBreaking(nThisCoastPoint, DBL_NODATA);
         m_VCoast[nCoast].SetBreakingDistance(nThisCoastPoint, DBL_NODATA);
      }

      return;
   }

   // OK, at least one of the two profiles is in the active zone
   // if (bFPIsEqual(dThisBreakingWaveHeight, DBL_NODATA, TOLERANCE))
   if (! bFPIsEqual(dNextBreakingWaveHeight, DBL_NODATA, TOLERANCE))
   {
      // The next profile must be in the active zone, so use values from the next profile
      for (int n = nThisCoastPoint; n < nNextCoastPoint; n++)
      {
         // Set the breaking wave height, breaking wave angle, and depth of breaking for this coast point TODO 056 This assert sometimes fails: why?
         assert(dNextBreakingWaveHeight >= 0);
         m_VCoast[nCoast].SetBreakingWaveHeight(n, dNextBreakingWaveHeight);
         m_VCoast[nCoast].SetBreakingWaveAngle(n, dNextBreakingWaveAngle);
         m_VCoast[nCoast].SetDepthOfBreaking(n, dNextBreakingDepth);
         m_VCoast[nCoast].SetBreakingDistance(n, nNextBreakingDist);
      }

      return;
   }

   // if (bFPIsEqual(dNextBreakingWaveHeight, DBL_NODATA, TOLERANCE))
   if (! bFPIsEqual(dThisBreakingWaveHeight, DBL_NODATA, TOLERANCE))
   {
      // This profile must be in the active zone, so use values from this profile
      for (int n = nThisCoastPoint + 1; n <= nNextCoastPoint; n++)
      {
         // Set the breaking wave height, breaking wave angle, and depth of breaking for this coast point TODO 056 This assert sometimes fails: why?
         assert(dThisBreakingWaveHeight >= 0);
         m_VCoast[nCoast].SetBreakingWaveHeight(n, dThisBreakingWaveHeight);
         m_VCoast[nCoast].SetBreakingWaveAngle(n, dThisBreakingWaveAngle);
         m_VCoast[nCoast].SetDepthOfBreaking(n, dThisBreakingDepth);
         m_VCoast[nCoast].SetBreakingDistance(n, nThisBreakingDist);
      }

      return;
   }

   // The values for both this profile point and the next profile point are fine, so do a weighted interpolation between this profile and the next profile
   for (int n = nThisCoastPoint + 1; n < nNextCoastPoint; n++)
   {
      int const nDist = n - nThisCoastPoint;

      double dBreakingWaveHeight = DBL_NODATA;
      double dBreakingWaveAngle = DBL_NODATA;
      double dBreakingDepth = DBL_NODATA;
      double dBreakingDist = DBL_NODATA;

      if ((dNextBreakingDepth > 0) && (dThisBreakingDepth > 0))
      {
         double const dThisWeight = (nDistBetween - nDist) / static_cast<double>(nDistBetween);
         double const dNextWeight = 1 - dThisWeight;

         dBreakingWaveHeight = (dThisWeight * dThisBreakingWaveHeight) + (dNextWeight * dNextBreakingWaveHeight),
         dBreakingWaveAngle = (dThisWeight * dThisBreakingWaveAngle) + (dNextWeight * dNextBreakingWaveAngle),
         dBreakingDepth = (dThisWeight * dThisBreakingDepth) + (dNextWeight * dNextBreakingDepth),
         dBreakingDist = (dThisWeight * nThisBreakingDist) + (dNextWeight * nNextBreakingDist);
      }
      else if (dNextBreakingDepth > 0)
      {
         dBreakingWaveHeight = dNextBreakingWaveHeight,
         dBreakingWaveAngle = dNextBreakingWaveAngle,
         dBreakingDepth = dNextBreakingDepth,
         dBreakingDist = nNextBreakingDist;
      }
      else if (dThisBreakingDepth > 0)
      {
         dBreakingWaveHeight = dThisBreakingWaveHeight,
         dBreakingWaveAngle = dThisBreakingWaveAngle,
         dBreakingDepth = dThisBreakingDepth,
         dBreakingDist = nThisBreakingDist;
      }

      // Set the breaking wave height, breaking wave angle, and depth of breaking for this coast point. Note that these could be DBL_NODATA
      m_VCoast[nCoast].SetBreakingWaveHeight(n, dBreakingWaveHeight);
      m_VCoast[nCoast].SetBreakingWaveAngle(n, dBreakingWaveAngle);
      m_VCoast[nCoast].SetDepthOfBreaking(n, dBreakingDepth);
      m_VCoast[nCoast].SetBreakingDistance(n, nRound(dBreakingDist));
   }
}

//===============================================================================================================================
//! Linearly interpolates wave properties from profiles to the coastline cells between two profiles
//===============================================================================================================================
void CSimulation::InterpolateWaveHeightToCoastPoints(int const nCoast)
{
   int const nCoastPoints = m_VCoast[nCoast].nGetCoastlineSize();

   // Initialize all vectors pairs (x,y) for each variable
   vector<int> nVCoastWaveHeightX;
   vector<double> dVCoastWaveHeightY;

   // Search all coast points for non NaN values and store them into temporary variables for later interporlation
   for (int n = 0; n < nCoastPoints; n++)
   {
      double const dCoastWaveHeight = m_VCoast[nCoast].dGetCoastWaveHeight(n);

      if (! bFPIsEqual(dCoastWaveHeight, DBL_NODATA, TOLERANCE))
      {
         nVCoastWaveHeightX.push_back(n);
         dVCoastWaveHeightY.push_back(dCoastWaveHeight);
      }
   }

   // Interpolate all coast points. Check first that x,y have more than 3 points and are both of equal size
   if ((nVCoastWaveHeightX.size() >= 3) && (nVCoastWaveHeightX.size() == dVCoastWaveHeightY.size()))
   {
      for (int n = 0; n < nCoastPoints; n++)
      {
         double const dInterpCoastWaveHeight = dGetInterpolatedValue(&nVCoastWaveHeightX, &dVCoastWaveHeightY, n, false);
         m_VCoast[nCoast].SetCoastWaveHeight(n, dInterpCoastWaveHeight);
      }
   }

   else
   {
      for (int n = 0; n < nCoastPoints; n++)
      {
         m_VCoast[nCoast].SetCoastWaveHeight(n, 0);
      }
   }
}

//===============================================================================================================================
//! Calculates tangents to a coastline: the tangent is assumed to be the orientation of energy/sediment flux along a coast. The tangent is specified as an angle (in degrees) measured clockwise from north. Based on a routine by Martin Hurst
//===============================================================================================================================
void CSimulation::CalcCoastTangents(int const nCoast)
{
   int const nCoastSize = m_VCoast[nCoast].nGetCoastlineSize();

   for (int nCoastPoint = 0; nCoastPoint < nCoastSize; nCoastPoint++)
   {
      double dXDiff;
      double dYDiff;

      if (nCoastPoint == 0)
      {
         // For the point at the start of the coastline: use the straight line from 'this' point to the next point
         CGeom2DPoint const PtThis = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nCoastPoint);      // In external CRS
         CGeom2DPoint const PtAfter = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nCoastPoint + 1); // In external CRS

         dXDiff = PtAfter.dGetX() - PtThis.dGetX();
         dYDiff = PtAfter.dGetY() - PtThis.dGetY();
      }

      else if (nCoastPoint == nCoastSize - 1)
      {
         // For the point at the end of the coastline: use the straight line from the point before to 'this' point
         CGeom2DPoint const PtBefore = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nCoastPoint - 1); // In external CRS
         CGeom2DPoint const PtThis = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nCoastPoint);       // In external CRS

         dXDiff = PtThis.dGetX() - PtBefore.dGetX();
         dYDiff = PtThis.dGetY() - PtBefore.dGetY();
      }

      else
      {
         // For coastline points not at the start or end of the coast: start with a straight line which links the coastline points before and after 'this' coastline point
         CGeom2DPoint const PtBefore = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nCoastPoint - 1); // In external CRS
         CGeom2DPoint const PtAfter = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nCoastPoint + 1);  // In external CRS

         dXDiff = PtAfter.dGetX() - PtBefore.dGetX();
         dYDiff = PtAfter.dGetY() - PtBefore.dGetY();
      }

      // Calculate angle between line and north point, measured clockwise (the azimuth)
      if (bFPIsEqual(dYDiff, 0.0, TOLERANCE))
      {
         // The linking line runs either W-E or E-W
         if (dXDiff > 0)
            // It runs W to E
            m_VCoast[nCoast].SetFluxOrientation(nCoastPoint, 90);

         else
            // It runs E to W
            m_VCoast[nCoast].SetFluxOrientation(nCoastPoint, 270);
      }

      else if (bFPIsEqual(dXDiff, 0.0, TOLERANCE))
      {
         // The linking line runs N-S or S-N
         if (dYDiff > 0)
            // It runs S to N
            m_VCoast[nCoast].SetFluxOrientation(nCoastPoint, 0);

         else
            // It runs N to S
            m_VCoast[nCoast].SetFluxOrientation(nCoastPoint, 180);
      }

      else
      {
         // The linking line runs neither W-E nor N-S so we have to work a bit harder to find the angle between it and the azimuth
         double dAzimuthAngle;

         if (dXDiff > 0)
            dAzimuthAngle = (180 / PI) * (PI * 0.5 - atan(dYDiff / dXDiff));

         else
            dAzimuthAngle = (180 / PI) * (PI * 1.5 - atan(dYDiff / dXDiff));

         m_VCoast[nCoast].SetFluxOrientation(nCoastPoint, dAzimuthAngle);
      }

      // LogStream << m_ulIter << ": nCoastPoint = " << nCoastPoint << " FluxOrientation = " << m_VCoast[nCoast].dGetFluxOrientation(nCoastPoint) << endl;
   }
}

//===============================================================================================================================
//! Calculates an average d50 for each polygon. Also fills in 'holes' in active zone and wave calcs i.e. orphan cells which should have been included in the active zone but which have been omitted because of rounding problems
//===============================================================================================================================
void CSimulation::CalcD50AndFillWaveCalcHoles(void)
{
   // Get the total number of polygons, all coasts
   int nTotPolygonAllCoasts = 0;
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
      nTotPolygonAllCoasts +=  m_VCoast[nCoast].nGetNumPolygons();

   // Vectors for D50 stuff
   vector<int> VnPolygonD50Count(nTotPolygonAllCoasts, 0);
   vector<double> VdPolygonD50(nTotPolygonAllCoasts, 0);

   for (int nX = 0; nX < m_nXGridSize; nX++)
   {
      for (int nY = 0; nY < m_nYGridSize; nY++)
      {
         if (m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea())
         {
            // This is a sea cell, first get polygon ID
            int const nPolyID = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonID();

            // Is it in the active zone?
            bool const bActive = m_pRasterGrid->m_Cell[nX][nY].bIsInActiveZone();

            if (bActive)
            {
               // It is in the active zone. Does it have unconsolidated sediment on it? Test this using the UnconD50 value: if dGetUnconsD50() returns DBL_NODATA, there is no unconsolidated sediment
               double const dTmpd50 = m_pRasterGrid->m_Cell[nX][nY].dGetUnconsD50();

               if (! bFPIsEqual(dTmpd50, DBL_NODATA, TOLERANCE))
               {
                  // It does have unconsolidated sediment, so which polygon is this cell in?
                  if (nPolyID != INT_NODATA)
                  {
                     VnPolygonD50Count[nPolyID]++;
                     VdPolygonD50[nPolyID] += dTmpd50;
                  }
               }
            }

            //             // Now fill in wave calc holes
            // if (m_pRasterGrid->m_Cell[nX][nY].dGetWaveHeight() == 0)
            // {
            // if (nID == INT_NODATA)
            // m_pRasterGrid->m_Cell[nX][nY].SetWaveHeight(m_dAllCellsDeepWaterWaveHeight);
            // }
            //
            // if (m_pRasterGrid->m_Cell[nX][nY].dGetWaveAngle() == 0)
            // {
            // if (nID == INT_NODATA)
            // m_pRasterGrid->m_Cell[nX][nY].SetWaveAngle(m_dAllCellsDeepWaterWaveAngle);
            // }

            // Next look at the cell's N-S and W-E neighbours
            int nXTmp;
            int nYTmp;
            int nActive = 0;
            int nShadow = 0;
            int nShadowNum = 0;
            int nDownDrift = 0;
            int nDownDriftNum = 0;
            int nCoast = 0;
            int nRead = 0;
            double dWaveHeight = 0;
            double dWaveAngle = 0;

            // North
            nXTmp = nX;
            nYTmp = nY - 1;

            if (bIsWithinValidGrid(nXTmp, nYTmp))
            {
               if (m_pRasterGrid->m_Cell[nXTmp][nYTmp].bIsInContiguousSea())
               {
                  nRead++;
                  dWaveHeight += m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetWaveHeight();
                  dWaveAngle += (m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetWaveAngle());

                  if (m_pRasterGrid->m_Cell[nXTmp][nYTmp].bIsInActiveZone())
                     nActive++;

                  int nTmp = m_pRasterGrid->m_Cell[nXTmp][nYTmp].nGetShadowZoneNumber();

                  if (nTmp != 0)
                  {
                     nShadow++;
                     nShadowNum = nTmp;
                  }

                  nTmp = m_pRasterGrid->m_Cell[nXTmp][nYTmp].nGetDownDriftZoneNumber();

                  if (nTmp > 0)
                  {
                     nDownDrift++;
                     nDownDriftNum = nTmp;
                  }
               }
               else
                  nCoast++;
            }

            // East
            nXTmp = nX + 1;
            nYTmp = nY;

            if (bIsWithinValidGrid(nXTmp, nYTmp))
            {
               if (m_pRasterGrid->m_Cell[nXTmp][nYTmp].bIsInContiguousSea())
               {
                  nRead++;
                  dWaveHeight += m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetWaveHeight();
                  dWaveAngle += (m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetWaveAngle());

                  if (m_pRasterGrid->m_Cell[nXTmp][nYTmp].bIsInActiveZone())
                     nActive++;

                  int nTmp = m_pRasterGrid->m_Cell[nXTmp][nYTmp].nGetShadowZoneNumber();

                  if (nTmp != 0)
                  {
                     nShadow++;
                     nShadowNum = nTmp;
                  }

                  nTmp = m_pRasterGrid->m_Cell[nXTmp][nYTmp].nGetDownDriftZoneNumber();

                  if (nTmp > 0)
                  {
                     nDownDrift++;
                     nDownDriftNum = nTmp;
                  }
               }
               else
                  nCoast++;
            }

            // South
            nXTmp = nX;
            nYTmp = nY + 1;

            if (bIsWithinValidGrid(nXTmp, nYTmp))
            {
               if (m_pRasterGrid->m_Cell[nXTmp][nYTmp].bIsInContiguousSea())
               {
                  nRead++;
                  dWaveHeight += m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetWaveHeight();
                  dWaveAngle += (m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetWaveAngle());

                  if (m_pRasterGrid->m_Cell[nXTmp][nYTmp].bIsInActiveZone())
                     nActive++;

                  int nTmp = m_pRasterGrid->m_Cell[nXTmp][nYTmp].nGetShadowZoneNumber();

                  if (nTmp != 0)
                  {
                     nShadow++;
                     nShadowNum = nTmp;
                  }

                  nTmp = m_pRasterGrid->m_Cell[nXTmp][nYTmp].nGetDownDriftZoneNumber();

                  if (nTmp > 0)
                  {
                     nDownDrift++;
                     nDownDriftNum = nTmp;
                  }
               }
               else
                  nCoast++;
            }

            // West
            nXTmp = nX - 1;
            nYTmp = nY;

            if (bIsWithinValidGrid(nXTmp, nYTmp))
            {
               if (m_pRasterGrid->m_Cell[nXTmp][nYTmp].bIsInContiguousSea())
               {
                  nRead++;
                  dWaveHeight += m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetWaveHeight();
                  dWaveAngle += (m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetWaveAngle());

                  if (m_pRasterGrid->m_Cell[nXTmp][nYTmp].bIsInActiveZone())
                     nActive++;

                  int nTmp = m_pRasterGrid->m_Cell[nXTmp][nYTmp].nGetShadowZoneNumber();

                  if (nTmp != 0)
                  {
                     nShadow++;
                     nShadowNum = nTmp;
                  }

                  nTmp = m_pRasterGrid->m_Cell[nXTmp][nYTmp].nGetDownDriftZoneNumber();

                  if (nTmp > 0)
                  {
                     nDownDrift++;
                     nDownDriftNum = nTmp;
                  }
               }
               else
                  nCoast++;
            }

            if (nRead > 0)
            {
               // Calculate the average of neighbours
               dWaveHeight /= nRead;
               dWaveAngle /= nRead;
               dWaveAngle = dKeepWithin360(dWaveAngle);

               // If this sea cell has four active-zone neighbours, then it must also be in the active zone: give it wave height and orientation which is the average of its neighbours
               if (nActive == 4)
               {
                  m_pRasterGrid->m_Cell[nX][nY].SetInActiveZone(true);
                  m_pRasterGrid->m_Cell[nX][nY].SetWaveHeight(dWaveHeight);
                  m_pRasterGrid->m_Cell[nX][nY].SetWaveAngle(dWaveAngle);
               }

               // If this sea cell has a wave height which is the same as its deep-water wave height, but its neighbours have a different average wave height, then give it the average of its neighbours
               double const dDeepWaterWaveHeight = m_pRasterGrid->m_Cell[nX][nY].dGetCellDeepWaterWaveHeight();

               if ((bFPIsEqual(m_pRasterGrid->m_Cell[nX][nY].dGetWaveHeight(), dDeepWaterWaveHeight, TOLERANCE)) && (! bFPIsEqual(dDeepWaterWaveHeight, dWaveHeight, TOLERANCE)))
               {
                  m_pRasterGrid->m_Cell[nX][nY].SetWaveHeight(dWaveHeight);
               }

               // If this sea cell has a wave orientation which is the same as its deep-water wave orientation, but its neighbours have a different average wave orientation, then give it the average of its neighbours
               double const dDeepWaterWaveAngle = m_pRasterGrid->m_Cell[nX][nY].dGetCellDeepWaterWaveAngle();

               if ((bFPIsEqual(m_pRasterGrid->m_Cell[nX][nY].dGetWaveAngle(), dDeepWaterWaveAngle, TOLERANCE)) && (! bFPIsEqual(dDeepWaterWaveAngle, dWaveAngle, TOLERANCE)))
               {
                  m_pRasterGrid->m_Cell[nX][nY].SetWaveAngle(dWaveAngle);
               }

               // Is this sea cell is not already marked as in a shadow zone (note could be marked as in a shadow zone but not yet processed: a -ve number)?
               int const nShadowZoneCode = m_pRasterGrid->m_Cell[nX][nY].nGetShadowZoneNumber();

               if (nShadowZoneCode <= 0)
               {
                  // If the cell has four neighbours which are all in a shadow zone, or four neighbours some of which are shadow zone and the remainder downdrift zone, or four neighbours some of which are shadow zone and the remainder coast; then it should also be in the shadow zone: give it the average of its neighbours
                  if ((nShadow == 4) || (nShadow + nDownDrift == 4) || (nShadow + nCoast == 4))
                  {
                     m_pRasterGrid->m_Cell[nX][nY].SetShadowZoneNumber(nShadowNum);
                     m_pRasterGrid->m_Cell[nX][nY].SetWaveHeight(dWaveHeight);
                     m_pRasterGrid->m_Cell[nX][nY].SetWaveAngle(dWaveAngle);
                  }
               }

               // If this sea cell is not marked as in a downdrift zone but has four neighbours which are in a downdrift zone, then it should also be in the downdrift zone: give it the average of its neighbours
               int const nDownDriftZoneCode = m_pRasterGrid->m_Cell[nX][nY].nGetDownDriftZoneNumber();

               if ((nDownDriftZoneCode == 0) && (nDownDrift == 4))
               {
                  m_pRasterGrid->m_Cell[nX][nY].SetDownDriftZoneNumber(nDownDriftNum);
                  m_pRasterGrid->m_Cell[nX][nY].SetWaveHeight(dWaveHeight);
                  m_pRasterGrid->m_Cell[nX][nY].SetWaveAngle(dWaveAngle);
               }
            }
         }
      }
   }

   // Calculate the average d50 for every polygon
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      for (int nPoly = 0; nPoly < m_VCoast[nCoast].nGetNumPolygons(); nPoly++)
      {
         CGeomCoastPolygon* pPolygon = m_VCoast[nCoast].pGetPolygon(nPoly);
         int const nPolyID = pPolygon->nGetPolygonCoastID();

         if (VnPolygonD50Count[nPolyID] > 0)
            VdPolygonD50[nPolyID] /= VnPolygonD50Count[nPolyID];

         pPolygon->SetAvgUnconsD50(VdPolygonD50[nPolyID]);
      }
   }
}
