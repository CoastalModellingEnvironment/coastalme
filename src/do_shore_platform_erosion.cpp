/*!
   \file do_shore_platform_erosion.cpp
   \brief Erodes the consolidated sediment of the shore platform. Eroded sediment from the shore platform becomes unconsolidated sediment stored in coastal polygons
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

#include <cmath>

#include <iostream>
using std::cerr;
using std::endl;
using std::ios;

#include <array>
using std::array;

#include <algorithm>
using std::shuffle;

#include "cme.h"
#include "hermite_cubic.h"
#include "simulation.h"
#include "coast.h"
#include "2di_point.h"

//===============================================================================================================================
//! Does platform erosion on all coastlines by first calculating platform erosion on coastline-normal profiles, then extrapolating this to cells between the profiles
//===============================================================================================================================
int CSimulation::nDoAllShorePlatFormErosion(void)
{
   if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
      LogStream << m_ulIter << ": Calculating shore platform erosion" << endl;

   // // DEBUG CODE ===========================================================================================================
   // string strOutFile = m_strOutPath + "01_polygon_shore_platform_test_";
   // strOutFile += to_string(m_ulIter);
   // strOutFile += ".tif";
   //
   // GDALDriver* pDriver = GetGDALDriverManager()->GetDriverByName("gtiff");
   // GDALDataset* pDataSet = pDriver->Create(strOutFile.c_str(), m_nXGridSize, m_nYGridSize, 1, GDT_Float64, m_papszGDALRasterOptions);
   // pDataSet->SetProjection(m_strGDALBasementDEMProjection.c_str());
   // pDataSet->SetGeoTransform(m_dGeoTransform);
   // double* pdRaster = new double[m_nXGridSize * m_nYGridSize];
   // int n = 0;
   // int nInPoly = 0;
   // int nNotInPoly = 0;
   //
   // for (int nY = 0; nY < m_nYGridSize; nY++)
   // {
   // for (int nX = 0; nX < m_nXGridSize; nX++)
   // {
   // int nID = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonID();
   // if (nID == INT_NODATA)
   // nNotInPoly++;
   // else
   // nInPoly++;
   //
   // pdRaster[n++] = nID;
   // }
   // }
   //
   // GDALRasterBand* pBand = pDataSet->GetRasterBand(1);
   // pBand->SetNoDataValue(m_dMissingValue);
   // int nRet = pBand->RasterIO(GF_Write, 0, 0, m_nXGridSize, m_nYGridSize, pdRaster, m_nXGridSize, m_nYGridSize, GDT_Float64, 0, 0, NULL);
   // if (nRet == CE_Failure)
   // LogStream << nRet << endl;
   //
   // GDALClose(pDataSet);
   // delete[] pdRaster;
   //
   // LogStream << m_ulIter << " Number of cells in a polygon = " << nInPoly << endl;
   // LogStream << m_ulIter << " Number of cells not in any polygon = " << nNotInPoly << endl;
   // // DEBUG CODE ===========================================================================================================

   // Set direction
   static bool bDownCoast = true;

   // Do this for each coast
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      int const nNumProfiles = m_VCoast[nCoast].nGetNumProfiles();

      // Calculate potential platform erosion along each coastline-normal profile, do in down-coast sequence
      for (int nn = 0; nn < nNumProfiles; nn++)
      {
         CGeomProfile *pProfile;

         if (bDownCoast)
            pProfile = m_VCoast[nCoast].pGetProfileWithDownCoastSeq(nn);
         else
            pProfile = m_VCoast[nCoast].pGetProfileWithUpCoastSeq(nn);

         int const nRet = nCalcPotentialPlatformErosionOnProfile(nCoast, pProfile);
         if (nRet != RTN_OK)
            return nRet;
      }

      // Calculate potential platform erosion between the coastline-normal profiles
      if (bDownCoast)
      {
         for (int nn = 0; nn < nNumProfiles - 1; nn++)
         {
            CGeomProfile *pProfile = m_VCoast[nCoast].pGetProfileWithDownCoastSeq(nn);

            // Calculate potential erosion for sea cells between this profile and the next profile (or up to the edge of the grid) on these cells
            int nRet = nCalcPotentialPlatformErosionBetweenProfiles(nCoast, pProfile, DIRECTION_DOWNCOAST);
            if (nRet != RTN_OK)
               return nRet;
         }
      }
      else
      {
         for (int nn = nNumProfiles - 1; nn > 0; nn--)
         {
            CGeomProfile *pProfile = m_VCoast[nCoast].pGetProfileWithDownCoastSeq(nn);

            // Calculate potential erosion for sea cells between this profile and the next profile (or up to the edge of the grid) on these cells
            int nRet = nCalcPotentialPlatformErosionBetweenProfiles(nCoast, pProfile, DIRECTION_UPCOAST);
            if (nRet != RTN_OK)
               return nRet;
         }
      }
   }

   // Swap direction for next timestep
   bDownCoast = ! bDownCoast;

   // Fill in 'holes' in the potential platform erosion i.e. orphan cells which get omitted because of rounding problems
   FillPotentialPlatformErosionHoles();

   // Do the same for beach protection
   FillInBeachProtectionHolesAndRemoveLegacyCliffs();

   // Finally calculate actual platform erosion on all sea cells (both on profiles, and between profiles)
   for (int nX = 0; nX < m_nXGridSize; nX++)
   {
      for (int nY = 0; nY < m_nYGridSize; nY++)
      {
         // If this cell has already had shore platform erosion during this timestep, then do nothing
         if (m_pRasterGrid->m_Cell[nX][nY].bGetPlatformErosionThisIter())
            continue;

         if (m_pRasterGrid->m_Cell[nX][nY].bPotentialPlatformErosion())
            // Calculate actual (supply-limited) shore platform erosion on each cell that has potential platform erosion, also add the eroded sand/coarse sediment to that cell's polygon, ready to be redistributed within the polygon during beach erosion/deposition
            DoActualPlatformErosionOnCell(nX, nY);
      }
   }

   if (m_nLogFileDetail >= LOG_FILE_ALL)
   {
      LogStream << m_ulIter << ":\t total potential shore platform erosion (m^3) = " << m_dThisIterPotentialPlatformErosion * m_dCellArea << " (on profiles = " << m_dTotPotentialPlatformErosionOnProfiles * m_dCellArea << ", between profiles = " << m_dTotPotentialPlatformErosionBetweenProfiles * m_dCellArea << ")" << endl;

      LogStream << m_ulIter << ":\t total actual shore platform erosion (m^3) = " << (m_dThisIterActualPlatformErosionFineCons + m_dThisIterActualPlatformErosionSandCons + m_dThisIterActualPlatformErosionCoarseCons) * m_dCellArea << " (fine = " << m_dThisIterActualPlatformErosionFineCons * m_dCellArea << ", sand = " << m_dThisIterActualPlatformErosionSandCons * m_dCellArea << ", coarse = " << m_dThisIterActualPlatformErosionCoarseCons * m_dCellArea << ")" << endl;
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Calculates potential (i.e. unconstrained by available sediment) erosional lowering of the shore platform for a single coastline-normal profile, due to wave action. This routine uses a behavioural rule to modify the original surface elevation profile geometry, in which erosion rate/slope = f(d/Db) based on Walkden & Hall (2005). Originally coded in Matlab by Andres Payo
//===============================================================================================================================
int CSimulation::nCalcPotentialPlatformErosionOnProfile(int const nCoast, CGeomProfile *pProfile)
{
   // Only work on this profile if it is problem-free TODO 024 Or if it has just hit dry land?
   if (! pProfile->bProfileOK())
      return RTN_OK;

   // Get the length of the profile (in cells) and the index of the coast point at which this profile starts
   int const nProfSize = pProfile->nGetNumCellsInProfile();
   int const nCoastPoint = pProfile->nGetCoastPoint();

   // Get the breaking depth for this profile from the coastline point
   double const dDepthOfBreaking = m_VCoast[nCoast].dGetDepthOfBreaking(nCoastPoint);

   if (bFPIsEqual(dDepthOfBreaking, DBL_NODATA, TOLERANCE))
      // This profile is not in the active zone, so no platform erosion here
      return RTN_OK;

   if (bFPIsEqual(dDepthOfBreaking, 0.0, TOLERANCE))
   {
      // Safety check, altho' this shouldn't happen
      if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
         LogStream << m_ulIter << ": depth of breaking is zero for profile " << pProfile->nGetProfileID() << " of coast " << nCoast << endl;

      return RTN_OK;
   }

   // LogStream << m_ulIter << ": ON PROFILE nProfile = " << nProfile << " dDepthOfBreaking = " << dDepthOfBreaking << endl;

   // Get the height of the associated breaking wave from the coast point: this height is used in beach protection calcs
   double const dBreakingWaveHeight = m_VCoast[nCoast].dGetBreakingWaveHeight(nCoastPoint);

   // Calculate the length of the profile in external CRS units
   int const nSegments = pProfile->nGetProfileSize() - 1;
   double dProfileLenXY = 0;

   for (int nSeg = 0; nSeg < nSegments; nSeg++)
   {
      // Do once for every line segment
      double const dSegStartX = pProfile->pPtGetPointInProfile(nSeg)->dGetX();
      double const dSegStartY = pProfile->pPtGetPointInProfile(nSeg)->dGetY();
      double const dSegEndX = pProfile->pPtGetPointInProfile(nSeg + 1)->dGetX(); // Is OK
      double const dSegEndY = pProfile->pPtGetPointInProfile(nSeg + 1)->dGetY();

      double const dSegLen = hypot(dSegStartX - dSegEndX, dSegStartY - dSegEndY);
      dProfileLenXY += dSegLen;
   }

   // Next calculate the average distance between profile points, again in external CRS units. Assume that the sample points are equally spaced along the profile (not quite true)
   double const dSpacingXY = dProfileLenXY / (nProfSize - 1);

   // Set up vectors for the coastline-normal profile elevations. The length of this vector line is given by the number of cells 'under' the profile. Thus each point on the vector relates to a single cell in the grid. This assumes that all points on the profile vector are equally spaced (not quite true, depends on the orientation of the line segments which comprise the profile). The elevation of each of these profile points is the elevation of the centroid of the cell that is 'under' the point. However we cannot always be confident that this is the 'true' elevation of the point on the vector since (unless the profile runs planview N-S or W-E) the vector does not always run exactly through the centroid of the cell
   vector<double> VdProfileZ(nProfSize, 0);                    // Initial (pre-erosion) elevation of both consolidated and unconsolidated sediment for cells 'under' the profile
   vector<double> VdProfileDistXY(nProfSize, 0);               // Along-profile distance measured from the coast, in external CRS units
   vector<double> dVConsProfileZ(nProfSize, 0);                // Initial (pre-erosion) elevation of consolidated sediment only for cells 'under' the profile
   vector<double> dVConsZDiff(nProfSize, 0);
   vector<double> dVConsSlope(nProfSize, 0);

   for (int i = 0; i < nProfSize; i++)
   {
      int const nX = pProfile->pPtiVGetCellsInProfile()->at(i).nGetX();
      int const nY = pProfile->pPtiVGetCellsInProfile()->at(i).nGetY();

      // Get the number of the highest layer with non-zero thickness
      int const nTopLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopNonZeroLayerAboveBasement();

      // Safety check
      if (nTopLayer == INT_NODATA)
         return RTN_ERR_NO_TOP_LAYER_DURING_PLATFORM_CALC;

      if (nTopLayer == NO_NONZERO_THICKNESS_LAYERS)
         // TODO 025 We are down to basement
         return RTN_OK;

      // Get the elevation for consolidated sediment only on this cell
      dVConsProfileZ[i] = m_pRasterGrid->m_Cell[nX][nY].dGetConsSedTopElevForLayerAboveBasement(nTopLayer);

      // Get the elevation for both consolidated and unconsolidated sediment on this cell (ignore any talus)
      VdProfileZ[i] = m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevOmitTalus();

      // And store the X-Y plane distance from the start of the profile
      VdProfileDistXY[i] = i * dSpacingXY;
   }

   for (int i = 0; i < nProfSize - 1; i++)
   {
      // For the consolidated-only profile, get the Z differences (already in external CRS units)
      dVConsZDiff[i] = dVConsProfileZ[i] - dVConsProfileZ[i + 1];

      // Calculate dZ/dXY, the Z slope (i.e. rise over run) in the XY direction. Note that we use the elevation difference on the seaward side of 'this' point
      dVConsSlope[i] = dVConsZDiff[i] / dSpacingXY;
   }

   // Sort out the final value
   dVConsSlope[nProfSize - 1] = dVConsSlope[nProfSize - 2];

   if (m_nProfileSmoothWindow > 0)
   {
      // Smooth the vector of slopes for the consolidated-only profile
      dVConsSlope = dVSmoothProfileSlope(&dVConsSlope);
   }

   vector<double> dVProfileDepthOverDB(nProfSize, 0);          // Depth over wave breaking depth at the coastline-normal sample points
   vector<double> dVProfileErosionPotential(nProfSize, 0);     // Erosion potential at the coastline-normal sample points

   // Calculate the erosion potential along this profile using the shape function
   double dTotalErosionPotential = 0;

   for (int i = 0; i < nProfSize; i++)
   {
      // Use the actual depth of water here (i.e. the depth to the top of the unconsolidated sediment, including the thickness of consolidated sediment beneath it)
      dVProfileDepthOverDB[i] = m_dThisIterSWL - VdProfileZ[i];
      dVProfileDepthOverDB[i] /= dDepthOfBreaking;

      // Constrain dDepthOverDB[i] to be between 0 (can get small -ve values due to rounding errors) and m_dDepthOverDBMax
      dVProfileDepthOverDB[i] = tMax(dVProfileDepthOverDB[i], 0.0);
      dVProfileDepthOverDB[i] = tMin(dVProfileDepthOverDB[i], m_dDepthOverDBMax);

      // And then use the look-up table to find the value of erosion potential at this point on the profile
      dVProfileErosionPotential[i] = dLookUpErosionPotential(dVProfileDepthOverDB[i]);

      // If erosion potential (a -ve value) is tiny, set it to zero
      if (dVProfileErosionPotential[i] > -SED_ELEV_TOLERANCE)
         dVProfileErosionPotential[i] = 0;

      // Keep track of the total erosion potential for this profile
      dTotalErosionPotential += dVProfileErosionPotential[i];
   }

   // Constrain erosion potential at every point on the profile, so that the integral of erosion potential on the whole profile is unity (Walkden and Hall 2005). Note that here, erosion potential is -ve so we must constrain to -1
   for (int i = 0; i < nProfSize; i++)
   {
      if (dTotalErosionPotential < 0)
         dVProfileErosionPotential[i] /= (-dTotalErosionPotential);
   }

   vector<double> dVRecessionXY(nProfSize, 0);
   vector<double> dVSCAPEXY(nProfSize, 0);

   // Calculate recession at every point on the coastline-normal profile
   for (int i = 0; i < nProfSize; i++)
   {
      // dRecession = dForce * (dBeachProtection / dR) * dErosionPotential * dSlope * dTime where:
      // dVRecession [m] is the landward migration distance defined in the profile relative (XY) CRS
      // dForce is given by Equation 4 in Walkden & Hall, 2005
      // dVBeachProtection [1] is beach protection factor [1, 0] = [no protection, fully protected]. (This is calculated later, see dCalcBeachProtectionFactor())
      // dVR  [m^(9/4)s^(2/3)] is the material strength and some hydrodynamic constant
      // dVProfileErosionPotential [?] is the erosion potential at each point along the profile
      // dVSlope [1] is the along-profile slope
      // m_dTimeStep * 3600 [s] is the time interval in seconds
      //
      // dRecession is horizontal recession (along the XY direction):
      //
      // dVRecessionXY[i] = (dForce * dVBeachProtection[i] * dVErosionPotentialFunc[i] * dVSlope[i] * m_dTimeStep * 3600) / dVR[i]
      //
      // XY recession must be -ve or zero. If it is +ve then it represents accretion not erosion, which must be described by a different set of equations. So we also need to constrain XY recession to be <= 0
      dVRecessionXY[i] = tMin(m_VCoast[nCoast].dGetWaveEnergyAtBreaking(nCoastPoint) * dVProfileErosionPotential[i] * dVConsSlope[i] / m_dR, 0.0);
      dVSCAPEXY[i] = VdProfileDistXY[i] - dVRecessionXY[i];
   }

   vector<double> dVChangeElevZ(nProfSize, 0);

   // We have calculated the XY-plane recession at every point on the profile, so now convert this to a change in Z-plane elevation at every inundated point on the profile (not the coast point). Again we use the elevation difference on the seaward side of 'this' point
   for (int i = 1; i < nProfSize - 1; i++)
   {
      // Vertical lowering dZ = dXY * tan(a), where tan(a) is the slope of the SCAPE profile in the XY direction
      double const dSCAPEHorizDist = dVSCAPEXY[i + 1] - dVSCAPEXY[i];
      double const dSCAPEVertDist = dVConsProfileZ[i] - dVConsProfileZ[i + 1];
      double const dSCAPESlope = dSCAPEVertDist / dSCAPEHorizDist;
      double dDeltaZ = dVRecessionXY[i] * dSCAPESlope;

      // Safety check: if thickness model has some jumps, dVConsProfileZ might be very high, limiting dSCAPESlope to 0 because all time erode a high fix quantity
      if (dSCAPESlope > 1)
         dDeltaZ = 0;

      int const nX = pProfile->pPtiVGetCellsInProfile()->at(i).nGetX();
      int const nY = pProfile->pPtiVGetCellsInProfile()->at(i).nGetY();

      // Store the local slope of the consolidated sediment, this is just for output display purposes
      m_pRasterGrid->m_Cell[nX][nY].SetLocalConsSedSlope(dVConsSlope[i]);

      // dDeltaZ is zero or -ve: if dDeltaZ is zero then do nothing, if -ve then remove some sediment from this cell
      if (dDeltaZ < 0)
      {
         // If there has already been potential erosion on this cell, then it must be a shared line segment (i.e. has co-incident profiles)
         double const dPrevPotentialErosion = -m_pRasterGrid->m_Cell[nX][nY].dGetPotentialPlatformErosion();
         if (dPrevPotentialErosion < 0)
         {
            // Average the two values
            LogStream << m_ulIter << ": [" << nX << "][" << nY << "] &&&& under profile " << pProfile->nGetProfileID() << " has previous potential platform erosion = " << dPrevPotentialErosion << endl;
            dDeltaZ = ((dDeltaZ + dPrevPotentialErosion) / 2);
         }

         // Constrain the lowering so we don't get negative slopes or +ve erosion amounts (dDeltaZ must be -ve), this is implicit in SCAPE
         dDeltaZ = tMax(dDeltaZ, -dVConsZDiff[i]);
         dDeltaZ = tMin(dDeltaZ, 0.0);
         dVChangeElevZ[i] = dDeltaZ;

         // Set the potential (unconstrained) erosion for this cell, is a +ve value
         m_pRasterGrid->m_Cell[nX][nY].SetPotentialPlatformErosion(-dDeltaZ);

         // Update this-timestep totals
         m_ulThisIterNumPotentialPlatformErosionCells++;
         m_dThisIterPotentialPlatformErosion -= dDeltaZ;       // Since dDeltaZ is a -ve value

         // Increment the check values
         m_ulTotPotentialPlatformErosionOnProfiles++;
         m_dTotPotentialPlatformErosionOnProfiles -= dDeltaZ;
      }

      // Finally, calculate the beach protection factor, this will be used in estimating actual (supply-limited) erosion
      double const dBeachProtectionFactor = dCalcBeachProtectionFactor(nX, nY, dBreakingWaveHeight);
      m_pRasterGrid->m_Cell[nX][nY].SetBeachProtectionFactor(dBeachProtectionFactor);
   }

   // If desired, save this coastline-normal consolidated-only profile data for checking purposes
   if (m_bOutputConsolidatedProfileData)
   {
      int const nRet = nSaveProfile(nCoast, pProfile, nProfSize, &VdProfileDistXY, &dVConsProfileZ, &dVProfileDepthOverDB, &dVProfileErosionPotential, &dVConsSlope, &dVRecessionXY, &dVChangeElevZ, pProfile->pPtiVGetCellsInProfile(), &dVSCAPEXY);
      if (nRet != RTN_OK)
         return nRet;
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Calculates potential platform erosion on cells to one side of a given coastline-normal profile, up to the next profile
//===============================================================================================================================
int CSimulation::nCalcPotentialPlatformErosionBetweenProfiles(int const nCoast, CGeomProfile *pProfile, int nDirection)
{
   // Only work on this profile if it is problem-free
   if (! pProfile->bProfileOK())
      return RTN_OK;

   int const nProfSize = pProfile->nGetNumCellsInProfile();
   int const nCoastProfileStart = pProfile->nGetCoastPoint();
   int const nProfileStartX = pProfile->pPtiVGetCellsInProfile()->at(0).nGetX();
   int const nProfileStartY = pProfile->pPtiVGetCellsInProfile()->at(0).nGetY();
   int const nCoastMax = m_VCoast[nCoast].nGetCoastlineSize();
   int nDistFromProfile = 0;
   int nParCoastXLast = nProfileStartX;
   int nParCoastYLast = nProfileStartY;

   // Start at the coast end of this coastline-normal profile, then move one cell (forward or backward) along the coast, then construct a parallel profile from this new coastline start cell. Calculate erosion along this parallel profile in the same way as above. Move another cell forward along the coastline, do the same. Keep going until both ends of the parallel profile have hit another profile
   bool bHitCoastEnd = false;
   bool bHitSeawardEnd = false;
   bool bEndProfileOnly = false;
   int nCellOK = 0;

   // Start the temporary parallel profile nDistFromProfile cells along the coastline from the coastline-normal profile, direction depending on nDirection
   int nThisPointOnCoast = nCoastProfileStart;
   while (true)
   {
      // Increment the distance from the start coast-normal profile
      nDistFromProfile++;

      if (nDirection == DIRECTION_DOWNCOAST)
         nThisPointOnCoast++;
      else
         nThisPointOnCoast--;

      // Is the coast end of the parallel profile at the beginning of the coast?
      if ((nDirection == DIRECTION_UPCOAST) && (nThisPointOnCoast < 0))
      {
         // TEST CODE ================
         pProfile = m_VCoast[nCoast].pGetProfileAtCoastPoint(0);
         bEndProfileOnly = true;
         nDirection = DIRECTION_DOWNCOAST;
         nDistFromProfile = 0;
         nThisPointOnCoast++;
         nParCoastXLast = INT_NODATA;
         nParCoastYLast = INT_NODATA;
         nCellOK = 0;

         LogStream << m_ulIter << ": CHANGED DIRECTION since temporary profile hit nThisPointOnCoast = " << nThisPointOnCoast << " while doing potential platform erosion " << (nDirection == DIRECTION_DOWNCOAST ? "down" : "up") << "-coast from profile = " << pProfile->nGetProfileID() << ", dist from profile = " <<  nDistFromProfile << endl;
         // TEST CODE ================

         // LogStream << m_ulIter << ": LEAVING LOOP since temporary profile hit nThisPointOnCoast = " << nThisPointOnCoast << " while doing potential platform erosion " << (nDirection == DIRECTION_DOWNCOAST ? "down" : "up") << "-coast from start profile = " << pProfile->nGetProfileID() << ", dist from profile = " <<  nDistFromProfile << endl;
         // break;
      }

      // Is the coast end of the parallel profile at the end of the coast?
      if ((nDirection == DIRECTION_DOWNCOAST) && (nThisPointOnCoast >= nCoastMax))
      {
         // TEST CODE ================
         pProfile = m_VCoast[nCoast].pGetProfileAtCoastPoint(nCoastMax-1);
         bEndProfileOnly = true;
         nDirection = DIRECTION_UPCOAST;
         nDistFromProfile = 0;
         nThisPointOnCoast--;
         nParCoastXLast = INT_NODATA;
         nParCoastYLast = INT_NODATA;
         nCellOK = 0;

         LogStream << m_ulIter << ": CHANGED DIRECTION since temporary profile hit nThisPointOnCoast = " << nThisPointOnCoast << " while doing potential platform erosion " << (nDirection == DIRECTION_DOWNCOAST ? "down" : "up") << "-coast from profile = " << pProfile->nGetProfileID() << ", dist from profile = " <<  nDistFromProfile << endl;
         // TEST CODE ================

         // LogStream << m_ulIter << ": LEAVING LOOP since temporary profile hit nThisPointOnCoast = " << nThisPointOnCoast << " while doing potential platform erosion " << (nDirection == DIRECTION_DOWNCOAST ? "down" : "up") << "-coast from profile = " << pProfile->nGetProfileID() << ", dist from profile = " <<  nDistFromProfile << endl;
         // break;
      }

      LogStream << m_ulIter << ": from start profile " << pProfile->nGetProfileID() << ", doing potential platform erosion " << (nDirection == DIRECTION_DOWNCOAST ? "down" : "up") << "-coast, dist from profile = " <<  nDistFromProfile << endl;

      double dDepthOfBreaking = m_VCoast[nCoast].dGetDepthOfBreaking(nThisPointOnCoast);

      if (bFPIsEqual(dDepthOfBreaking, DBL_NODATA, TOLERANCE))
      {
         // This parallel profile is not in the active zone, so no platform erosion here. Move on to the next point along the coastline in this direction
         LogStream << m_ulIter << ": temporary profile not in active zone at coastline " << nCoast << " coast point " << nThisPointOnCoast << " when constructing parallel profile for potential platform erosion. Working from profile " << pProfile->nGetProfileID() << ", " << (nDirection == DIRECTION_DOWNCOAST ? "down" : "up") << "-coast, dist from profile = " << nDistFromProfile << endl;
         continue;
      }

      LogStream << m_ulIter << ": between profiles " << (nDirection == DIRECTION_DOWNCOAST ? "down" : "up") << "-coast from start profile " << pProfile->nGetProfileID() << " nThisPointOnCoast = " << nThisPointOnCoast << " dDepthOfBreaking = " << dDepthOfBreaking << endl;

      // All is OK, so get the grid coordinates of this point, which is the coastline start point for the parallel profile
      int const nParCoastX = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nThisPointOnCoast)->nGetX();
      int const nParCoastY = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nThisPointOnCoast)->nGetY();

      if ((nParCoastX == nParCoastXLast) && (nParCoastY == nParCoastYLast))
      {
         // Should not happen, but can due to rounding errors
         LogStream << WARN << m_ulIter << ": rounding problem on coast " << nCoast << " profile " << pProfile->nGetProfileID() << " at [" << nParCoastX << "][" << nParCoastY << "]" << endl;

         // So move on to the next point along the coastline in this direction
         continue;
      }

      // Is this coastline start point the start point of an adjacent coastline-normal vector?
      if ((nDistFromProfile > 1) && m_pRasterGrid->m_Cell[nParCoastX][nParCoastY].bIsProfile())
      {
         LogStream << m_ulIter << ": coast " << nCoast << ", temporary profile parallel to start profile " << pProfile->nGetProfileID() << " hit another profile at coast end, nThisPointOnCoast = " << nThisPointOnCoast << " while doing potential platform erosion " << (nDirection == DIRECTION_DOWNCOAST ? "down" : "up") << "-coast from profile = " << pProfile->nGetProfileID() << ", dist from start profile = " <<  nDistFromProfile << endl;
         bHitCoastEnd = true;
      }

      // Get the height of the associated breaking wave from the coast point: this height is used in beach protection calcs. Note that it will be DBL_NODATA if not in active zone
      double const dBreakingWaveHeight = m_VCoast[nCoast].dGetBreakingWaveHeight(nThisPointOnCoast);

      // OK, now construct a parallel profile
      vector<CGeom2DIPoint> PtiVGridParProfile; // Integer coords (grid CRS) of cells under the parallel profile
      vector<CGeom2DPoint> PtVExtCRSParProfile; // coordinates (external CRS) of cells under the parallel profile

      ConstructParallelProfile(nProfileStartX, nProfileStartY, nParCoastX, nParCoastY, nProfSize, pProfile->pPtiVGetCellsInProfile(), &PtiVGridParProfile, &PtVExtCRSParProfile);

      int const nParProfSize = static_cast<int>(PtiVGridParProfile.size());

      // We have a parallel profile which starts at the coast, but is it long enough to be useful? May have been cut short because it extended outside the grid, or we hit an adjacent profile
      if (nParProfSize < 3)
      {
         // We cannot use this parallel profile, it is too short to calculate along-profile slope, so abandon it and move on to the next parallel profile in this direction
         nParCoastXLast = nParCoastX;
         nParCoastYLast = nParCoastY;
         LogStream << m_ulIter << ": temporary profile from start profile " << pProfile->nGetProfileID() << " abandoned since too short, starts at [" << nParCoastX << "][" << nParCoastY << "] coastline point " << nThisPointOnCoast << ", length = " << nParProfSize << endl;
         continue;
      }

      // Now see whether the seaward end of the temporary profile has hit another profile
      int const nXParEnd = PtiVGridParProfile[nParProfSize-1].nGetX();
      int const nYParEnd = PtiVGridParProfile[nParProfSize-1].nGetY();
      if (m_pRasterGrid->m_Cell[nXParEnd][nYParEnd].bIsProfile())
      {
         LogStream << m_ulIter << ": coast " << nCoast << ", temporary profile parallel to start profile " << pProfile->nGetProfileID() << " hit another profile at seaward end while doing potential platform erosion " << (nDirection == DIRECTION_DOWNCOAST ? "down" : "up") << "-coast from profile = " << pProfile->nGetProfileID() << ", dist from start profile = " <<  nDistFromProfile << endl;
         bHitSeawardEnd = true;
      }

      // Calculate potential erosion along the parallel profile. First calculate the length of the parallel profile in external CRS units
      double const dParProfileLenXY = dGetDistanceBetween(&PtVExtCRSParProfile[0], &PtVExtCRSParProfile[nParProfSize - 1]);

      // Next calculate the distance between profile points, again in external CRS units. Assume that the sample points are equally spaced along the parallel profile (not quite true)
      double dParSpacingXY = dParProfileLenXY / (nParProfSize - 1);

      // Safety check
      if (bFPIsEqual(dParSpacingXY, 0.0, TOLERANCE))
         dParSpacingXY = TOLERANCE;

      // LogStream << "dParSpacingXY = " << dParSpacingXY << endl;

      vector<double> dVParProfileZ(nParProfSize, 0);      // Initial (pre-erosion) elevation of both consolidated and unconsolidated sediment for cells 'under' the parallel profile
      vector<double> dVParProfileDistXY(nParProfSize, 0); // Along-profile distance measured from the coast, in external CRS units
      vector<double> dVParConsProfileZ(nParProfSize, 0);  // Initial (pre-erosion) elevation of consolidated sediment only for cells 'under' the parallel profile
      vector<double> dVParConsZDiff(nParProfSize, 0);
      vector<double> dVParConsSlope(nParProfSize, 0);

      for (int i = 0; i < nParProfSize; i++)
      {
         int const nXPar = PtiVGridParProfile[i].nGetX();
         int const nYPar = PtiVGridParProfile[i].nGetY();

//          // If this cell has already had platform erosion this timestep, do nothing
//          if (m_pRasterGrid->m_Cell[nXPar][nYPar].bGetPlatformErosionThisIter())
//             continue;
//
//          // Is this isn't a sea cell, do nothing
//          if (! m_pRasterGrid->m_Cell[nXPar][nYPar].bIsInundated())
//          {
//             // It isn't so move along, nothing to do here
//             // LogStream << m_ulIter << " : [" << nXPar << "][" << nYPar << "] is not inundated" << endl;
//             continue;
//          }
//
//          // Is this cell in a polygon?
//          int const nPolyID = m_pRasterGrid->m_Cell[nXPar][nYPar].nGetPolygonID();
//          if (nPolyID == INT_NODATA)
//          {
//             // It isn't. This can happen at the seaward end of polygons TODO 026 Is it a problem?
//             // LogStream << m_ulIter << " : [" << nXPar << "][" << nYPar << "] = {" << dGridCentroidXToExtCRSX(nXPar) << ", " << dGridCentroidYToExtCRSY(nYPar) << "} is not in a polygon" << endl;
//             continue;
//          }

         // Get the number of the highest layer with non-zero thickness
         int const nTopLayer = m_pRasterGrid->m_Cell[nXPar][nYPar].nGetTopNonZeroLayerAboveBasement();

         // Safety check
         if (nTopLayer == INT_NODATA)
            return RTN_ERR_NO_TOP_LAYER_DURING_PLATFORM_CALC;

         if (nTopLayer == NO_NONZERO_THICKNESS_LAYERS)
            // TODO 025 We are down to basement
            return RTN_OK;

         // Get the elevation for consolidated sediment only on this cell
         dVParConsProfileZ[i] = m_pRasterGrid->m_Cell[nXPar][nYPar].dGetConsSedTopElevForLayerAboveBasement(nTopLayer);

         // Get the elevation for both consolidated and unconsolidated sediment on this cell (ignore any talus)
         dVParProfileZ[i] = m_pRasterGrid->m_Cell[nXPar][nYPar].dGetAllSedTopElevOmitTalus();

         // And store the X-Y plane distance from the start of the profile
         dVParProfileDistXY[i] = i * dParSpacingXY;
      }

      for (int i = 0; i < nParProfSize - 1; i++)
      {
         // For the consolidated-only profile, get the Z differences (already in external CRS units)
         dVParConsZDiff[i] = dVParConsProfileZ[i] - dVParConsProfileZ[i + 1];

         // Calculate dZ/dXY, the Z slope (i.e. rise over run) in the XY direction. Note that we use the elevation difference on the seaward side of 'this' point
         dVParConsSlope[i] = dVParConsZDiff[i] / dParSpacingXY;
      }

      // Sort out the final slope value
      dVParConsSlope[nParProfSize - 1] = dVParConsSlope[nParProfSize - 2];

      if (m_nProfileSmoothWindow > 0)
      {
         // Smooth the vector of slopes for the consolidated-only profile
         dVParConsSlope = dVSmoothProfileSlope(&dVParConsSlope);
      }

      // Initialize the parallel profile vector with depth / m_dWaveBreakingDepth
      vector<double> dVParProfileDepthOverDB(nParProfSize, 0);      // Depth / wave breaking depth at the parallel profile sample points
      vector<double> dVParProfileErosionPotential(nParProfSize, 0); // Erosion potential at the parallel profile sample points

      // Calculate the erosion potential along this profile using the shape function which we read in earlier
      double dTotalErosionPotential = 0;

      // Safety check TODO 061 Is this safety check to depth of breaking a reasonable thing to do?
      if (dDepthOfBreaking <= 0.0)
         dDepthOfBreaking = 1e-10;

      for (int i = 0; i < nParProfSize; i++)
      {
         // Use the actual depth of water here (i.e. the depth to the top of the unconsolidated sediment, including the thickness of consolidated sediment beneath it)
         dVParProfileDepthOverDB[i] = m_dThisIterSWL - dVParProfileZ[i];
         dVParProfileDepthOverDB[i] /= dDepthOfBreaking;

         // Constrain dDepthOverDB[i] to be between 0 (can get small -ve due to rounding errors) and m_dDepthOverDBMax
         dVParProfileDepthOverDB[i] = tMax(dVParProfileDepthOverDB[i], 0.0);
         dVParProfileDepthOverDB[i] = tMin(dVParProfileDepthOverDB[i], m_dDepthOverDBMax);

         // And then use the look-up table to find the value of erosion potential at this point on the profile
         dVParProfileErosionPotential[i] = dLookUpErosionPotential(dVParProfileDepthOverDB[i]);

         // If erosion potential (a -ve value) is tiny, set it to zero
         if (dVParProfileErosionPotential[i] > -SED_ELEV_TOLERANCE)
            dVParProfileErosionPotential[i] = 0;

         // Keep track of the total erosion potential for this profile
         dTotalErosionPotential += dVParProfileErosionPotential[i];
      }

      // Constrain erosion potential at every point on the profile, so that the integral of erosion potential on the whole profile is unity (Walkden and Hall 2005). Note that here, erosion potential is -ve so we must constrain to -1
      for (int i = 0; i < nParProfSize; i++)
      {
         if (dTotalErosionPotential < 0)
            dVParProfileErosionPotential[i] /= (-dTotalErosionPotential);
      }

      vector<double> dVParRecessionXY(nParProfSize, 0);
      vector<double> dVParSCAPEXY(nParProfSize, 0);

      // Calculate recession at every point on the parallel profile
      for (int i = 0; i < nParProfSize; i++)
      {
         // dRecession = dForce * (dBeachProtection / dR) * dErosionPotential * dSlope * dTime
         // where:
         // dVRecession [m] is the landward migration distance defined in the profile relative (XY) CRS
         // dForce is given by Equation 4 in Walkden & Hall, 2005
         // dVBeachProtection [1] is beach protection factor [1, 0] = [no protection, fully protected] (This is calculated later, see dCalcBeachProtectionFactor())
         // dVR  [m^(9/4)s^(2/3)] is the material strength and some hydrodynamic constant
         // dVProfileErosionPotential [?] is the erosion potential at each point along the profile
         // dVSlope [1] is the along-profile slope
         // m_dTimeStep * 3600 [s] is the time interval in seconds
         //
         // dRecession is horizontal recession (along the XY direction):
         //
         // dVRecessionXY[i] = (dForce * dVBeachProtection[i] * dVErosionPotentialFunc[i] * dVSlope[i] * m_dTimeStep * 3600) / dVR[i]
         //
         //
         // XY recession must be -ve or zero. If it is +ve then it represents accretion not erosion, which must be described by a different set of equations. So we also need to constrain XY recession to be <= 0
         dVParRecessionXY[i] = tMin(m_VCoast[nCoast].dGetWaveEnergyAtBreaking(nThisPointOnCoast) * dVParProfileErosionPotential[i] * dVParConsSlope[i] / m_dR, 0.0);
         dVParSCAPEXY[i] = dVParProfileDistXY[i] - dVParRecessionXY[i];

         // LogStream << m_ulIter << ": [" << nXPar << "][" << nYPar << "] = {" << dGridCentroidXToExtCRSX(nXPar) << ", " <<  dGridCentroidYToExtCRSY(nYPar) << "} wave energy = " << m_VCoast[nCoast].dGetWaveEnergyAtBreaking(nThisPointOnCoast) << " erosion potential = " << dVParProfileErosionPotential[i] << " slope = " << dVParProfileSlope[i] << " dVParZDiff[i] = " << dVParZDiff[i] << " nParProfSize = " << nParProfSize << endl;
      }

      vector<double> dVParDeltaZ(nParProfSize, 0);

      // We have calculated the XY-plane recession at every point on the profile, so now convert this to a change in Z-plane elevation at every inundated point on the profile (not the coast point). Again we use the elevation difference on the seaward side of 'this' point
      for (int i = 1; i < nParProfSize - 1; i++)
      {
         int const nXPar = PtiVGridParProfile[i].nGetX();
         int const nYPar = PtiVGridParProfile[i].nGetY();

         // If this cell has already had platform erosion this timestep, do nothing
         if (m_pRasterGrid->m_Cell[nXPar][nYPar].bGetPlatformErosionThisIter())
         {
            continue;
         }

         // Is this isn't a sea cell, do nothing
         if (! m_pRasterGrid->m_Cell[nXPar][nYPar].bIsInundated())
         {
            // LogStream << m_ulIter << " : [" << nXPar << "][" << nYPar << "] is not inundated" << endl;
            continue;
         }

         // Is this cell in a polygon?
         int const nPolyID = m_pRasterGrid->m_Cell[nXPar][nYPar].nGetPolygonID();
         if (nPolyID == INT_NODATA)
         {
            // It isn't. This can happen at the seaward end of polygons TODO 026 Is it a problem?
            // LogStream << m_ulIter << " : [" << nXPar << "][" << nYPar << "] = {" << dGridCentroidXToExtCRSX(nXPar) << ", " << dGridCentroidYToExtCRSY(nYPar) << "} is not in a polygon" << endl;
            continue;
         }

         // This cell is OK, so we will calculate potential platform erosion here
         nCellOK++;

         // Vertical lowering dZ = dXY * tan(a), where tan(a) is the slope of the SCAPE profile in the XY direction
         double const dSCAPEHorizDist = dVParSCAPEXY[i + 1] - dVParSCAPEXY[i];

         // Safety check
         if (bFPIsEqual(dSCAPEHorizDist, 0.0, TOLERANCE))
            continue;

         double const dSCAPEVertDist = dVParConsProfileZ[i] - dVParConsProfileZ[i + 1];
         double const dSCAPESlope = dSCAPEVertDist / dSCAPEHorizDist;
         double dDeltaZ = dVParRecessionXY[i] * dSCAPESlope;

         // Safety check: if thickness model has some jumps, dVConsProfileZ might be very high, limiting dSCAPESlope to 0
         if (dSCAPESlope > 1)
            dDeltaZ = 0;

         // Store the local slope of the consolidated sediment, this is just for output display purposes
         LogStream << m_ulIter << ": [" << nXPar << "][" << nYPar << "] local slope = " << dVParConsSlope[i] << endl;
         m_pRasterGrid->m_Cell[nXPar][nYPar].SetLocalConsSedSlope(dVParConsSlope[i]);

         // dDeltaZ is zero or -ve: if dDeltaZ is zero then do nothing, if -ve then remove some sediment from this cell
         if (dDeltaZ < 0)
         {
            // // Has this cell already had platform erosion during this timestep?
            // if (m_pRasterGrid->m_Cell[nXPar][nYPar].bPotentialPlatformErosion())
            // {
            //    // It has
            //    double const dPrevPotentialErosion = -m_pRasterGrid->m_Cell[nXPar][nYPar].dGetPotentialPlatformErosion();
            //
            //    // LogStream << m_ulIter << ": [" << nXPar << "][" << nYPar << "] parallel profile " << nDistFromProfile << " coast points " << (nDirection == DIRECTION_DOWNCOAST ? "down" : "up") << "-coast from profile " << nProfile << " has previous potential platform erosion = " << dPrevPotentialErosion << ", current potential platform erosion = " << dDeltaZ << ", max value = " << tMin(dPrevPotentialErosion, dDeltaZ) << endl;
            //
            //    // Use the larger of the two -ve values
            //    dDeltaZ = tMin(dPrevPotentialErosion, dDeltaZ);
            //
            //    // Adjust this-timestep totals, since this cell has already been eroded
            //    m_ulThisIterNumPotentialPlatformErosionCells--;
            //    m_dThisIterPotentialPlatformErosion += dPrevPotentialErosion; // Since dPrevPotentialErosion is +ve
            //    // assert(isfinite(m_dThisIterPotentialPlatformErosion));
            //    // assert(m_dThisIterPotentialPlatformErosion >= 0);
            //
            //    // And also adjust the check values
            //    m_ulTotPotentialPlatformErosionBetweenProfiles--;
            //    m_dTotPotentialPlatformErosionBetweenProfiles += dPrevPotentialErosion; // Since -ve
            // }

            // Constrain the lowering so we don't get negative slopes or +ve erosion amounts (dDeltaZ must be -ve), this is implicit in SCAPE
            dDeltaZ = tMax(dDeltaZ, -dVParConsZDiff[i]);
            dDeltaZ = tMin(dDeltaZ, 0.0);
            dVParDeltaZ[i] = dDeltaZ;

            // Set the potential (unconstrained) erosion for this cell, it is a +ve value
            m_pRasterGrid->m_Cell[nXPar][nYPar].SetPotentialPlatformErosion(-dDeltaZ);
            // LogStream << "[" << nXPar << "][" << nYPar << "] = {" << dGridCentroidXToExtCRSX(nXPar) << ", " <<  dGridCentroidYToExtCRSY(nYPar) << "} has potential platform erosion = " << -dDeltaZ << endl;

            // Update this-timestep totals
            m_ulThisIterNumPotentialPlatformErosionCells++;
            m_dThisIterPotentialPlatformErosion -= dDeltaZ; // Since dDeltaZ is a -ve value
            // assert(isfinite(m_dThisIterPotentialPlatformErosion));
            // assert(m_dThisIterPotentialPlatformErosion >= 0);

            // Increment the check values
            m_ulTotPotentialPlatformErosionBetweenProfiles++;
            m_dTotPotentialPlatformErosionBetweenProfiles -= dDeltaZ; // Since -ve
         }

         // Finally, calculate the beach protection factor, this will be used in estimating actual (supply-limited) erosion
         double const dBeachProtectionFactor = dCalcBeachProtectionFactor(nXPar, nYPar, dBreakingWaveHeight);
         m_pRasterGrid->m_Cell[nXPar][nYPar].SetBeachProtectionFactor(dBeachProtectionFactor);
      }

      // If desired, save this parallel coastline-normal profile for checking purposes
      if (m_bOutputParallelProfileData)
      {
         int const nRet = nSaveParProfile(nCoast, pProfile, nParProfSize, nDirection, nDistFromProfile, &dVParProfileDistXY, &dVParConsProfileZ, &dVParProfileDepthOverDB, &dVParProfileErosionPotential, &dVParConsSlope, &dVParRecessionXY, &dVParDeltaZ, pProfile->pPtiVGetCellsInProfile(), &dVParSCAPEXY);

         if (nRet != RTN_OK)
            return nRet;
      }

      // If this temporary profile, or a previous temporary profile, has hit another profile at both the coast and seaward end, then quit
      if (bHitCoastEnd && bHitSeawardEnd)
      {
         LogStream << m_ulIter << ": coast " << nCoast << ", temporary parallel profiles from start profile " << pProfile->nGetProfileID() << " hit another profile at both coast and seaward end while doing potential platform erosion " << (nDirection == DIRECTION_DOWNCOAST ? "down" : "up") << "-coast from profile = " << pProfile->nGetProfileID() << ", dist from start profile = " <<  nDistFromProfile << endl;
         break;
      }

      // If we did not calculate potential platform on any cells while processing this temporary profile, then quit
      if (nCellOK == 0)
      {
         LogStream << m_ulIter << ": coast " << nCoast << ", no cells have potential platform erosion while processing temporary parallel profile from start profile " << pProfile->nGetProfileID() << " " << (nDirection == DIRECTION_DOWNCOAST ? "down" : "up") << "-coast from profile = " << pProfile->nGetProfileID() << ", dist from start profile = " <<  nDistFromProfile << endl;
         break;
      }

      // Update for next time round the loop
      nParCoastXLast = nParCoastX;
      nParCoastYLast = nParCoastY;
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Calculates actual (constrained by available sediment) erosion of the consolidated shore platform on a single sea cell
//===============================================================================================================================
void CSimulation::DoActualPlatformErosionOnCell(int const nX, int const nY)
{
   // LogStream << m_ulIter << ": doing platform erosion on cell [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;

   // Get the beach protection factor, which quantifies the extent to which unconsolidated sediment on the shore platform (beach) protects the shore platform
   double const dBeachProtectionFactor = m_pRasterGrid->m_Cell[nX][nY].dGetBeachProtectionFactor();

   if (bFPIsEqual(dBeachProtectionFactor, 0.0, TOLERANCE))
      // The beach is sufficiently thick to prevent any platform erosion (or we are down to basement)
      return;

   // Get the potential depth of potential erosion, considering beach protection
   double const dThisPotentialErosion = m_pRasterGrid->m_Cell[nX][nY].dGetPotentialPlatformErosion() * dBeachProtectionFactor;

   // We will be eroding the topmost layer that has non-zero thickness
   int const nThisLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopNonZeroLayerAboveBasement();

   // Safety check
   if (nThisLayer == INT_NODATA)
   {
      cerr << ERR << "no sediment layer while doing actual shore platform erosion" << endl;
      return;
   }

   if (nThisLayer == NO_NONZERO_THICKNESS_LAYERS)
      // No layer with non-zero thickness left, we are down to basement
      return;

   // OK, we have a layer that can be eroded so find out how much consolidated sediment we have available on this cell
   double const dExistingAvailableFine = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetConsolidatedSediment()->dGetFineDepth();
   double const dExistingAvailableSand = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetConsolidatedSediment()->dGetSandDepth();
   double const dExistingAvailableCoarse = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetConsolidatedSediment()->dGetCoarseDepth();

   // Now partition the total lowering for this cell between the three size fractions: do this by relative erodibility
   int const nFineWeight = (dExistingAvailableFine > 0 ? 1 : 0);
   int const nSandWeight = (dExistingAvailableSand > 0 ? 1 : 0);
   int const nCoarseWeight = (dExistingAvailableCoarse > 0 ? 1 : 0);

   double const dTotErodibility = (nFineWeight * m_dFineErodibilityNormalized) + (nSandWeight * m_dSandErodibilityNormalized) + (nCoarseWeight * m_dCoarseErodibilityNormalized);
   double dTotActualErosion = 0;
   double dSandEroded = 0;
   double dCoarseEroded = 0;

   if (nFineWeight)
   {
      // Erode some fine-sized consolidated sediment
      double const dFineLowering = (m_dFineErodibilityNormalized * dThisPotentialErosion) / dTotErodibility;

      // Make sure we don't get -ve amounts left on the cell
      double const dFineEroded = tMin(dExistingAvailableFine, dFineLowering);
      double const dRemaining = dExistingAvailableFine - dFineEroded;

      dTotActualErosion += dFineEroded;

      // Set the value for this layer
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetConsolidatedSediment()->SetFineDepth(dRemaining);

      // And set the changed-this-timestep switch
      m_bConsSedChangedThisIter[nThisLayer] = true;

      // And increment the per-timestep total, also add to the suspended sediment load
      m_dThisIterActualPlatformErosionFineCons += dFineEroded;
      m_dThisIterFineSedimentToSuspension += dFineEroded;
   }

   if (nSandWeight)
   {
      // Erode some sand-sized consolidated sediment
      double const dSandLowering = (m_dSandErodibilityNormalized * dThisPotentialErosion) / dTotErodibility;

      // Make sure we don't get -ve amounts left on the source cell
      dSandEroded = tMin(dExistingAvailableSand, dSandLowering);
      double const dRemaining = dExistingAvailableSand - dSandEroded;

      dTotActualErosion += dSandEroded;

      // Set the new value of sand consolidated sediment depth for this layer
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetConsolidatedSediment()->SetSandDepth(dRemaining);

      // And add this to the depth of sand unconsolidated sediment for this layer
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetUnconsolidatedSediment()->AddSandDepth(dSandEroded);

      // Set the changed-this-timestep switch
      m_bConsSedChangedThisIter[nThisLayer] = true;

      // And increment the per-timestep total
      m_dThisIterActualPlatformErosionSandCons += dSandEroded;
   }

   if (nCoarseWeight)
   {
      // Erode some coarse-sized consolidated sediment
      double const dCoarseLowering = (m_dCoarseErodibilityNormalized * dThisPotentialErosion) / dTotErodibility;

      // Make sure we don't get -ve amounts left on the source cell
      dCoarseEroded = tMin(dExistingAvailableCoarse, dCoarseLowering);
      double const dRemaining = dExistingAvailableCoarse - dCoarseEroded;

      dTotActualErosion += dCoarseEroded;

      // Set the new value of coarse consolidated sediment depth for this layer
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetConsolidatedSediment()->SetCoarseDepth(dRemaining);

      // And add this to the depth of coarse unconsolidated sediment for this layer
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetUnconsolidatedSediment()->AddCoarseDepth(dCoarseEroded);

      // Set the changed-this-timestep switch
      m_bConsSedChangedThisIter[nThisLayer] = true;

      // And increment the per-timestep total
      m_dThisIterActualPlatformErosionCoarseCons += dCoarseEroded;
   }

   // Did we erode anything?
   if (dTotActualErosion > 0)
   {
      // We did, so set the actual erosion value for this cell
      m_pRasterGrid->m_Cell[nX][nY].SetActualPlatformErosion(dTotActualErosion);

      // Recalculate the elevation of every layer
      m_pRasterGrid->m_Cell[nX][nY].CalcAllLayerElevsAndD50();

      // And update the cell's sea depth
      m_pRasterGrid->m_Cell[nX][nY].SetSeaDepth();

      // Update per-timestep totals
      m_ulThisIterNumActualPlatformErosionCells++;

      // Add eroded sand/coarse sediment for this cell to the polygon that contains the cell, ready for redistribution during beach erosion/deposition (fine sediment has already been dealt with)
      int nPolyID = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonID();
      int nPolyCoastID = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonCoastID();

      if (nPolyID == INT_NODATA)
      {
         // Can get occasional problems with polygon rasterization near the coastline, so also search the eight adjacent cells
         array<int, 8> nDirection = {NORTH, NORTH_EAST, EAST, SOUTH_EAST, SOUTH, SOUTH_WEST, WEST, NORTH_WEST};
         shuffle(nDirection.begin(), nDirection.end(), m_Rand[0]);

         for (int n = 0; n < 8; n++)
         {
            int nXAdj;
            int nYAdj;

            if (nDirection[n] == NORTH)
            {
               nXAdj = nX;
               nYAdj = nY - 1;

               if (nYAdj >= 0)
               {
                  nPolyID = m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetPolygonID();

                  if (nPolyID != INT_NODATA)
                  {
                     nPolyCoastID = m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetPolygonCoastID();
                     break;
                  }
               }
            }
            else if (nDirection[n] == NORTH_EAST)
            {
               nXAdj = nX + 1;
               nYAdj = nY - 1;

               if ((nXAdj < m_nXGridSize) && (nYAdj >= 0))
               {
                  nPolyID = m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetPolygonID();

                  if (nPolyID != INT_NODATA)
                  {
                     nPolyCoastID = m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetPolygonCoastID();
                     break;
                  }
               }
            }
            else if (nDirection[n] == EAST)
            {
               nXAdj = nX + 1;
               nYAdj = nY;

               if (nXAdj < m_nXGridSize)
               {
                  nPolyID = m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetPolygonID();

                  if (nPolyID != INT_NODATA)
                  {
                     nPolyCoastID = m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetPolygonCoastID();
                     break;
                  }
               }
            }
            else if (nDirection[n] == SOUTH_EAST)
            {
               nXAdj = nX + 1;
               nYAdj = nY + 1;

               if ((nXAdj < m_nXGridSize) && (nYAdj < m_nYGridSize))
               {
                  nPolyID = m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetPolygonID();

                  if (nPolyID != INT_NODATA)
                  {
                     nPolyCoastID = m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetPolygonCoastID();
                     break;
                  }
               }
            }
            else if (nDirection[n] == SOUTH)
            {
               nXAdj = nX;
               nYAdj = nY + 1;

               if (nYAdj < m_nYGridSize)
               {
                  nPolyID = m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetPolygonID();

                  if (nPolyID != INT_NODATA)
                  {
                     nPolyCoastID = m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetPolygonCoastID();
                     break;
                  }
               }
            }
            else if (nDirection[n] == SOUTH_WEST)
            {
               nXAdj = nX - 1;
               nYAdj = nY + 1;

               if ((nXAdj >= 0) && (nXAdj < m_nXGridSize))
               {
                  nPolyID = m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetPolygonID();

                  if (nPolyID != INT_NODATA)
                  {
                     nPolyCoastID = m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetPolygonCoastID();
                     break;
                  }
               }
            }
            else if (nDirection[n] == WEST)
            {
               nXAdj = nX - 1;
               nYAdj = nY;

               if (nXAdj >= 0)
               {
                  nPolyID = m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetPolygonID();

                  if (nPolyID != INT_NODATA)
                  {
                     nPolyCoastID = m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetPolygonCoastID();
                     break;
                  }
               }
            }
            else if (nDirection[n] == NORTH_WEST)
            {
               nXAdj = nX - 1;
               nYAdj = nY - 1;

               if ((nXAdj >= 0) && (nYAdj >= 0))
               {
                  nPolyID = m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetPolygonID();

                  if (nPolyID != INT_NODATA)
                  {
                     nPolyCoastID = m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetPolygonCoastID();
                     break;
                  }
               }
            }
         }
      }

      // assert(nPolyID < m_VCoast[0].nGetNumPolygons());

      // Safety check
      if (nPolyID == INT_NODATA)
      {
         // Uh-oh, we have a problem
         if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
            LogStream << m_ulIter << ": " << WARN << "platform erosion on cell [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "} but this is not in a polygon" << endl;

         // m_dDepositionSandDiff and m_dDepositionCoarseDiff are both +ve
         m_dDepositionSandDiff += dSandEroded;
         m_dDepositionCoarseDiff += dCoarseEroded;

         return;
      }

      // All OK, so add this to the polygon's total of unconsolidated sand/coarse sediment, to be deposited or moved later. These values are +ve (deposition)
      m_VCoast[nPolyCoastID].pGetPolygon(nPolyID)->AddPlatformErosionUnconsSand(dSandEroded);
      m_VCoast[nPolyCoastID].pGetPolygon(nPolyID)->AddPlatformErosionUnconsCoarse(dCoarseEroded);
   }
}

//===============================================================================================================================
//! Creates a look-up table for erosion potential, given depth over DB
//===============================================================================================================================
bool CSimulation::bCreateErosionPotentialLookUp(vector<double> *VdDepthOverDBIn, vector<double> *VdErosionPotentialIn, vector<double> *VdErosionPotentialFirstDerivIn)
{
   // Set up a temporary vector to hold the incremental DepthOverDB values
   vector<double> VdDepthOverDB;
   double dTempDOverDB = 0;

   while (dTempDOverDB <= 1.1) // Arbitrary max value, we will adjust this later
   {
      VdDepthOverDB.push_back(dTempDOverDB); // These are the incremental sample values of DepthOverDB
      dTempDOverDB += DEPTH_OVER_DB_INCREMENT;

      m_VdErosionPotential.push_back(0); // This will hold the corresponding value of erosion potential for each sample point
   }

   int const nSize = static_cast<int>(VdDepthOverDB.size());
   vector<double> VdDeriv(nSize, 0);   // First derivative at the sample points: calculated by the spline function but not subsequently used
   vector<double> VdDeriv2(nSize, 0.); // Second derivative at the sample points, ditto
   vector<double> VdDeriv3(nSize, 0.); // Third derivative at the sample points, ditto

   // Calculate the value of erosion potential (is a -ve value) for each of the sample values of DepthOverDB, and store it for use in the look-up function
   hermite_cubic_spline_value(static_cast<int>(VdDepthOverDBIn->size()), &(VdDepthOverDBIn->at(0)), &(VdErosionPotentialIn->at(0)), &(VdErosionPotentialFirstDerivIn->at(0)), nSize, &(VdDepthOverDB[0]), &(m_VdErosionPotential[0]), &(VdDeriv[0]), &(VdDeriv2[0]), &(VdDeriv3[0]));

   // Tidy the erosion potential look-up data: cut off values (after the first) for which erosion potential is no longer -ve
   int nLastVal = -1;

   for (int n = 1; n < nSize - 1; n++)
   {
      if (m_VdErosionPotential[n] > 0)
      {
         nLastVal = n;
         break;
      }
   }

   if (nLastVal > 0)
   {
      // Erosion potential is no longer -ve at this value of DepthOverDB, so set the maximum value of DepthOverDB that will be used in the simulation (any DepthOverDB value greater than this produces zero erosion potential)
      m_dDepthOverDBMax = VdDepthOverDB[nLastVal];
      m_VdErosionPotential.erase(m_VdErosionPotential.begin() + nLastVal + 1, m_VdErosionPotential.end());
      m_VdErosionPotential.back() = 0;
   }

   else
      // Erosion potential is unbounded, i.e. it is still -ve when we have reached the end of the look-up vector
      return false;

   // All OK
   return true;
}

//===============================================================================================================================
//! The erosion potential lookup: it returns a value for erosion potential given a value of Depth Over DB
//===============================================================================================================================
double CSimulation::dLookUpErosionPotential(double const dDepthOverDB)
{
   // If dDepthOverDB exceeds the maximum value which we calculated earlier, erosion potential is zero
   if (dDepthOverDB > m_dDepthOverDBMax)
      return 0;

   // OK, dDepthOverDB is less than the maximum so look up a corresponding value for erosion potential. The look-up index is dDepthOverDB divided by (the Depth Over DB increment used when creating the look-up vector). But since this look-up index may not be an integer, split the look-up index into integer and fractional parts and deal with each separately
   double const dErosionPotential = dGetInterpolatedValue(&m_VdDepthOverDB, &m_VdErosionPotential, dDepthOverDB, false);

   return dErosionPotential;
}

//===============================================================================================================================
//! Calculates the (inverse) beach protection factor as in SCAPE: 0 is fully protected, 1 = no protection
//===============================================================================================================================
double CSimulation::dCalcBeachProtectionFactor(int const nX, int const nY, double const dBreakingWaveHeight)
{
   // Safety check
   if (bFPIsEqual(dBreakingWaveHeight, DBL_NODATA, TOLERANCE))
      return 0;

   // We are considering the unconsolidated sediment (beach) of the topmost layer that has non-zero thickness
   int const nThisLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopNonZeroLayerAboveBasement();

   // Safety check
   if (nThisLayer == INT_NODATA)
   {
      cerr << ERR << "no sediment layer in dCalcBeachProtectionFactor()" << endl;
      return 0;
   }

   if (nThisLayer == NO_NONZERO_THICKNESS_LAYERS)
      // There are no layers with non-zero thickness left (i.e. we are down to basement) so no beach protection
      return 0;

   // In SCAPE, 0.23 * the significant breaking wave height is assumed to be the maximum depth of beach that waves can penetrate to erode a platform. For depths less than this, the beach protective ability is assumed to vary linearly
   double const dBeachDepth = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->dGetAllUnconsDepth();
   double const dMaxPenetrationDepth = BEACH_PROTECTION_HB_RATIO * dBreakingWaveHeight;
   double dFactor = 0;

   if (dMaxPenetrationDepth > 0)
      dFactor = tMax(1 - (dBeachDepth / dMaxPenetrationDepth), 0.0);

   // LogStream << m_ulIter << ": cell[" << nX << "][" << nY << "] has beach protection factor = " << dFactor << endl;

   return dFactor;
}

//===============================================================================================================================
//! Fills in 'holes' in the beach protection i.e. orphan cells which get omitted because of rounding problems. Also removes 'legacy' cliff notches
//===============================================================================================================================
void CSimulation::FillInBeachProtectionHolesAndRemoveLegacyCliffs(void)
{
   for (int nX = 0; nX < m_nXGridSize; nX++)
   {
      for (int nY = 0; nY < m_nYGridSize; nY++)
      {
         // Find any 'legacy' ciff cells: cells with an erosional notch apex elevation which is now - due to shore platform erosion - above the top of the consolidated sediment
         double const dNotchApexElev = m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->dGetCliffNotchApexElev();
         if (! bFPIsEqual(dNotchApexElev, DBL_NODATA, TOLERANCE))
         {
            // This cell has an erosional notch
            double const dSedTopElevNoTalus = m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevOmitTalus();
            if (dNotchApexElev >= dSedTopElevNoTalus)
            {
               // The apex elevation of the notch is above the top of the consolidated sediment, so this notch has been removed
               m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetCliffNotchApexElev(DBL_NODATA);
               m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetCliffNotchIncisionDepth(DBL_NODATA);

               // Now determine the landform category
               int const nTopLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopNonZeroLayerAboveBasement();
               CRWCellLayer* pTopLayer = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer);

               if (pTopLayer->bHasTalus())
               {
                  // There is talus here
                  m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetLFCategory(LF_DRIFT_TALUS);
               }
               else if (pTopLayer->bHasUncons())
               {
                  // There is some unconsolidated sediment here
                  m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetLFCategory(LF_DRIFT_BEACH);
               }
               else
               {
                  // Set as hinterland
                  m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetLFCategory(LF_HINTERLAND);
               }
            }
         }

         // Now look at beach protection
         if ((m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea()) && (bFPIsEqual(m_pRasterGrid->m_Cell[nX][nY].dGetBeachProtectionFactor(), DBL_NODATA, TOLERANCE)))
         {
            // This is a sea cell, and it has an initialised beach protection value. So look at its eight neighbours
            int nXTmp;
            int nYTmp;
            int nAdjacent = 0;
            double dBeachProtection = 0;

            // North
            nXTmp = nX;
            nYTmp = nY - 1;

            if ((bIsWithinValidGrid(nXTmp, nYTmp)) && (! bFPIsEqual(m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetBeachProtectionFactor(), DBL_NODATA, TOLERANCE)))
            {
               nAdjacent++;
               dBeachProtection += m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetBeachProtectionFactor();
            }

            // North-east
            nXTmp = nX + 1;
            nYTmp = nY - 1;

            if ((bIsWithinValidGrid(nXTmp, nYTmp)) && (! bFPIsEqual(m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetBeachProtectionFactor(), DBL_NODATA, TOLERANCE)))
            {
               nAdjacent++;
               dBeachProtection += m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetBeachProtectionFactor();
            }

            // East
            nXTmp = nX + 1;
            nYTmp = nY;

            if ((bIsWithinValidGrid(nXTmp, nYTmp)) && (! bFPIsEqual(m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetBeachProtectionFactor(), DBL_NODATA, TOLERANCE)))
            {
               nAdjacent++;
               dBeachProtection += m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetBeachProtectionFactor();
            }

            // South-east
            nXTmp = nX + 1;
            nYTmp = nY + 1;

            if ((bIsWithinValidGrid(nXTmp, nYTmp)) && (! bFPIsEqual(m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetBeachProtectionFactor(), DBL_NODATA, TOLERANCE)))
            {
               nAdjacent++;
               dBeachProtection += m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetBeachProtectionFactor();
            }

            // South
            nXTmp = nX;
            nYTmp = nY + 1;

            if ((bIsWithinValidGrid(nXTmp, nYTmp)) && (! bFPIsEqual(m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetBeachProtectionFactor(), DBL_NODATA, TOLERANCE)))
            {
               nAdjacent++;
               dBeachProtection += m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetBeachProtectionFactor();
            }

            // South-west
            nXTmp = nX - 1;
            nYTmp = nY + 1;

            if ((bIsWithinValidGrid(nXTmp, nYTmp)) && (! bFPIsEqual(m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetBeachProtectionFactor(), DBL_NODATA, TOLERANCE)))
            {
               nAdjacent++;
               dBeachProtection += m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetBeachProtectionFactor();
            }

            // West
            nXTmp = nX - 1;
            nYTmp = nY;

            if ((bIsWithinValidGrid(nXTmp, nYTmp)) && (! bFPIsEqual(m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetBeachProtectionFactor(), DBL_NODATA, TOLERANCE)))
            {
               nAdjacent++;
               dBeachProtection += m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetBeachProtectionFactor();
            }

            // North-west
            nXTmp = nX - 1;
            nYTmp = nY - 1;

            if ((bIsWithinValidGrid(nXTmp, nYTmp)) && (! bFPIsEqual(m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetBeachProtectionFactor(), DBL_NODATA, TOLERANCE)))
            {
               nAdjacent++;
               dBeachProtection += m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetBeachProtectionFactor();
            }

            // If this sea cell has eight neighbours with initialised beach protection values, then assume that it should not have an uninitialised beach protection value. Set it to the average of its neighbours
            if (nAdjacent == 8)
            {
               m_pRasterGrid->m_Cell[nX][nY].SetBeachProtectionFactor(dBeachProtection / 8);
            }
         }
      }
   }
}

//===============================================================================================================================
//! Fills in 'holes' in the potential platform erosion i.e. orphan cells which get omitted because of rounding problems
//===============================================================================================================================
void CSimulation::FillPotentialPlatformErosionHoles(void)
{
   for (int nX = 0; nX < m_nXGridSize; nX++)
   {
      for (int nY = 0; nY < m_nYGridSize; nY++)
      {
         if ((m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea()) && (bFPIsEqual(m_pRasterGrid->m_Cell[nX][nY].dGetPotentialPlatformErosion(), 0.0, TOLERANCE)))
         {
            // This is a sea cell, it has a zero potential platform erosion value. So look at its N-S and W-E neighbours
            int nXTmp;
            int nYTmp;
            int nAdjacent = 0;
            double dPotentialPlatformErosion = 0;

            // North
            nXTmp = nX;
            nYTmp = nY - 1;

            if ((bIsWithinValidGrid(nXTmp, nYTmp)) && (! bFPIsEqual(m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetPotentialPlatformErosion(), 0.0, TOLERANCE)))
            {
               nAdjacent++;
               dPotentialPlatformErosion += m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetPotentialPlatformErosion();
            }

            // East
            nXTmp = nX + 1;
            nYTmp = nY;

            if ((bIsWithinValidGrid(nXTmp, nYTmp)) && (! bFPIsEqual(m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetPotentialPlatformErosion(), 0.0, TOLERANCE)))
            {
               nAdjacent++;
               dPotentialPlatformErosion += m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetPotentialPlatformErosion();
            }

            // South
            nXTmp = nX;
            nYTmp = nY + 1;

            if ((bIsWithinValidGrid(nXTmp, nYTmp)) && (! bFPIsEqual(m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetPotentialPlatformErosion(), 0.0, TOLERANCE)))
            {
               nAdjacent++;
               dPotentialPlatformErosion += m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetPotentialPlatformErosion();
            }

            // West
            nXTmp = nX - 1;
            nYTmp = nY;

            if ((bIsWithinValidGrid(nXTmp, nYTmp)) && (! bFPIsEqual(m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetPotentialPlatformErosion(), 0.0, TOLERANCE)))
            {
               nAdjacent++;
               dPotentialPlatformErosion += m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetPotentialPlatformErosion();
            }

            // If this sea cell has four neighbours with non-zero potential platform erosion values, then assume that it should not have a zero potential platform erosion value. Set it to the average of its neighbours
            if (nAdjacent == 4)
            {
               double const dThisPotentialPlatformErosion = dPotentialPlatformErosion / 4;

               m_pRasterGrid->m_Cell[nX][nY].SetPotentialPlatformErosion(dThisPotentialPlatformErosion);

               // Update this-timestep totals
               m_ulThisIterNumPotentialPlatformErosionCells++;
               m_dThisIterPotentialPlatformErosion += dThisPotentialPlatformErosion;
               // assert(isfinite(m_dThisIterPotentialPlatformErosion));

               // Increment the check values
               m_ulTotPotentialPlatformErosionBetweenProfiles++;
               m_dTotPotentialPlatformErosionBetweenProfiles += dThisPotentialPlatformErosion;
            }
         }
      }
   }
}

//===============================================================================================================================
//! Constructs a temporary coast profile which is parallel to the start profile
//===============================================================================================================================
void CSimulation::ConstructParallelProfile(int const nProfileStartX, int const nProfileStartY, int const nParCoastX, int const nParCoastY, int const nProfSize, vector<CGeom2DIPoint> *const pPtViGridProfile, vector<CGeom2DIPoint> *pPtiVGridParProfile, vector<CGeom2DPoint> *pPtVExtCRSParProfile)
{
   // OK, we have the coastline start point for the parallel profile. Now construct a temporary profile, parallel to the coastline-normal profile, starting from this point
   int const nXOffset = nParCoastX - nProfileStartX;
   int const nYOffset = nParCoastY - nProfileStartY;

   // Append co-ord values for the temporary profile
   for (int nProfileStartPoint = 0; nProfileStartPoint < nProfSize; nProfileStartPoint++)
   {
      // Get the grid coordinates of the cell which is this distance seaward, from the coastline-normal profile
      int const nXProf = pPtViGridProfile->at(nProfileStartPoint).nGetX();
      int const nYProf = pPtViGridProfile->at(nProfileStartPoint).nGetY();

      // Now calculate the grid coordinates of this cell, which is potentially in the parallel profile
      int const nXPar = nXProf + nXOffset;
      int const nYPar = nYProf + nYOffset;

      // Is this cell within the grid? If not, cut short the profile
      if (! bIsWithinValidGrid(nXPar, nYPar))
      {
         LogStream << "NOT WITHIN GRID [" << nXPar << "][" << nYPar << "]" << endl;
         return;
      }

      // // Have we hit an adjacent coastline-normal profile? If so, cut short
      // if (m_pRasterGrid->m_Cell[nXPar][nYPar].bIsProfile())
      // {
      //    // LogStream << "HIT PROFILE " << m_pRasterGrid->m_Cell[nXPar][nYPar].nGetProfileID() << " at [" << nXPar << "][" << nYPar << "] = {" << dGridCentroidXToExtCRSX(nXPar) << ", " <<  dGridCentroidYToExtCRSY(nYPar) << "}" << endl;
      //    return;
      // }

      // OK, append the cell details
      pPtiVGridParProfile->push_back(CGeom2DIPoint(nXPar, nYPar));
      pPtVExtCRSParProfile->push_back(CGeom2DPoint(dGridCentroidXToExtCRSX(nXPar), dGridCentroidYToExtCRSY(nYPar)));
   }
}
