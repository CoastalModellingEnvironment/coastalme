/*!
   \file do_beach_within_polygon.cpp
   \brief Does within-polygon actual erosion and distribution of transported beach sediment
   \details TODO 001 A more detailed description of these routines.
   \author David Favis-Mortlock
   \author Andres Payo
   \author Wilf Chun
   \date 2025
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
using std::endl;

#include <algorithm>
using std::shuffle;

#include "cme.h"
#include "simulation.h"
#include "coast.h"
#include "2d_point.h"
#include "2di_point.h"

//===============================================================================================================================
//! Erodes unconsolidated beach sediment of one texture class on the cells within a polygon. This is done by working down the coastline and constructing profiles which are parallel to the up-coast polygon boundary; then reversing direction and going up-coast, constructing profiles parallel to the down-coast boundary. Then iteratively fit a Dean equilibrium profile until the normal's share of the change in total depth of unconsolidated sediment is accommodated under the revised profile. For erosion, this reduces the beach volume
//===============================================================================================================================
int CSimulation::nDoUnconsErosionOnPolygon(int const nCoast, CGeomCoastPolygon* pPolygon, int const nTexture, double const dErosionTargetOnPolygon, double& dEroded)
{
   string strTexture;

   if (nTexture == TEXTURE_FINE)
      strTexture = "fine";

   else if (nTexture == TEXTURE_SAND)
      strTexture = "sand";

   else if (nTexture == TEXTURE_COARSE)
      strTexture = "coarse";

   // Intialise
   double dStillToErodeOnPolygon = dErosionTargetOnPolygon;
   dEroded = 0;

   // Get the up-coast and down-coast boundary details
   int const nUpCoastProfile = pPolygon->nGetUpCoastProfile();
   CGeomProfile* pUpCoastProfile = m_VCoast[nCoast].pGetProfile(nUpCoastProfile);

   int const nDownCoastProfile = pPolygon->nGetDownCoastProfile();
   CGeomProfile const* pDownCoastProfile = m_VCoast[nCoast].pGetProfile(nDownCoastProfile);

   // We will use only part of the up-coast boundary profile, seaward as far as the depth of closure. First find the seaward end point of this up-coast part-profile. Note that this does not change as the landwards offset changes
   int const nIndex = pUpCoastProfile->nGetCellGivenDepth(m_pRasterGrid, m_dDepthOfClosure);

   if (nIndex == INT_NODATA)
   {
      LogStream << m_ulIter << ": " << ERR << "while eroding unconsolidated " + strTexture + " sediment on polygon " << pPolygon->nGetPolygonCoastID() << ", could not find the seaward end point of the up-coast profile (" << nUpCoastProfile << ") for depth of closure = " << m_dDepthOfClosure << endl;

      return RTN_ERR_NO_SEAWARD_END_OF_PROFILE_BEACH_EROSION;
   }

   // The part-profile length is one greater than nIndex, since nGetCellGivenDepth() returns the index of the cell at the depth of closure
   int const nUpCoastPartProfileLen = nIndex + 1;

   // assert(bIsWithinValidGrid(&PtiUpCoastPartProfileSeawardEnd));

   // Store the cell coordinates of the boundary part-profile in reverse (sea to coast) order so we can append to the coastward end as we move inland (i.e. as nInlandOffset increases)
   vector<CGeom2DIPoint> PtiVUpCoastPartProfileCell;

   // for (int n = 0; n < nUpCoastPartProfileLen; n++)
   // PtiVUpCoastPartProfileCell.push_back(*pUpCoastProfile->pPtiGetCellInProfile(nUpCoastPartProfileLen - n - 1));
   for (int n = nUpCoastPartProfileLen - 1; n >= 0; n--)
   {
      CGeom2DIPoint const Pti = *pUpCoastProfile->pPtiGetCellInProfile(n);
      PtiVUpCoastPartProfileCell.push_back(Pti);
   }

   int const nUpCoastProfileCoastPoint = pUpCoastProfile->nGetCoastPoint();
   int const nDownCoastProfileCoastPoint = pDownCoastProfile->nGetCoastPoint();
   int const nXUpCoastProfileExistingCoastPoint = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nUpCoastProfileCoastPoint)->nGetX();
   int const nYUpCoastProfileExistingCoastPoint = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nUpCoastProfileCoastPoint)->nGetY();
   int nCoastSegLen;

   // Store the coast point numbers for this polygon so that we can shuffle them
   vector<int> nVCoastPoint;

   if (nDownCoastProfileCoastPoint == m_VCoast[nCoast].nGetCoastlineSize() - 1)
   {
      // This is the final down-coast polygon, so also include the down-coast polygon boundary
      nCoastSegLen = nDownCoastProfileCoastPoint - nUpCoastProfileCoastPoint + 1;

      for (int nCoastPoint = nUpCoastProfileCoastPoint; nCoastPoint <= nDownCoastProfileCoastPoint; nCoastPoint++)
         nVCoastPoint.push_back(nCoastPoint);
   }
   else
   {
      // This is not the final down-coast polygon, so do not include the polygon's down-coast boundary
      nCoastSegLen = nDownCoastProfileCoastPoint - nUpCoastProfileCoastPoint;

      for (int nCoastPoint = nUpCoastProfileCoastPoint; nCoastPoint < nDownCoastProfileCoastPoint; nCoastPoint++)
         nVCoastPoint.push_back(nCoastPoint);
   }

   // Estimate the volume of sediment which is to be eroded from each parallel profile
   double const dAllTargetPerProfile = dErosionTargetOnPolygon / nCoastSegLen;

   // Shuffle the coast points, this is necessary so that leaving the loop does not create sequence-related artefacts
   shuffle(nVCoastPoint.begin(), nVCoastPoint.begin() + nCoastSegLen, m_Rand[1]);

   // Traverse the polygon's existing coastline in a DOWN-COAST (i.e. increasing coastpoint indices) sequence, at each coast point fitting a Dean profile which is parallel to the up-coast polygon boundary
   for (int n = 0; n < nCoastSegLen; n++)
   {
      // Get the coast point
      int const nCoastPoint = nVCoastPoint[n];

      CGeom2DIPoint const PtiCoastPoint = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint);
      int const nCoastX = PtiCoastPoint.nGetX();
      int const nCoastY = PtiCoastPoint.nGetY();

      // LogStream << m_ulIter << ": nCoastX = " << nCoastX << ", nCoastY = " << nCoastY << ", this is {" << dGridCentroidXToExtCRSX(nCoastX) << ", " <<  dGridCentroidYToExtCRSY(nCoastY) << "}" << endl;

      // Is the coast cell an intervention structure?
      if (bIsInterventionCell(nCoastX, nCoastY))
      {
         // No erosion possible on this parallel profile, so move on
         // LogStream << m_ulIter << ": intervention structure at coast point [" << nCoastX << "][" << nCoastY << "] = {" << dGridCentroidXToExtCRSX(nCoastX) << ", " <<  dGridCentroidYToExtCRSY(nCoastY) << "}, cannot erode this parallel profile" << endl;

         continue;
      }

      // LogStream << m_ulIter << ": PtiVUpCoastPartProfileCell.back().nGetX() = " << PtiVUpCoastPartProfileCell.back().nGetX() << ", PtiVUpCoastPartProfileCell.back().nGetY() = " << PtiVUpCoastPartProfileCell.back().nGetY() << ", this is {" << dGridCentroidXToExtCRSX(PtiVUpCoastPartProfileCell.back().nGetX()) << ", " <<  dGridCentroidYToExtCRSY(PtiVUpCoastPartProfileCell.back().nGetY()) << "}" << endl;

      // Not an intervention structure, so calculate the x-y offset between this coast point, and the coast point of the up-coast normal
      int const nXOffset = nCoastX - PtiVUpCoastPartProfileCell.back().nGetX();
      int const nYOffset = nCoastY - PtiVUpCoastPartProfileCell.back().nGetY();

      // Get the x-y coords of a profile starting from this coast point and parallel to the up-coast polygon boundary profile (these are in reverse sequence, like the boundary part-profile)
      vector<CGeom2DIPoint> VPtiParProfile;

      for (int m = 0; m < nUpCoastPartProfileLen; m++)
      {
         // TODO 017 Check that each point is within valid grid, do same for other similar places in rest of model
         CGeom2DIPoint const PtiTmp(PtiVUpCoastPartProfileCell[m].nGetX() + nXOffset, PtiVUpCoastPartProfileCell[m].nGetY() + nYOffset);
         VPtiParProfile.push_back(PtiTmp);
      }

      // Get the elevations of the start and end points of the parallel profiles (as we extend the profile inland, the elevation of the new coast point of the Dean profile is set to the elevation of the original coast point)
      int nParProfEndX = VPtiParProfile[0].nGetX();
      int nParProfEndY = VPtiParProfile[0].nGetY();

      // Safety check
      if (! bIsWithinValidGrid(nParProfEndX, nParProfEndY))
      {
         // if (m_nLogFileDetail >= LOG_FILE_ALL)
         // LogStream << WARN << "while eroding unconsolidated " + strTexture + " sediment for coast " << nCoast << " polygon " << pPolygon->nGetPolygonCoastID() << ", hit edge of grid at [" << nParProfEndX << "][" << nParProfEndY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its seaward end" << endl;

         KeepWithinValidGrid(nCoastX, nCoastY, nParProfEndX, nParProfEndY);
         VPtiParProfile[0].SetX(nParProfEndX);
         VPtiParProfile[0].SetY(nParProfEndY);
      }

      bool bHitEdge = false;
      bool bEndProfile = false;
      bool bZeroGradient = false;
      bool bEnoughEroded = false;

      int nParProfLen;
      int nInlandOffset = -1;

      // Note that we must not consider the extra elevation due to to any talus here, since talus is removed quickly compared with the time required to create a Dean profile
      double const dParProfCoastElev = m_pRasterGrid->m_Cell[nCoastX][nCoastY].dGetAllSedTopElevOmitTalus();
      double const dParProfEndElev = m_pRasterGrid->m_Cell[nParProfEndX][nParProfEndY].dGetAllSedTopElevOmitTalus();

      vector<double> VdParProfileDeanElev;

      // These are for saving values for each offset
      vector<int> VnParProfLenEachOffset;
      vector<double> VdAmountEachOffset;
      vector<vector<CGeom2DIPoint>> VVPtiParProfileEachOffset;
      vector<vector<double>> VVdParProfileDeanElevEachOffset;

      // OK, loop either until we can erode sufficient unconsolidated sediment, or until the landwards-moving parallel profile hits the grid edge
      while (true)
      {
         // Move inland by one cell
         nInlandOffset++;

         if (nInlandOffset > 0)
         {
            if (nInlandOffset > (pUpCoastProfile->nGetNumCellsInProfile() - 1))
            {
               if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
                  LogStream << m_ulIter << ":\t coast " << nCoast << " reached end of up-coast profile " << nUpCoastProfile << " during down-coast erosion of unconsolidated " + strTexture + " sediment for coast " << nCoast << " polygon " << pPolygon->nGetPolygonCoastID() << " (nCoastPoint = " << nCoastPoint << " nInlandOffset = " << nInlandOffset << ")" << endl;

               bEndProfile = true;
               break;
            }

            // Extend the parallel profile by one cell in the coastward direction. First get the coords of the cell that is nInlandOffset cells seaward from the existing up-coast part-profile start point
            CGeom2DIPoint const PtiUpCoastTmp = *pUpCoastProfile->pPtiGetCellInProfile(nInlandOffset);

            // Then get the offset between this PtiUpCoastTmp cell and the existing up-coast part-profile start point, and use the reverse of this offset to get the coordinates of the cell that extends the existing up-coast part-profile landwards
            int const nXUpCoastStartOffset = PtiUpCoastTmp.nGetX() - nXUpCoastProfileExistingCoastPoint;
            int const nYUpCoastStartOffset = PtiUpCoastTmp.nGetY() - nYUpCoastProfileExistingCoastPoint;
            int const nXUpCoastThisStart = nCoastX - nXUpCoastStartOffset;
            int const nYUpCoastThisStart = nCoastY - nYUpCoastStartOffset;

            // Is the new landwards point within the raster grid?
            if (! bIsWithinValidGrid(nXUpCoastThisStart, nYUpCoastThisStart))
            {
               // It isn't
               // LogStream << WARN << "reached edge of grid at [" << nXUpCoastThisStart << "][" << nYUpCoastThisStart << "] = {" << dGridCentroidXToExtCRSX(nXUpCoastThisStart) << ", " << dGridCentroidYToExtCRSY(nYUpCoastThisStart) << "} during DOWN-COAST beach erosion for coast " << nCoast << " polygon " << nPoly << ", nCoastPoint = " << nCoastPoint << ", this is [" << nCoastX << "][" << nCoastY << "] = {" << dGridCentroidXToExtCRSX(nCoastX) << ", " <<  dGridCentroidYToExtCRSY(nCoastY) << "}, nInlandOffset = " << nInlandOffset << ")" << endl << endl;

               // TODO 018 Need to improve this: at present we just abandon erosion on this coast point and move to another coast point
               bHitEdge = true;
               break;
            }

            // CGeom2DIPoint PtiThisUpCoastStart(nXUpCoastThisStart, nYUpCoastThisStart);

            // Calculate the coordinates of a possible new landwards cell for the parallel profile
            int nXParNew = nXUpCoastThisStart + nXOffset;
            int nYParNew = nYUpCoastThisStart + nYOffset;

            // Safety check
            if (! bIsWithinValidGrid(nXParNew, nYParNew))
            {
               // LogStream << WARN << "while eroding beach on coast " << nCoast << " polygon " << nPoly << " (nInlandOffset = " << nInlandOffset << "), outside valid grid at [" << nXParNew << "][" << nYParNew << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its landward end, is now [";

               KeepWithinValidGrid(nCoastX, nCoastY, nXParNew, nYParNew);

               // LogStream << "[" << nXParNew << "][" << nYParNew << "] = {" << dGridCentroidXToExtCRSX(nXParNew) << ", " <<  dGridCentroidYToExtCRSY(nYParNew) << "}" << endl;

               // Is this cell already in the parallel profile?
               if ((VPtiParProfile.back().nGetX() != nXParNew) || (VPtiParProfile.back().nGetY() != nYParNew))
               {
                  // It isn't, so append it to the parallel profile
                  CGeom2DIPoint const PtiTmp(nXParNew, nYParNew);
                  VPtiParProfile.push_back(PtiTmp);
               }
            }

            else
            {
               // No problem, so just append this to the parallel profile
               CGeom2DIPoint const PtiTmp(nXParNew, nYParNew);
               VPtiParProfile.push_back(PtiTmp);
            }
         }

         nParProfLen = static_cast<int>(VPtiParProfile.size());

         if (nParProfLen < MIN_PARALLEL_PROFILE_SIZE)
         {
            // Can't have a meaningful parallel profile with very few points
            // LogStream << m_ulIter << ": only " << nParProfLen << " points in parallel profile, min is " << MIN_PARALLEL_PROFILE_SIZE << ", abandoning" << endl;

            continue;
         }

         // for (int m = 0; m < static_cast<int>(VPtiParProfile.size()); m++)
         // LogStream << "[" << VPtiParProfile[m].nGetX() << "][" << VPtiParProfile[m].nGetY() << "] ";
         // LogStream << endl;

         // Get the distance between the start and end of the parallel profile, in external CRS units. Note that the parallel profile coordinates are in reverse sequence
         CGeom2DPoint const PtStart = PtGridCentroidToExt(&VPtiParProfile.back());
         CGeom2DPoint const PtEnd = PtGridCentroidToExt(&VPtiParProfile[0]);

         // Calculate the length of the parallel profile
         double dParProfileLen = dGetDistanceBetween(&PtStart, &PtEnd);

         // Calculate the elevation difference between the start and end of the parallel profile
         double const dElevDiff = dParProfCoastElev - dParProfEndElev;

         if (bFPIsEqual(dElevDiff, 0.0, TOLERANCE))
         {
            // Can't have a meaningful Dean profile with a near-zero elevation difference
            // TODO 019 Need to improve this: at present we just abandon erosion on this coast point and move to another coast point
            // LogStream << m_ulIter << ": zero gradient on parallel profile, abandoning" << endl;

            bZeroGradient = true;
            break;
         }

         // Safety check
         if (dParProfileLen <= 0)
            dParProfileLen = 1;

         // Solve for dA so that the existing elevations at the end of the parallel profile, and at the end of a Dean equilibrium profile on that part-normal, are the same
         double const dParProfA = dElevDiff / pow(dParProfileLen, DEAN_POWER);
         VdParProfileDeanElev.resize(nParProfLen, 0);
         double const dInc = dParProfileLen / (nParProfLen - 1);

         // For the parallel profile, calculate the Dean equilibrium profile of the unconsolidated sediment h(y) = A * y^(2/3) where h(y) is the distance below the highest point in the profile at a distance y from the landward start of the profile
         CalcDeanProfile(&VdParProfileDeanElev, dInc, dParProfCoastElev, dParProfA, false, 0, 0);

         vector<double> dVParProfileNow(nParProfLen, 0);
         vector<bool> bVProfileValid(nParProfLen, true);

         for (int m = 0; m < nParProfLen; m++)
         {
            int const nX = VPtiParProfile[nParProfLen - m - 1].nGetX();
            int const nY = VPtiParProfile[nParProfLen - m - 1].nGetY();

            // Safety check
            if (! bIsWithinValidGrid(nX, nY))
            {
               // if (m_nLogFileDetail >= LOG_FILE_ALL)
               // LogStream << WARN << "while constructing parallel profile for erosion of unconsolidated " + strTexture + " sediment on coast " << nCoast << " polygon " << nPoly << ", hit edge of grid at [" << nX << "][" << nY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its seaward end" << endl;

               bVProfileValid[m] = false;
               continue;
            }

            // Don't erode intervention cells
            if (bIsInterventionCell(nX, nY))
               bVProfileValid[m] = false;

            dVParProfileNow[m] = m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevOmitTalus();
         }

         // Get the total difference in elevation (present profile - Dean profile)
         double const dParProfTotDiff = dSubtractProfiles(&dVParProfileNow, &VdParProfileDeanElev, &bVProfileValid);

         // // DEBUG CODE -----------------------------------------------------
         // LogStream << m_ulIter<< ": eroding polygon " << nPoly << " at coast point " << nCoastPoint << ", parallel profile with nInlandOffset = " << nInlandOffset << ", from [" << VPtiParProfile.back().nGetX() << "][" << VPtiParProfile.back().nGetY() << "] = {" << dGridCentroidXToExtCRSX(VPtiParProfile.back().nGetX()) << ", " <<  dGridCentroidYToExtCRSY(VPtiParProfile.back().nGetY()) << "} to [" << VPtiParProfile[0].nGetX() << "][" << VPtiParProfile[0].nGetY() << "] = {" << dGridCentroidXToExtCRSX(VPtiParProfile[0].nGetX()) << ", " <<  dGridCentroidYToExtCRSY(VPtiParProfile[0].nGetY()) << "}, nParProfLen = " << nParProfLen << " dParProfileLen = " << dParProfileLen << " dParProfCoastElev = " << dParProfCoastElev << " dParProfEndElev = " << dParProfEndElev << " dParProfA = " << dParProfA << endl;

         // LogStream << "Profile now:" << endl;
         // for (int n = 0; n < nParProfLen; n++)
         // {
         // if (bVProfileValid[n])
         // LogStream << dVParProfileNow[n] << " ";
         // else
         // LogStream << "XXX ";
         // }
         // LogStream << endl << endl;;
         //
         // LogStream << "Parallel Dean profile for erosion:" << endl;
         // for (int n = 0; n < nParProfLen; n++)
         // {
         // if (bVProfileValid[n])
         // LogStream << VdParProfileDeanElev[n] << " ";
         // else
         // LogStream << "XXX ";
         // }
         // LogStream << endl << endl;
         //
         // LogStream << "Difference (present profile minus Dean profile):" << endl;
         // for (int n = 0; n < nParProfLen; n++)
         // {
         // if (bVProfileValid[n])
         // LogStream << dVParProfileNow[n] - VdParProfileDeanElev[n] << " ";
         // else
         // LogStream << "XXX ";
         // }
         // LogStream << endl << endl;
         //
         // LogStream << "dParProfTotDiff = " << dParProfTotDiff << " dAllTargetPerProfile = " << dAllTargetPerProfile << endl;
         // // DEBUG CODE -----------------------------------------------------

         // So will we be able to erode as much as is needed?
         if (dParProfTotDiff > dAllTargetPerProfile)
         {
            // LogStream << m_ulIter << ": eroding polygon " << nPoly << " at coast point " << nCoastPoint << ", on parallel profile with nInlandOffset = " << nInlandOffset << ", can meet erosion target: dParProfTotDiff = " << dParProfTotDiff << " dAllTargetPerProfile = " << dAllTargetPerProfile << endl;

            bEnoughEroded = true;
            break;
         }

         // We have not been able to reach the erosion target for this parallel profile. Now, if we have moved at least MIN_INLAND_OFFSET_UNCONS_EROSION cells inland, and dParProfTotDiff is zero, break out of the loop
         if ((nInlandOffset >= MIN_INLAND_OFFSET_UNCONS_EROSION) && (bFPIsEqual(dParProfTotDiff, 0.0, TOLERANCE)))
         {
            if (m_nLogFileDetail >= LOG_FILE_ALL)
               LogStream << m_ulIter << ": leaving loop because nInlandOffset (" << nInlandOffset << ") >= MIN_INLAND_OFFSET_UNCONS_EROSION) and dParProfTotDiff = " << dParProfTotDiff << endl;

            break;
         }

         // Save the amount which can be eroded for this offset
         VdAmountEachOffset.push_back(dParProfTotDiff);
         VnParProfLenEachOffset.push_back(nParProfLen);
         VVPtiParProfileEachOffset.push_back(VPtiParProfile);
         VVdParProfileDeanElevEachOffset.push_back(VdParProfileDeanElev);
      }

      // If we hit the edge of the grid, or have a zero gradient on the profile, or this is an end profile, then abandon this profile and do the next parallel profile
      // TODO 019 TODO 018 Improve this, see above
      if (bHitEdge || bEndProfile || bZeroGradient)
         continue;

      // OK we will do some erosion on this parallel profile, set the target for this profile to be the full target amount
      double dStillToErodeOnProfile = dAllTargetPerProfile;

      if (! bEnoughEroded)
      {
         // We have not been able to reach the target for erosion on this parallel profile. So find the offset that gives us the largest erosion amount
         int nOffsetForLargestPossible = -1;
         double dLargestPossibleErosion = 0;

         for (unsigned int nn = 0; nn < VdAmountEachOffset.size(); nn++)
         {
            if (VdAmountEachOffset[nn] > dLargestPossibleErosion)
            {
               dLargestPossibleErosion = VdAmountEachOffset[nn];
               nOffsetForLargestPossible = nn;
            }
         }

         // If every offset gave zero erosion then abandon this profile and do the next parallel profile
         if (nOffsetForLargestPossible < 0)
            continue;

         // OK, we have an offset which gives us the largest possible erosion (but less than the full target amount), continue with this
         nInlandOffset = nOffsetForLargestPossible;
         dStillToErodeOnProfile = dLargestPossibleErosion;
         nParProfLen = VnParProfLenEachOffset[nInlandOffset];
         VPtiParProfile = VVPtiParProfileEachOffset[nInlandOffset];
         VdParProfileDeanElev = VVdParProfileDeanElevEachOffset[nInlandOffset];

         // LogStream << m_ulIter << ": eroding polygon " << nPoly << " at coast point " << nCoastPoint << ", for parallel profile with nInlandOffset = " << nInlandOffset << ", could not meet erosion target dAllTargetPerProfile = " << dAllTargetPerProfile << ", instead using best possible: nInlandOffset = " << nInlandOffset << " which gives dStillToErodeOnProfile = " << dStillToErodeOnProfile << endl;
      }

      // This value of nInlandOffset gives us some (tho' maybe not enough) erosion. So do the erosion of this sediment size class, by working along the parallel profile from the landward end (which is inland from the existing coast, if nInlandOffset > 0). Note that dStillToErodeOnProfile and dStillToErodeOnPolygon are changed within nDoParallelProfileUnconsErosion()
      int const nRet = nDoParallelProfileUnconsErosion(pPolygon, nCoastPoint, nCoastX, nCoastY, nTexture, nInlandOffset, nParProfLen, &VPtiParProfile, &VdParProfileDeanElev, dStillToErodeOnProfile, dStillToErodeOnPolygon, dEroded);

      if (nRet != RTN_OK)
         return nRet;

      // LogStream << m_ulIter << ": eroding polygon " << nPoly << ": finished at coast point " << nCoastPoint << " dStillToErodeOnProfile = " << dStillToErodeOnProfile << " dStillToErodeOnPolygon = " << dStillToErodeOnPolygon << " dAllTargetPerProfile = " << dAllTargetPerProfile << endl;

      // if (dStillToErodeOnProfile > 0)
      // {
      //    //          LogStream << "                     X Still to erode on profile = " << dStillToErodeOnProfile << " dStillToErodeOnPolygon WAS = " << dStillToErodeOnPolygon;
      // dStillToErodeOnPolygon += dStillToErodeOnProfile;
      // dAllTargetPerProfile += (dStillToErodeOnProfile / (nCoastSegLen - n));
      //    //          LogStream << " dStillToErodeOnPolygon NOW = " << dStillToErodeOnPolygon << endl;
      // }

      if (dStillToErodeOnPolygon <= 0)
      {
         // LogStream << m_ulIter << ": YYYYYYYYYYYYYYY in polygon " << nPoly << ", leaving loop because dStillToErodeOnPolygon = " << dStillToErodeOnPolygon << endl;
         break;
      }
   }

   // How much have we been able to erode?
   // LogStream << m_ulIter << ": nPoly = " << nPoly << " dEroded = " << dEroded << endl;

   // // Is the depth eroded within TOLERANCE of the target depth-equivalent?
   // if (bFPIsEqual(dEroded, dErosionTargetOnPolygon, TOLERANCE))
   // {
   // if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
   // LogStream << m_ulIter << ": polygon " << nPoly << " actual erosion of " + strTexture + " sediment is approximately equal to estimate, actual = " << dEroded * m_dCellArea << " estimate = " << dErosionTargetOnPolygon * m_dCellArea << endl;
   // }
   // else
   // {
   // if (dEroded < dErosionTargetOnPolygon)
   // {
   // if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
   // LogStream << m_ulIter << ": polygon " << nPoly << " actual erosion of " + strTexture + " sediment is less than estimate, actual = " << dEroded * m_dCellArea << " estimate = " << dErosionTargetOnPolygon * m_dCellArea << " dStillToErodeOnPolygon = " << dStillToErodeOnPolygon * m_dCellArea << ". This reduces the volume of " + strTexture + " sediment exported from this polygon to " << (dErosionTargetOnPolygon - dStillToErodeOnPolygon) * m_dCellArea << endl;
   // }
   // else
   // {
   // if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
   // LogStream << ERR << "polygon " << nPoly << " actual erosion of " + strTexture + " sediment is greater than estimate, actual = " << dEroded * m_dCellArea << " estimate = " << dErosionTargetOnPolygon * m_dCellArea << " dStillToErodeOnPolygon = " << dStillToErodeOnPolygon * m_dCellArea << endl;
   // }
   // }

   return RTN_OK;
}

//===============================================================================================================================
//! This routine erodes unconsolidated beach sediment (either fine, sand, or coarse) on a parallel profile
//===============================================================================================================================
int CSimulation::nDoParallelProfileUnconsErosion(CGeomCoastPolygon* pPolygon, int const nCoastPoint, int const nCoastX, int const nCoastY, int const nTexture, int const nInlandOffset, int const nParProfLen, vector<CGeom2DIPoint> const* pVPtiParProfile, vector<double> const* pVdParProfileDeanElev, double& dStillToErodeOnProfile, double& dStillToErodeOnPolygon, double& dTotEroded)
{
   for (int nDistSeawardFromNewCoast = 0; nDistSeawardFromNewCoast < nParProfLen; nDistSeawardFromNewCoast++)
   {
      // Leave the loop if we have eroded enough for this polygon
      if (dStillToErodeOnPolygon <= 0)
      {
         if (m_nLogFileDetail >= LOG_FILE_ALL)
            LogStream << m_ulIter << ": AAA in polygon " << pPolygon->nGetPolygonCoastID() << " at coast point " << nCoastPoint << " nInlandOffset = " << nInlandOffset << ", leaving loop because enough erosion for polygon, dStillToErodeOnPolygon = " << dStillToErodeOnPolygon << " dStillToErodeOnProfile = " << dStillToErodeOnProfile << endl;

         break;
      }

      // Leave the loop if we have eroded enough for this parallel profile
      if (dStillToErodeOnProfile <= 0)
      {
         if (m_nLogFileDetail >= LOG_FILE_ALL)
            LogStream << m_ulIter << ": BBB in polygon " << pPolygon->nGetPolygonCoastID() << " at coast point " << nCoastPoint << " nInlandOffset = " << nInlandOffset << ", leaving loop because enough erosion for profile, dStillToErodeOnPolygon = " << dStillToErodeOnPolygon << " dStillToErodeOnProfile = " << dStillToErodeOnProfile << endl;

         break;
      }

      CGeom2DIPoint PtiTmp = pVPtiParProfile->at(nParProfLen - nDistSeawardFromNewCoast - 1);
      int
          nX = PtiTmp.nGetX(),
          nY = PtiTmp.nGetY();

      // Safety check
      if (! bIsWithinValidGrid(nX, nY))
      {
         // LogStream << WARN << "04 @@@@ while eroding polygon " << nPoly << ", hit edge of grid at [" << nX << "][" << nY << "] for parallel profile from coast point " << nCoastPoint << ". Constraining this parallel profile at its seaward end" << endl;

         KeepWithinValidGrid(nCoastX, nCoastY, nX, nY);
         PtiTmp.SetX(nX);
         PtiTmp.SetY(nY);
      }

      // Don't do anything to intervention cells
      if (bIsInterventionCell(nX, nY))
         continue;

      // Don't do cells twice
      if (! m_pRasterGrid->m_Cell[nX][nY].bBeachErosionOrDepositionThisIter())
      {
         // Get this cell's current elevation
         double const dThisElevNow = m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevOmitTalus();

         // LogStream << "\tnPoly = " << nPoly << ", [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " <<  dGridCentroidYToExtCRSY(nY) << "}  nCoastPoint = " << nCoastPoint << " nDistSeawardFromNewCoast = " << nDistSeawardFromNewCoast << " dThisElevNow = " << dThisElevNow << " Dean Elev = " << VdParProfileDeanElev[nDistSeawardFromNewCoast] << endl;

         // Subtract the two elevations
         double const dElevDiff = dThisElevNow - pVdParProfileDeanElev->at(nDistSeawardFromNewCoast);

         if ((dElevDiff > 0) && (m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea()))
         {
            // The current elevation is higher than the Dean elevation, so we have possible beach erosion (i.e. if not constrained by availability of unconsolidated sediment) here
            // LogStream << "\tnPoly = " << nPoly << " doing DOWN-COAST, possible beach erosion = " << dElevDiff << " at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " <<  dGridCentroidYToExtCRSY(nY) << "} nCoastPoint = " << nCoastPoint << " nDistSeawardFromNewCoast = " << nDistSeawardFromNewCoast << " dThisElevNow = " << dThisElevNow << " Dean Elev = " << pVdParProfileDeanElev->at(nDistSeawardFromNewCoast) << endl;

            // Now get the number of the highest layer with non-zero thickness
            int const nThisLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopNonZeroLayerAboveBasement();

            // Safety check
            if (nThisLayer == INT_NODATA)
               return RTN_ERR_NO_TOP_LAYER_DURING_BEACH_CALC;

            if (nThisLayer != NO_NONZERO_THICKNESS_LAYERS)
            {
               // We still have at least one layer left with non-zero thickness (i.e. we are not down to basement), and the cell's current elevation is higher than the Dean equilibrium profile elevation. So do some beach erosion
               double const dToErode = tMin(dElevDiff, dStillToErodeOnProfile, dStillToErodeOnPolygon);
               double dRemoved = 0;

               // assert(dToErode > 0);

               // Erode this sediment size class
               ErodeCellBeachSedimentSupplyLimited(nX, nY, nThisLayer, nTexture, dToErode, dRemoved);

               if (dRemoved > 0)
               {
                  // Update totals for the polygon and the profile
                  dTotEroded += dRemoved;
                  dStillToErodeOnProfile -= dRemoved;
                  dStillToErodeOnPolygon -= dRemoved;

                  // Update this-timestep totals
                  m_ulThisIterNumActualBeachErosionCells++;

                  // LogStream << m_ulIter << ": in polygon " << nPoly << ", actual beach erosion = " << dTmpTot << " at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " <<  dGridCentroidYToExtCRSY(nY) << "} nCoastPoint = " << nCoastPoint << " nDistSeawardFromNewCoast = " << nDistSeawardFromNewCoast << endl;

                  // // DEBUG CODE ==================================================================================================
                  // if (nTexture == TEXTURE_SAND)
                  // {
                  // LogStream << m_ulIter << ": $$$$$$$$$ BEACH 2 Dean profile lower than existing profile, SAND depth eroded on [" << nX << "][" << nY << "] in Poly " << nPoly << " is " << dRemoved * m_dCellArea << endl;
                  // }
                  // else
                  // {
                  // LogStream << m_ulIter << ": $$$$$$$$$ BEACH 2 Dean profile lower than existing profile, COARSE depth eroded on [" << nX << "][" << nY << "] in Poly " << nPoly << " is " << dRemoved * m_dCellArea << endl;
                  // }
                  // // DEBUG CODE ==================================================================================================

                  // Store the this-polygon depth of sediment eroded during Dean profile deposition of beach unconsolidated sediment
                  if (nTexture == TEXTURE_SAND)
                     pPolygon->AddBeachSandErodedDeanProfile(dRemoved);

                  else
                     pPolygon->AddBeachCoarseErodedDeanProfile(dRemoved);
               }
            }
         }

         else if ((dElevDiff < 0) && (m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea()))
         {
            if ((nTexture == TEXTURE_SAND) || (nTexture == TEXTURE_COARSE))
            {
               // The current elevation is below the Dean elevation, so we have can have unconsolidated sediment deposition here provided that we have some previously-eroded unconsolidated sediment (sand and coarse only) to deposit
               if (dTotEroded > 0)
               {
                  double dTotToDeposit = tMin(-dElevDiff, dTotEroded);

                  int const nTopLayer = m_pRasterGrid->m_Cell[nX][nY].nGetNumOfTopLayerAboveBasement();

                  // Safety check
                  if (nTopLayer == INT_NODATA)
                     return RTN_ERR_NO_TOP_LAYER_DURING_BEACH_CALC;

                  if (dTotToDeposit > 0)
                  {
                     dTotToDeposit = tMin(dTotToDeposit, dTotEroded);

                     if (nTexture == TEXTURE_SAND)
                     {
                        double const dSandNow = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->dGetSandDepth();
                        m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->SetSandDepth(dSandNow + dTotToDeposit);

                        // Set the changed-this-timestep switch
                        m_bUnconsChangedThisIter[nTopLayer] = true;

                        dTotEroded -= dTotToDeposit;

                        dStillToErodeOnProfile += dTotToDeposit;
                        dStillToErodeOnPolygon += dTotToDeposit;

                        // Update per-timestep totals
                        m_ulThisIterNumBeachDepositionCells++;
                        // m_dThisIterBeachDepositionSand += dTotToDeposit;
                     }

                     if (nTexture == TEXTURE_COARSE)
                     {
                        double const dCoarseNow = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->dGetCoarseDepth();
                        m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->SetCoarseDepth(dCoarseNow + dTotToDeposit);

                        // Set the changed-this-timestep switch
                        m_bUnconsChangedThisIter[nTopLayer] = true;

                        dTotEroded -= dTotToDeposit;

                        dStillToErodeOnProfile += dTotToDeposit;
                        dStillToErodeOnPolygon += dTotToDeposit;

                        // Update per-timestep totals
                        m_ulThisIterNumBeachDepositionCells++;
                        m_dThisIterBeachDepositionCoarse += dTotToDeposit;
                     }
                  }

                  // Now update the cell's layer elevations
                  m_pRasterGrid->m_Cell[nX][nY].CalcAllLayerElevsAndD50();

                  // Update the cell's sea depth
                  m_pRasterGrid->m_Cell[nX][nY].SetSeaDepth();

                  // Update the cell's beach deposition, and total beach deposition, values
                  m_pRasterGrid->m_Cell[nX][nY].IncrBeachDeposition(dTotToDeposit);

                  // And set the landform category
                  CRWCellLandform* pLandform = m_pRasterGrid->m_Cell[nX][nY].pGetLandform();
                  int const nCat = pLandform->nGetLFCategory();

                  if ((nCat != LF_SEDIMENT_INPUT_UNCONSOLIDATED) && (nCat != LF_SEDIMENT_INPUT_CONSOLIDATED))
                     pLandform->SetLFCategory(LF_DRIFT_BEACH);

                  // LogStream << m_ulIter << ": nPoly = " << nPoly << ", beach deposition = " << dTotToDeposit << " at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " <<  dGridCentroidYToExtCRSY(nY) << "} nCoastPoint = " << nCoastPoint << " nDistSeawardFromNewCoast = " << nDistSeawardFromNewCoast << endl;
               }
            }
         }
      }
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Erodes the unconsolidated beach sediment on a single cell, for a single size class, and returns the depth-equivalents of sediment removed
//===============================================================================================================================
void CSimulation::ErodeCellBeachSedimentSupplyLimited(int const nX, int const nY, int const nThisLayer, int const nTexture, double const dMaxToErode, double& dRemoved)
{
   // Find out how much unconsolidated sediment of this size class we have available on this cell
   double dExistingAvailable = 0;
   double dErodibility = 0;

   if (nTexture == TEXTURE_FINE)
   {
      dExistingAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetUnconsolidatedSediment()->dGetFineDepth();
      dErodibility = m_dFineErodibilityNormalized;
   }
   else if (nTexture == TEXTURE_SAND)
   {
      dExistingAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetUnconsolidatedSediment()->dGetSandDepth();
      dErodibility = m_dSandErodibilityNormalized;
   }
   else if (nTexture == TEXTURE_COARSE)
   {
      dExistingAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetUnconsolidatedSediment()->dGetCoarseDepth();
      dErodibility = m_dCoarseErodibilityNormalized;
   }

   // Is there any unconsolidated sediment on this cell?
   if (dExistingAvailable <= 0)
   {
      return;
   }

   // Erode some unconsolidated sediment
   double const dLowering = dErodibility * dMaxToErode;

   // Make sure we don't get -ve amounts left on the cell
   dRemoved = tMin(dExistingAvailable, dLowering);
   double const dRemaining = dExistingAvailable - dRemoved;

   if (nTexture == TEXTURE_FINE)
   {
      // Set the value for this layer
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetUnconsolidatedSediment()->SetFineDepth(dRemaining);
   }
   else if (nTexture == TEXTURE_SAND)
   {
      // Set the value for this layer
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetUnconsolidatedSediment()->SetSandDepth(dRemaining);
   }
   else if (nTexture == TEXTURE_COARSE)
   {
      // Set the value for this layer
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetUnconsolidatedSediment()->SetCoarseDepth(dRemaining);
   }

   // And set the changed-this-timestep switch
   m_bUnconsChangedThisIter[nThisLayer] = true;

   // Set the actual erosion value for this cell
   m_pRasterGrid->m_Cell[nX][nY].SetActualBeachErosion(dRemoved);

   if (dRemoved > 0)
   {
      // Recalculate the elevation of every layer
      m_pRasterGrid->m_Cell[nX][nY].CalcAllLayerElevsAndD50();

      // And update the cell's sea depth
      m_pRasterGrid->m_Cell[nX][nY].SetSeaDepth();
   }
}

//===============================================================================================================================
//! Deposits unconsolidated beach sediment (sand or coarse) on the cells within a polygon. This is done by working down the coastline and constructing profiles which are parallel to the up-coast polygon boundary; then reversing direction and going up-coast, constructing profiles parallel to the down-coast boundary. Then iteratively fit a Dean equilibrium profile until the normal's share of the change in total depth of unconsolidated sediment is accommodated under the revised profile. For deposition, this adds to the beach volume
//===============================================================================================================================
int CSimulation::nDoUnconsDepositionOnPolygon(int const nCoast, CGeomCoastPolygon* pPolygon, int const nTexture, double dTargetToDepositOnPoly, double& dDepositedOnPoly)
{
   // // DEBUG CODE #####################
   // if (m_ulIter == 1)
   // {
   // int nPoly = pPolygon->nGetPolygonCoastID();
   // LogStream << m_ulIter << ": entered nDoUnconsDepositionOnPolygon() nCoast = " << nCoast << " nPoly = " << nPoly << " dTargetToDepositOnPoly = " << dTargetToDepositOnPoly * m_dCellArea << " dDepositedOnPoly = " << dDepositedOnPoly * m_dCellArea << endl;
   // }
   // // DEBUG CODE #####################

   string strTexture;

   if (nTexture == TEXTURE_SAND)
      strTexture = "sand";

   else if (nTexture == TEXTURE_COARSE)
      strTexture = "coarse";

   // Don't bother with tiny amounts of deposition
   if (dTargetToDepositOnPoly < SED_ELEV_TOLERANCE)
      return RTN_OK;

   // Get the grid cell coordinates of this polygon's up-coast and down-coast profiles
   int const nUpCoastProfile = pPolygon->nGetUpCoastProfile();
   int const nDownCoastProfile = pPolygon->nGetDownCoastProfile();

   CGeomProfile* pUpCoastProfile = m_VCoast[nCoast].pGetProfile(nUpCoastProfile);
   CGeomProfile* pDownCoastProfile = m_VCoast[nCoast].pGetProfile(nDownCoastProfile);

   // We are using only part of each profile, seaward as far as the depth of closure. First find the seaward end point of the up-coast part-profile
   // CGeom2DIPoint PtiUpCoastPartProfileSeawardEnd;
   // int nIndex =  pUpCoastProfile->nGetCellGivenDepth(m_pRasterGrid, m_dDepthOfClosure, &PtiUpCoastPartProfileSeawardEnd);
   int const nIndex = pUpCoastProfile->nGetCellGivenDepth(m_pRasterGrid, m_dDepthOfClosure);

   if (nIndex == INT_NODATA)
   {
      LogStream << m_ulIter << ": " << ERR << "while depositing " + strTexture + " unconsolidated sediment for coast " << nCoast << " polygon " << pPolygon->nGetPolygonCoastID() << ", could not find the seaward end point of the up-coast profile (" << nUpCoastProfile << ") for depth of closure = " << m_dDepthOfClosure << endl;

      return RTN_ERR_NO_SEAWARD_END_OF_PROFILE_UPCOAST_BEACH_DEPOSITION;
   }

   // The part-profile length is one greater than nIndex, since nGetCellGivenDepth() returns the index of the cell at depth of closure. This will be the number of cells in the Dean profile portion of every parallel profile
   int const nUpCoastDeanLen = nIndex + 1;

   // assert(bIsWithinValidGrid(&PtiUpCoastPartProfileSeawardEnd));

   // Get the distance between the start and end of the part-profile (the Dean length), in external CRS units
   // CGeom2DPoint
   // PtUpCoastProfileStart = *pUpCoastProfile->pPtGetPointInProfile(0),
   // PtUpCoastProfileEnd = PtGridCentroidToExt(&PtiUpCoastPartProfileSeawardEnd);

   // double dUpCoastDeanLen = dGetDistanceBetween(&PtUpCoastProfileStart, &PtUpCoastProfileEnd);

   int const nUpCoastProfileCoastPoint = pUpCoastProfile->nGetCoastPoint();
   int const nDownCoastProfileCoastPoint = pDownCoastProfile->nGetCoastPoint();
   int const nXUpCoastProfileExistingCoastPoint = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nUpCoastProfileCoastPoint)->nGetX();
   int const nYUpCoastProfileExistingCoastPoint = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nUpCoastProfileCoastPoint)->nGetY();
   int nCoastSegLen;

   // Store the coast point numbers for this polygon so that we can shuffle them
   vector<int> nVCoastPoint;

   if (nDownCoastProfileCoastPoint == m_VCoast[nCoast].nGetCoastlineSize() - 1)
   {
      // This is the final down-coast polygon, so also include the down-coast polygon boundary
      nCoastSegLen = nDownCoastProfileCoastPoint - nUpCoastProfileCoastPoint + 1;

      for (int nCoastPoint = nUpCoastProfileCoastPoint; nCoastPoint <= nDownCoastProfileCoastPoint; nCoastPoint++)
         nVCoastPoint.push_back(nCoastPoint);
   }
   else
   {
      // This is not the final down-coast polygon, so do not include the polygon's down-coast boundary
      nCoastSegLen = nDownCoastProfileCoastPoint - nUpCoastProfileCoastPoint;

      for (int nCoastPoint = nUpCoastProfileCoastPoint; nCoastPoint < nDownCoastProfileCoastPoint; nCoastPoint++)
         nVCoastPoint.push_back(nCoastPoint);
   }

   double dStillToDepositOnPoly = dTargetToDepositOnPoly;
   double dTargetToDepositOnProfile = dTargetToDepositOnPoly / nCoastSegLen;
   double dStillToDepositOnProfile; // = dTargetToDepositOnProfile;

   // Shuffle the coast points, this is necessary so that leaving the loop does not create sequence-related artefacts
   shuffle(nVCoastPoint.begin(), nVCoastPoint.begin() + nCoastSegLen, m_Rand[1]);

   //    // Get the volume of sediment which is to be deposited on the polygon and on each parallel profile. Note that if dSandToMoveOnPoly is -ve, then don't do any sand deposition. Similarly if dCoarseToMoveOnPoly is -ve then don't do any coarse deposition
   // double
   // dSandToMoveOnPoly = pPolygon->dGetToDoBeachDepositionUnconsSand(),          // +ve is deposition, -ve is erosion
   // dCoarseToMoveOnPoly = pPolygon->dGetToDoBeachDepositionUnconsCoarse();      // +ve is deposition, -ve is erosion
   //
   // double
   // dTmpLower = 0,
   // dSandRatio = 0,
   // dCoarseRatio = 0,
   // dSandTargetToDepositOnPoly = 0,
   // dCoarseTargetToDepositOnPoly = 0;
   // if (dSandToMoveOnPoly > 0)
   // dTmpLower += dSandToMoveOnPoly;
   // if (dCoarseToMoveOnPoly > 0)
   // dTmpLower += dCoarseToMoveOnPoly;
   //
   // if ((dSandToMoveOnPoly > 0) && (dTmpLower > 0))
   // dSandRatio = dSandToMoveOnPoly / dTmpLower;
   // if (dCoarseToMoveOnPoly > 0)
   // dCoarseRatio = 1 - dSandRatio;
   //
   //    // LogStream << "####### nPoly = " << nPoly << " dSandRatio = " << dSandRatio << " dCoarseRatio = " << dCoarseRatio << endl;
   //
   // dSandTargetToDepositOnPoly = dAllTargetToDepositOnPoly * dSandRatio,
   // dCoarseTargetToDepositOnPoly = dAllTargetToDepositOnPoly * dCoarseRatio;
   //
   // double
   // dAllTargetPerProfile = dTargetToDepositOnPoly / nCoastSegLen,
   // dTargetStillToDepositOnProfile = dTargetToDepositOnPoly / nCoastSegLen,
   // dSandDepositedOnPoly = 0,                                              // Grand total for the whole polygon
   // dCoarseDepositedOnPoly = 0;                                            // Grand total for the whole polygon
   //
   //    // LogStream << "####### nPoly = " << nPoly << " dSandToMoveOnPoly = " << dSandToMoveOnPoly << " dCoarseToMoveOnPoly = " << dCoarseToMoveOnPoly << endl;
   //
   //    // LogStream << "####### nPoly = " << nPoly << " dSandTargetToDepositOnPoly = " << dSandTargetToDepositOnPoly << " dCoarseTargetToDepositOnPoly = " << dCoarseTargetToDepositOnPoly << " dAllTargetToDepositOnPoly = " << dAllTargetToDepositOnPoly << endl;

   // Now traverse the polygon's existing coastline in a random (but broadly DOWN-COAST i.e. increasing coast point indices) sequence, fitting a Dean profile at each coast point
   // DOWN-COAST (increasing coastpoint indices) ================================================================================
   // LogStream << "####### DOWN-COAST nPoly = " << nPoly << " nCoastSegLen = " << nCoastSegLen << endl;
   for (int n = 0; n < nCoastSegLen; n++)
   {
      // Pick a random coast point
      int const nCoastPoint = nVCoastPoint[n];
      CGeom2DIPoint const PtiCoastPoint = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint);
      int const nCoastX = PtiCoastPoint.nGetX();
      int const nCoastY = PtiCoastPoint.nGetY();

      // Calculate the x-y offset between this coast point, and the coast point of the up-coast normal
      int const nXOffset = nCoastX - nXUpCoastProfileExistingCoastPoint;
      int const nYOffset = nCoastY - nYUpCoastProfileExistingCoastPoint;
      int nSeawardOffset = -1;
      // nParProfLen;
      vector<CGeom2DIPoint> PtiVParProfile;
      vector<double> VdParProfileDeanElev;

      // OK, loop until we can deposit sufficient unconsolidated sediment on the parallel profile starting at this coast point
      while (true)
      {
         // Move seaward by one cell
         nSeawardOffset++;

         // And lengthen the parallel profile
         int nParProfLen = nUpCoastDeanLen + nSeawardOffset;

         if (nParProfLen > (pUpCoastProfile->nGetNumCellsInProfile()))
         {
            // We've reached the seaward end of the up-coast profile, need to quit
            // if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
            // LogStream << m_ulIter << ": " << WARN << "reached seaward end of up-coast profile during DOWN-COAST deposition of unconsolidated sediment for coast " << nCoast << " polygon " << pPolygon->nGetPolygonCoastID() << " (nCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << ")" << endl;

            break;
         }

         // Get the x-y coords of a profile starting from this coast point and parallel to the up-coast polygon boundary profile (these are in natural sequence, like the boundary part-profile)
         PtiVParProfile.resize(0);

         for (int m = 0; m < nParProfLen; m++)
         {
            CGeom2DIPoint const PtiProf = *pUpCoastProfile->pPtiGetCellInProfile(m);
            CGeom2DIPoint const PtiTmp(PtiProf.nGetX() + nXOffset, PtiProf.nGetY() + nYOffset);
            PtiVParProfile.push_back(PtiTmp);
         }

         // Get the existing elevation of the seaward end of the parallel profile
         int nSeaEndX = PtiVParProfile.back().nGetX();
         int nSeaEndY = PtiVParProfile.back().nGetY();

         // Safety check
         if (! bIsWithinValidGrid(nSeaEndX, nSeaEndY))
         {
            // if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
            // LogStream << WARN << "09 @@@@ while doing DOWN-COAST deposition on polygon " << nPoly << ", hit edge of grid at [" << nSeaEndX << "][" << nSeaEndY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its seaward end" << endl;

            KeepWithinValidGrid(nCoastX, nCoastY, nSeaEndX, nSeaEndY);
            PtiVParProfile.back().SetX(nSeaEndX);
            PtiVParProfile.back().SetY(nSeaEndY);
         }

         double const dParProfEndElev = m_pRasterGrid->m_Cell[nSeaEndX][nSeaEndY].dGetAllSedTopElevOmitTalus();

         // Set the start elevation for the Dean profile just a bit above mean SWL for this timestep (i.e. so that it is a Bruun profile)
         double const dParProfStartElev = m_dThisIterMeanSWL + m_dDeanProfileStartAboveSWL;

         // Calculate the total length of the parallel profile, including any seaward offset
         double const dParProfLen = dGetDistanceBetween(&PtiVParProfile.front(), &PtiVParProfile.back());

         // Now calculate the length of the Dean profile-only part i.e. without any seaward offset. The approach used here is approximate but probably OK
         double dParProfDeanLen = dParProfLen - (nSeawardOffset * m_dCellSide);

         // Safety check
         if (dParProfDeanLen <= 0)
            dParProfDeanLen = 1;

         // Solve for dA so that the existing elevations at the end of the parallel profile, and at the end of a Dean equilibrium profile on that part-normal, are the same
         double const dParProfA = (dParProfStartElev - dParProfEndElev) / pow(dParProfDeanLen, DEAN_POWER);

         nParProfLen = static_cast<int>(PtiVParProfile.size());
         VdParProfileDeanElev.resize(nParProfLen, 0);

         // for (int m = 0; m < static_cast<int>(PtiVParProfile.size()); m++)
         // LogStream << "[" << PtiVParProfile[m].nGetX() << "][" << PtiVParProfile[m].nGetY() << "] ";
         // LogStream << endl;

         double const dInc = dParProfDeanLen / (nParProfLen - nSeawardOffset - 2);

         // The elevation of the coast point in the Dean profile is the same as the elevation of the current coast point, ignoring any talus (since talus is assumed to be removed quickly) TODO 020 Is this correct? Should it be dParProfStartElev?
         double const dCoastElev = m_pRasterGrid->m_Cell[nCoastX][nCoastY].dGetAllSedTopElevOmitTalus();

         // For this depositing parallel profile, calculate the Dean equilibrium profile of the unconsolidated sediment h(y) = A * y^(2/3) where h(y) is the distance below the highest point in the profile at a distance y from the landward start of the profile
         CalcDeanProfile(&VdParProfileDeanElev, dInc, dParProfStartElev, dParProfA, true, nSeawardOffset, dCoastElev);

         double dParProfTotDiff = 0;

         for (int m = 0; m < nParProfLen; m++)
         {
            CGeom2DIPoint PtiTmp = PtiVParProfile[m];
            int nX = PtiTmp.nGetX();
            int nY = PtiTmp.nGetY();

            // Safety check
            if (! bIsWithinValidGrid(nX, nY))
            {
               // LogStream << WARN << "10 @@@@ while doing DOWN-COAST deposition on polygon " << nPoly << ", hit edge of grid at [" << nX << "][" << nY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its seaward end" << endl;

               KeepWithinValidGrid(nCoastX, nCoastY, nX, nY);
               PtiTmp.SetX(nX);
               PtiTmp.SetY(nY);
            }

            // Don't do anything to intervention cells
            if (bIsInterventionCell(nX, nY))
               continue;

            // Don't do cells twice
            if (! m_pRasterGrid->m_Cell[nX][nY].bBeachErosionOrDepositionThisIter())
            {
               double const dTmpElev = m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevOmitTalus();
               double const dDiff = VdParProfileDeanElev[m] - dTmpElev;

               dParProfTotDiff += dDiff;
            }
         }

         //          // DEBUG CODE -----------------------------------------------------
         // LogStream << endl << "\tFor polygon " << nPoly << " doing DOWN-COAST deposition, nSeawardOffset = " << nSeawardOffset << ", parallel profile from [" << PtiVParProfile[0].nGetX() << "][" << PtiVParProfile[0].nGetY() << "] to [" << PtiVParProfile.back().nGetX() << "][" << PtiVParProfile.back().nGetY() << "], nUpCoastDeanLen = " << nUpCoastDeanLen << " dUpCoastDeanLen = " << dUpCoastDeanLen << " nParProfLen = " << nParProfLen << " dParProfDeanLen = " << dParProfDeanLen << " dInc = " << dInc << " dParProfStartElev = " << dParProfStartElev << " dParProfEndElev = " << dParProfEndElev << " dParProfA = " << dParProfA << endl;
         //
         // LogStream << "\tExisting profile for deposition = ";
         // for (int n = 0; n < nParProfLen; n++)
         // {
         // CGeom2DIPoint PtiTmp = PtiVParProfile[n];
         // int
         // nX = PtiTmp.nGetX(),
         // nY = PtiTmp.nGetY();
         //
         //             // Safety check
         // if (! bIsWithinValidGrid(nX, nY))
         // KeepWithinValidGrid(nX, nY);
         //
         // LogStream << m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevOmitTalus() << " ";
         // }
         // LogStream << endl;
         // LogStream << "\tParallel Dean equilibrium profile for deposition = ";
         // for (int n = 0; n < nParProfLen; n++)
         // {
         // LogStream << VdParProfileDeanElev[n] << " ";
         // }
         // LogStream << endl;
         //
         // LogStream << "\tnCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << " difference = ";
         // for (int n = 0; n < nParProfLen; n++)
         // {
         // CGeom2DIPoint PtiTmp = PtiVParProfile[n];
         // int
         // nX = PtiTmp.nGetX(),
         // nY = PtiTmp.nGetY();
         //
         //             // Safety check
         // if (! bIsWithinValidGrid(nX, nY))
         // KeepWithinValidGrid(nX, nY);
         //
         // double
         // dTmpElev = m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevOmitTalus(),
         // dDiff = VdParProfileDeanElev[n] - dTmpElev;
         //
         // LogStream << dDiff << " ";
         // }
         // LogStream << endl;
         //          // END DEBUG CODE -----------------------------------------------------

         // So will we be able to deposit as much as is needed?
         if (dParProfTotDiff >= dTargetToDepositOnProfile)
         {
            // LogStream << "        DOWN-COAST nPoly = " << pPolygon->nGetPolygonCoastID() << " nCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << " dParProfTotDiff = " << dParProfTotDiff << " dTargetToDepositOnProfile = " << dTargetToDepositOnProfile << " WILL BE ABLE TO DEPOSIT ENOUGH" << endl;

            break;
         }
      }

      // assert(dParProfTotDiff > 0);

      // OK, this value of nSeawardOffset gives us enough deposition. So start depositing on the parallel profile from the coastward end
      double dDepositedOnProfile = 0; // Total for this parallel profile

      // LogStream << "        DOWN-COAST nPoly = " << nPoly << " nCoastPoint = " << nCoastPoint << " doing deposition on parallel profile, dSandTargetStillToDepositOnProfile = " << dSandTargetStillToDepositOnProfile << " dCoarseTargetToDepositOnProfile = " << dCoarseTargetToDepositOnProfile << endl;

      dStillToDepositOnProfile = dTargetToDepositOnProfile;

      for (unsigned int nSeawardFromCoast = 0; nSeawardFromCoast < PtiVParProfile.size(); nSeawardFromCoast++)
      {
         // Move along this parallel profile starting from the coast. Leave the loop if we have deposited enough on this polygon
         if (dStillToDepositOnPoly < SED_ELEV_TOLERANCE)
         {
            // LogStream << m_ulIter << ": nPoly = " << nPoly << " texture = " << strTexture << " DOWN-COAST nCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << " leaving loop because enough deposited on polygon, dTargetToDepositOnPoly = " << dTargetToDepositOnPoly << " dDepositedOnPoly = " << dDepositedOnPoly << endl;

            break;
         }

         // Leave the loop if we have deposited enough on this parallel profile
         if (dStillToDepositOnProfile < SED_ELEV_TOLERANCE)
         {
            // LogStream << m_ulIter << ": nPoly = " << nPoly << " texture = " << strTexture << " DOWN-COAST nCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << " leaving loop because enough deposited on parallel profile, dTargetStillToDepositOnProfile = " << dTargetStillToDepositOnProfile << " dDepositedOnProfile = " << dDepositedOnProfile << endl;

            break;
         }

         // assert(nSeawardFromCoast < PtiVParProfile.size());
         CGeom2DIPoint PtiTmp = PtiVParProfile[nSeawardFromCoast];
         int nX = PtiTmp.nGetX();
         int nY = PtiTmp.nGetY();

         // Safety check
         if (! bIsWithinValidGrid(nX, nY))
         {
            // LogStream << WARN << "11 @@@@ while doing DOWN-COAST deposition on polygon " << nPoly << ", hit edge of grid at [" << nX << "][" << nY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its seaward end" << endl;

            KeepWithinValidGrid(nCoastX, nCoastY, nX, nY);
            PtiTmp.SetX(nX);
            PtiTmp.SetY(nY);
         }

         // Don't do anything to intervention cells
         if (bIsInterventionCell(nX, nY))
            continue;

         // Don't do cells twice
         if (! m_pRasterGrid->m_Cell[nX][nY].bBeachErosionOrDepositionThisIter())
         {
            // Get this cell's current elevation
            double const dThisElevNow = m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevOmitTalus();

            // LogStream << "\tnPoly = " << nPoly << ", [" << nX << "][" << nY << "] nCoastPoint = " << nCoastPoint << " nSeawardFromCoast = " << nSeawardFromCoast << " dThisElevNow = " << dThisElevNow << " Dean Elev = " << VdParProfileDeanElev[nSeawardFromCoast] << endl;

            // Subtract the two elevations
            double const dElevDiff = VdParProfileDeanElev[nSeawardFromCoast] - dThisElevNow;
            CRWCellLandform* pLandform = m_pRasterGrid->m_Cell[nX][nY].pGetLandform();

            if (dElevDiff > SED_ELEV_TOLERANCE)
            {
               bool bDeposited = false;
               double dToDepositHere = 0;

               // The current elevation is below the Dean elevation, so we have can have beach deposition here
               if (dStillToDepositOnProfile > SED_ELEV_TOLERANCE)
               {
                  dToDepositHere = tMin(dElevDiff, dStillToDepositOnProfile, dStillToDepositOnPoly);

                  // LogStream << "        DOWN-COAST nPoly = " << nPoly << " nCoastPoint = " << nCoastPoint << " texture = " << strTexture << " dToDepositHere = " << dToDepositHere << endl;

                  int const nTopLayer = m_pRasterGrid->m_Cell[nX][nY].nGetNumOfTopLayerAboveBasement();

                  // Safety check
                  if (nTopLayer == INT_NODATA)
                     return RTN_ERR_NO_TOP_LAYER_DURING_BEACH_CALC;

                  if (dToDepositHere > SED_ELEV_TOLERANCE)
                  {
                     bDeposited = true;

                     if (nTexture == TEXTURE_SAND)
                     {
                        double const dSandNow = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->dGetSandDepth();

                        m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->SetSandDepth(dSandNow + dToDepositHere);

                        // Set the changed-this-timestep switch
                        m_bUnconsChangedThisIter[nTopLayer] = true;

                        dDepositedOnProfile += dToDepositHere;
                        dDepositedOnPoly += dToDepositHere;

                        // LogStream << "XXXXXXX texture = " << strTexture << " dDepositedOnProfile = " << dDepositedOnProfile << " dDepositedOnPoly = " << dDepositedOnPoly << endl;

                        // Update the cell's beach deposition, and total beach deposition, values
                        m_pRasterGrid->m_Cell[nX][nY].IncrBeachDeposition(dToDepositHere);

                        dStillToDepositOnPoly -= dToDepositHere;
                        dStillToDepositOnProfile -= dToDepositHere;

                        // if (dDepositedOnPoly > dTargetToDepositOnPoly)
                        // {
                        //    // LogStream << "££££££ 1 nPoly = " << nPoly << " nCoastPoint = " << nCoastPoint << " dDepositedOnProfile = " << dDepositedOnProfile << " dDepositedOnPoly = " << dDepositedOnPoly << " dTargetToDepositOnPoly = " << dTargetToDepositOnPoly << endl;
                        // }

                        // LogStream << m_ulIter << ": nPoly = " << nPoly << " texture = " << strTexture << " dToDepositHere = " << dToDepositHere << " at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " <<  dGridCentroidYToExtCRSY(nY) << "} nCoastPoint = " << nCoastPoint << " nSeawardFromCoast = " << nSeawardFromCoast << endl;
                     }
                     else if (nTexture == TEXTURE_COARSE)
                     {
                        double const dCoarseNow = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->dGetCoarseDepth();

                        m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->SetCoarseDepth(dCoarseNow + dToDepositHere);

                        // Set the changed-this-timestep switch
                        m_bUnconsChangedThisIter[nTopLayer] = true;

                        // if (dDepositedOnPoly > dTargetToDepositOnPoly)
                        // {
                        //    // LogStream << "££££££ 2 nPoly = " << nPoly << " dDepositedOnPoly = " << dDepositedOnPoly << " dTargetToDepositOnPoly = " << dTargetToDepositOnPoly << endl;
                        // }

                        dDepositedOnProfile += dToDepositHere;
                        dDepositedOnPoly += dToDepositHere;

                        // LogStream << "XXXXXXX texture = " << strTexture << " dDepositedOnProfile = " << dDepositedOnProfile << " dDepositedOnPoly = " << dDepositedOnPoly << endl;

                        // Update the cell's beach deposition, and total beach deposition, values
                        m_pRasterGrid->m_Cell[nX][nY].IncrBeachDeposition(dToDepositHere);

                        dStillToDepositOnPoly -= dToDepositHere;
                        dStillToDepositOnProfile -= dToDepositHere;

                        // LogStream << m_ulIter << ": nPoly = " << nPoly << " texture = " << strTexture << " dToDepositHere = " << dToDepositHere << " at [" << nX << "][" << nY << "] nCoastPoint = " << nCoastPoint << " nSeawardFromCoast = " << nSeawardFromCoast << endl;
                     }
                  }
               }

               // if (dDepositedOnPoly > dTargetToDepositOnPoly)
               // {
               //    // LogStream << "££££££ 3 nPoly = " << nPoly << " texture = " << strTexture << " dDepositedOnPoly = " << dDepositedOnPoly << " dTargetToDepositOnPoly = " << dTargetToDepositOnPoly << endl;
               // }

               if (bDeposited)
               {
                  // Update the cell's layer elevations
                  m_pRasterGrid->m_Cell[nX][nY].CalcAllLayerElevsAndD50();

                  // Update the cell's sea depth
                  m_pRasterGrid->m_Cell[nX][nY].SetSeaDepth();

                  // And set the landform category
                  int const nCat = pLandform->nGetLFCategory();

                  if ((nCat != LF_SEDIMENT_INPUT_CONSOLIDATED) && (nCat != LF_SEDIMENT_INPUT_UNCONSOLIDATED))
                     pLandform->SetLFCategory(LF_DRIFT_BEACH);

                  // Update this-timestep totals
                  m_ulThisIterNumBeachDepositionCells++;

                  if (nTexture == TEXTURE_SAND)
                     m_dThisIterBeachDepositionSand += dToDepositHere;

                  else if (nTexture == TEXTURE_COARSE)
                     m_dThisIterBeachDepositionCoarse += dToDepositHere;
               }
            }
            else if (bElevAboveDeanElev(nX, nY, dElevDiff, pLandform))
            {
               // The current elevation is higher than the Dean elevation, so we have potential beach erosion (i.e. not constrained by availability of unconsolidated sediment) here
               m_ulThisIterNumPotentialBeachErosionCells++;

               m_pRasterGrid->m_Cell[nX][nY].SetPotentialBeachErosion(-dElevDiff);

               // LogStream << m_ulIter << ": nPoly = " << nPoly << " texture = " << strTexture << " going DOWN-COAST, potential beach erosion = " << -dElevDiff << " at [" << nX << "][" << nY << "] nCoastPoint = " << nCoastPoint << " nSeawardFromCoast = " << nSeawardFromCoast << " dThisElevNow = " << dThisElevNow << " Dean Elev = " << VdParProfileDeanElev[nSeawardFromCoast] << endl;

               // Now get the number of the highest layer with non-zero thickness
               int const nThisLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopNonZeroLayerAboveBasement();

               // Safety check
               if (nThisLayer == INT_NODATA)
                  return RTN_ERR_NO_TOP_LAYER_DURING_BEACH_CALC;

               if (nThisLayer != NO_NONZERO_THICKNESS_LAYERS)
               {
                  // We still have at least one layer left with non-zero thickness (i.e. we are not down to basement), and the cell's current elevation is higher than the Dean equilibrium profile elevation. So do some beach erosion on this cell (Note: we are in nDoUnconsDepositionOnPolygon() still)
                  if (nTexture == TEXTURE_SAND)
                  {
                     double dSandRemoved = 0;
                     ErodeCellBeachSedimentSupplyLimited(nX, nY, nThisLayer, TEXTURE_SAND, -dElevDiff, dSandRemoved);

                     // Update totals for this parallel profile
                     dDepositedOnProfile -= dSandRemoved;
                     dStillToDepositOnProfile += dSandRemoved;

                     // Update totals for the polygon (note that these may become slightly -ve if erosion occurs early during this routine, but this is not a problem since the values will eventually become +ve again
                     dDepositedOnPoly -= dSandRemoved;
                     dStillToDepositOnPoly += dSandRemoved;

                     // if (dDepositedOnPoly < 0)
                     // LogStream << m_ulIter << ": BBB LESS THAN ZERO dDepositedOnPoly = " << dDepositedOnPoly << endl;

                     // Update this-timestep totals
                     m_ulThisIterNumActualBeachErosionCells++;

                     // // DEBUG CODE ==================================================================================================
                     // LogStream << m_ulIter << ": $$$$$$$$$ BEACH 1 Dean profile lower than existing profile, SAND depth eroded on [" << nX << "][" << nY << "] in Poly " << nPoly << " is " << dSandRemoved * m_dCellArea << endl;
                     // // DEBUG CODE ==================================================================================================

                     // Store the this-polygon depth of sand sediment eroded during Dean profile deposition of beach unconsolidated sediment
                     pPolygon->AddBeachSandErodedDeanProfile(dSandRemoved);
                  }
                  else if (nTexture == TEXTURE_COARSE)
                  {
                     double dCoarseRemoved = 0;
                     ErodeCellBeachSedimentSupplyLimited(nX, nY, nThisLayer, TEXTURE_COARSE, -dElevDiff, dCoarseRemoved);

                     // Update totals for this parallel profile
                     dDepositedOnProfile -= dCoarseRemoved;
                     dStillToDepositOnProfile += dCoarseRemoved;

                     // Update totals for the polygon (note that these may become slightly -ve if erosion occurs early during this routine, but this is not a problem since the values will eventually become +ve again
                     dDepositedOnPoly -= dCoarseRemoved;
                     dStillToDepositOnPoly += dCoarseRemoved;

                     // Update this-timestep totals
                     m_ulThisIterNumActualBeachErosionCells++;

                     // // DEBUG CODE ==================================================================================================
                     // LogStream << m_ulIter << ": $$$$$$$$$ BEACH 1 Dean profile lower than existing profile, COARSE depth eroded on [" << nX << "][" << nY << "] in Poly " << nPoly << " is " << dCoarseRemoved * m_dCellArea << endl;
                     // // DEBUG CODE ==================================================================================================

                     // Store the this-polygon depth of coarse sediment eroded during Dean profile deposition of beach unconsolidated sediment
                     pPolygon->AddBeachCoarseErodedDeanProfile(dCoarseRemoved);
                  }
               }

               // if (dDepositedOnPoly > dTargetToDepositOnPoly)
               // {
               // LogStream << "££££££ 4 nPoly = " << nPoly << " texture = " << strTexture << " dDepositedOnPoly = " << dDepositedOnPoly << " dTargetToDepositOnPoly = " << dTargetToDepositOnPoly << endl;
               // }

               // LogStream << m_ulIter << ": in polygon " << nPoly << " have net deposition for " << strTexture << ", actual beach erosion = " << dSand + dCoarse << " at [" << nX << "][" << nY << "]" << endl;
            }
         }
      }

      // LogStream << m_ulIter << ": nPoly = " << nPoly << " texture = " << strTexture << " nCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << " dTargetStillToDepositOnProfile = " << dTargetStillToDepositOnProfile << endl;
   }

   // Have we deposited enough?
   if (dTargetToDepositOnPoly > 0)
   {
      // LogStream << "##### texture = " << strTexture << " dDepositedOnPoly = " << dDepositedOnPoly << " dTargetToDepositOnPoly = " << dTargetToDepositOnPoly << " Not enough deposited, so now going UP-COAST" << endl;
      // No, so do the same in an UP-COAST (i.e. decreasing coastpoint indices) direction
      // UP-COAST (decreasing coastpoint indices) ==================================================================================

      // Start by finding the seaward end point of the down-coast part-profile, as before we are using only part of each profile, seaward as far as the depth of closure
      // CGeom2DIPoint PtiDownCoastPartProfileSeawardEnd;
      // int nIndex1 = pDownCoastProfile->nGetCellGivenDepth(m_pRasterGrid, m_dDepthOfClosure, &PtiDownCoastPartProfileSeawardEnd);
      int const nIndex1 = pDownCoastProfile->nGetCellGivenDepth(m_pRasterGrid, m_dDepthOfClosure);

      if (nIndex1 == INT_NODATA)
      {
         LogStream << m_ulIter << ": " << ERR << "while depositing beach for coast " << nCoast << " polygon " << pPolygon->nGetPolygonCoastID() << ", could not find the seaward end point of the up-coast profile (" << nUpCoastProfile << ") for depth of closure = " << m_dDepthOfClosure << endl;

         return RTN_ERR_NO_SEAWARD_END_OF_PROFILE_DOWNCOAST_BEACH_DEPOSITION;
      }

      // The part-profile length is one greater than nIndex1, since nGetCellGivenDepth() returns the index of the cell at depth of closure. This will be the number of cells in the Dean profile portion of every parallel profile
      int const nDownCoastDeanLen = nIndex1 + 1;

      // assert(bIsWithinValidGrid(&PtiDownCoastPartProfileSeawardEnd));

      // Get the distance between the start and end of the part-profile (the Dean length), in external CRS units
      // CGeom2DPoint
      // PtDownCoastProfileStart = *pDownCoastProfile->pPtGetPointInProfile(0),
      // PtDownCoastProfileEnd = PtGridCentroidToExt(&PtiDownCoastPartProfileSeawardEnd);

      // double dDownCoastDeanLen = dGetDistanceBetween(&PtDownCoastProfileStart, &PtDownCoastProfileEnd);

      int const nXDownCoastProfileExistingCoastPoint = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nDownCoastProfileCoastPoint)->nGetX();
      int const nYDownCoastProfileExistingCoastPoint = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nDownCoastProfileCoastPoint)->nGetY();

      // int nCoastSegLen;

      // Store the coast point numbers for this polygon so that we can shuffle them
      nVCoastPoint.resize(0);

      if (nUpCoastProfileCoastPoint == 0)
      {
         // This is the final up-coast polygon, so also include the up-coast polygon boundary
         nCoastSegLen = nDownCoastProfileCoastPoint - nUpCoastProfileCoastPoint + 1;

         for (int nCoastPoint = nUpCoastProfileCoastPoint; nCoastPoint <= nDownCoastProfileCoastPoint; nCoastPoint++)
            nVCoastPoint.push_back(nCoastPoint);
      }
      else
      {
         // This is not the final down-coast polygon, so do not include the polygon's down-coast boundary
         nCoastSegLen = nDownCoastProfileCoastPoint - nUpCoastProfileCoastPoint;

         for (int nCoastPoint = nUpCoastProfileCoastPoint; nCoastPoint < nDownCoastProfileCoastPoint; nCoastPoint++)
            nVCoastPoint.push_back(nCoastPoint);
      }

      // Shuffle the coast points, this is necessary so that leaving the loop does not create sequence-related artefacts
      shuffle(nVCoastPoint.begin(), nVCoastPoint.begin() + nCoastSegLen, m_Rand[1]);

      // Recalc the targets for deposition per profile
      dTargetToDepositOnProfile = dStillToDepositOnPoly / nCoastSegLen;

      // Set a switch
      bool bEnoughDepositedOnPolygon = false;

      // Now traverse the polygon's existing coastline in a random (but broadly UP-COAST, i.e. decreasing coastpoint indices) sequence, fitting a Dean profile at each coast point
      for (int n = 0; n < nCoastSegLen; n++)
      {
         // Check the switch
         if (bEnoughDepositedOnPolygon)
            break;

         // Pick a random coast point
         int const nCoastPoint = nVCoastPoint[n];
         CGeom2DIPoint const PtiCoastPoint = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint);
         int const nCoastX = PtiCoastPoint.nGetX();
         int const nCoastY = PtiCoastPoint.nGetY();

         // Calculate the x-y offset between this coast point, and the coast point of the down-coast normal
         int const nXOffset = nCoastX - nXDownCoastProfileExistingCoastPoint;
         int const nYOffset = nCoastY - nYDownCoastProfileExistingCoastPoint;
         int nSeawardOffset = -1;
         // nParProfLen;
         vector<CGeom2DIPoint> PtiVParProfile;
         vector<double> VdParProfileDeanElev;

         // OK, loop until we can deposit sufficient unconsolidated sediment
         while (true)
         {
            // Move seaward by one cell
            nSeawardOffset++;

            // And lengthen the parallel profile
            int nParProfLen = nDownCoastDeanLen + nSeawardOffset;

            if (nParProfLen > (pDownCoastProfile->nGetNumCellsInProfile()))
            {
               // We've reached the seaward end of the down-coast profile, need to quit
               // if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
               // LogStream << m_ulIter << ": " << WARN << "reached seaward end of down-coast profile during UP-COAST deposition of unconsolidated sediment for coast " << nCoast << " polygon " << nPoly << " (nCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << ")" << endl;

               break;
            }

            // Get the x-y coords of a profile starting from this coast point and parallel to the down-coast polygon boundary profile (these are in natural sequence, like the boundary part-profile)
            PtiVParProfile.resize(0);

            for (int m = 0; m < nParProfLen; m++)
            {
               CGeom2DIPoint const PtiProf = *pDownCoastProfile->pPtiGetCellInProfile(m);
               CGeom2DIPoint const PtiTmp(PtiProf.nGetX() + nXOffset, PtiProf.nGetY() + nYOffset);
               PtiVParProfile.push_back(PtiTmp);
            }

            // Get the existing elevation of the seaward end of the parallel profile
            int nSeaEndX = PtiVParProfile.back().nGetX();
            int nSeaEndY = PtiVParProfile.back().nGetY();

            // Safety check
            if (! bIsWithinValidGrid(nSeaEndX, nSeaEndY))
            {
               // LogStream << WARN << "12 @@@@ while doing UP-COAST deposition on polygon " << nPoly << ", hit edge of grid at [" << nSeaEndX << "][" << nSeaEndY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its seaward end" << endl;

               KeepWithinValidGrid(nCoastX, nCoastY, nSeaEndX, nSeaEndY);
               PtiVParProfile.back().SetX(nSeaEndX);
               PtiVParProfile.back().SetY(nSeaEndY);
            }

            double const dParProfEndElev = m_pRasterGrid->m_Cell[nSeaEndX][nSeaEndY].dGetAllSedTopElevOmitTalus();

            // Set the start elevation for the Dean profile just a bit above mean SWL for this timestep (i.e. so that it is a Bruun profile)
            double const dParProfStartElev = m_dThisIterMeanSWL + m_dDeanProfileStartAboveSWL;

            // Calculate the total length of the parallel profile, including any seaward offset
            double const dParProfLen = dGetDistanceBetween(&PtiVParProfile.front(), &PtiVParProfile.back());

            // Now calculate the length of the Dean profile-only part i.e. without any seaward offset. The approach used here is approximate but probably OK
            double dParProfDeanLen = dParProfLen - (nSeawardOffset * m_dCellSide);

            // Safety check
            if (dParProfDeanLen <= 0)
               dParProfDeanLen = 1;

            // Solve for dA so that the existing elevations at the end of the parallel profile, and at the end of a Dean equilibrium profile on that part-normal, are the same
            double const dParProfA = (dParProfStartElev - dParProfEndElev) / pow(dParProfDeanLen, DEAN_POWER);

            nParProfLen = static_cast<int>(PtiVParProfile.size());
            VdParProfileDeanElev.resize(nParProfLen, 0);

            double const dInc = dParProfDeanLen / (nParProfLen - nSeawardOffset - 2);

            // The elevation of the coast point in the Dean profile is the same as the elevation of the current coast point TODO 020 Is this correct? Should it be dParProfStartElev?
            double const dCoastElev = m_pRasterGrid->m_Cell[nCoastX][nCoastY].dGetAllSedTopElevOmitTalus();

            // For this depositing parallel profile, calculate the Dean equilibrium profile of the unconsolidated sediment h(y) = A * y^(2/3) where h(y) is the distance below the highest point in the profile at a distance y from the landward start of the profile
            CalcDeanProfile(&VdParProfileDeanElev, dInc, dParProfStartElev, dParProfA, true, nSeawardOffset, dCoastElev);

            double dParProfTotDiff = 0;

            for (int m = 0; m < nParProfLen; m++)
            {
               CGeom2DIPoint PtiTmp = PtiVParProfile[m];
               int nX = PtiTmp.nGetX();
               int nY = PtiTmp.nGetY();

               // Safety check
               if (! bIsWithinValidGrid(nX, nY))
               {
                  // LogStream << WARN << "13 @@@@ while constructing parallel profile for UP-COAST deposition on polygon " << nPoly << ", hit edge of grid at [" << nX << "][" << nY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its seaward end" << endl;

                  KeepWithinValidGrid(nCoastX, nCoastY, nX, nY);
                  PtiTmp.SetX(nX);
                  PtiTmp.SetY(nY);
               }

               // Don't do anything to intervention cells
               if (bIsInterventionCell(nX, nY))
                  continue;

               // Don't do cells twice
               if (! m_pRasterGrid->m_Cell[nX][nY].bBeachErosionOrDepositionThisIter())
               {
                  double const dTmpElev = m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevOmitTalus();
                  double const dDiff = VdParProfileDeanElev[m] - dTmpElev;

                  dParProfTotDiff += dDiff;
               }
            }

            //             // DEBUG CODE -----------------------------------------------------
            // LogStream << endl << "\tFor polygon " << nPoly << " doing UP-COAST deposition, nSeawardOffset = " << nSeawardOffset << ", parallel profile from [" << PtiVParProfile[0].nGetX() << "][" << PtiVParProfile[0].nGetY() << "] to [" << PtiVParProfile.back().nGetX() << "][" << PtiVParProfile.back().nGetY() << "], nDownCoastDeanLen = " << nDownCoastDeanLen << " dDownCoastDeanLen = " << dDownCoastDeanLen << " nParProfLen = " << nParProfLen << " dParProfDeanLen = " << dParProfDeanLen << " dInc = " << dInc << " dParProfStartElev = " << dParProfStartElev << " dParProfEndElev = " << dParProfEndElev << " dParProfA = " << dParProfA << endl;
            //
            // LogStream << "\tExisting profile for deposition = ";
            // for (int n = 0; n < nParProfLen; n++)
            // {
            // CGeom2DIPoint PtiTmp = PtiVParProfile[n];
            // int
            // nX = PtiTmp.nGetX(),
            // nY = PtiTmp.nGetY();
            //
            //                // Safety check
            // if (! bIsWithinValidGrid(nX, nY))
            // KeepWithinValidGrid(nX, nY);
            //
            // LogStream << m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevOmitTalus() << " ";
            // }
            // LogStream << endl;
            // LogStream << "\tParallel Dean equilibrium profile for deposition = ";
            // for (int n = 0; n < nParProfLen; n++)
            // {
            // LogStream << VdParProfileDeanElev[n] << " ";
            // }
            // LogStream << endl;
            //
            // LogStream << "\tnCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << " difference = ";
            // for (int n = 0; n < nParProfLen; n++)
            // {
            // CGeom2DIPoint PtiTmp = PtiVParProfile[n];
            // int
            // nX = PtiTmp.nGetX(),
            // nY = PtiTmp.nGetY();
            //
            //                // Safety check
            // if (! bIsWithinValidGrid(nX, nY))
            // KeepWithinValidGrid(nX, nY);
            //
            // double
            // dTmpElev = m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevOmitTalus(),
            // dDiff = VdParProfileDeanElev[n] - dTmpElev;
            //
            // LogStream << dDiff << " ";
            // }
            // LogStream << endl;
            //             // END DEBUG CODE -----------------------------------------------------

            // So will we be able to deposit as much as is needed?
            if (dParProfTotDiff >= dTargetToDepositOnProfile)
            {
               // LogStream << "        UP-COAST nCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << " dParProfTotDiff = " << dParProfTotDiff << " dTargetToDepositOnProfile = " << dTargetToDepositOnProfile << " WILL BE ABLE TO DEPOSIT ENOUGH" << endl;

               break;
            }
         }

         // assert(dParProfTotDiff > 0);

         // OK, this value of nSeawardOffset gives us enough deposition. So start depositing on the parallel profile from the coastward end
         double dDepositedOnProfile = 0; // Total for this parallel profile

         dStillToDepositOnProfile = dTargetToDepositOnProfile;

         for (unsigned int nSeawardFromCoast = 0; nSeawardFromCoast < PtiVParProfile.size(); nSeawardFromCoast++)
         {
            // Move along this parallel profile starting from the coast. Leave the loop if we have deposited enough on this polygon
            if (dStillToDepositOnPoly < SED_ELEV_TOLERANCE)
            {
               // LogStream << m_ulIter << ": nPoly = " << nPoly << " texture = " << strTexture << " UP-COAST nCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << " leaving loop because enough deposited on polygon, dTargetToDepositOnPoly = " << dTargetToDepositOnPoly << " dDepositedOnPoly = " << dDepositedOnPoly << endl;

               bEnoughDepositedOnPolygon = true;

               break;
            }

            // Leave the loop if we have deposited enough on this parallel profile
            if (dStillToDepositOnProfile < SED_ELEV_TOLERANCE)
            {
               // LogStream << m_ulIter << ": nPoly = " << nPoly << " texture = " << strTexture << " UP-COAST nCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << " leaving loop because enough deposited on parallel profile, dTargetStillToDepositOnProfile = " << dTargetStillToDepositOnProfile << " dDepositedOnProfile = " << dDepositedOnProfile << endl;

               break;
            }

            // assert(nSeawardFromCoast < PtiVParProfile.size());
            CGeom2DIPoint PtiTmp = PtiVParProfile[nSeawardFromCoast];
            int nX = PtiTmp.nGetX();
            int nY = PtiTmp.nGetY();

            // Safety check
            if (! bIsWithinValidGrid(nX, nY))
            {
               // LogStream << WARN << "14 @@@@ while constructing parallel profile for UP-COAST deposition on polygon " << nPoly << ", hit edge of grid at [" << nX << "][" << nY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its seaward end" << endl;

               KeepWithinValidGrid(nCoastX, nCoastY, nX, nY);
               PtiTmp.SetX(nX);
               PtiTmp.SetY(nY);
            }

            // Don't do anything to intervention cells
            if (bIsInterventionCell(nX, nY))
               continue;

            // Don't do cells twice
            if (! m_pRasterGrid->m_Cell[nX][nY].bBeachErosionOrDepositionThisIter())
            {
               // Get this cell's current elevation
               double const dThisElevNow = m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevOmitTalus();

               // LogStream << "\tnPoly = " << nPoly << " going UP-COAST, [" << nX << "][" << nY << "] nCoastPoint = " << nCoastPoint << " nSeawardFromCoast = " << nSeawardFromCoast << " dThisElevNow = " << dThisElevNow << " Dean Elev = " << VdParProfileDeanElev[nSeawardFromCoast] << endl;

               // Subtract the two elevations
               double const dElevDiff = VdParProfileDeanElev[nSeawardFromCoast] - dThisElevNow;
               CRWCellLandform* pLandform = m_pRasterGrid->m_Cell[nX][nY].pGetLandform();

               if (dElevDiff > SED_ELEV_TOLERANCE)
               {
                  bool bDeposited = false;
                  double dToDepositHere = 0;

                  // The current elevation is below the Dean elevation, so we have can have beach deposition here
                  if (dStillToDepositOnProfile > SED_ELEV_TOLERANCE)
                  {
                     dToDepositHere = tMin(dElevDiff, dStillToDepositOnProfile, dStillToDepositOnPoly);

                     // LogStream << "          UP-COAST nPoly = " << nPoly << " nCoastPoint = " << nCoastPoint << " texture = " << strTexture << " dToDepositHere = " << dToDepositHere << endl;

                     int const nTopLayer = m_pRasterGrid->m_Cell[nX][nY].nGetNumOfTopLayerAboveBasement();

                     // Safety check
                     if (nTopLayer == INT_NODATA)
                        return RTN_ERR_NO_TOP_LAYER_DURING_BEACH_CALC;

                     if (dToDepositHere > SED_ELEV_TOLERANCE)
                     {
                        bDeposited = true;

                        if (nTexture == TEXTURE_SAND)
                        {
                           double const dSandNow = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->dGetSandDepth();

                           m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->SetSandDepth(dSandNow + dToDepositHere);

                           // Set the changed-this-timestep switch
                           m_bUnconsChangedThisIter[nTopLayer] = true;

                           dDepositedOnProfile += dToDepositHere;
                           dDepositedOnPoly += dToDepositHere;

                           // LogStream << "XXXXXXX texture = " << strTexture << " dDepositedOnProfile = " << dDepositedOnProfile << " dDepositedOnPoly = " << dDepositedOnPoly << endl;

                           // Update the cell's beach deposition, and total beach deposition, values
                           m_pRasterGrid->m_Cell[nX][nY].IncrBeachDeposition(dToDepositHere);

                           dStillToDepositOnPoly -= dToDepositHere;
                           dStillToDepositOnProfile -= dToDepositHere;

                           // if (dDepositedOnPoly > dTargetToDepositOnPoly)
                           // {
                           //    // LogStream << "££££££ 5 nPoly = " << nPoly << " nCoastPoint = " << nCoastPoint << " dDepositedOnProfile = " << dDepositedOnProfile << " dDepositedOnPoly = " << dDepositedOnPoly << " dTargetToDepositOnPoly = " << dTargetToDepositOnPoly << endl;
                           // }

                           // LogStream << m_ulIter << ": nPoly = " << nPoly << " texture = " << strTexture << " dToDepositHere = " << dToDepositHere << " at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " <<  dGridCentroidYToExtCRSY(nY) << "} nCoastPoint = " << nCoastPoint << " nSeawardFromCoast = " << nSeawardFromCoast << endl;
                        }
                        else if (nTexture == TEXTURE_COARSE)
                        {
                           double const dCoarseNow = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->dGetCoarseDepth();

                           m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->SetCoarseDepth(dCoarseNow + dToDepositHere);

                           // Set the changed-this-timestep switch
                           m_bUnconsChangedThisIter[nTopLayer] = true;

                           if (dDepositedOnPoly > dTargetToDepositOnPoly)
                           {
                              // LogStream << "££££££ 6 nPoly = " << nPoly << " dDepositedOnPoly = " << dDepositedOnPoly << " dTargetToDepositOnPoly = " << dTargetToDepositOnPoly << endl;
                           }

                           dDepositedOnProfile += dToDepositHere;
                           dDepositedOnPoly += dToDepositHere;

                           // LogStream << "XXXXXXX texture = " << strTexture << " dDepositedOnProfile = " << dDepositedOnProfile << " dDepositedOnPoly = " << dDepositedOnPoly << endl;

                           // Update the cell's beach deposition, and total beach deposition, values
                           m_pRasterGrid->m_Cell[nX][nY].IncrBeachDeposition(dToDepositHere);

                           dStillToDepositOnPoly -= dToDepositHere;
                           dStillToDepositOnProfile -= dToDepositHere;

                           // LogStream << m_ulIter << ": nPoly = " << nPoly << " texture = " << strTexture << " dToDepositHere = " << dToDepositHere << " at [" << nX << "][" << nY << "] nCoastPoint = " << nCoastPoint << " nSeawardFromCoast = " << nSeawardFromCoast << endl;
                        }
                     }
                  }

                  if (bDeposited)
                  {
                     // Now update the cell's layer elevations
                     m_pRasterGrid->m_Cell[nX][nY].CalcAllLayerElevsAndD50();

                     // Update the cell's sea depth
                     m_pRasterGrid->m_Cell[nX][nY].SetSeaDepth();

                     int const nCat = pLandform->nGetLFCategory();

                     if ((nCat != LF_SEDIMENT_INPUT_CONSOLIDATED) && (nCat != LF_SEDIMENT_INPUT_UNCONSOLIDATED))
                        pLandform->SetLFCategory(LF_DRIFT_BEACH);

                     // Update this-timestep totals
                     m_ulThisIterNumBeachDepositionCells++;

                     if (nTexture == TEXTURE_SAND)
                        m_dThisIterBeachDepositionSand += dToDepositHere;

                     else if (nTexture == TEXTURE_COARSE)
                        m_dThisIterBeachDepositionCoarse += dToDepositHere;
                  }
               }
               else if (bElevAboveDeanElev(nX, nY, dElevDiff, pLandform))
               {
                  // The current elevation is higher than the Dean elevation, so we could have beach erosion here
                  m_ulThisIterNumPotentialBeachErosionCells++;

                  m_pRasterGrid->m_Cell[nX][nY].SetPotentialBeachErosion(-dElevDiff);

                  // LogStream << m_ulIter << ": nPoly = " << nPoly << " texture = " << strTexture << " going UP-COAST, potential beach erosion = " << -dElevDiff << " at [" << nX << "][" << nY << "] nCoastPoint = " << nCoastPoint << " nSeawardFromCoast = " << nSeawardFromCoast << " dThisElevNow = " << dThisElevNow << " Dean Elev = " << VdParProfileDeanElev[nSeawardFromCoast] << endl;

                  // Now get the number of the highest layer with non-zero thickness
                  int const nThisLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopNonZeroLayerAboveBasement();

                  // Safety check
                  if (nThisLayer == INT_NODATA)
                     return RTN_ERR_NO_TOP_LAYER_DURING_BEACH_CALC;

                  if (nThisLayer != NO_NONZERO_THICKNESS_LAYERS)
                  {
                     // We still have at least one layer left with non-zero thickness (i.e. we are not down to basement), and the cell's current elevation is higher than the Dean equilibrium profile elevation. So do some beach erosion on this cell (Note: we are in nDoUnconsDepositionOnPolygon() still)
                     if (nTexture == TEXTURE_SAND)
                     {
                        double dSandRemoved = 0;
                        ErodeCellBeachSedimentSupplyLimited(nX, nY, nThisLayer, TEXTURE_SAND, -dElevDiff, dSandRemoved);

                        // Update total for this parallel profile
                        dDepositedOnProfile -= dSandRemoved;
                        dStillToDepositOnProfile += dSandRemoved;

                        // Update totals for the polygon (note that these may become slightly -ve if erosion occurs early during this routine, but this is not a problem since the values will eventually become +ve again
                        dDepositedOnPoly -= dSandRemoved;
                        dStillToDepositOnPoly += dSandRemoved;

                        // if (dDepositedOnPoly < 0)
                        // LogStream << m_ulIter << ": AAA LESS THAN ZERO dDepositedOnPoly = " << dDepositedOnPoly << endl;

                        // Update this-timestep totals
                        m_ulThisIterNumActualBeachErosionCells++;

                        // // DEBUG CODE ==============================================================================================
                        // LogStream << m_ulIter << ": $$$$$$$$$ BEACH 3 Dean profile lower than existing profile, SAND depth eroded on [" << nX << "][" << nY << "] in Poly " << nPoly << " is " << dSandRemoved * m_dCellArea << endl;
                        // // DEBUG CODE ==============================================================================================

                        // Store the this-polygon depth of sand sediment eroded during Dean profile deposition of beach unconsolidated sediment
                        pPolygon->AddBeachSandErodedDeanProfile(dSandRemoved);
                     }
                     else if (nTexture == TEXTURE_COARSE)
                     {
                        double dCoarseRemoved = 0;
                        ErodeCellBeachSedimentSupplyLimited(nX, nY, nThisLayer, TEXTURE_COARSE, -dElevDiff, dCoarseRemoved);

                        // Update total for this parallel profile
                        dDepositedOnProfile -= dCoarseRemoved;
                        dStillToDepositOnProfile += dCoarseRemoved;

                        // Update totals for the polygon (note that these may become slightly -ve if erosion occurs early during this routine, but this is not a problem since the values will eventually become +ve again
                        dDepositedOnPoly -= dCoarseRemoved;
                        dStillToDepositOnPoly += dCoarseRemoved;

                        // Update this-timestep totals
                        m_ulThisIterNumActualBeachErosionCells++;

                        // // DEBUG CODE ==============================================================================================
                        // LogStream << m_ulIter << ": $$$$$$$$$ BEACH 3 Dean profile lower than existing profile, COARSE depth eroded on [" << nX << "][" << nY << "] in Poly " << nPoly << " is " << dCoarseRemoved * m_dCellArea << endl;
                        // // DEBUG CODE ==============================================================================================

                        // Store the this-polygon depth of coarse sediment eroded during Dean profile deposition of beach unconsolidated sediment
                        pPolygon->AddBeachCoarseErodedDeanProfile(dCoarseRemoved);
                     }
                  }

                  // if (dDepositedOnPoly > dTargetToDepositOnPoly)
                  // {
                  //    // LogStream << "££££££ 4 nPoly = " << nPoly << " texture = " << strTexture << " dDepositedOnPoly = " << dDepositedOnPoly << " dTargetToDepositOnPoly = " << dTargetToDepositOnPoly << endl;
                  // }

                  // LogStream << m_ulIter << ": in polygon " << nPoly << " have net deposition for " << strTexture << ", actual beach erosion = " << dSand + dCoarse << " at [" << nX << "][" << nY << "]" << endl;
               }
            }
         }
      }
   }

   // Store amounts deposited on this polygon
   if (nTexture == TEXTURE_SAND)
   {
      pPolygon->SetBeachDepositionUnconsSand(dDepositedOnPoly);

      // Check mass balance for sand deposited
      if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
         if (! bFPIsEqual(pPolygon->dGetToDoBeachDepositionUnconsSand(), dDepositedOnPoly, MASS_BALANCE_TOLERANCE))
            LogStream << m_ulIter << ":\t coast " << nCoast << " polygon " << pPolygon->nGetPolygonCoastID() << " NOT equal SAND dGetToDoBeachDepositionUnconsSand() = " << pPolygon->dGetToDoBeachDepositionUnconsSand() << " dDepositedOnPoly = " << dDepositedOnPoly << endl;

      pPolygon->SetZeroToDoDepositionUnconsSand();
   }
   else if (nTexture == TEXTURE_COARSE)
   {
      pPolygon->SetBeachDepositionUnconsCoarse(dDepositedOnPoly);

      // Check mass balance for coarse deposited
      if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
         if (! bFPIsEqual(pPolygon->dGetToDoBeachDepositionUnconsCoarse(), dDepositedOnPoly, MASS_BALANCE_TOLERANCE))
            LogStream << m_ulIter << ":\t coast " << nCoast << " polygon " << pPolygon->nGetPolygonCoastID() << " NOT equal COARSE dGetToDoBeachDepositionUnconsCoarse() = " << pPolygon->dGetToDoBeachDepositionUnconsCoarse() << " dDepositedOnPoly = " << dDepositedOnPoly << endl;

      pPolygon->SetZeroToDoDepositionUnconsCoarse();
   }

   // // DEBUG CODE #####################
   // if (m_ulIter == 5)
   // {
   // LogStream << m_ulIter << ": leaving nDoUnconsDepositionOnPolygon() nCoast = " << nCoast << " nPoly = " << nPoly << " strTexture = " << strTexture << " dTargetToDepositOnPoly = " << dTargetToDepositOnPoly * m_dCellArea << " dDepositedOnPoly = " << dDepositedOnPoly * m_dCellArea << endl;
   //
   // double dTmpSandUncons = 0;
   // for (int nX1 = 0; nX1 < m_nXGridSize; nX1++)
   // {
   // for (int nY1 = 0; nY1 < m_nYGridSize; nY1++)
   // {
   // dTmpSandUncons += m_pRasterGrid->m_Cell[nX1][nY1].dGetUnconsSandDepthAllLayers();
   // }
   // }
   //
   // LogStream << m_ulIter << ": leaving nDoUnconsDepositionOnPolygon() for nPoly = " << nPoly << " TOTAL UNCONSOLIDATED SAND ON ALL CELLS = " << dTmpSandUncons * m_dCellArea << endl;
   // }
   // // DEBUG CODE #####################

   return RTN_OK;
}

//===============================================================================================================================
//! Return true if the given elevation is higher than the Dean elevation (and other land use conditions are met), which means that we could have beach erosion on this cell
//===============================================================================================================================
bool CSimulation::bElevAboveDeanElev(int const nX, int const nY, double const dElevDiff, CRWCellLandform const* pLandform)
{
   // TODO 075 What if it is bedrock that sticks above Dean profile?
   if (dElevDiff <= 0)
   {
      if (m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea())
         return true;

      int const nCat = pLandform->nGetLFCategory();

      if ((nCat == LF_DRIFT_BEACH) || (nCat == LF_DRIFT_DUNES) || (nCat == LF_SEDIMENT_INPUT_UNCONSOLIDATED))
         return true;
   }

   return false;
}
