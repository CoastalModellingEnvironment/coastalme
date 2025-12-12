/*!
   \file multiple_coastlines.cpp
   \brief Routines relating to ukltiple coastlines
   \details TODO 001 A more detailed description of these routines.
   \author David Favis-Mortlock
   \author Andres Payo
   \author Wilf Chun
   \date 2025
   \copyright GNU General Public License
*/

/* ===============================================================================================================================

   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

===============================================================================================================================*/
#include <assert.h>

#include <iostream>
using std::endl;

#include <algorithm>
using std::shuffle;

#include "cme.h"
#include "simulation.h"
#include "coast.h"
#include "2di_point.h"
#include "2d_point.h"
#include "multi_line.h"

class CRWCoast;      // Forward declaration

//===============================================================================================================================
//! Checks all profiles on all coasts for intersections between profiles belonging to different coasts
//===============================================================================================================================
int CSimulation::nDoMultipleCoastlines(void)
{
   int const nCoastSize = static_cast<int>(m_VCoast.size());

   if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
      LogStream << endl << m_ulIter << ": " << nCoastSize << " coastlines found, checking coast-coast intersections" << endl;

   // Create a vector of coast IDs and randomly shuffle it to avoid sequence-related artefacts
   vector<int> VnCoastID(nCoastSize);
   for (int nn = 0; nn < nCoastSize; nn++)
      VnCoastID[nn] = nn;
   shuffle(VnCoastID.begin(), VnCoastID.end(), m_Rand[1]);

   // Check all coastlines
   for (int nn = 0; nn < nCoastSize; nn++)
   {
      // Use the shuffled coast ID
      int const nCoast = VnCoastID[nn];

      // Create a vector of this coast's profile IDs and randomly shuffle it to avoid sequence-related artefacts
      vector<int> VnProfileID(m_VCoast[nCoast].nGetNumProfiles());
      for (int m = 0; m < m_VCoast[nCoast].nGetNumProfiles(); m++)
         VnProfileID[m] = m;
      shuffle(VnProfileID.begin(), VnProfileID.end(), m_Rand[1]);

      // Check all profiles
      for (int n = 0; n < m_VCoast[nCoast].nGetNumProfiles(); n++)
      {
         // Use the shuffled profile ID
         CGeomProfile* pProfile = m_VCoast[nCoast].pGetProfileWithDownCoastSeq(VnProfileID[n]);
         int const nProfile = pProfile->nGetProfileID();

         // Check every cell that is 'under' this profile, start one cell seaward of coastline
         for (int nCell = 1; nCell < pProfile->nGetNumCellsInProfile(); nCell++)
         {
            CGeom2DIPoint const* pCell = pProfile->pPtiGetCellInProfile(nCell);
            int const nX = pCell->nGetX();
            int const nY = pCell->nGetY();

            // Have we hit a cell which is 'under' another coastline?
            if (m_pRasterGrid->m_Cell[nX][nY].bIsCoastline())
            {
               // Yes this is a coastline cell, as well as a coast-normal profile cell
               int const nHitCoast = m_pRasterGrid->m_Cell[nX][nY].nGetCoastline();
               if (nHitCoast != nCoast)
               {
                  // We have hit a different coastline, so truncate this profile
                  int const nRtn = nTruncateProfileHitDifferentCoast(nCoast, nProfile, nX, nY);
                  if (nRtn != RTN_OK)
                     return nRtn;
               }
            }
            else
            {
               // We also need to check an adjacent cell, doesn't matter which one
               int nYTmp = nY+1;
               if (nY+1 >= m_nYGridSize)
                  nYTmp = nY-1;

               if (m_pRasterGrid->m_Cell[nX][nYTmp].bIsCoastline())
               {
                  // Yes this is a coastline cell, as well as a coast-normal profile cell
                  int const nHitCoast = m_pRasterGrid->m_Cell[nX][nYTmp].nGetCoastline();
                  if (nHitCoast != nCoast)
                  {
                     // We have hit a different coastline, so truncate this profile
                     int const nRtn = nTruncateProfileHitDifferentCoast(nCoast, nProfile, nX, nY);
                     if (nRtn != RTN_OK)
                        return nRtn;
                  }
               }
            }

            // For grid-edge cells, don't check for intersection with other-coast profiles, instead wait until this profile hits a different coast
            if (! pProfile->bIsGridEdge())
            {
               // Have we hit a cell which is 'under' a coast-normal profile belonging to another coast? NOTE Is a problem if get more than two coast normals passing through this cell
               int nHitProfileCoast = m_pRasterGrid->m_Cell[nX][nY].nGetProfileCoastID();
               if ((nHitProfileCoast != INT_NODATA) && (nHitProfileCoast != nCoast))
               {
                  // Yes, we have hit a profile which belongs to a different coast
                  int const nHitProfile = m_pRasterGrid->m_Cell[nX][nY].nGetProfileID();

                  // Safety check
                  if (nHitProfile != INT_NODATA)
                  {
                     // Truncate both profiles
                     int const nRtn = nTruncateProfilesDifferentCoasts(nCoast, nProfile, nHitProfileCoast, nHitProfile, nX, nY);
                     if (nRtn != RTN_OK)
                        return nRtn;
                  }
               }
               else
               {
                  // Try again with an adjacent point, does not matter wehich one NOTE Is a problem if get more than two coast normals passing through this cell
                  int nYTmp = nY+1;
                  if (nY+1 >= m_nYGridSize)
                     nYTmp = nY-1;

                  nHitProfileCoast = m_pRasterGrid->m_Cell[nX][nYTmp].nGetProfileCoastID();
                  if ((nHitProfileCoast != INT_NODATA) && (nHitProfileCoast != nCoast))
                  {
                     // Yes, we have hit a profile which belongs to a different coast
                     int const nHitProfile = m_pRasterGrid->m_Cell[nX][nYTmp].nGetProfileID();

                     // Safety check
                     if (nHitProfile != INT_NODATA)
                     {
                        // Truncate both profiles
                        int const nRtn = nTruncateProfilesDifferentCoasts(nCoast, nProfile, nHitProfileCoast, nHitProfile, nX, nYTmp);
                        if (nRtn != RTN_OK)
                           return nRtn;
                     }
                  }
               }
            }
         }
      }
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
      // // DEBUG CODE ================
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Truncates two intersecting coast-normal profile belonging to different coasts
//===============================================================================================================================
int CSimulation::nTruncateProfilesDifferentCoasts(int const nThisProfileCoast, int const nThisProfile, int const nHitProfileCoast, int const nHitProfile, int const nXIntersect, int const nYIntersect)
{
   // LogStream << m_ulIter << ": coast " << nThisProfileCoast << " profile " << nThisProfile << " and coast " << nHitProfileCoast << " profile "<< nHitProfile << " intersect at [" << nXIntersect << "][" << nYIntersect << "] = {" << dGridCentroidXToExtCRSX(nXIntersect) << ", " << dGridCentroidYToExtCRSY(nYIntersect) << "}" << endl;

   // OK, get pointers to 'this' profile and to the hit profile
   CGeomProfile* pThisProfile = m_VCoast[nThisProfileCoast].pGetProfile(nThisProfile);
   CGeomProfile* pHitProfile = m_VCoast[nHitProfileCoast].pGetProfile(nHitProfile);

   // Get the length of each profile (grid CRS)
   int const nThisProfileLen = pThisProfile->nGetNumCellsInProfile();
   int const nHitProfileLen = pHitProfile->nGetNumCellsInProfile();

   // Get the start points of both profiles
   CGeom2DIPoint const* pPtiThisProfileStart = pThisProfile->pPtiGetFirstCellInProfile();
   CGeom2DIPoint const* pPtiHitProfileStart = pHitProfile->pPtiGetFirstCellInProfile();

   double const dXThisProfileStart = pPtiThisProfileStart->nGetX();
   double const dYThisProfileStart = pPtiThisProfileStart->nGetY();
   double const dXHitProfileStart = pPtiHitProfileStart->nGetX();
   double const dYHitProfileStart = pPtiHitProfileStart->nGetY();
   double dClosestX;
   double dClosestY;

   // On a line joining the start points of both profiles, find the point on that line which is closest to the intersection of both profiles
   FindClosestPointOnStraightLine(dXThisProfileStart, dYThisProfileStart, dXHitProfileStart, dYHitProfileStart, nXIntersect, nYIntersect, dClosestX, dClosestY);

   // Get the distance between this closest point and each profile start point
   double const dThisDist = dGetDistanceBetween(pPtiThisProfileStart->nGetX(), pPtiThisProfileStart->nGetY(), dClosestX, dClosestY);
   double const dHitDist = dGetDistanceBetween(pPtiHitProfileStart->nGetX(), pPtiHitProfileStart->nGetY(), dClosestX, dClosestY);
   double const dTotDist = dThisDist + dHitDist;

   // Calculate the proportion of each profile that is to be retained
   double const dThisPropToRetain = dThisDist / dTotDist;
   double const dHitPropToRetain = dHitDist / dTotDist;

   // Now calculate the indices of the new endpoints for each profile
   int const nThisProfileEndpointIndex = tMax(nRound(nThisProfileLen * dThisPropToRetain) - GAP_BETWEEN_DIFFERENT_COAST_PROFILES, MIN_PROFILE_SIZE);
   int const nHitProfileEndpointIndex = tMax(nRound(nHitProfileLen * dHitPropToRetain) - GAP_BETWEEN_DIFFERENT_COAST_PROFILES, MIN_PROFILE_SIZE);

   // Get pointers to the cells in the two profiles
   vector<CGeom2DIPoint>* pVThisProfileCells = pThisProfile->pPtiVGetCellsInProfile();
   vector<CGeom2DIPoint>* pVHitProfileCells = pHitProfile->pPtiVGetCellsInProfile();

   // New endpoints of the two profiles (grid CRS)
   int const nXThisProfileEndPoint = pVThisProfileCells->at(nThisProfileEndpointIndex).nGetX();
   int const nYThisProfileEndPoint = pVThisProfileCells->at(nThisProfileEndpointIndex).nGetY();
   int const nXHitProfileEndPoint = pVHitProfileCells->at(nHitProfileEndpointIndex).nGetX();
   int const nYHitProfileEndPoint = pVHitProfileCells->at(nHitProfileEndpointIndex).nGetY();

   // New endpoints of the two profiles (external CRS)
   double const dXThisProfileEndPoint = dGridCentroidXToExtCRSX(nXThisProfileEndPoint);
   double const dYThisProfileEndPoint = dGridCentroidYToExtCRSY(nYThisProfileEndPoint);
   double const dXHitProfileEndPoint = dGridCentroidXToExtCRSX(nXHitProfileEndPoint);
   double const dYHitProfileEndPoint = dGridCentroidYToExtCRSY(nYHitProfileEndPoint);

   // LogStream << m_ulIter << ": new end point of coast " << nThisProfileCoast << " profile " << nThisProfile << " is at [" << nXThisProfileEndPoint << "][" << nYThisProfileEndPoint << "] = {" << dXThisProfileEndPoint << ", " << dYThisProfileEndPoint << "}, new end point of coast " << nHitProfileCoast << " profile " << nHitProfile << " is at [" << nXHitProfileEndPoint << "][" << nYHitProfileEndPoint << "] = {" << dXHitProfileEndPoint << ", " << dYHitProfileEndPoint << "}" << endl;

   // Next unmark the cells that will no longer be 'under' this profile
   for (int nn = nThisProfileEndpointIndex+1; nn < nThisProfileLen; nn++)
   {
      int const nXTmp = pVThisProfileCells->at(nn).nGetX();
      int const nYTmp = pVThisProfileCells->at(nn).nGetY();

      if ((m_pRasterGrid->m_Cell[nXTmp][nYTmp].nGetProfileID() == nThisProfile) && (m_pRasterGrid->m_Cell[nXTmp][nYTmp].nGetProfileCoastID() == nThisProfileCoast))
         m_pRasterGrid->m_Cell[nXTmp][nYTmp].SetCoastAndProfileID(INT_NODATA, INT_NODATA);

      // LogStream << "For coast " << nThisProfileCoast << " profile " << nThisProfile << " unmarking [" << nXTmp << "][" << nYTmp << "]" << endl;
   }

   // And unmark the cells that will no longer be 'under' the hit profile
   for (int nn = nHitProfileEndpointIndex+1; nn < nHitProfileLen; nn++)
   {
      int const nXTmp = pVHitProfileCells->at(nn).nGetX();
      int const nYTmp = pVHitProfileCells->at(nn).nGetY();

      if ((m_pRasterGrid->m_Cell[nXTmp][nYTmp].nGetProfileID() == nHitProfile) && (m_pRasterGrid->m_Cell[nXTmp][nYTmp].nGetProfileCoastID() == nHitProfileCoast))
         m_pRasterGrid->m_Cell[nXTmp][nYTmp].SetCoastAndProfileID(INT_NODATA, INT_NODATA);

      // LogStream << "For coast " << nHitProfileCoast << " profile " << nHitProfile << " unmarking [" << nXTmp << "][" << nYTmp << "]" << endl;
   }

   // Truncate this profile's CGeomMultiLine (external CRS)
   int nRtn = nTruncateProfileMultiLineDifferentCoasts(pThisProfile, dXThisProfileEndPoint, dYThisProfileEndPoint);
   if (nRtn != RTN_OK)
      return nRtn;

   // Truncate the hit profile's CGeomMultiLine (external CRS)
   nRtn = nTruncateProfileMultiLineDifferentCoasts(pHitProfile, dXHitProfileEndPoint, dYHitProfileEndPoint);
   if (nRtn != RTN_OK)
      return nRtn;

   // Truncate the list of cells in this profile
   pVThisProfileCells->resize(nThisProfileEndpointIndex);
   pThisProfile->SetCellsInProfile(pVThisProfileCells);

   // Truncate the list of cells in the hit profile
   pVHitProfileCells->resize(nHitProfileEndpointIndex);
   pHitProfile->SetCellsInProfile(pVHitProfileCells);

   // Flag this profile as truncated
   pThisProfile->SetTruncatedDifferentCoast(true);

   // Flag the hit profile as truncated
   pHitProfile->SetTruncatedDifferentCoast(true);

   if (m_nLogFileDetail >= LOG_FILE_ALL)
   {
      string strTmp;
      if (pThisProfile->bIsGridEdge())
      {
         strTmp += " grid-edge ";

         if (pThisProfile->bStartOfCoast())
            strTmp += "start";

         if (pThisProfile->bEndOfCoast())
            strTmp += "end";
      }

      LogStream << m_ulIter << ": coast " << nThisProfileCoast << strTmp << " profile " << nThisProfile << " hit by coast " << nHitProfileCoast << " profile " << nHitProfile << " at [" << nXIntersect << "][" << nYIntersect << "] = {" << dGridCentroidXToExtCRSX(nXIntersect) << ", " << dGridCentroidYToExtCRSY(nYIntersect) << "}. Profile truncated, new endpoint is [" << nXThisProfileEndPoint << "][" << nYThisProfileEndPoint << "] = {" << dGridCentroidXToExtCRSX(nXThisProfileEndPoint) << ", " << dGridCentroidYToExtCRSY(nYThisProfileEndPoint) << "}. Length was " << nThisProfileLen << " cells, length is now " << nThisProfileEndpointIndex << " cells (" << 100 * nThisProfileEndpointIndex / nThisProfileLen << "%)" << endl;

      if (pHitProfile->bIsGridEdge())
      {
         strTmp = " grid-edge ";

         if (pHitProfile->bStartOfCoast())
            strTmp += "start";

         if (pHitProfile->bEndOfCoast())
            strTmp += "end";
      }

     LogStream << m_ulIter << ": coast " << nHitProfileCoast << strTmp << " profile " << nHitProfile << " hit by coast " << nThisProfileCoast << " profile " << nThisProfile << " at [" << nXIntersect << "][" << nYIntersect << "] = {" << dGridCentroidXToExtCRSX(nXIntersect) << ", " << dGridCentroidYToExtCRSY(nYIntersect) << "}. Profile truncated, new endpoint is [" << nXHitProfileEndPoint << "][" << nYHitProfileEndPoint << "] = {" << dGridCentroidXToExtCRSX(nXHitProfileEndPoint) << ", " << dGridCentroidYToExtCRSY(nYHitProfileEndPoint) << "}. Length was " << nHitProfileLen << " cells, length is now " << nHitProfileEndpointIndex << " cells (" << 100 * nHitProfileEndpointIndex / nHitProfileLen << "%)" << endl;
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Truncates a profile which has hit a different coast
//===============================================================================================================================
int CSimulation::nTruncateProfileHitDifferentCoast(int const nCoast, int const nProfile, int const nXHitCoast, int const nYHitCoast)
{
   // OK, get a pointer to 'this' profile
   CGeomProfile* pProfile = m_VCoast[nCoast].pGetProfile(nProfile);

   // And get the start point of this profile
   CGeom2DIPoint const* pPtiProfileStart = pProfile->pPtiGetFirstCellInProfile();

   // Get the distance between the point where the profile has hit the different coast, and the profile start point
   double const dThisDist = dGetDistanceBetween(pPtiProfileStart->nGetX(), pPtiProfileStart->nGetY(), nXHitCoast, nYHitCoast);

   // Then halve this distance: the result is the new length of the profile
   int const nNewLenProfile = tMax(nRound(dThisDist / 2) - GAP_BETWEEN_DIFFERENT_COAST_PROFILES, MIN_PROFILE_SIZE);

   // Now get pointers to the cells in this profile
   vector<CGeom2DIPoint>* pVProfileCells = pProfile->pPtiVGetCellsInProfile();
   int const nProfileLen = static_cast<int>(pVProfileCells->size());

   // We have the new profile length so next unmark the cells that will no longer be 'under' the profile
   for (int nn = nNewLenProfile-1; nn < nProfileLen; nn++)
   {
      int const nXThis = pVProfileCells->at(nn).nGetX();
      int const nYThis = pVProfileCells->at(nn).nGetY();

      if ((m_pRasterGrid->m_Cell[nXThis][nYThis].nGetProfileID() == nProfile) && (m_pRasterGrid->m_Cell[nXThis][nYThis].nGetProfileCoastID() == nCoast))
         m_pRasterGrid->m_Cell[nXThis][nYThis].SetCoastAndProfileID(INT_NODATA, INT_NODATA);
   }

   // Truncate the list of cells in the profile, and then update
   pVProfileCells->resize(nNewLenProfile);
   pProfile->SetCellsInProfile(pVProfileCells);

   // Get the endpoint of the truncated profile (external CRS)
   double const dXEndPoint = dGridCentroidXToExtCRSX(pVProfileCells->back().nGetX());
   double const dYEndPoint = dGridCentroidYToExtCRSY(pVProfileCells->back().nGetY());

   // Now truncate the profile's CGeomMultiLine (external CRS)
   int const nRtn = nTruncateProfileMultiLineDifferentCoasts(pProfile, dXEndPoint, dYEndPoint);
   if (nRtn != RTN_OK)
      return nRtn;

   // And flag as truncated
   pProfile->SetTruncatedDifferentCoast(true);

   if (m_nLogFileDetail >= LOG_FILE_ALL)
   {
      string strTmp;
      if (pProfile->bStartOfCoast() || pProfile->bEndOfCoast())
      {
         strTmp += " grid-edge ";

         if (pProfile->bStartOfCoast())
            strTmp += "start";

         if (pProfile->bEndOfCoast())
            strTmp += "end";
      }

      LogStream << m_ulIter << ": coast " << nCoast << strTmp << " profile " << nProfile << " hit another coast at [" << nXHitCoast << "][" << nYHitCoast << "] = {" << dGridCentroidXToExtCRSX(nXHitCoast) << ", " << dGridCentroidYToExtCRSY(nYHitCoast) << "}. Profile truncated, length of profile " << nProfile << " was " << nProfileLen << " cells, is now " << nNewLenProfile << " cells (" << 100 * nNewLenProfile / nProfileLen << "%)" << endl;
   }

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
   // if (! bWriteRasterGISFile(RASTER_PLOT_POLYGON, &RASTER_PLOT_POLYGON_TITLE))
   //    return false;
   // // DEBUG CODE ================

   return RTN_OK;
}

//===============================================================================================================================
//! Truncates the CGeomMultiLine (external CRS) of a profile which has hit a different coast, or a coast-normal profile belonging to a different coast
//===============================================================================================================================
int CSimulation::nTruncateProfileMultiLineDifferentCoasts(CGeomProfile* pProfile, double const dX, double const dY)
{
   // Find which multiline line segment 'contains' this ext CRS end point, then truncate at this point
   bool bFound = false;

   int const nProfileLineSegments = pProfile->CGeomMultiLine::nGetNumLineSegments();
   for (int nSeg = 0; nSeg < nProfileLineSegments; nSeg++)
   {
      // Search each segment
      int const nNumCoinc = pProfile->CGeomMultiLine::nGetNumCoincidentProfilesInLineSegment(nSeg);
      for (int nCoinc = 0; nCoinc < nNumCoinc; nCoinc++)
      {
         vector<CGeom2DPoint>& pVPt = pProfile->CGeomMultiLine::pGetPoints();

         for (int nLin = 0; nLin < static_cast<int>(pVPt.size())-1; nLin++)
         {
            double const dX1 = pVPt[nLin].dGetX();
            double const dY1 = pVPt[nLin].dGetY();
            double const dX2 = pVPt[nLin+1].dGetX();
            double const dY2 = pVPt[nLin+1].dGetY();

            double const dXMin = tMin(dX1, dX2);
            double const dXMax = tMax(dX1, dX2);
            double const dYMin = tMin(dY1, dY2);
            double const dYMax = tMax(dY1, dY2);

            // LogStream << m_ulIter << ": coast " << pProfile->nGetCoastID() << " profile " << pProfile->nGetProfileID() << " min is [" << dXMin << ", " << dYMin << "], search point is [" << dX << ", " << dY << "], max is [" << dXMax << ", " << dYMax << "]" << endl;

            if ((dX >= dXMin) && (dX <= dXMax) && (dY >= dYMin) && (dY <= dYMax))
            {
               bFound = true;
               // LogStream << "FOUND" << endl;

               pProfile->TruncateProfile(nSeg+1);
               pProfile->AppendPointInProfile(dX, dY);

               break;
            }
         }
         if (bFound)
             break;
      }
      if (bFound)
         break;
   }

   if (! bFound)
   {
      // LogStream << "NOT FOUND" << endl;
      return RTN_ERR_POINT_NOT_FOUND_IN_MULTILINE_DIFFERENT_COASTS;
   }

   return RTN_OK;
}
