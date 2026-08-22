/*!
   \file calc_shadow_zones.cpp
   \brief Locates shadow zones, is part of wave propagation calculations
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
using std::endl;

#include <stack>
using std::stack;

#include <deque>
using std::deque;

#include <numeric>
using std::accumulate;

#include "cme.h"
#include "coast.h"
#include "simulation.h"
#include "raster_grid.h"
#include "2d_point.h"
#include "2di_point.h"
#include "line.h"
#include "i_line.h"

namespace
{
//===============================================================================================================================
//! Helper function used when sorting unsigned coastline curvature values, to locate start points of normal profiles. If the first argument must be ordered before the second, return true
//===============================================================================================================================
bool bCurvaturePairCompareDescending(const pair<int, double>& prLeft, const pair<int, double>& prRight)
{
   // Sort curvature (low values are straight, high values are curved)
   return prLeft.second > prRight.second;
}
} // namespace

//===============================================================================================================================
//! Determines whether the wave orientation at this point on a coast is on-shore or off-shore, and up-coast (i.e. along the coast in the direction of decreasing coastline point numbers) or down-coast (i.e. along the coast in the direction of increasing coastline point numbers). Note that wave orientation is the oceanographic convention i.e. direction TOWARDS which the waves move (in degrees clockwise from north)
//===============================================================================================================================
bool CSimulation::bOnOrOffShoreAndUpOrDownCoast(double const dCoastAngle, double const dWaveAngle, int const nSeaHand, bool& bDownCoast)
{
   bool bOnShore;
   double const dWaveToCoastAngle = fmod((dWaveAngle - dCoastAngle + 360), 360);

   bDownCoast = ((dWaveToCoastAngle > 270) || (dWaveToCoastAngle < 90)) ? true : false;

   if (nSeaHand == RIGHT_HANDED)
      // The sea is on the RHS when travelling down-coast (i.e. along the coast in the direction of increasing coastline point numbers)
      bOnShore = dWaveToCoastAngle > 180 ? true : false;
   else
      // The sea is on the LHS when travelling down-coast (i.e. along the coast in the direction of increasing coastline point numbers)
      bOnShore = dWaveToCoastAngle > 180 ? false : true;

   return bOnShore;
}

//===============================================================================================================================
//! Given a cell and a wave orientation, finds the 'upwave' cell. Note that wave orientation is the oceanographic convention i.e. direction TOWARDS which the waves move (in degrees clockwise from north)
//===============================================================================================================================
CGeom2DIPoint CSimulation::PtiFollowWaveAngle(CGeom2DIPoint const* pPtiLast, double const dWaveAngleIn, double& dCorrection)
{
   int const nXLast = pPtiLast->nGetX();
   int const nYLast = pPtiLast->nGetY();
   int nXNext = nXLast;
   int nYNext = nYLast;

   double const dWaveAngle = dWaveAngleIn - dCorrection;

   if (dWaveAngle < 22.5)
   {
      nYNext--;
      dCorrection = 22.5 - dWaveAngle;
   }
   else if (dWaveAngle < 67.5)
   {
      nYNext--;
      nXNext++;
      dCorrection = 67.5 - dWaveAngle;
   }
   else if (dWaveAngle < 112.5)
   {
      nXNext++;
      dCorrection = 112.5 - dWaveAngle;
   }
   else if (dWaveAngle < 157.5)
   {
      nXNext++;
      nYNext++;
      dCorrection = 157.5 - dWaveAngle;
   }
   else if (dWaveAngle < 202.5)
   {
      nYNext++;
      dCorrection = 202.5 - dWaveAngle;
   }
   else if (dWaveAngle < 247.5)
   {
      nXNext--;
      nYNext++;
      dCorrection = 247.5 - dWaveAngle;
   }
   else if (dWaveAngle < 292.5)
   {
      nXNext--;
      dCorrection = 292.5 - dWaveAngle;
   }
   else if (dWaveAngle < 337.5)
   {
      nXNext--;
      nYNext--;
      dCorrection = 337.5 - dWaveAngle;
   }
   else
   {
      nYNext--;
      dCorrection = 22.5 - dWaveAngle;
   }

   dCorrection = dKeepWithin360(dCorrection);

   return CGeom2DIPoint(nXNext, nYNext);
}

//===============================================================================================================================
//! Finds wave shadow zones and modifies waves in and near them. Note that where up-coast and down-coast shadow zones overlap, the effects on wave values in the overlap area is an additive decrease in wave energy. Changes to wave energy in any down-drift increased-energy zones are also additive.
//===============================================================================================================================
int CSimulation::nDoAllShadowZones(void)
{
   if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
      LogStream << endl << m_ulIter << ": Finding shadow zones" << endl;

   int nZone = 0;

   // Do this once for each coastline
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      int const nSeaHand = m_VCoast[nCoast].nGetSeaHandedness();
      int const nCoastSize = m_VCoast[nCoast].nGetCoastlineSize();

      // Now create a vector of pairs: the first value of the pair is the coastline point, the second is the coastline's curvature at that point
      vector<pair<int, double>> prVCurvature;

      for (int nCoastPoint = 0; nCoastPoint < nCoastSize; nCoastPoint++)
      {
         double dCurvature;

         int const nCat = m_VCoast[nCoast].pGetCoastLandform(nCoastPoint)->nGetLandFormCategory();
         if ((nCat == LF_INTERVENTION_STRUCT) || (nCat == LF_INTERVENTION_NON_STRUCT))
         {
            // This is an intervention coast point, which is likely to have some sharp angles. So store the detailed curvature
            dCurvature = m_VCoast[nCoast].dGetDetailedCurvature(nCoastPoint);
         }
         else
         {
            // Not an intervention coast point, so store the smoothed curvature
            dCurvature = m_VCoast[nCoast].dGetSmoothCurvature(nCoastPoint);
         }

         // Store the coast point and curvature
         prVCurvature.push_back(make_pair(nCoastPoint, dCurvature));
      }

      // Sort this pair vector in descending order, so that the most convex points are first
      sort(prVCurvature.begin(), prVCurvature.end(), bCurvaturePairCompareDescending);

      // // DEBUG CODE =======================================================================================================================
      // for (int n = 0; n < static_cast<int>(prVCurvature.size()); n++)
      // {
      //    if (prVCurvature[n].second > 0)
      //    {
      //       CGeom2DPoint PtTmp = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(prVCurvature[n].first);
      //       LogStream << prVCurvature[n].first << "\t" << "{" << PtTmp.dGetX() << ", " << PtTmp.dGetY() << "}" << "\t" << prVCurvature[n].second << endl;
      //    }
      // }
      // LogStream << endl;
      // // DEBUG CODE =======================================================================================================================

      // Now process each coastpoint, starting with the most convex
      for (int nSortedPoint = 0; nSortedPoint < static_cast<int>(prVCurvature.size()); nSortedPoint++)
      {
         double const dCurvature = prVCurvature[nSortedPoint].second ;

         // Quit when we get to non-convex points
         if (dCurvature <= 0)
            break;

         int const nCoastPoint = prVCurvature[nSortedPoint].first;

         // // DEBUG CODE ======================
         // CGeom2DIPoint PtiTmp2 = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint);
         // CGeom2DPoint PtTmp2 = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nCoastPoint);
         // LogStream << "Processing nSortedPoint = " << nSortedPoint << " nCoastPoint = " << nCoastPoint << " at [" << PtiTmp2.nGetX() << "][" << PtiTmp2.nGetY() << "] = {" << PtTmp2.dGetX() << ", " << PtTmp2.dGetY() << "} dCurvature = " << dCurvature << endl;
         // // DEBUG CODE ======================

         // OK, the coast is convex (+ve) here, now get the flux orientation (a tangent to the coastline)
         double const dFluxOrientation = m_VCoast[nCoast].dGetFluxOrientation(nCoastPoint);

         // If this coast point is in the active zone, use the breaking wave orientation, otherwise use the deep water wave orientation
         double dWaveAngle;

         double dDepthOfBreaking = m_VCoast[nCoast].dGetDepthOfBreaking(nCoastPoint);

         if (bFPIsEqual(dDepthOfBreaking, DBL_NODATA, TOLERANCE))
            // Not in active zone
            dWaveAngle = m_VCoast[nCoast].dGetCoastDeepWaterWaveAngle(nCoastPoint);
         else
            // In active zone
            dWaveAngle = m_VCoast[nCoast].dGetBreakingWaveAngle(nCoastPoint);

         // At this point on the coast, are waves on- or off-shore, and up- or down-coast?
         bool bDownCoast = false;
         bool const bOnShore = bOnOrOffShoreAndUpOrDownCoast(dFluxOrientation, dWaveAngle, nSeaHand, bDownCoast);
         m_VCoast[nCoast].SetWavesOnShore(nCoastPoint, bOnShore);
         m_VCoast[nCoast].SetWavesDownCoast(nCoastPoint, bDownCoast);

         // if (m_nLogFileDetail >= LOG_FILE_ALL)
         // {
         //    CGeom2DIPoint PtiTmp1 = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint);
         //    CGeom2DPoint PtTmp1 = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nCoastPoint);
         //    LogStream << m_ulIter << ": coast " << nCoast << " coast point " << nCoastPoint << " at [" << PtiTmp1.nGetX() << "][" << PtiTmp1.nGetY() << "] = {" << PtTmp1.dGetX() << ", " << PtTmp1.dGetY() << "} has " << (bDownCoast ? "down-coast " : "up-coast ") << (bOnShore ? "on-shore" : "off-shore") << " waves, dWaveAngle = " << dWaveAngle << " dFluxOrientation = " << dFluxOrientation << endl;
         // }

         // Find the "previous" coast point i.e. the coast point that the waves get to before they get to this point
         int nPrevCoastPoint;
         if (bDownCoast)
            nPrevCoastPoint = nCoastPoint - 1;
         else
            nPrevCoastPoint = nCoastPoint + 1;

         // Check that previous coast point is not beyond the coast endpoints
         if ((nPrevCoastPoint < 0)  || (nPrevCoastPoint > (nCoastSize-1)))
            // It is, so give up for this coast point
            continue;

         // Get the flux orientation (a tangent to the coastline) of the previous coast point
         double const dPrevFluxOrientation = m_VCoast[nCoast].dGetFluxOrientation(nPrevCoastPoint);

         // If the previous coast point is in the active zone, use the breaking wave orientation, otherwise use the deep water wave orientation
         double dPrevWaveAngle;
         double const dPrevDepthOfBreaking = m_VCoast[nCoast].dGetDepthOfBreaking(nPrevCoastPoint);

         if (bFPIsEqual(dPrevDepthOfBreaking, DBL_NODATA, TOLERANCE))
            // Not in active zone
            dPrevWaveAngle = m_VCoast[nCoast].dGetCoastDeepWaterWaveAngle(nPrevCoastPoint);
         else
            // In active zone
            dPrevWaveAngle = m_VCoast[nCoast].dGetBreakingWaveAngle(nPrevCoastPoint);

         // At the previous coast point, are waves on- or off-shore, and up- or down-coast?
         bool bPrevDownCoast = false;
         bool const bPrevOnShore = bOnOrOffShoreAndUpOrDownCoast(dPrevFluxOrientation, dPrevWaveAngle, nSeaHand, bPrevDownCoast);

         // if (m_nLogFileDetail >= LOG_FILE_ALL)
         // {
         //    CGeom2DIPoint PtiTmp1 = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nPrevCoastPoint);
         //    CGeom2DPoint PtTmp1 = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nPrevCoastPoint);
         //    LogStream << m_ulIter << ": coast " << nCoast << " previous coast point " << nPrevCoastPoint << " at [" << PtiTmp1.nGetX() << "][" << PtiTmp1.nGetY() << "] = {" << PtTmp1.dGetX() << ", " << PtTmp1.dGetY() << "} has " << (bPrevDownCoast ? "down-coast " : "up-coast ") << (bPrevOnShore ? "on-shore" : "off-shore") << " waves, dPrevWaveAngle = " << dPrevWaveAngle << " dPrevFluxOrientation = " << dPrevFluxOrientation << endl;
         // }

         // // DEBUG CODE ====================================
         // LogStream << m_ulIter << ": bPrevDownCoast = " << (bPrevDownCoast ? "true" : "false") << " bDownCoast = " << (bDownCoast ? "true" : "false") << " bPrevOnShore = " << (bPrevOnShore ? "true" : "false") << " bOnShore = " << (bOnShore ? "true" : "false") << endl;
         //
         // if ((bPrevDownCoast == bDownCoast) && (bPrevOnShore && (! bOnShore)))
         //    LogStream << "YES" << endl;
         // else
         //    LogStream << "NO" << endl << endl;
         // // DEBUG CODE ====================================

         // If the waves are in the same direction for this coast point and the previous coast point, and the waves are onshore for the previous point but not onshore for this coast point, then this could be the start of a shadow boundary
         if ((bPrevDownCoast == bDownCoast) && (bPrevOnShore && (! bOnShore)))
         {
            // if (m_nLogFileDetail >= LOG_FILE_ALL)
            // {
            //    CGeom2DIPoint PtiTmp1 = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nPrevCoastPoint);
            //    CGeom2DPoint PtTmp1 = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nCoastPoint);
            //    LogStream << m_ulIter << ": coast " << nCoast << " possible shadow boundary start at coast point " << nCoastPoint << ", dCurvature " << dCurvature << ", at [" << PtiTmp1.nGetX() << "][" << PtiTmp1.nGetY() << "] = {" << PtTmp1.dGetX() << ", " << PtTmp1.dGetY() << "}" << endl;
            // }

            // OK this coast point could be the start of a shadow boundary. But have we already marked it as a shadow zone boundary?
            CGeom2DIPoint PtiStart = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint);
            if (m_pRasterGrid->m_Cell[PtiStart.nGetX()][PtiStart.nGetY()].bIsShadowZoneBoundary())
               // Yes, we have already been here
               continue;

            // OK, the next stage is to trace this boundary
            bool bHitEdge = false;
            bool bHitCoast = false;
            bool bHitSea = false;
            bool bStillInland = false;
            int nShadowBoundaryEndCoastPoint = -1;

            // For the shadow zone boundary
            CGeomILine ILShadowBoundary;
            ILShadowBoundary.Append(&PtiStart);

            // TEST TODO ===================================
            bool bFollowUpwave = false;

            if (bFollowUpwave)
            {
               // Trace the shadow zone boundary by following waves in the up-wave direction. The shadow zone boundary may be curved
               int nRtn = nFindShadowZoneBoundaryUpWave(nCoast, nCoastPoint, nShadowBoundaryEndCoastPoint, bHitEdge, bHitCoast, bHitSea, bStillInland, &PtiStart, dPrevWaveAngle, &ILShadowBoundary);
               if (nRtn != RTN_OK)
                  // Do next coast point
                  continue;
            }
            else
            {
               // Trace the shadow zone boundary as a straight line
               int nRtn = nFindShadowZoneBoundaryLine(nCoast, nCoastPoint, nShadowBoundaryEndCoastPoint, bHitEdge, bHitCoast, bHitSea, bStillInland, &PtiStart, dPrevWaveAngle, &ILShadowBoundary);
               if (nRtn != RTN_OK)
                  // Do next coast point
                  continue;
            }

            // The number of the coast point at which the shadow boundary ends
            int nShadowEndCoastPoint;

            // Is this shadow boundary valid?
            if (bHitEdge)
            {
               if (CREATE_SHADOW_ZONE_IF_HITS_GRID_EDGE)
               {
                  // User choice is to create shadow zones if we hit the grid edge. OK, but is the shadow zone line trivially short?
                  double const dShadowLen = dGetDistanceBetween(&ILShadowBoundary[0], &ILShadowBoundary.Back()) * m_dCellSide;

                  if (dShadowLen < MIN_LENGTH_OF_SHADOW_ZONE_LINE)
                  {
                     // Too short, so forget about it
                     if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
                        LogStream << m_ulIter << ":\t coast " << nCoast << " Possible shadow boundary start point " << nCoastPoint << " too short. Length " << dShadowLen << " m minimum length is " << MIN_LENGTH_OF_SHADOW_ZONE_LINE << " m. Starts at [" << ILShadowBoundary[0].nGetX() << "][" << ILShadowBoundary[0].nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILShadowBoundary[0].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILShadowBoundary[0].nGetY()) << "} hits grid edge at [" << ILShadowBoundary.Back().nGetX() << "][" << ILShadowBoundary.Back().nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILShadowBoundary.Back().nGetX()) << ", " << dGridCentroidYToExtCRSY(ILShadowBoundary.Back().nGetY()) << "}" << endl;

                     // Shadow boundary is too short, so go to the next coast point
                     continue;
                  }

                  // We've found a valid grid-edge shadow zone, but we need a distance (in cells) between the shadow boundary start and the 'virtual' shadow boundary end: this is the off-grid point where the shadow boundary would have intersected the coastline, if the grid were big enough. This is of course unknowable. So as a best guess, we choose the shorter of the two distances between the point where the shadow boundary hits the valid edge of the grid, and the start or end of the coast
                  // int nCoastSize = m_VCoast[nCoast].nGetCoastlineSize();
                  CGeom2DIPoint const PtiCoastStart = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(0);
                  CGeom2DIPoint const PtiCoastEnd = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastSize - 1);

                  if (dGetDistanceBetween(&ILShadowBoundary.Back(), &PtiCoastStart) < dGetDistanceBetween(&ILShadowBoundary.Back(), &PtiCoastEnd))
                     nShadowEndCoastPoint = 0;
                  else
                     nShadowEndCoastPoint = nCoastSize-1;

                  // if (m_nLogFileDetail >= LOG_FILE_ALL)
                  //    LogStream << m_ulIter << ":\t coast " << nCoast << " possible shadow boundary from start point " << nCoastPoint << " defines a valid shadow zone. Start point [" << ILShadowBoundary[0].nGetX() << "][" << ILShadowBoundary[0].nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILShadowBoundary[0].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILShadowBoundary[0].nGetY()) << "}, hit grid edge at [" << ILShadowBoundary.Back().nGetX() << "][" << ILShadowBoundary.Back().nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILShadowBoundary.Back().nGetX()) << ", " << dGridCentroidYToExtCRSY(ILShadowBoundary.Back().nGetY()) << "}" << endl;
               }
               else
               {
                  // User choice is to not create shadow zones if we hit the grid edge
                  if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
                     LogStream << m_ulIter << ":\t coast " << nCoast << " possible shadow boundary from start point " << nCoastPoint << " hits a grid edge: ignored. Starts at [" << ILShadowBoundary[0].nGetX() << "][" << ILShadowBoundary[0].nGetY() << "]" << endl;

                  // Go to next coast point
                  continue;
               }
            }
            else
            {
               // The shodow boundary does not hit a grid edge
               CGeom2DIPoint PtiEndTmp = ILShadowBoundary.Back();
               nShadowEndCoastPoint = m_VCoast[nCoast].nGetCoastPointGivenCell(&PtiEndTmp);

               // Safety check
               if (nShadowEndCoastPoint == INT_NODATA)
               {
                  // Couldn't find the shadow boundary endpoint in the coast
                  LogStream << ERR << "could not find shadow zone end point [" << PtiEndTmp.nGetX() << "][" << PtiEndTmp.nGetY() << "] in coast points";

                  continue;
               }

               // Did nGetCoastPointGivenCell() find an adjacent cell to the shadow zone endpoint, rather than the endpoint?
               if (PtiEndTmp != ILShadowBoundary.Back())
                  ILShadowBoundary.Append(&PtiEndTmp);

               if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
                  LogStream << m_ulIter << ":\t coast " << nCoast << " coast " << nCoast << " shadow boundary start point " << nCoastPoint << " is valid shadow zone. Start [" << ILShadowBoundary[0].nGetX() << "][" << ILShadowBoundary[0].nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILShadowBoundary[0].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILShadowBoundary[0].nGetY()) << "} hits coast at [" << ILShadowBoundary.Back().nGetX() << "][" << ILShadowBoundary.Back().nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILShadowBoundary.Back().nGetX()) << ", " << dGridCentroidYToExtCRSY(ILShadowBoundary.Back().nGetY()) << "} coast point " << nShadowBoundaryEndCoastPoint << endl;

            }

            // OK, this is a valid shadow zone
            if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
               LogStream << m_ulIter << ": coast " << nCoast << " valid shadow boundary from coast point " << nCoastPoint << " at [" << ILShadowBoundary[0].nGetX() << "][" << ILShadowBoundary[0].nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILShadowBoundary[0].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILShadowBoundary[0].nGetY()) << "}" << endl;

            nZone++;
            int nShadowLineLen = ILShadowBoundary.nGetSize();
            CGeom2DIPoint const PtiEnd = ILShadowBoundary[nShadowLineLen - 1];

            // The vector shadow boundary (external CRS)
            CGeomLine LShadowBoundary;

            for (int nn = 0; nn < nShadowLineLen; nn++)
            {
               int const nTmpX = ILShadowBoundary[nn].nGetX();
               int const nTmpY = ILShadowBoundary[nn].nGetY();

               // Mark the cells as shadow zone boundary
               m_pRasterGrid->m_Cell[nTmpX][nTmpY].SetShadowZoneBoundary();

               // If this is a sea cell, mark the shadow zone boundary cell as being in the shadow zone, but not yet processed (a -ve number)
               if (m_pRasterGrid->m_Cell[nTmpX][nTmpY].bIsInContiguousSea())
                  m_pRasterGrid->m_Cell[nTmpX][nTmpY].SetShadowZoneNumber(-(nZone + 1));

               // If not already there, append this value to the two shadow boundary vectors
               LShadowBoundary.AppendIfNotPrevious(dGridCentroidXToExtCRSX(nTmpX), dGridCentroidYToExtCRSY(nTmpY));
               ILShadowBoundary.AppendIfNotPrevious(nTmpX, nTmpY);

               if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
                  LogStream << m_ulIter << ":\t coast " << nCoast << " shadow zone " << nZone << ", which starts at [" << ILShadowBoundary[0].nGetX() << "][" << ILShadowBoundary[0].nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILShadowBoundary[0].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILShadowBoundary[0].nGetY()) << "} has cell [" << nTmpX << "][" << nTmpY << "] = {" << dGridCentroidXToExtCRSX(nTmpX) << ", " << dGridCentroidYToExtCRSY(nTmpY) << "}marked as shadow zone boundary" << endl;
            }

            // Put the ext CRS vector shadow boundary into reverse sequence (i.e. start point is last)
            LShadowBoundary.Reverse();

            // Store the reversed ext-CRS shadow zone boundary
            m_VCoast[nCoast].AppendShadowBoundary(&LShadowBoundary);

            // Next, store the coastline part of the whole shadow zone in ILShadowBoundary
            int nEndCoastPoint = m_VCoast[nCoast].nGetCoastPointGivenCell(&ILShadowBoundary[nShadowLineLen - 1]);

            // Safety check
            if (nEndCoastPoint == INT_NODATA)
            {
               // Couldn't find the shadow boundary endpoint in the coast
               LogStream << ERR << "could not find shadow zone end point [" << ILShadowBoundary[nShadowLineLen - 1].nGetX() << "][" << ILShadowBoundary[nShadowLineLen - 1].nGetY() << "] in coast points";

               continue;
            }

            if (nEndCoastPoint > nCoastPoint)
            {
               for (int n = nEndCoastPoint; n > nCoastPoint; n--)
               {
                  CGeom2DIPoint PtiTmp = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(n);
                  ILShadowBoundary.Append(&PtiTmp);

                  // Mark the cell as shadow zone boundary
                  m_pRasterGrid->m_Cell[PtiTmp.nGetX()][PtiTmp.nGetY()].SetShadowZoneBoundary();
               }
            }
            else
            {
               for (int n = nEndCoastPoint; n < nCoastPoint; n++)
               {
                  CGeom2DIPoint PtiTmp = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(n);
                  ILShadowBoundary.Append(&PtiTmp);

                  // Mark the cell as shadow zone boundary
                  m_pRasterGrid->m_Cell[PtiTmp.nGetX()][PtiTmp.nGetY()].SetShadowZoneBoundary();
               }
            }

            // DEBUG CODE =======================================================================================================================
            LogStream << endl;
            for (int k = 0; k < ILShadowBoundary.nGetSize(); k++)
            {
               CGeom2DIPoint PtiTmp = *ILShadowBoundary.pPtiGetAt(k);
               LogStream << k << " [" << PtiTmp.nGetX() << "][" << PtiTmp.nGetY() << "]" << endl;
            }
            LogStream << endl;
            // DEBUG CODE =======================================================================================================================

            // Calculate the centroid
            CGeom2DIPoint const PtiCentroid = PtiPolygonCentroid(ILShadowBoundary.pPtiVGetPoints());

            // Safety check
            if (! bIsWithinValidGrid(&PtiCentroid))
            {
               LogStream << ERR << "coast " << nCoast << " start point for shadow zone cell-by-cell fill [" << PtiCentroid.nGetX() << "][" << PtiCentroid.nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiCentroid.nGetX()) << ", " << dGridCentroidYToExtCRSY(PtiCentroid.nGetY()) << "} is outside grid" << endl;

               continue;
            }

            if (m_nLogFileDetail >= LOG_FILE_ALL)
               LogStream << m_ulIter << ":\t coast " << nCoast << " start point for shadow zone cell-by-cell fill is [" << PtiCentroid.nGetX() << "][" << PtiCentroid.nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiCentroid.nGetX()) << ", " << dGridCentroidYToExtCRSY(PtiCentroid.nGetY()) << "}" << endl;

            // Cell-by-cell fill the shadow zone
            int const nRet = nCellByCellFillShadowZone(nCoast, nZone, &PtiCentroid, &PtiStart, &PtiEnd);
            if (nRet != RTN_OK)
            {
               // Could not find start point for cell-by-cell fill. How serious we judge this to be depends on the length of the shadow zone line
               if (nShadowLineLen < MAX_LEN_SHADOW_LINE_TO_IGNORE)
               {
                  if (m_nLogFileDetail >= LOG_FILE_ALL)
                     LogStream << m_ulIter << ":\t" << WARN << "could not find start point for cell-by-cell fill of shadow zone " << nZone << " but continuing because this is a small shadow zone (shadow line length = " << nShadowLineLen << " cells)" << endl;

                  continue;
               }
               else
               {
                  LogStream << m_ulIter << ":\t " << ERR << "could not find start point for cell-by-cell fill of shadow zone " << nZone << " (shadow line length = " << nShadowLineLen << " cells)" << endl;

                  // TODO Improve this
                  continue;
               }
            }

            // Finally sweep the shadow zone, changing wave orientation and height
            ModifyWavesInShadowZoneAndDownDriftZone(nCoast, nZone, nCoastPoint, nEndCoastPoint);

         }
      }
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Does a cell-by-cell fill of a shadow zone, starting from the centroid
//===============================================================================================================================
int CSimulation::nCellByCellFillShadowZone(int const nCoast, int const nZone, CGeom2DIPoint const* pPtiCentroid, CGeom2DIPoint const* pPtiShadowBoundaryStart, CGeom2DIPoint const* pPtiShadowBoundaryEnd)
{
   // Is the centroid a sea cell?
   bool bStartPointOK = true;
   bool bAllPointNotSea = true;
   CGeom2DIPoint PtiFloodFillStart = *pPtiCentroid;

   if (! m_pRasterGrid->m_Cell[PtiFloodFillStart.nGetX()][PtiFloodFillStart.nGetY()].bIsInContiguousSea())
   {
      // No it isn't: so try to find a cell that is
      bStartPointOK = false;
      double dWeight = 0.05;

      while ((! bStartPointOK) && (dWeight < 1))
      {
         // Find a start point for the Cell-by-cell fill. Because shadow zones are generally triangular, start by choosing a low weighting so that the start point is close to the centroid, but a bit towards the coast. If this doesn't work, go further coastwards
         PtiFloodFillStart = PtiWeightedAverage(pPtiShadowBoundaryEnd, pPtiCentroid, dWeight);

         // Safety check
         if (PtiFloodFillStart == *pPtiCentroid)
         {
            dWeight += 0.05;
            continue;
         }

         // Safety check
         if (! bIsWithinValidGrid(&PtiFloodFillStart))
         {
            if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
               LogStream << m_ulIter << ": " << ERR << "start point [" << PtiFloodFillStart.nGetX() << "][" << PtiFloodFillStart.nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiFloodFillStart.nGetX()) << ", " << dGridCentroidYToExtCRSY(PtiFloodFillStart.nGetY()) << "} for cell-by-cell fill of shadow zone is outside grid" << endl;

            return RTN_ERR_SHADOW_ZONE_FLOOD_FILL_NOGRID;
         }

         if (m_pRasterGrid->m_Cell[PtiFloodFillStart.nGetX()][PtiFloodFillStart.nGetY()].bIsInContiguousSea())
         {
            // Start point is a sea cell, all OK
            bStartPointOK = true;
            bAllPointNotSea = false;
         }
         else
         {
            // Start point is not a sea cell
            if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
               LogStream << m_ulIter << ":\t coast " << nCoast << " shadow zone cell-by-cell fill start point [" << PtiFloodFillStart.nGetX() << "][" << PtiFloodFillStart.nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiFloodFillStart.nGetX()) << ", " << dGridCentroidYToExtCRSY(PtiFloodFillStart.nGetY()) << "} is NOT a sea cell for shadow boundary from cape point [" << pPtiShadowBoundaryStart->nGetX() << "][" << pPtiShadowBoundaryStart->nGetY() << "] = {" << dGridCentroidXToExtCRSX(pPtiShadowBoundaryStart->nGetX()) << ", " << dGridCentroidYToExtCRSY(pPtiShadowBoundaryStart->nGetY()) << "} to [" << pPtiShadowBoundaryEnd->nGetX() << "][" << pPtiShadowBoundaryEnd->nGetY() << "] = {" << dGridCentroidXToExtCRSX(pPtiShadowBoundaryEnd->nGetX()) << ", " << dGridCentroidYToExtCRSY(pPtiShadowBoundaryEnd->nGetY()) << "}, dWeight = " << dWeight << endl;

            dWeight += 0.05;
         }
      }
   }

   if ((! bStartPointOK) && (! bAllPointNotSea))
   {
      if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
         LogStream << m_ulIter << ":\t coast " << nCoast << " " << ERR << "could not find shadow zone cell-by-cell fill start point" << endl;

      return RTN_ERR_SHADOW_ZONE_FLOOD_START_POINT;
   }

   if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
      LogStream << m_ulIter << ":\t coast " << nCoast << " shadow zone cell-by-cell fill start point [" << PtiFloodFillStart.nGetX() << "][" << PtiFloodFillStart.nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiFloodFillStart.nGetX()) << ", " << dGridCentroidYToExtCRSY(PtiFloodFillStart.nGetY()) << "} OK for shadow boundary from [" << pPtiShadowBoundaryStart->nGetX() << "][" << pPtiShadowBoundaryStart->nGetY() << "] = {" << dGridCentroidXToExtCRSX(pPtiShadowBoundaryStart->nGetX()) << ", " << dGridCentroidYToExtCRSY(pPtiShadowBoundaryStart->nGetY()) << "} to [" << pPtiShadowBoundaryEnd->nGetX() << "][" << pPtiShadowBoundaryEnd->nGetY() << "] = {" << dGridCentroidXToExtCRSX(pPtiShadowBoundaryEnd->nGetX()) << ", " << dGridCentroidYToExtCRSY(pPtiShadowBoundaryEnd->nGetY()) << "}" << endl;

   // All OK, so create an empty stack
   stack<CGeom2DIPoint> PtiStack;

   // We have a cell-by-cell fill start point so push this point onto the stack
   PtiStack.push(PtiFloodFillStart);

   // Then do the cell-by-cell fill: loop until there are no more cell coordinates on the stack
   while (! PtiStack.empty())
   {
      CGeom2DIPoint const Pti = PtiStack.top();
      PtiStack.pop();

      int nX = Pti.nGetX();
      int nY = Pti.nGetY();

      while ((nX >= 0) && m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea() && (! m_pRasterGrid->m_Cell[nX][nY].bIsinThisShadowZone(-nZone - 1)) && (! m_pRasterGrid->m_Cell[nX][nY].bIsShadowZoneBoundary()) && (! m_pRasterGrid->m_Cell[nX][nY].bIsCoastline()))
         nX--;

      nX++;

      bool bSpanAbove = false;
      bool bSpanBelow = false;

      while ((nX < m_nXGridSize) && m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea() && (! m_pRasterGrid->m_Cell[nX][nY].bIsinThisShadowZone(-nZone - 1)) && (! m_pRasterGrid->m_Cell[nX][nY].bIsShadowZoneBoundary()) && (! m_pRasterGrid->m_Cell[nX][nY].bIsCoastline()))
      {
         // Mark the cell as being in the shadow zone but not yet processed (a -ve number, with -1 being zone 1)
         m_pRasterGrid->m_Cell[nX][nY].SetShadowZoneNumber(-nZone - 1);

         // LogStream << m_ulIter << ": [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "} marked as shadow zone" << endl;

         if ((! bSpanAbove) && (nY > 0) && m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea() && (! m_pRasterGrid->m_Cell[nX][nY - 1].bIsinThisShadowZone(-nZone - 1)) && (! m_pRasterGrid->m_Cell[nX][nY - 1].bIsShadowZoneBoundary()) && (! m_pRasterGrid->m_Cell[nX][nY - 1].bIsCoastline()))
         {
            PtiStack.push(CGeom2DIPoint(nX, nY - 1));
            bSpanAbove = true;
         }
         else if (bSpanAbove && (nY > 0) && ((!m_pRasterGrid->m_Cell[nX][nY - 1].bIsInContiguousSea()) || m_pRasterGrid->m_Cell[nX][nY - 1].bIsinThisShadowZone(-nZone - 1) || m_pRasterGrid->m_Cell[nX][nY - 1].bIsShadowZoneBoundary() || m_pRasterGrid->m_Cell[nX][nY - 1].bIsCoastline()))
         {
            bSpanAbove = false;
         }

         if ((! bSpanBelow) && m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea() && (nY < m_nYGridSize - 1) && (! m_pRasterGrid->m_Cell[nX][nY + 1].bIsinThisShadowZone(-nZone - 1)) && (!m_pRasterGrid->m_Cell[nX][nY + 1].bIsShadowZoneBoundary()) && (! m_pRasterGrid->m_Cell[nX][nY + 1].bIsCoastline()))
         {
            PtiStack.push(CGeom2DIPoint(nX, nY + 1));
            bSpanBelow = true;
         }
         else if (bSpanBelow && (nY < m_nYGridSize - 1) && ((! m_pRasterGrid->m_Cell[nX][nY + 1].bIsInContiguousSea()) || m_pRasterGrid->m_Cell[nX][nY + 1].bIsinThisShadowZone(-nZone - 1) || m_pRasterGrid->m_Cell[nX][nY + 1].bIsShadowZoneBoundary() || m_pRasterGrid->m_Cell[nX][nY + 1].bIsCoastline()))
         {
            bSpanBelow = false;
         }

         nX++;
      }
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Traverse the shadow zone, changing wave orientation and height, and the down-drift zone, changing only wave height. Do this by following the coast between the shadow boundary start point and end point, and following the downdrift boundary between the same points. At each step, trace a linking line, then move along this line and change wave properties
//===============================================================================================================================
void CSimulation::ModifyWavesInShadowZoneAndDownDriftZone(int const nCoast, int const nZone, int const nStartPoint, int const nShadowBoundaryEndPoint)
{
   int const nCoastSeaHand = m_VCoast[nCoast].nGetSeaHandedness();
   int nShadowZoneCoastToCapeSeaHand;

   if (nCoastSeaHand == LEFT_HANDED)
      nShadowZoneCoastToCapeSeaHand = RIGHT_HANDED;
   else
      nShadowZoneCoastToCapeSeaHand = LEFT_HANDED;

   // We will traverse the coastline from the start point of the shadow zone line, going toward the end point. Which direction is this?
   bool bSweepDownCoast = true;

   if (nShadowBoundaryEndPoint < nStartPoint)
      bSweepDownCoast = false;

   // Get the distance (in cells) from the shadow boundary start point to the shadow boundary end point, going along the coast
   int const nAlongCoastDistanceToShadowEndpoint = tAbs(nShadowBoundaryEndPoint - nStartPoint - 1);

   // // DEBUG CODE =========================
   // CGeom2DIPoint PtiDownDriftStarttointTmp = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nStartPoint);
   // LogStream << "nStartPoint = " << nStartPoint << " [" << PtiDownDriftStarttointTmp.nGetX() << "][" << PtiDownDriftStarttointTmp.nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiDownDriftStarttointTmp.nGetX()) << ", " << dGridCentroidYToExtCRSY(PtiDownDriftStarttointTmp.nGetY()) << "}" << endl;
   // LogStream << "nAlongCoastDistanceToShadowEndpoint = " << nAlongCoastDistanceToShadowEndpoint << endl;
   // // DEBUG CODE =========================

   // Calculate the point on the coastline which is 2 * nAlongCoastDistanceToShadowEndpoint from the shadow boundary start point, this will be the end point of the downdrift zone. This point may be beyond the end of the coastline in either direction
   int nDownDriftEndPoint;
   int const nTotAlongCoastDistanceToDownDriftEndpoint = 2 * nAlongCoastDistanceToShadowEndpoint;

   if (bSweepDownCoast)
      nDownDriftEndPoint = nStartPoint + nTotAlongCoastDistanceToDownDriftEndpoint;
   else
      nDownDriftEndPoint = nStartPoint - nTotAlongCoastDistanceToDownDriftEndpoint;

   // // DEBUG CODE =========================
   // CGeom2DIPoint PtiDownDriftEndPointTmp = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nDownDriftEndPoint);
   // LogStream << "nShadowBoundaryEndPoint = " << nStartPoint << " [" << PtiDownDriftEndPointTmp.nGetX() << "][" << PtiDownDriftEndPointTmp.nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiDownDriftEndPointTmp.nGetX()) << ", " << dGridCentroidYToExtCRSY(PtiDownDriftEndPointTmp.nGetY()) << "}" << endl;
   // // DEBUG CODE =========================

   // Next find the actual (i.e. within-grid) end of the downdrift line
   CGeom2DIPoint PtiDownDriftEndPoint;

   // Is the downdrift end point beyond the start or end of the coastline?
   if (nDownDriftEndPoint < 0)
   {
      // Is beyond the start of the coastline
      int const nStartEdge = m_VCoast[nCoast].nGetStartEdge();

      if (nStartEdge == NORTH)
      {
         PtiDownDriftEndPoint.SetX(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(0)->nGetX());
         PtiDownDriftEndPoint.SetY(nDownDriftEndPoint);
      }
      else if (nStartEdge == SOUTH)
      {
         PtiDownDriftEndPoint.SetX(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(0)->nGetX());
         PtiDownDriftEndPoint.SetY(m_nYGridSize - nDownDriftEndPoint - 1);
      }
      else if (nStartEdge == WEST)
      {
         PtiDownDriftEndPoint.SetX(nDownDriftEndPoint);
         PtiDownDriftEndPoint.SetY(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(0)->nGetY());
      }
      else if (nStartEdge == EAST)
      {
         PtiDownDriftEndPoint.SetX(m_nXGridSize - nDownDriftEndPoint - 1);
         PtiDownDriftEndPoint.SetY(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(0)->nGetY());
      }
   }
   else if (nDownDriftEndPoint >= m_VCoast[nCoast].nGetCoastlineSize())
   {
      // Is beyond the end of the coastline
      int const nEndEdge = m_VCoast[nCoast].nGetEndEdge();
      int const nCoastSize = m_VCoast[nCoast].nGetCoastlineSize();

      if (nEndEdge == NORTH)
      {
         PtiDownDriftEndPoint.SetX(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastSize - 1)->nGetX());
         PtiDownDriftEndPoint.SetY(-nDownDriftEndPoint);
      }
      else if (nEndEdge == SOUTH)
      {
         PtiDownDriftEndPoint.SetX(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastSize - 1)->nGetX());
         PtiDownDriftEndPoint.SetY(m_nYGridSize + nDownDriftEndPoint);
      }
      else if (nEndEdge == WEST)
      {
         PtiDownDriftEndPoint.SetX(-nDownDriftEndPoint);
         PtiDownDriftEndPoint.SetY(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastSize - 1)->nGetY());
      }
      else if (nEndEdge == EAST)
      {
         PtiDownDriftEndPoint.SetX(m_nXGridSize + nDownDriftEndPoint);
         PtiDownDriftEndPoint.SetY(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastSize - 1)->nGetY());
      }
   }
   else
   {
      // Is on the coastline, so get the location (grid CRS)
      PtiDownDriftEndPoint = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nDownDriftEndPoint);
   }

   // // DEBUG CODE =========================
   // LogStream << "After within-grid check, PtiDownDriftEndPoint = [" << PtiDownDriftEndPoint.nGetX() << "][" << PtiDownDriftEndPoint.nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiDownDriftEndPoint.nGetX()) << ", " << dGridCentroidYToExtCRSY(PtiDownDriftEndPoint.nGetY()) << "}" << endl;
   // // DEBUG CODE =========================

   // Get the location (grid CRS) of the shadow boundary start point: this is also the start point of the downdrift boundary
   CGeom2DIPoint const* pPtiDownDriftBoundaryStartPoint = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nStartPoint);

   // Now trace the down-drift boundary line: interpolate between cells by a simple DDA line algorithm, see http://en.wikipedia.org/wiki/Digital_differential_analyzer_(graphics_algorithm) Note that Bresenham's algorithm gave occasional gaps
   int const nXStart = pPtiDownDriftBoundaryStartPoint->nGetX();
   int const nYStart = pPtiDownDriftBoundaryStartPoint->nGetY();
   int const nXEnd = PtiDownDriftEndPoint.nGetX();
   int const nYEnd = PtiDownDriftEndPoint.nGetY();

   // Safety check
   if ((nXStart == nXEnd) && (nYStart == nYEnd))
      return;

   double const dXStart = dGridCentroidXToExtCRSX(nXStart);
   double const dYStart = dGridCentroidYToExtCRSY(nYStart);
   double const dXEnd = dGridCentroidXToExtCRSX(nXEnd);
   double const dYEnd = dGridCentroidYToExtCRSY(nYEnd);
   double dXInc = dXEnd - dXStart;
   double dYInc = dYEnd - dYStart;
   double const dLength = tMax(tAbs(dXInc), tAbs(dYInc));

   dXInc /= dLength;
   dYInc /= dLength;

   int nTotDownDriftBoundaryDistance = 0;
   double dX = dXStart;
   double dY = dYStart;

   CGeomLine LDownDriftBoundary;

   // Process each interpolated point
   for (int m = 0; m <= nRound(dLength); m++)
   {
      int const nX = nRound(dExtCRSXToGridX(dX));
      int const nY = nRound(dExtCRSYToGridY(dY));

      if (! bIsWithinValidGrid(nX, nY))
      {
         // Safety check
         break;
      }

      // OK, this is part of the downdrift boundary so store this coordinate and mark the cell
      CGeom2DPoint const PtThis(dX, dY);

      // Make sure we have not already stored this coordinate (can happen, due to rounding)
      if ((LDownDriftBoundary.nGetSize() == 0) || (PtThis != LDownDriftBoundary.pPtBack()))
      {
         // Store this coordinate
         LDownDriftBoundary.Append(&PtThis);

         // Mark the cell (a +ve number, same as the associated shadow zone number i.e. starting from 1)
         m_pRasterGrid->m_Cell[nX][nY].SetDownDriftZoneNumber(nZone + 1);

         // Increment the boundary length
         nTotDownDriftBoundaryDistance++;

         // LogStream << "DownDrift boundary [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;

         dX += dXInc;
         dY += dYInc;
      }
   }

   // Store the downdrift boundary (external CRS), with the start point first
   m_VCoast[nCoast].AppendShadowDowndriftBoundary(&LDownDriftBoundary);

   // Compare the lengths of the along-coast and the along-downdrift boundaries. The increment will be 1 for the smaller of the two, will be > 1 for the larger of the two
   int nMaxDistance;
   double dAlongCoastIncrement = 1;
   double dDownDriftBoundaryIncrement = 1;

   if (nTotAlongCoastDistanceToDownDriftEndpoint < nTotDownDriftBoundaryDistance)
   {
      // The downdrift boundary distance is the larger, so change it
      dDownDriftBoundaryIncrement = static_cast<double>(nTotDownDriftBoundaryDistance) / nTotAlongCoastDistanceToDownDriftEndpoint;
      nMaxDistance = nTotDownDriftBoundaryDistance;
   }
   else
   {
      // The along-coast distance is the larger, so change it
      dAlongCoastIncrement = static_cast<double>(nTotAlongCoastDistanceToDownDriftEndpoint) / nTotDownDriftBoundaryDistance;
      nMaxDistance = nTotAlongCoastDistanceToDownDriftEndpoint;
   }

   double dCoastDistSoFar = 0;
   double dDownDriftBoundaryDistSoFar = 0;

   // Now traverse the along-coast line and the down-drift boundary line, but with different increments for each
   for (int n = 1; n < nMaxDistance - 1; n++)
   {
      dCoastDistSoFar += dAlongCoastIncrement;
      dDownDriftBoundaryDistSoFar += dDownDriftBoundaryIncrement;

      // // DEBUG CODE =========================
      // LogStream << "dDownDriftBoundaryDistSoFar = " << dDownDriftBoundaryDistSoFar << " nTotDownDriftBoundaryDistance = " << nTotDownDriftBoundaryDistance << endl;
      // // DEBUG CODE =========================

      if ((dCoastDistSoFar >= nTotAlongCoastDistanceToDownDriftEndpoint) || (dDownDriftBoundaryDistSoFar >= nTotDownDriftBoundaryDistance))
         break;

      bool bPastShadowEnd = false;
      int nAlongCoast;

      if (bSweepDownCoast)
      {
         nAlongCoast = nStartPoint + nRound(dCoastDistSoFar);

         if (nAlongCoast >= m_VCoast[nCoast].nGetCoastlineSize())
            break;

         if (nAlongCoast >= nShadowBoundaryEndPoint)
            bPastShadowEnd = true;
      }
      else
      {
         nAlongCoast = nStartPoint - nRound(dCoastDistSoFar);

         if (nAlongCoast < 0)
            break;

         if (nAlongCoast <= nShadowBoundaryEndPoint)
            bPastShadowEnd = true;
      }

      int const nAlongDownDriftBoundary = nRound(dDownDriftBoundaryDistSoFar);

      // // DEBUG CODE =========================
      // LogStream << m_ulIter << ":\t dCoastDistSoFar = " << dCoastDistSoFar << " (nTotAlongCoastDistanceToDownDriftEndpoint = " << nTotAlongCoastDistanceToDownDriftEndpoint << ") dDownDriftBoundaryDistSoFar = " << dDownDriftBoundaryDistSoFar << " (nTotDownDriftBoundaryDistance = " << nTotDownDriftBoundaryDistance << ") nAlongCoast = " << nAlongCoast << ", nShadowBoundaryEndPoint = " << nShadowBoundaryEndPoint << ",  nAlongDownDriftBoundary = " << nAlongDownDriftBoundary << endl;
      // // DEBUG CODE =========================

      // Get the two endpoints of the linking line
      CGeom2DIPoint const* pPtiCoast = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nAlongCoast);
      int const nCoastX = pPtiCoast->nGetX();
      int const nCoastY = pPtiCoast->nGetY();

      int const nDownDriftX = nRound(dExtCRSXToGridX(LDownDriftBoundary[nAlongDownDriftBoundary].dGetX()));

      // Safety check
      if (nCoastX >= m_nXGridSize)
         continue;

      int const nDownDriftY = nRound(dExtCRSYToGridY(LDownDriftBoundary[nAlongDownDriftBoundary].dGetY()));

      // Safety check
      if (nCoastY >= m_nYGridSize)
         continue;

      // Safety check, in case the two points are identical (can happen due to rounding)
      if ((nCoastX == nDownDriftX) && (nCoastY == nDownDriftY))
      {
         if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
            LogStream << m_ulIter << ":\t coast " << nCoast << " coast point and downdrift boundary point [" << nCoastX << "][" << nCoastY << "] = {" << dGridCentroidXToExtCRSX(nCoastX) << ", " << dGridCentroidYToExtCRSY(nCoastY) << "} are identical, ignoring" << endl;

         continue;
      }

      // Traverse the linking line, interpolating between cells by a simple DDA line algorithm, see http://en.wikipedia.org/wiki/Digital_differential_analyzer_(graphics_algorithm)
      dXInc = nDownDriftX - nCoastX;
      dYInc = nDownDriftY - nCoastY;
      double const dLinkingLineLength = tMax(tAbs(dXInc), tAbs(dYInc));

      dXInc /= dLinkingLineLength;
      dYInc /= dLinkingLineLength;

      dX = nCoastX,
      dY = nCoastY;

      // Process each interpolated point along the linking line
      int nXLast = -1;
      int nYLast = -1;
      int nShadowZoneLength = 0;
      vector<int> VnShadowCellX, VnShadowCellY;

      for (int m = 0; m < dLinkingLineLength; m++)
      {
         int const nX = nRound(dX);
         int const nY = nRound(dY);

         // Check to see if we just processed this point, can happen due to rounding
         if ((nX == nXLast) && (nY == nYLast))
         {
            // // DEBUG CODE =========================
            // LogStream << m_ulIter << ":\t n = " << n << ", m = " << m << ", dLinkingLineLength = " << dLinkingLineLength << ", dCoastDistSoFar = " << dCoastDistSoFar << " (nTotAlongCoastDistanceToDownDriftEndpoint = " << nTotAlongCoastDistanceToDownDriftEndpoint << "), dDownDriftBoundaryDistSoFar = " << dDownDriftBoundaryDistSoFar << " (nTotDownDriftBoundaryDistance = " << nTotDownDriftBoundaryDistance << ") same as last point at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;
            // // DEBUG CODE =========================

            // Set for next time
            nXLast = nX;
            nYLast = nY;
            dX += dXInc;
            dY += dYInc;

            continue;
         }

         // Outside valid grid?
         if (! bIsWithinValidGrid(nX, nY))
         {
            // // DEBUG CODE =========================
            // LogStream << m_ulIter << ":\t n = " << n << ", m = " << m << ", dLinkingLineLength = " << dLinkingLineLength << ", dCoastDistSoFar = " << dCoastDistSoFar << " (nTotAlongCoastDistanceToDownDriftEndpoint = " << nTotAlongCoastDistanceToDownDriftEndpoint << "), dDownDriftBoundaryDistSoFar = " << dDownDriftBoundaryDistSoFar << " (nTotDownDriftBoundaryDistance = " << nTotDownDriftBoundaryDistance << ") outside valid grid at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;
            // // DEBUG CODE =========================

            // Set for next time
            nXLast = nX;
            nYLast = nY;
            dX += dXInc;
            dY += dYInc;

            continue;
         }

         // Not a sea cell?
         if (! m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea())
         {
            // Not a sea cell
            // // DEBUG CODE =========================
            // LogStream << m_ulIter << ":\t n = " << n << ", m = " << m << ", dLinkingLineLength = " << dLinkingLineLength << ", dCoastDistSoFar = " << dCoastDistSoFar << " (nTotAlongCoastDistanceToDownDriftEndpoint = " << nTotAlongCoastDistanceToDownDriftEndpoint << "), dDownDriftBoundaryDistSoFar = " << dDownDriftBoundaryDistSoFar << " (nTotDownDriftBoundaryDistance = " << nTotDownDriftBoundaryDistance << ") not a sea cell at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;
            // // DEBUG CODE =========================

            // Set for next time
            nXLast = nX;
            nYLast = nY;
            dX += dXInc;
            dY += dYInc;

            continue;
         }

         // Have we gone past the point where the shadow boundary meets the coast (i.e. the shadow boundary end point)?
         if (! bPastShadowEnd)
         {
            // We have not, so the linking line has two parts: one between the coast and the shadow boundary, one between the shadow boundary and the downdrift boundary
            bool bInShadowZone = true;

            if (! m_pRasterGrid->m_Cell[nX][nY].bIsinAnyShadowZone())
            {
               // We have left the shadow zone
               // // DEBUG CODE =========================
               // LogStream << "[" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "} LEFT SHADOW ZONE" << endl;
               // // DEBUG CODE =========================

               bInShadowZone = false;

               // Go back over stored cell coords and set their wave properties
               for (unsigned int mm = 0; mm < VnShadowCellX.size(); mm++)
               {
                  // Process this shadow zone cell
                  ProcessShadowZoneCell(VnShadowCellX[mm], VnShadowCellY[mm], nShadowZoneCoastToCapeSeaHand, pPtiCoast, VnShadowCellX.back(), VnShadowCellY.back(), nZone);

                  // Also process adjacent cells
                  if (mm > 0)
                  {
                     CGeom2DIPoint const PtiLeft = PtiGetPerpendicular(VnShadowCellX[mm], VnShadowCellY[mm], VnShadowCellX[mm - 1], VnShadowCellY[mm - 1], 1, RIGHT_HANDED);
                     CGeom2DIPoint const PtiRight = PtiGetPerpendicular(VnShadowCellX[mm], VnShadowCellY[mm], VnShadowCellX[mm - 1], VnShadowCellY[mm - 1], 1, LEFT_HANDED);

                     if ((PtiLeft.nGetX() != INT_NODATA) && bIsWithinValidGrid(&PtiLeft))
                        ProcessShadowZoneCell(PtiLeft.nGetX(), PtiLeft.nGetY(), nShadowZoneCoastToCapeSeaHand, pPtiCoast, VnShadowCellX.back(), VnShadowCellY.back(), nZone);

                     if ((PtiRight.nGetX() != INT_NODATA) && bIsWithinValidGrid(&PtiRight))
                        ProcessShadowZoneCell(PtiRight.nGetX(), PtiRight.nGetY(), nShadowZoneCoastToCapeSeaHand, pPtiCoast, VnShadowCellX.back(), VnShadowCellY.back(), nZone);
                  }
               }
            }

            if (bInShadowZone)
            {
               // Save coords for later
               VnShadowCellX.push_back(nX);
               VnShadowCellY.push_back(nY);

               nShadowZoneLength++;
            }
            else
            {
               // In downdrift zone
            //    // DEBUG CODE =========================
            //    LogStream << m_ulIter << ":\t n = " << n << ", m = " << m << ", dLinkingLineLength = " << dLinkingLineLength << ", dCoastDistSoFar = " << dCoastDistSoFar << " (nTotAlongCoastDistanceToDownDriftEndpoint = " << nTotAlongCoastDistanceToDownDriftEndpoint << "), dDownDriftBoundaryDistSoFar = " << dDownDriftBoundaryDistSoFar << " (nTotDownDriftBoundaryDistance = " << nTotDownDriftBoundaryDistance << ") has [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "} in downdrift zone" << endl;
            // // DEBUG CODE =========================

               // Process this downdrift cell
               ProcessDownDriftCell(nX, nY, (m - nShadowZoneLength), (dLinkingLineLength - nShadowZoneLength), nZone);

               // Also process adjacent cells
               CGeom2DIPoint const PtiLeft = PtiGetPerpendicular(nX, nY, nXLast, nYLast, 1, RIGHT_HANDED);
               CGeom2DIPoint const PtiRight = PtiGetPerpendicular(nX, nY, nXLast, nYLast, 1, LEFT_HANDED);

               if ((PtiLeft.nGetX() != INT_NODATA) && bIsWithinValidGrid(&PtiLeft))
                  ProcessDownDriftCell(PtiLeft.nGetX(), PtiLeft.nGetY(), (m - nShadowZoneLength), (dLinkingLineLength - nShadowZoneLength), nZone);

               if ((PtiRight.nGetX() != INT_NODATA) && bIsWithinValidGrid(&PtiRight))
                  ProcessDownDriftCell(PtiRight.nGetX(), PtiRight.nGetY(), (m - nShadowZoneLength), (dLinkingLineLength - nShadowZoneLength), nZone);
            }
         }
         else
         {
            // We have, so the linking line has only one part: between the coast and the downdrift boundary
            // // DEBUG CODE =========================
            // LogStream << m_ulIter << ":\t n = " << n << ", m = " << m << ", dLinkingLineLength = " << dLinkingLineLength << ", dCoastDistSoFar = " << dCoastDistSoFar << " (nTotAlongCoastDistanceToDownDriftEndpoint = " << nTotAlongCoastDistanceToDownDriftEndpoint << "), dDownDriftBoundaryDistSoFar = " << dDownDriftBoundaryDistSoFar << " (nTotDownDriftBoundaryDistance = " << nTotDownDriftBoundaryDistance << ") has [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "} in downdrift zone" << endl;
            // // DEBUG CODE =========================

            // Process this downdrift cell
            ProcessDownDriftCell(nX, nY, m, dLinkingLineLength, nZone);

            // Also process adjacent cells
            if ((nXLast != -1) && (nYLast != -1))
            {
               CGeom2DIPoint const PtiLeft = PtiGetPerpendicular(nX, nY, nXLast, nYLast, 1, RIGHT_HANDED);
               CGeom2DIPoint const PtiRight = PtiGetPerpendicular(nX, nY, nXLast, nYLast, 1, LEFT_HANDED);

               if ((PtiLeft.nGetX() != INT_NODATA) && bIsWithinValidGrid(&PtiLeft))
                  ProcessDownDriftCell(PtiLeft.nGetX(), PtiLeft.nGetY(), m, dLinkingLineLength, nZone);

               if ((PtiRight.nGetX() != INT_NODATA) && bIsWithinValidGrid(&PtiRight))
                  ProcessDownDriftCell(PtiRight.nGetX(), PtiRight.nGetY(), m, dLinkingLineLength, nZone);
            }
         }

         // Set for next time
         nXLast = nX;
         nYLast = nY;
         dX += dXInc;
         dY += dYInc;
      }
   }
}

//===============================================================================================================================
//! Process a single cell which is in the downdrift zone, changing its wave height
//===============================================================================================================================
void CSimulation::ProcessDownDriftCell(int const nX, int const nY, int const nTraversed, double const dTotalToTraverse, int const nZone)
{
   // Get the pre-existing (i.e. shore-parallel) wave height
   double const dWaveHeight = m_pRasterGrid->m_Cell[nX][nY].dGetWaveHeight();

   if (bFPIsEqual(dWaveHeight, DBL_NODATA, TOLERANCE))
   {
      // Is not a sea cell
      // LogStream << m_ulIter << ": [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "} ignored, not a sea cell" << endl;

      return;
   }

   int nZoneCode = m_pRasterGrid->m_Cell[nX][nY].nGetShadowZoneNumber();

   if (nZoneCode == (nZone + 1))
   {
      // This cell is in the associated shadow zone, so don't change it
      // LogStream << m_ulIter << ": [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "} ignored, is in associated shadow zone (" << nZone+1 << ")" << endl;

      return;
   }

   if (nZoneCode < 0)
   {
      // This cell is in a shadow zone but is not yet processed, so don't change it
      // LogStream << m_ulIter << ": [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "} ignored, is an unprocessed cell in shadow zone " << nZoneCode << endl;

      return;
   }

   nZoneCode = m_pRasterGrid->m_Cell[nX][nY].nGetDownDriftZoneNumber();

   if (nZoneCode == (nZone + 1))
   {
      // We have already processed this cell for this downdrift zone
      // LogStream << m_ulIter << ": [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "} ignored, already done for this down-drift zone " << nZone+1 << endl;

      return;
   }

   // OK, we are downdrift of the shadow zone area and have not yet processed this cell for this zone, so mark it
   m_pRasterGrid->m_Cell[nX][nY].SetDownDriftZoneNumber(nZone + 1);

   // Equation 14 from Hurst et al. TODO 056 Check this! Could not get this to work (typo in paper?), so used the equation below instead
   // double dKp = 0.5 * (1.0 - sin((PI * 90.0 * nSweep) / (180.0 * nSweepLength)));
   double const dKp = 0.5 + (0.5 * sin((PI * nTraversed) / (2.0 * dTotalToTraverse)));

   // Set the modified wave height
   m_pRasterGrid->m_Cell[nX][nY].SetWaveHeight(dKp * dWaveHeight);

   // LogStream << m_ulIter << ": [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}, nTraversed = " << nTraversed << " dTotalToTraverse = " << dTotalToTraverse << " fraction traversed = " << nTraversed / dTotalToTraverse << endl << "m_pRasterGrid->m_Cell[" << nX << "][" << nY << "].dGetCellDeepWaterWaveHeight() = " << m_pRasterGrid->m_Cell[nX][nY].dGetCellDeepWaterWaveHeight() << " m, original dWaveHeight = " << dWaveHeight << " m, dKp = " << dKp << ", modified wave height = " << dKp * dWaveHeight << " m" << endl << endl;
}

//===============================================================================================================================
//! Process a single cell which is in the shadow zone, changing its wave height and orientation
//===============================================================================================================================
void CSimulation::ProcessShadowZoneCell(int const nX, int const nY, int const nShadowZoneCoastToCapeSeaHand, CGeom2DIPoint const* pPtiCoast, int const nShadowEndX, int const nShadowEndY, int const nZone)
{
   int const nZoneCode = m_pRasterGrid->m_Cell[nX][nY].nGetShadowZoneNumber();

   if (nZoneCode == (-nZone - 1))
   {
      // OK, we are in the shadow zone and have not already processed this cell, so mark it (a +ve number, starting from 1)
      m_pRasterGrid->m_Cell[nX][nY].SetShadowZoneNumber(nZone + 1);

      // Next calculate wave angle here: first calculate dOmega, the signed angle subtended between this end point and the start point, and this end point and the end of the shadow boundary
      CGeom2DIPoint const PtiThis(nX, nY);
      CGeom2DIPoint const PtiShadowBoundary(nShadowEndX, nShadowEndY);
      double const dOmega = 180 * dAngleSubtended(pPtiCoast, &PtiThis, &PtiShadowBoundary) / PI;

      // If dOmega is 90 degrees or more in either direction, set both wave angle and wave height to zero
      if (tAbs(dOmega) >= 90)
      {
         m_pRasterGrid->m_Cell[nX][nY].SetWaveAngle(0);
         m_pRasterGrid->m_Cell[nX][nY].SetWaveHeight(0);

         // LogStream << m_ulIter << ": on shadow linking line with coast end [" << pPtiCoast->nGetX() << "][" << pPtiCoast->nGetY() << "] = {" << dGridCentroidXToExtCRSX(pPtiCoast->nGetX()) << ", " << dGridCentroidYToExtCRSY(pPtiCoast->nGetY()) << "} and shadow boundary end [" << PtiShadowBoundary.nGetX() << "][" << PtiShadowBoundary.nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiShadowBoundary.nGetX()) << ", " << dGridCentroidYToExtCRSY(PtiShadowBoundary.nGetY()) << "}, this point [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl << "angle subtended = " << dOmega << " degrees, m_pRasterGrid->m_Cell[" << nX << "][" << nY << "].dGetCellDeepWaterWaveHeight() = " << m_pRasterGrid->m_Cell[nX][nY].dGetCellDeepWaterWaveHeight() << " degrees, wave orientation = 0 degrees, wave height = 0 m" << endl;
      }
      else
      {
         // Adapted from equation 12 in Hurst et al.
         double const dDeltaShadowWaveAngle = 1.5 * dOmega;

         // Get the pre-existing (i.e. shore-parallel) wave orientation
         double const dWaveAngle = m_pRasterGrid->m_Cell[nX][nY].dGetWaveAngle();

         double dShadowWaveAngle;

         if (nShadowZoneCoastToCapeSeaHand == LEFT_HANDED)
            dShadowWaveAngle = dWaveAngle + dDeltaShadowWaveAngle;

         else
            dShadowWaveAngle = dWaveAngle - dDeltaShadowWaveAngle;

         // Set the shadow zone wave orientation
         m_pRasterGrid->m_Cell[nX][nY].SetWaveAngle(dKeepWithin360(dShadowWaveAngle));

         // Now calculate wave height within the shadow zone, use equation 13 from Hurst et al.
         double const dKp = 0.5 * cos(dOmega * PI / 180);

         // Get the pre-existing (i.e. shore-parallel) wave height
         double const dWaveHeight = m_pRasterGrid->m_Cell[nX][nY].dGetWaveHeight();

         // Set the shadow zone wave height
         m_pRasterGrid->m_Cell[nX][nY].SetWaveHeight(dKp * dWaveHeight);

         // LogStream << m_ulIter << ": on shadow linking line with coast end [" << pPtiCoast->nGetX() << "][" << pPtiCoast->nGetY() << "] = {" << dGridCentroidXToExtCRSX(pPtiCoast->nGetX()) << ", " << dGridCentroidYToExtCRSY(pPtiCoast->nGetY()) << "} and shadow boundary end [" << PtiShadowBoundary.nGetX() << "][" << PtiShadowBoundary.nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiShadowBoundary.nGetX()) << ", " << dGridCentroidYToExtCRSY(PtiShadowBoundary.nGetY()) << "}, this point [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}, angle subtended = " << dOmega << " degrees, m_pRasterGrid->m_Cell[" << nX << "][" << nY << "].dGetCellDeepWaterWaveHeight() = " << m_pRasterGrid->m_Cell[nX][nY].dGetCellDeepWaterWaveHeight() << " m, dDeltaShadowWaveAngle = " << dDeltaShadowWaveAngle << " degrees, dWaveAngle = " << dWaveAngle << " degrees, dShadowWaveAngle = " << dShadowWaveAngle << " degrees, dWaveHeight = " << dWaveHeight << " m, dKp = " << dKp << ", shadow zone wave height = " << dKp * dWaveHeight << " m" << endl;
      }
   }
}

//===============================================================================================================================
//! Trace a shadow zone boundary by following waves in the upwave direction (the boundary may be curved)
//===============================================================================================================================
int CSimulation::nFindShadowZoneBoundaryUpWave(int const nCoast, int const nStartPoint, int& nEndPoint, bool& bHitEdge, bool& bHitCoast, bool& bHitSea, bool& bStillInland, CGeom2DIPoint const* pPtiStart, double dPrevWaveAngle, CGeomILine* pILShadowBoundary)
{
   int nDist = 0;
   double dCorrection = 0;
   CGeom2DIPoint PtiPrev = *pPtiStart;
   deque<double> DQdPrevOrientations;

   while ((! bHitEdge) && (! bHitCoast))
   {
      if (nDist > 0)
      {
         int const nXPrev = PtiPrev.nGetX();
         int const nYPrev = PtiPrev.nGetY();

         if (! m_pRasterGrid->m_Cell[nXPrev][nYPrev].bIsInActiveZone())
         {
            // The previous cell was outside the active zone, so use its wave orientation value
            dPrevWaveAngle = m_pRasterGrid->m_Cell[nXPrev][nYPrev].dGetWaveAngle();
         }
         else
         {
            // The previous cell was in the active zone
            if (bHitSea)
            {
               // If this shadow boundary has already hit sea, then we must be getting near a coast: use the average-so-far wave orientation
               double const dAvgOrientationSoFar = accumulate(DQdPrevOrientations.begin(), DQdPrevOrientations.end(), 0.0) / static_cast<double>(DQdPrevOrientations.size());

               dPrevWaveAngle = dAvgOrientationSoFar;
            }
            else
            {
               // This shadow boundary has not already hit sea, just use the wave orientation from the previous cell
               dPrevWaveAngle = m_pRasterGrid->m_Cell[nXPrev][nYPrev].dGetWaveAngle();

               // LogStream << m_ulIter << ": not already hit sea, using previous cell's wave orientation for cell [" << nXPrev << "][" << nYPrev << "] = {" << dGridCentroidXToExtCRSX(nXPrev) << ", " << dGridCentroidYToExtCRSY(nYPrev) << "}" << endl;
            }
         }

         if (bFPIsEqual(dPrevWaveAngle, DBL_NODATA, TOLERANCE))
         {
            // LogStream << m_ulIter << ": dPrevWaveAngle == DBL_NODATA for cell [" << nXPrev << "][" << nYPrev << "] = {" << dGridCentroidXToExtCRSX(nXPrev) << ", " << dGridCentroidYToExtCRSY(nYPrev) << "}" << endl;

            if (! m_pRasterGrid->m_Cell[nXPrev][nYPrev].bIsInContiguousSea())
            {
               // The previous cell was an inland cell, so use the deep water wave orientation
               dPrevWaveAngle = m_pRasterGrid->m_Cell[nXPrev][nYPrev].dGetCellDeepWaterWaveAngle();
            }
            else
            {
               double const dAvgOrientationSoFar = accumulate(DQdPrevOrientations.begin(), DQdPrevOrientations.end(), 0.0) / static_cast<double>(DQdPrevOrientations.size());

               dPrevWaveAngle = dAvgOrientationSoFar;
            }
         }

         if (DQdPrevOrientations.size() == MAX_NUM_PREV_ORIENTATION_VALUES)
            DQdPrevOrientations.pop_front();

         DQdPrevOrientations.push_back(dPrevWaveAngle);
      }

      // Go upwave along the previous cell's wave orientation to find the new boundary cell
      CGeom2DIPoint PtiNew = PtiFollowWaveAngle(&PtiPrev, dPrevWaveAngle, dCorrection);

      // Get the coordinates of 'this' cell
      int const nX = PtiNew.nGetX();
      int const nY = PtiNew.nGetY();

      // Have we hit the edge of the valid part of the grid?
      if ((! bIsWithinValidGrid(&PtiNew)) || (m_pRasterGrid->m_Cell[nX][nY].bIsBoundingBoxEdge()))
      {
         // Yes we have
         bHitEdge = true;

         if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
            LogStream << m_ulIter << ":\t coast " << nCoast << " shadow boundary " << nStartPoint << " hit edge cell at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;

         continue;
      }

      // OK so far. Have we hit a sea cell yet?
      if ((nDist > MAX_LAND_LENGTH_OF_SHADOW_ZONE_LINE) && (! bHitSea))
      {
         if (m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea())
            bHitSea = true;
         else
         {
            if (nDist >= MAX_LAND_LENGTH_OF_SHADOW_ZONE_LINE)
            {
               // If we have travelled MAX_LAND_LENGTH_OF_SHADOW_ZONE_LINE cells without hitting sea, then abandon this shadow boundary
               bStillInland = true;
               break;
            }
         }
      }

      // Store the coordinates of every cell which we cross
      pILShadowBoundary->Append(&PtiNew);

      // LogStream << m_ulIter << ": at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;

      // Having hit sea, have we now hit we hit a coast point? Note that two diagonal(ish) raster lines can cross each other without any intersection, so must also test an adjacent cell for intersection (does not matter which adjacent cell)
      if (bHitSea)
      {
         if (m_pRasterGrid->m_Cell[nX][nY].bIsCoastline() || (bIsWithinValidGrid(nX, nY + 1) && m_pRasterGrid->m_Cell[nX][nY + 1].bIsCoastline()))
         {
            bHitCoast = true;

            if (m_nLogFileDetail >= LOG_FILE_ALL)
               LogStream << m_ulIter << ":\t coast " << nCoast << " possible shadow boundary from start point " << nStartPoint << " hit the coast at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;
         }
      }

      // For next time
      PtiPrev = PtiNew;
      nDist++;
   }

   if (bStillInland)
   {
      // Shadow line is still inland after crossing MAX_LAND_LENGTH_OF_SHADOW_ZONE_LINE calls
      if (m_nLogFileDetail >= LOG_FILE_ALL)
         LogStream << m_ulIter << ":\t coast " << nCoast << " possible shadow boundary from start point " << nStartPoint << " is still inland after crossing " << MAX_LAND_LENGTH_OF_SHADOW_ZONE_LINE << " cells, abandoning. Starts at [" << pILShadowBoundary->at(0).nGetX() << "][" << pILShadowBoundary->at(0).nGetY() << "] = {" << dGridCentroidXToExtCRSX(pILShadowBoundary->at(0).nGetX()) << ", " << dGridCentroidYToExtCRSY(pILShadowBoundary->at(0).nGetY()) << "} abandoned at [" << pILShadowBoundary->Back().nGetX() << "][" << pILShadowBoundary->Back().nGetY() << "] = {" << dGridCentroidXToExtCRSX(pILShadowBoundary->Back().nGetX()) << ", " << dGridCentroidYToExtCRSY(pILShadowBoundary->Back().nGetY()) << "}" << endl;

      pILShadowBoundary->Clear();

      return RTN_ERR_SHADOW_BOUNDARY_NOGOOD;
   }

   if (bHitCoast)
   {
      // The shadow zone boundary has hit a coast, but is the shadow zone line trivially short?
      CGeom2DIPoint PtiEnd = pILShadowBoundary->Back();

      double const dShadowLen = dGetDistanceBetween(pPtiStart, &PtiEnd) * m_dCellSide;

      if (dShadowLen < MIN_LENGTH_OF_SHADOW_ZONE_LINE)
      {
         // Too short, so forget about it
         if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
            LogStream << m_ulIter << ":\t coast " << nCoast << ", possible shadow boundary from start point " << nStartPoint << " is too short. Length " << dShadowLen << " m minimum length " << MIN_LENGTH_OF_SHADOW_ZONE_LINE << " m. Starts at [" << pILShadowBoundary->at(0).nGetX() << "][" << pILShadowBoundary->at(0).nGetY() << "] = {" << dGridCentroidXToExtCRSX(pILShadowBoundary->at(0).nGetX()) << ", " << dGridCentroidYToExtCRSY(pILShadowBoundary->at(0).nGetY()) << "} hits coast at [" << pILShadowBoundary->Back().nGetX() << "][" << pILShadowBoundary->Back().nGetY() << "] = {" << dGridCentroidXToExtCRSX(pILShadowBoundary->Back().nGetX()) << ", " << dGridCentroidYToExtCRSY(pILShadowBoundary->Back().nGetY()) << "}" << endl;

         pILShadowBoundary->Clear();

         return RTN_ERR_SHADOW_BOUNDARY_NOGOOD;
      }

      // We've found a valid shadow zone. Check the last point in the shadow boundary. Note that occasionally this last cell is not 'above' a cell but is above one of its neighbouring cells is: in which case, replace the last point in the shadow boundary with the coordinates of this neighbouring cell
      nEndPoint = m_VCoast[nCoast].nGetCoastPointGivenCell(&PtiEnd);

      if (nEndPoint == INT_NODATA)
      {
         // Could not find a neighbouring cell which is 'under' the coastline
         if (m_nLogFileDetail >= LOG_FILE_ALL)
            LogStream << m_ulIter << ":\t coast " << nCoast << ", no coast point under {" << dGridCentroidXToExtCRSX(pILShadowBoundary->Back().nGetX()) << ", " << dGridCentroidYToExtCRSY(pILShadowBoundary->Back().nGetY()) << "}" << endl;

         // TODO 004 Need to fix this, for the moment just abandon this shadow zone and carry on
         pILShadowBoundary->Clear();

         return RTN_ERR_SHADOW_BOUNDARY_NOGOOD;
         // return RTN_ERR_NO_CELL_UNDER_COASTLINE;
      }
   }

   return RTN_OK;
}


//===============================================================================================================================
//! Trace a shadow zone boundary as a straight line
//===============================================================================================================================
int CSimulation::nFindShadowZoneBoundaryLine(int const nCoast, int const nStartPoint, int& nEndPoint, bool& bHitEdge, bool& bHitCoast, bool& bHitSea, bool& bStillInland, CGeom2DIPoint const* pPtiStart, double dWaveAngle, CGeomILine* pILShadowBoundary)
{
   int const nXStart = pPtiStart->nGetX();
   int const nYStart = pPtiStart->nGetY();

   if (bFPIsEqual(dWaveAngle, DBL_NODATA, TOLERANCE))
   {
      LogStream << m_ulIter << ": dWaveAngle == DBL_NODATA for cell [" << nXStart << "][" << nYStart << "] = {" << dGridCentroidXToExtCRSX(nXStart) << ", " << dGridCentroidYToExtCRSY(nYStart) << "}" << endl;

      return RTN_ERR_SHADOW_BOUNDARY_NOGOOD;
   }

   int nDist = 0;
   int nX = nXStart;
   int nY = nYStart;

   do
   {
      nDist++;

      double dDeltaX = sin(dWaveAngle * (PI/180));
      double dDeltaY = -cos(dWaveAngle * (PI/180));

      nX += nRound(dDeltaX);
      nY += nRound(dDeltaY);

      // LogStream << " [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "} dDeltaX = " << dDeltaX << " dDeltaY = " << dDeltaY << endl;

      // Have we hit the edge of the valid part of the grid?
      if (! bIsWithinValidGrid(nX, nY))
      {
         // Yes we have
         bHitEdge = true;

         if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
            LogStream << m_ulIter << "\t outside valid grid at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "} " << endl;

         if (CREATE_SHADOW_ZONE_IF_HITS_GRID_EDGE)
            // The shadow boundary hits the grid edge but accept it anyway
            break;
         else
         {
            pILShadowBoundary->Clear();

            LogStream << m_ulIter << "\t abandoning shadow boundary which starts at [" << nXStart << "][" << nYStart << "] = {" << dGridCentroidXToExtCRSX(nXStart) << ", " << dGridCentroidYToExtCRSY(nYStart) << "}" << endl;

            return RTN_ERR_SHADOW_BOUNDARY_NOGOOD;
         }
      }

      // // Is this the same as the previous cell? Can get this because Have we been to this cell before?
      // if (pILShadowBoundary->bIsPresent(nX, nY))
      // {
      //    // We have, so we are in a loop. Abandon this shadow line
      //    bInLoop = true;
      //    break;
      // }

      // OK so far. Have we hit a sea cell yet?
      if ((nDist > MAX_LAND_LENGTH_OF_SHADOW_ZONE_LINE) && (! bHitSea))
      {
         if (m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea())
            bHitSea = true;
         else
         {
            if (nDist >= MAX_LAND_LENGTH_OF_SHADOW_ZONE_LINE)
            {
               // If we have travelled MAX_LAND_LENGTH_OF_SHADOW_ZONE_LINE cells without hitting sea, then abandon this shadow boundary
               bStillInland = true;
               break;
            }
         }
      }

      // Store the coordinates of every cell which we cross
      pILShadowBoundary->Append(nX, nY);

      // LogStream << m_ulIter << ": at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;

      // Having hit sea, have we now hit we hit a coast point? Note that two diagonal(ish) raster lines can cross each other without any intersection, so must also test an adjacent cell for intersection (does not matter which adjacent cell)
      if (bHitSea)
      {
         if (m_pRasterGrid->m_Cell[nX][nY].bIsCoastline() || (bIsWithinValidGrid(nX, nY + 1) && m_pRasterGrid->m_Cell[nX][nY + 1].bIsCoastline()))
         {
            bHitCoast = true;

            if (m_nLogFileDetail >= LOG_FILE_ALL)
               LogStream << m_ulIter << ":\t coast " << nCoast << " possible shadow boundary from start point " << nStartPoint << " hit the coast at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;
         }
      }

      if (m_pRasterGrid->m_Cell[nX][nY].bIsCoastline())
      {
         // LogStream << m_ulIter << "\t hit coast at [" << nX << "][" << nY << "]" << endl;

         break;
      }
   } while (true);

   if (bStillInland)
   {
      // Shadow line is still inland after crossing MAX_LAND_LENGTH_OF_SHADOW_ZONE_LINE calls
      if (m_nLogFileDetail >= LOG_FILE_ALL)
         LogStream << m_ulIter << ":\t coast " << nCoast << " possible shadow boundary from start point " << nStartPoint << " is still inland after crossing " << MAX_LAND_LENGTH_OF_SHADOW_ZONE_LINE << " cells, abandoning. Starts at [" << pILShadowBoundary->at(0).nGetX() << "][" << pILShadowBoundary->at(0).nGetY() << "] = {" << dGridCentroidXToExtCRSX(pILShadowBoundary->at(0).nGetX()) << ", " << dGridCentroidYToExtCRSY(pILShadowBoundary->at(0).nGetY()) << "} abandoned at [" << pILShadowBoundary->Back().nGetX() << "][" << pILShadowBoundary->Back().nGetY() << "] = {" << dGridCentroidXToExtCRSX(pILShadowBoundary->Back().nGetX()) << ", " << dGridCentroidYToExtCRSY(pILShadowBoundary->Back().nGetY()) << "}" << endl;

      pILShadowBoundary->Clear();

      return RTN_ERR_SHADOW_BOUNDARY_NOGOOD;
   }

   if (bHitCoast)
   {
      // The shadow zone boundary has hit a coast, but is the shadow zone line trivially short?
      CGeom2DIPoint PtiEnd = pILShadowBoundary->Back();

      double const dShadowLen = dGetDistanceBetween(pPtiStart, &PtiEnd) * m_dCellSide;

      if (dShadowLen < MIN_LENGTH_OF_SHADOW_ZONE_LINE)
      {
         // Too short, so forget about it
         if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
            LogStream << m_ulIter << ":\t coast " << nCoast << ", possible shadow boundary from start point " << nStartPoint << " is too short. Length " << dShadowLen << " m minimum length " << MIN_LENGTH_OF_SHADOW_ZONE_LINE << " m. Starts at [" << pILShadowBoundary->at(0).nGetX() << "][" << pILShadowBoundary->at(0).nGetY() << "] = {" << dGridCentroidXToExtCRSX(pILShadowBoundary->at(0).nGetX()) << ", " << dGridCentroidYToExtCRSY(pILShadowBoundary->at(0).nGetY()) << "} hits coast at [" << pILShadowBoundary->Back().nGetX() << "][" << pILShadowBoundary->Back().nGetY() << "] = {" << dGridCentroidXToExtCRSX(pILShadowBoundary->Back().nGetX()) << ", " << dGridCentroidYToExtCRSY(pILShadowBoundary->Back().nGetY()) << "}" << endl;

         pILShadowBoundary->Clear();

         return RTN_ERR_SHADOW_BOUNDARY_NOGOOD;
      }

      // We've found a valid shadow zone. Check the last point in the shadow boundary. Note that occasionally this last cell is not 'above' a cell but is above one of its neighbouring cells is: in which case, replace the last point in the shadow boundary with the coordinates of this neighbouring cell
      nEndPoint = m_VCoast[nCoast].nGetCoastPointGivenCell(&PtiEnd);

      if (nEndPoint == INT_NODATA)
      {
         // Could not find a neighbouring cell which is 'under' the coastline
         if (m_nLogFileDetail >= LOG_FILE_ALL)
            LogStream << m_ulIter << ":\t coast " << nCoast << ", no coast point under {" << dGridCentroidXToExtCRSX(pILShadowBoundary->Back().nGetX()) << ", " << dGridCentroidYToExtCRSY(pILShadowBoundary->Back().nGetY()) << "}" << endl;

         // TODO 004 Need to fix this, for the moment just abandon this shadow zone and carry on
         pILShadowBoundary->Clear();
         return RTN_ERR_SHADOW_BOUNDARY_NOGOOD;
         // return RTN_ERR_NO_CELL_UNDER_COASTLINE;
      }
   }

   return RTN_OK;
}
