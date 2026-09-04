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

#include <utility>
using std::make_pair;

#include <algorithm>
using std::sort;

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
//! Helper function used when sorting unsigned coastline curvature values, to locate possible start points of wave shadow zones. If the first argument must be ordered before the second, return true
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
//! Given a cell and a wave orientation, finds the next cell to which the wave travels. Note that wave orientation is the oceanographic convention i.e. direction TOWARDS which the waves move (in degrees clockwise from north)
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

      for (int nStartCoastPoint = 0; nStartCoastPoint < nCoastSize; nStartCoastPoint++)
      {
         double dCurvature;

         int const nCat = m_VCoast[nCoast].pGetCoastLandform(nStartCoastPoint)->nGetLandFormCategory();
         if ((nCat == LF_INTERVENTION_STRUCT) || (nCat == LF_INTERVENTION_NON_STRUCT))
         {
            // This is an intervention coast point, which is likely to have some sharp angles. So store the detailed curvature
            dCurvature = m_VCoast[nCoast].dGetDetailedCurvature(nStartCoastPoint);
         }
         else
         {
            // Not an intervention coast point, so store the smoothed curvature
            dCurvature = m_VCoast[nCoast].dGetSmoothCurvature(nStartCoastPoint);
         }

         // Store the coast point and curvature
         prVCurvature.push_back(make_pair(nStartCoastPoint, dCurvature));
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

      // Now process each coast point, starting with the most convex
      for (int nSortedPoint = 0; nSortedPoint < static_cast<int>(prVCurvature.size()); nSortedPoint++)
      {
         double const dCurvature = prVCurvature[nSortedPoint].second ;

         // Quit when we get to non-convex points
         if (dCurvature <= 0)
            break;

         int const nStartCoastPoint = prVCurvature[nSortedPoint].first;

         // // DEBUG CODE ======================
         // CGeom2DIPoint PtiTmp2 = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nStartCoastPoint);
         // CGeom2DPoint PtTmp2 = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nStartCoastPoint);
         // LogStream << "Processing nSortedPoint = " << nSortedPoint << " nStartCoastPoint = " << nStartCoastPoint << " at [" << PtiTmp2.nGetX() << "][" << PtiTmp2.nGetY() << "] = {" << PtTmp2.dGetX() << ", " << PtTmp2.dGetY() << "} dCurvature = " << dCurvature << endl;
         // // DEBUG CODE ======================

         // OK, the coast is convex (+ve) here, now get the flux orientation (a tangent to the coastline)
         double const dFluxOrientation = m_VCoast[nCoast].dGetFluxOrientation(nStartCoastPoint);

         // If this coast point is in the active zone, use the breaking wave orientation, otherwise use the deep water wave orientation
         double dWaveAngle;

         double const dDepthOfBreaking = m_VCoast[nCoast].dGetDepthOfBreaking(nStartCoastPoint);

         if (bFPIsEqual(dDepthOfBreaking, DBL_NODATA, TOLERANCE))
            // Not in active zone
            dWaveAngle = m_VCoast[nCoast].dGetCoastDeepWaterWaveAngle(nStartCoastPoint);
         else
            // In active zone
            dWaveAngle = m_VCoast[nCoast].dGetBreakingWaveAngle(nStartCoastPoint);

         // At this point on the coast, are waves on- or off-shore, and up- or down-coast?
         bool bDownCoast = false;
         bool const bOnShore = bOnOrOffShoreAndUpOrDownCoast(dFluxOrientation, dWaveAngle, nSeaHand, bDownCoast);
         m_VCoast[nCoast].SetWavesOnShore(nStartCoastPoint, bOnShore);
         m_VCoast[nCoast].SetWavesDownCoast(nStartCoastPoint, bDownCoast);

         // if (m_nLogFileDetail >= LOG_FILE_ALL)
         // {
         //    CGeom2DIPoint PtiTmp1 = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nStartCoastPoint);
         //    CGeom2DPoint PtTmp1 = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nStartCoastPoint);
         //    LogStream << m_ulIter << ": coast " << nCoast << " coast point " << nStartCoastPoint << " at [" << PtiTmp1.nGetX() << "][" << PtiTmp1.nGetY() << "] = {" << PtTmp1.dGetX() << ", " << PtTmp1.dGetY() << "} has " << (bDownCoast ? "down-coast " : "up-coast ") << (bOnShore ? "on-shore" : "off-shore") << " waves, dWaveAngle = " << dWaveAngle << " dFluxOrientation = " << dFluxOrientation << endl;
         // }

         // Find the "previous" coast point i.e. the coast point that the waves got to before they get to this point
         int nPrevCoastPoint;
         if (bDownCoast)
            nPrevCoastPoint = nStartCoastPoint - 1;
         else
            nPrevCoastPoint = nStartCoastPoint + 1;

         // Check that previous coast point is not beyond the coast endpoints
         if ((nPrevCoastPoint < 0)  || (nPrevCoastPoint > (nCoastSize-1)))
            // It is, so abandon this coast point
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
            //    CGeom2DPoint PtTmp1 = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nStartCoastPoint);
            //    LogStream << m_ulIter << ": coast " << nCoast << " possible shadow boundary start at coast point " << nStartCoastPoint << ", dCurvature " << dCurvature << ", at [" << PtiTmp1.nGetX() << "][" << PtiTmp1.nGetY() << "] = {" << PtTmp1.dGetX() << ", " << PtTmp1.dGetY() << "}" << endl;
            // }

            // OK this coast point could be the start of a shadow boundary. But have we already marked it as a shadow zone boundary?
            CGeom2DIPoint const PtiStart = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nStartCoastPoint);
            if (m_pRasterGrid->m_Cell[PtiStart.nGetX()][PtiStart.nGetY()].bIsShadowZone())
               // OK, we have already been here
               continue;

            // The next stage is to trace this boundary
            bool bHitEdge = false;
            bool bHitCoast = false;
            bool bHitSea = false;
            bool bStillInland = false;
            int nShadowBoundaryEndCoastPoint = -1;

            // For the shadow zone boundary
            CGeomILine ILShadowBoundary;

            if (m_bShadowFollowWaveDirection)
            {
               // Trace the shadow zone boundary by following waves. The shadow zone boundary may be curved
               int const nRtn = nFindShadowZoneBoundaryFollowWave(nCoast, nStartCoastPoint, nShadowBoundaryEndCoastPoint, bHitEdge, bHitCoast, bHitSea, bStillInland, &PtiStart, dPrevWaveAngle, &ILShadowBoundary);
               if (nRtn != RTN_OK)
                  // Do next coast point
                  continue;
            }
            else
            {
               // Trace the shadow zone boundary as a straight line
               int const nRtn = nFindShadowZoneBoundaryLine(nCoast, nStartCoastPoint, nShadowBoundaryEndCoastPoint, bHitEdge, bHitCoast, bHitSea, bStillInland, &PtiStart, dPrevWaveAngle, &ILShadowBoundary);
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
                  if (ILShadowBoundary.pPtiVGetPoints()->empty())
                  {
                     if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
                        LogStream << m_ulIter << ":\t coast " << nCoast << " Possible shadow boundary start point " << nStartCoastPoint << " has zero points" << endl;

                     continue;
                  }

                  double const dShadowLen = dGetDistanceBetween(&ILShadowBoundary[0], &ILShadowBoundary.Back()) * m_dCellSide;

                  if (dShadowLen < MIN_LENGTH_OF_SHADOW_ZONE_LINE)
                  {
                     // Too short, so forget about it
                     if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
                        LogStream << m_ulIter << ":\t coast " << nCoast << " Possible shadow boundary start point " << nStartCoastPoint << " too short. Length " << dShadowLen << " m minimum length is " << MIN_LENGTH_OF_SHADOW_ZONE_LINE << " m. Starts at [" << ILShadowBoundary[0].nGetX() << "][" << ILShadowBoundary[0].nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILShadowBoundary[0].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILShadowBoundary[0].nGetY()) << "} hits grid edge at [" << ILShadowBoundary.Back().nGetX() << "][" << ILShadowBoundary.Back().nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILShadowBoundary.Back().nGetX()) << ", " << dGridCentroidYToExtCRSY(ILShadowBoundary.Back().nGetY()) << "}" << endl;

                     // Shadow boundary is too short, so go to the next coast point
                     continue;
                  }

                  // We've found a valid grid-edge shadow zone, but we need a distance (in cells) between the shadow boundary start and the 'virtual' shadow boundary end: this is the off-grid point where the shadow boundary would have intersected the coastline, if the grid were big enough. This is of course unknowable. So as a best guess, we choose the shorter of the two distances between the point where the shadow boundary hits the valid edge of the grid, and the start or end of the coast
                  // CGeom2DIPoint const PtiCoastStart = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(0);
                  // CGeom2DIPoint const PtiCoastEnd = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastSize - 1);

                  // if (dGetDistanceBetween(&ILShadowBoundary.Back(), &PtiCoastStart) < dGetDistanceBetween(&ILShadowBoundary.Back(), &PtiCoastEnd))
                  //    nShadowEndCoastPoint = 0;
                  // else
                  //    nShadowEndCoastPoint = nCoastSize-1;

                  // if (m_nLogFileDetail >= LOG_FILE_ALL)
                  //    LogStream << m_ulIter << ":\t coast " << nCoast << " possible shadow boundary from start point " << nStartCoastPoint << " [" << m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nStartCoastPoint)->nGetX() << "][" << m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nStartCoastPoint)->nGetY() <<"][ = {" << dGridCentroidXToExtCRSX(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nStartCoastPoint)->nGetX()) << ", " << dGridCentroidYToExtCRSY(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nStartCoastPoint)->nGetY()) << "} defines a valid shadow zone. Start point [" << ILShadowBoundary[0].nGetX() << "][" << ILShadowBoundary[0].nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILShadowBoundary[0].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILShadowBoundary[0].nGetY()) << "}, hit grid edge at [" << ILShadowBoundary.Back().nGetX() << "][" << ILShadowBoundary.Back().nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILShadowBoundary.Back().nGetX()) << ", " << dGridCentroidYToExtCRSY(ILShadowBoundary.Back().nGetY()) << "}" << endl;
               }
               else
               {
                  // User choice is to not create shadow zones if we hit the grid edge
                  if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
                     LogStream << m_ulIter << ":\t coast " << nCoast << " possible shadow boundary from start point " << nStartCoastPoint << " [" << m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nStartCoastPoint)->nGetX() << "][" << m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nStartCoastPoint)->nGetY() <<"][ = {" << dGridCentroidXToExtCRSX(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nStartCoastPoint)->nGetX()) << ", " << dGridCentroidYToExtCRSY(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nStartCoastPoint)->nGetY()) << "} hits a grid edge: ignored. Starts at [" << ILShadowBoundary[0].nGetX() << "][" << ILShadowBoundary[0].nGetY() << "]" << endl;

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
                  LogStream << m_ulIter << ":\t coast " << nCoast << " coast " << nCoast << " shadow boundary start point " << nStartCoastPoint << " is valid shadow zone. Start [" << ILShadowBoundary[0].nGetX() << "][" << ILShadowBoundary[0].nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILShadowBoundary[0].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILShadowBoundary[0].nGetY()) << "} hits coast at [" << ILShadowBoundary.Back().nGetX() << "][" << ILShadowBoundary.Back().nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILShadowBoundary.Back().nGetX()) << ", " << dGridCentroidYToExtCRSY(ILShadowBoundary.Back().nGetY()) << "} coast point " << nShadowBoundaryEndCoastPoint << endl;

            }

            // OK, this is a valid shadow zone
            // if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
            //    LogStream << m_ulIter << ": coast " << nCoast << " valid shadow boundary from coast point " << nStartCoastPoint << " at [" << ILShadowBoundary[0].nGetX() << "][" << ILShadowBoundary[0].nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILShadowBoundary[0].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILShadowBoundary[0].nGetY()) << "}" << endl;

            // Create a copy of the shadow line only (i.e. without the coast, which we are just about to add), this will be used later to find the downdrift zone centroid
            CGeomILine ILTmp = ILShadowBoundary;

            nZone++;
            int const nShadowLineLen = ILShadowBoundary.nGetSize();
            CGeom2DIPoint const PtiShadowEnd = ILShadowBoundary.Back();

            // The vector shadow boundary (external CRS)
            CGeomLine LShadowBoundary;

            // The shadow zone cells
            vector<CGeom2DIPoint> VPtiShadowCells;

            // Mark the cells under the grid-CRS shadow boundary, also create an ext-CRS boundary
            for (int nn = 0; nn < nShadowLineLen; nn++)
            {
               int const nTmpX = ILShadowBoundary[nn].nGetX();
               int const nTmpY = ILShadowBoundary[nn].nGetY();

               // If this shadow line cell is a sea cell, mark it as wave shadow zone
               if (m_pRasterGrid->m_Cell[nTmpX][nTmpY].bIsInContiguousSea())
               {
                  m_pRasterGrid->m_Cell[nTmpX][nTmpY].SetShadowZoneNumber(nZone);

                  // Save cell, for later wave modification
                  VPtiShadowCells.push_back(CGeom2DIPoint(nTmpX, nTmpY));

                  // LogStream << m_ulIter << ": [" << nTmpX << "][" << nTmpY << "] = {" << dGridCentroidXToExtCRSX(nTmpX) << ", " << dGridCentroidYToExtCRSY(nTmpY) << "} marked as shadow zone " << nZone << endl;
               }

               // If not already there, append this value to the ext-CRS shadow boundary vector
               LShadowBoundary.AppendIfNotPrevious(dGridCentroidXToExtCRSX(nTmpX), dGridCentroidYToExtCRSY(nTmpY));

               // if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
               //    LogStream << m_ulIter << ":\t coast " << nCoast << " shadow zone " << nZone << " has cell [" << nTmpX << "][" << nTmpY << "] = {" << dGridCentroidXToExtCRSX(nTmpX) << ", " << dGridCentroidYToExtCRSY(nTmpY) << "} stored as raster and vector shadow zone boundary" << endl;
            }

            // Put the ext CRS vector shadow boundary into reverse sequence (i.e. start point is last)
            LShadowBoundary.Reverse();

            // Store the reversed ext-CRS shadow zone boundary
            m_VCoast[nCoast].AppendShadowBoundary(&LShadowBoundary);

            // // Next, store the coastline part of the whole shadow zone in ILShadowBoundary
            int const nEndSZCoastPoint = m_VCoast[nCoast].nGetCoastPointGivenCell(&ILShadowBoundary[nShadowLineLen - 1]);

            // Safety check
            if (nEndSZCoastPoint == INT_NODATA)
            {
               // Couldn't find the shadow boundary endpoint in the coast
               LogStream << ERR << "could not find shadow zone end point [" << ILShadowBoundary[nShadowLineLen - 1].nGetX() << "][" << ILShadowBoundary[nShadowLineLen - 1].nGetY() << "] in coast points";

               continue;
            }

            if (nEndSZCoastPoint > nStartCoastPoint)
            {
               for (int n = nEndSZCoastPoint; n > nStartCoastPoint; n--)
               {
                  CGeom2DIPoint const PtiTmp = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(n);
                  ILShadowBoundary.Append(&PtiTmp);

                  // // Mark this coastline cell as wave shadow zone. Set it as a -ve number to indicate that it is in the shadow zone, but wave values have not yet been changed
                  // int nTmpX = PtiTmp.nGetX();
                  // int nTmpY = PtiTmp.nGetY();
                  // // m_pRasterGrid->m_Cell[nTmpX][nTmpY].SetShadowZoneNumber(nZone);
                  //
                  // // Save cell, for later wave modification
                  // VPtiShadowCells.push_back(CGeom2DIPoint(nTmpX, nTmpY));
                  //
                  // LogStream << m_ulIter << ": BBB [" << nTmpX << "][" << nTmpY << "] = {" << dGridCentroidXToExtCRSX(nTmpX) << ", " << dGridCentroidYToExtCRSY(nTmpY) << "} marked as shadow zone " << nZone << endl;
               }
            }
            else
            {
               for (int n = nEndSZCoastPoint; n < nStartCoastPoint; n++)
               {
                  CGeom2DIPoint const PtiTmp = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(n);
                  ILShadowBoundary.Append(&PtiTmp);

                  // // Mark this coastline cell as wave shadow zone. Set it as a -ve number to indicate that it is in the shadow zone, but wave values have not yet been changed
                  // int nTmpX = PtiTmp.nGetX();
                  // int nTmpY = PtiTmp.nGetY();
                  // // m_pRasterGrid->m_Cell[nTmpX][nTmpY].SetShadowZoneNumber(nZone);
                  //
                  // // Save cell, for later wave modification
                  // VPtiShadowCells.push_back(CGeom2DIPoint(nTmpX, nTmpY));
                  //
                  // LogStream << m_ulIter << ": CCC [" << nTmpX << "][" << nTmpY << "] = {" << dGridCentroidXToExtCRSX(nTmpX) << ", " << dGridCentroidYToExtCRSY(nTmpY) << "} marked as shadow zone " << nZone << endl;
               }
            }

            // // DEBUG CODE =======================================================================================================================
            // LogStream << endl;
            // for (int k = 0; k < ILShadowBoundary.nGetSize(); k++)
            // {
            //    CGeom2DIPoint PtiTmp = *ILShadowBoundary.pPtiGetAt(k);
            //    LogStream << k << " [" << PtiTmp.nGetX() << "][" << PtiTmp.nGetY() << "]" << endl;
            // }
            // LogStream << endl;
            // // DEBUG CODE =======================================================================================================================

            // Calculate the centroid (grid CRS) of the shadow zone
            CGeom2DIPoint PtiCentroid = PtiPolygonCentroid(ILShadowBoundary.pPtiVGetPoints());

            // Safety check
            if (! bIsWithinValidGrid(&PtiCentroid))
            {
               LogStream << ERR << "coast " << nCoast << " start point for shadow zone cell-by-cell fill [" << PtiCentroid.nGetX() << "][" << PtiCentroid.nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiCentroid.nGetX()) << ", " << dGridCentroidYToExtCRSY(PtiCentroid.nGetY()) << "} is outside grid" << endl;

               continue;
            }

            // if (m_nLogFileDetail >= LOG_FILE_ALL)
            //    LogStream << m_ulIter << ":\t coast " << nCoast << " start point for shadow zone cell-by-cell fill is [" << PtiCentroid.nGetX() << "][" << PtiCentroid.nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiCentroid.nGetX()) << ", " << dGridCentroidYToExtCRSY(PtiCentroid.nGetY()) << "}" << endl;

            // Cell-by-cell fill the shadow zone
            int nRet = nCellByCellFillShadowZone(nCoast, nZone, &PtiCentroid, &PtiStart, &PtiShadowEnd, &VPtiShadowCells);
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

                  return RTN_ERR_SHADOW_ZONE_FLOOD_START_POINT;
               }
            }

            // For the downdrift zone boundary
            CGeomILine ILDownDriftBoundary;

            // The downdrift zone cells
            vector<CGeom2DIPoint> VPtiDownDriftCells;

            // Locate the downdrift boundary
            int nAlongCoastDistanceDownDriftStartToEnd;
            nRet = nFindDownDriftBoundaryLine(nCoast, nZone, &PtiStart, nStartCoastPoint, nEndSZCoastPoint, &ILDownDriftBoundary, &VPtiDownDriftCells, nAlongCoastDistanceDownDriftStartToEnd);
            if (nRet != RTN_OK)
               // Do next coast point
               continue;

            // Calculate the centroid (grid CRS) of the downdrift zone
            ILTmp.Reverse();
            ILTmp.Append(&ILDownDriftBoundary);
            PtiCentroid = PtiPolygonCentroid(ILTmp.pPtiVGetPoints());

            // if (m_nLogFileDetail >= LOG_FILE_ALL)
            //    LogStream << m_ulIter << ":\t coast " << nCoast << " start point for downdrift zone cell-by-cell fill is [" << PtiCentroid.nGetX() << "][" << PtiCentroid.nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiCentroid.nGetX()) << ", " << dGridCentroidYToExtCRSY(PtiCentroid.nGetY()) << "}" << endl;

            // Safety check
            if (! bIsWithinValidGrid(&PtiCentroid))
            {
               // LogStream << ERR << "coast " << nCoast << " start point for downdrift zone cell-by-cell fill [" << PtiCentroid.nGetX() << "][" << PtiCentroid.nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiCentroid.nGetX()) << ", " << dGridCentroidYToExtCRSY(PtiCentroid.nGetY()) << "} is outside grid" << endl;

               continue;
            }

            // if (m_nLogFileDetail >= LOG_FILE_ALL)
            //    LogStream << m_ulIter << ":\t coast " << nCoast << " start point for downdrift zone cell-by-cell fill is [" << PtiCentroid.nGetX() << "][" << PtiCentroid.nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiCentroid.nGetX()) << ", " << dGridCentroidYToExtCRSY(PtiCentroid.nGetY()) << "}" << endl;

            // Cell-by-cell fill the downdrift zone
            CGeom2DIPoint const PtiDownDriftEnd = ILDownDriftBoundary.Back();

            nRet = nCellByCellFillDownDriftZone(nCoast, nZone, &PtiCentroid, &PtiStart, &PtiDownDriftEnd, &VPtiDownDriftCells);
            if (nRet != RTN_OK)
            {
               // Could not find start point for cell-by-cell fill. How serious we judge this to be depends on the length of the downdrift zone line
               if (nShadowLineLen < (3 * MAX_LEN_SHADOW_LINE_TO_IGNORE))
               {
                  if (m_nLogFileDetail >= LOG_FILE_ALL)
                     LogStream << m_ulIter << ":\t" << WARN << "could not find start point for cell-by-cell fill of downdrift zone " << nZone << " but continuing because this is a small shadow zone (downdrift line length = " << nShadowLineLen << " cells)" << endl;

                  continue;
               }
               else
               {
                  LogStream << m_ulIter << ":\t " << ERR << "could not find start point for cell-by-cell fill of downdrift zone " << nZone << " (shadow line length = " << nShadowLineLen << " cells)" << endl;

                  return RTN_ERR_DOWNDRIFT_ZONE_FLOOD_START_POINT;
               }
            }

            // Modify waves in the shadow zone (change wave orientation and height)
            int const nCoastSeaHand = m_VCoast[nCoast].nGetSeaHandedness();
            int nShadowZoneCoastToCapeSeaHand;

            if (nCoastSeaHand == LEFT_HANDED)
               nShadowZoneCoastToCapeSeaHand = RIGHT_HANDED;
            else
               nShadowZoneCoastToCapeSeaHand = LEFT_HANDED;

            // The end of the shadow zone boundary
            CGeom2DIPoint const PtiEnd = ILTmp.Back();

            for (int n = 0; n < static_cast<int>(VPtiShadowCells.size()); n++)
            {
               int const nX = VPtiShadowCells[n].nGetX();
               int const nY = VPtiShadowCells[n].nGetY();

               ModifyWavesOnShadowZoneCell(nX, nY, nShadowZoneCoastToCapeSeaHand, &PtiStart, &PtiEnd);
            }

            // Modify waves in the downdrift zone (change only wave height)
            CGeom2DIPoint const PtiSweepStart = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nStartCoastPoint);
            if (nStartCoastPoint < nEndSZCoastPoint)
            {
               // Going down-coast
               int const nEndDDCoastPoint = nEndSZCoastPoint + nAlongCoastDistanceDownDriftStartToEnd;

               int nSweepStart = -1;
               int const nSweepLength = nEndSZCoastPoint - nEndDDCoastPoint;

               for (int n = nEndSZCoastPoint; n <= nEndDDCoastPoint; n++)
               {
                  nSweepStart++;

                  // Is the coast-point end of this sweep outside the grid?
                  if (n >= nCoastSize)
                  {
                     // TODO
                     LogStream << m_ulIter << ":  n = " << n << " coast size = " << nCoastSize << endl;
                     continue;
                  }

                  // Is the coast-point end of this sweep outside the grid?
                  if (n < 0)
                  {
                     // TODO
                     LogStream << m_ulIter << ":  n = " << n << endl;
                     continue;
                  }

                  // Create a line
                  CGeom2DIPoint const PtiSweepEnd = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(n);

                  double const dXStart = dGridCentroidXToExtCRSX(PtiSweepStart.nGetX());
                  double const dYStart = dGridCentroidYToExtCRSY(PtiSweepStart.nGetY());
                  double dX = dXStart;
                  double dY = dYStart;
                  double const dXEnd = dGridCentroidXToExtCRSX(PtiSweepEnd.nGetX());
                  double const dYEnd = dGridCentroidYToExtCRSY(PtiSweepEnd.nGetY());

                  // Interpolate between start and end cells by a simple DDA line algorithm, see http://en.wikipedia.org/wiki/Digital_differential_analyzer_(graphics_algorithm) Note that Bresenham's algorithm gave occasional gaps
                  double dXInc = dXEnd - dXStart;
                  double dYInc = dYEnd - dYStart;
                  double const dLength = tMax(tAbs(dXInc), tAbs(dYInc));

                  dXInc /= dLength;
                  dYInc /= dLength;

                  // Process each interpolated point along the sweep line
                  for (int m = 0; m <= nRound(dLength); m++)
                  {
                     int const nX = nRound(dExtCRSXToGridX(dX));
                     int const nY = nRound(dExtCRSYToGridY(dY));

                     if (! bIsWithinValidGrid(nX, nY))
                     {
                        // Safety check
                        // LogStream << m_ulIter << ":  outside grid [" << nX << "][" << nY << "] = {" << dX << ", " << dY << "}" << endl;

                        dX -= dXInc;
                        dY -= dYInc;
                        break;
                     }

                     if (! m_pRasterGrid->m_Cell[nX][nY].bIsInThisDownDriftZone(nZone))
                     {
                        // Safety check
                        // LogStream << m_ulIter << ":  not in downdrift zone " << nZone << " [" << nX << "][" << nY << "] = {" << dX << ", " << dY << "}" << endl;

                        dX -= dXInc;
                        dY -= dYInc;
                        continue;
                     }

                     ModifyWavesOnDownDriftCell(nX, nY, nSweepStart, nSweepLength);

                     // LogStream << m_ulIter << ":  processing downdrift cell on sweep line [" << nX << "][" << nY << "] = {" << dX << ", " << dY << "}" << endl;

                     dX -= dXInc;
                     dY -= dYInc;
                  }
               }
            }
            else
            {
               // Going up-coast
               int const nEndDDCoastPoint = nEndSZCoastPoint - nAlongCoastDistanceDownDriftStartToEnd;

               int nSweepStart = -1;
               int const nSweepLength = nEndSZCoastPoint - nEndDDCoastPoint;

               for (int n = nEndSZCoastPoint; n >= nEndDDCoastPoint; n--)
               {
                  nSweepStart++;

                  // Is the coast-point end of this sweep outside the grid?
                  if (n >= nCoastSize)
                  {
                     // TODO
                     LogStream << m_ulIter << ":  n = " << n << " coast size = " << nCoastSize << endl;
                     continue;
                  }

                  // Is the coast-point end of this sweep outside the grid?
                  if (n < 0)
                  {
                     // TODO
                     LogStream << m_ulIter << ":  n = " << n << endl;
                     continue;
                  }

                  // Create a line
                  CGeom2DIPoint const PtiSweepEnd = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(n);

                  double const dXStart = dGridCentroidXToExtCRSX(PtiSweepStart.nGetX());
                  double const dYStart = dGridCentroidYToExtCRSY(PtiSweepStart.nGetY());
                  double dX = dXStart;
                  double dY = dYStart;
                  double const dXEnd = dGridCentroidXToExtCRSX(PtiSweepEnd.nGetX());
                  double const dYEnd = dGridCentroidYToExtCRSY(PtiSweepEnd.nGetY());

                  // Interpolate between start and end cells by a simple DDA line algorithm, see http://en.wikipedia.org/wiki/Digital_differential_analyzer_(graphics_algorithm) Note that Bresenham's algorithm gave occasional gaps
                  double dXInc = dXEnd - dXStart;
                  double dYInc = dYEnd - dYStart;
                  double const dLength = tMax(tAbs(dXInc), tAbs(dYInc));

                  dXInc /= dLength;
                  dYInc /= dLength;

                  // Process each interpolated point along the sweep line
                  for (int m = 0; m <= nRound(dLength); m++)
                  {
                     int const nX = nRound(dExtCRSXToGridX(dX));
                     int const nY = nRound(dExtCRSYToGridY(dY));

                     if (! bIsWithinValidGrid(nX, nY))
                     {
                        // Safety check
                        // LogStream << m_ulIter << ":  outside grid [" << nX << "][" << nY << "] = {" << dX << ", " << dY << "}" << endl;

                        dX += dXInc;
                        dY += dYInc;
                        break;
                     }

                     if (! m_pRasterGrid->m_Cell[nX][nY].bIsInThisDownDriftZone(nZone))
                     {
                        // Safety check
                        // LogStream << m_ulIter << ":  not in downdrift zone " << nZone << " [" << nX << "][" << nY << "] = {" << dX << ", " << dY << "}" << endl;

                        dX += dXInc;
                        dY += dYInc;
                        continue;
                     }

                     ModifyWavesOnDownDriftCell(nX, nY, nSweepStart, nSweepLength);

                     // LogStream << m_ulIter << ":  processing downdrift cell on sweep line [" << nX << "][" << nY << "] = {" << dX << ", " << dY << "}" << endl;

                     dX += dXInc;
                     dY += dYInc;
                  }
               }
            }
         }
      }
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Does a cell-by-cell fill of a wave shadow zone, starting from the centroid
//===============================================================================================================================
int CSimulation::nCellByCellFillShadowZone(int const nCoast, int const nZone, CGeom2DIPoint const* pPtiCentroid, CGeom2DIPoint const* pPtiShadowBoundaryStart, CGeom2DIPoint const* pPtiShadowBoundaryEnd, vector<CGeom2DIPoint>* pVPtiShadowCells)
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
         // Find a start point for the cell-by-cell fill. Because shadow zones are generally triangular, start by choosing a low weighting so that the start point is close to the centroid, but a bit towards the coast. If this doesn't work, go further coastwards
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
      // if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
      //    LogStream << m_ulIter << ":\t coast " << nCoast << " " << ERR << "could not find shadow zone cell-by-cell fill start point" << endl;

      return RTN_ERR_SHADOW_ZONE_FLOOD_START_POINT;
   }

   // if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
   //    LogStream << m_ulIter << ":\t coast " << nCoast << " shadow zone cell-by-cell fill start point [" << PtiFloodFillStart.nGetX() << "][" << PtiFloodFillStart.nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiFloodFillStart.nGetX()) << ", " << dGridCentroidYToExtCRSY(PtiFloodFillStart.nGetY()) << "} OK for shadow boundary from [" << pPtiShadowBoundaryStart->nGetX() << "][" << pPtiShadowBoundaryStart->nGetY() << "] = {" << dGridCentroidXToExtCRSX(pPtiShadowBoundaryStart->nGetX()) << ", " << dGridCentroidYToExtCRSY(pPtiShadowBoundaryStart->nGetY()) << "} to [" << pPtiShadowBoundaryEnd->nGetX() << "][" << pPtiShadowBoundaryEnd->nGetY() << "] = {" << dGridCentroidXToExtCRSX(pPtiShadowBoundaryEnd->nGetX()) << ", " << dGridCentroidYToExtCRSY(pPtiShadowBoundaryEnd->nGetY()) << "}" << endl;

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
      int const nY = Pti.nGetY();

      while ((nX >= 0) && m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea() && (! m_pRasterGrid->m_Cell[nX][nY].bIsInThisShadowZone(nZone)) && (! m_pRasterGrid->m_Cell[nX][nY].bIsShadowZone()) && (! m_pRasterGrid->m_Cell[nX][nY].bIsCoastline()))
         nX--;

      nX++;

      bool bSpanAbove = false;
      bool bSpanBelow = false;

      while ((nX < m_nXGridSize) && m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea() && (! m_pRasterGrid->m_Cell[nX][nY].bIsInThisShadowZone(nZone)) && (! m_pRasterGrid->m_Cell[nX][nY].bIsShadowZone()) && (! m_pRasterGrid->m_Cell[nX][nY].bIsCoastline()))
      {
         // Mark this sea cell as wave shadow zone
         m_pRasterGrid->m_Cell[nX][nY].SetShadowZoneNumber(nZone);

         // Save cell, for later wave modification
         pVPtiShadowCells->push_back(CGeom2DIPoint(nX, nY));

         // LogStream << m_ulIter << ": DDD [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "} marked as shadow zone " << nZone << endl;

         if ((! bSpanAbove) && (nY > 0) && m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea() && (! m_pRasterGrid->m_Cell[nX][nY - 1].bIsInThisShadowZone(nZone)) && (! m_pRasterGrid->m_Cell[nX][nY - 1].bIsShadowZone()) && (! m_pRasterGrid->m_Cell[nX][nY - 1].bIsCoastline()))
         {
            PtiStack.push(CGeom2DIPoint(nX, nY - 1));
            bSpanAbove = true;
         }
         else if (bSpanAbove && (nY > 0) && ((! m_pRasterGrid->m_Cell[nX][nY - 1].bIsInContiguousSea()) || m_pRasterGrid->m_Cell[nX][nY - 1].bIsInThisShadowZone(nZone) || m_pRasterGrid->m_Cell[nX][nY - 1].bIsShadowZone() || m_pRasterGrid->m_Cell[nX][nY - 1].bIsCoastline()))
         {
            bSpanAbove = false;
         }

         if ((! bSpanBelow) && m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea() && (nY < m_nYGridSize - 1) && (! m_pRasterGrid->m_Cell[nX][nY + 1].bIsInThisShadowZone(nZone)) && (!m_pRasterGrid->m_Cell[nX][nY + 1].bIsShadowZone()) && (! m_pRasterGrid->m_Cell[nX][nY + 1].bIsCoastline()))
         {
            PtiStack.push(CGeom2DIPoint(nX, nY + 1));
            bSpanBelow = true;
         }
         else if (bSpanBelow && (nY < m_nYGridSize - 1) && ((! m_pRasterGrid->m_Cell[nX][nY + 1].bIsInContiguousSea()) || m_pRasterGrid->m_Cell[nX][nY + 1].bIsInThisShadowZone(nZone - 1) || m_pRasterGrid->m_Cell[nX][nY + 1].bIsShadowZone() || m_pRasterGrid->m_Cell[nX][nY + 1].bIsCoastline()))
         {
            bSpanBelow = false;
         }

         nX++;
      }
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Does a cell-by-cell fill of a downdrift zone, starting from the centroid
//===============================================================================================================================
int CSimulation::nCellByCellFillDownDriftZone(int const nCoast, int const nZone, CGeom2DIPoint const* pPtiCentroid, CGeom2DIPoint const* pPtiShadowBoundaryStart, CGeom2DIPoint const* pPtiShadowBoundaryEnd, vector<CGeom2DIPoint>* pVPtiDownDriftCells)
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
         // Find a start point for the cell-by-cell fill. Because downdrift zones are generally triangular, start by choosing a low weighting so that the start point is close to the centroid, but a bit towards the coast. If this doesn't work, go further coastwards
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
               LogStream << m_ulIter << ": " << ERR << "start point [" << PtiFloodFillStart.nGetX() << "][" << PtiFloodFillStart.nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiFloodFillStart.nGetX()) << ", " << dGridCentroidYToExtCRSY(PtiFloodFillStart.nGetY()) << "} for cell-by-cell fill of downdrift zone is outside grid" << endl;

            return RTN_ERR_DOWNDRIFT_ZONE_FLOOD_FILL_NOGRID;
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
               LogStream << m_ulIter << ":\t coast " << nCoast << " downdrift zone cell-by-cell fill start point [" << PtiFloodFillStart.nGetX() << "][" << PtiFloodFillStart.nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiFloodFillStart.nGetX()) << ", " << dGridCentroidYToExtCRSY(PtiFloodFillStart.nGetY()) << "} is NOT a sea cell for shadow boundary from cape point [" << pPtiShadowBoundaryStart->nGetX() << "][" << pPtiShadowBoundaryStart->nGetY() << "] = {" << dGridCentroidXToExtCRSX(pPtiShadowBoundaryStart->nGetX()) << ", " << dGridCentroidYToExtCRSY(pPtiShadowBoundaryStart->nGetY()) << "} to [" << pPtiShadowBoundaryEnd->nGetX() << "][" << pPtiShadowBoundaryEnd->nGetY() << "] = {" << dGridCentroidXToExtCRSX(pPtiShadowBoundaryEnd->nGetX()) << ", " << dGridCentroidYToExtCRSY(pPtiShadowBoundaryEnd->nGetY()) << "}, dWeight = " << dWeight << endl;

            dWeight += 0.05;
         }
      }
   }

   if ((! bStartPointOK) && (! bAllPointNotSea))
   {
      if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
         LogStream << m_ulIter << ":\t coast " << nCoast << " " << ERR << "could not find downdrift zone cell-by-cell fill start point" << endl;

      return RTN_ERR_DOWNDRIFT_ZONE_FLOOD_START_POINT;
   }

   if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
      LogStream << m_ulIter << ":\t coast " << nCoast << " downdrift zone cell-by-cell fill start point [" << PtiFloodFillStart.nGetX() << "][" << PtiFloodFillStart.nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiFloodFillStart.nGetX()) << ", " << dGridCentroidYToExtCRSY(PtiFloodFillStart.nGetY()) << "} OK for downdrift boundary from [" << pPtiShadowBoundaryStart->nGetX() << "][" << pPtiShadowBoundaryStart->nGetY() << "] = {" << dGridCentroidXToExtCRSX(pPtiShadowBoundaryStart->nGetX()) << ", " << dGridCentroidYToExtCRSY(pPtiShadowBoundaryStart->nGetY()) << "} to [" << pPtiShadowBoundaryEnd->nGetX() << "][" << pPtiShadowBoundaryEnd->nGetY() << "] = {" << dGridCentroidXToExtCRSX(pPtiShadowBoundaryEnd->nGetX()) << ", " << dGridCentroidYToExtCRSY(pPtiShadowBoundaryEnd->nGetY()) << "}" << endl;

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
      int const nY = Pti.nGetY();

      while ((nX >= 0) && m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea() && (! m_pRasterGrid->m_Cell[nX][nY].bIsInThisDownDriftZone(nZone)) && (! m_pRasterGrid->m_Cell[nX][nY].bIsDownDriftZone()) && (! m_pRasterGrid->m_Cell[nX][nY].bIsShadowZone()) && (! m_pRasterGrid->m_Cell[nX][nY].bIsCoastline()))
         nX--;

      nX++;

      bool bSpanAbove = false;
      bool bSpanBelow = false;

      while ((nX < m_nXGridSize) && m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea() && (! m_pRasterGrid->m_Cell[nX][nY].bIsInThisDownDriftZone(nZone)) && (! m_pRasterGrid->m_Cell[nX][nY].bIsShadowZone()) && (! m_pRasterGrid->m_Cell[nX][nY].bIsCoastline()))
      {
         // Mark this sea cell as wave downdrift zone
         m_pRasterGrid->m_Cell[nX][nY].SetDownDriftZoneNumber(nZone);

         // Save cell, for later wave modification
         pVPtiDownDriftCells->push_back(CGeom2DIPoint(nX, nY));

         // LogStream << m_ulIter << ": DDD [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "} marked as downdrift zone " << nZone << endl;

         if ((! bSpanAbove) && (nY > 0) && m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea() && (! m_pRasterGrid->m_Cell[nX][nY - 1].bIsInThisDownDriftZone(nZone)) && (! m_pRasterGrid->m_Cell[nX][nY - 1].bIsDownDriftZone()) && (! m_pRasterGrid->m_Cell[nX][nY].bIsShadowZone()) && (! m_pRasterGrid->m_Cell[nX][nY - 1].bIsCoastline()))
         {
            PtiStack.push(CGeom2DIPoint(nX, nY - 1));
            bSpanAbove = true;
         }
         else if (bSpanAbove && (nY > 0) && ((! m_pRasterGrid->m_Cell[nX][nY - 1].bIsInContiguousSea()) || m_pRasterGrid->m_Cell[nX][nY - 1].bIsInThisDownDriftZone(nZone) || m_pRasterGrid->m_Cell[nX][nY - 1].bIsShadowZone() || m_pRasterGrid->m_Cell[nX][nY - 1].bIsCoastline()))
         {
            bSpanAbove = false;
         }

         if ((! bSpanBelow) && m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea() && (nY < m_nYGridSize - 1) && (! m_pRasterGrid->m_Cell[nX][nY + 1].bIsInThisDownDriftZone(nZone)) && (!m_pRasterGrid->m_Cell[nX][nY + 1].bIsShadowZone()) && (! m_pRasterGrid->m_Cell[nX][nY + 1].bIsCoastline()))
         {
            PtiStack.push(CGeom2DIPoint(nX, nY + 1));
            bSpanBelow = true;
         }
         else if (bSpanBelow && (nY < m_nYGridSize - 1) && ((! m_pRasterGrid->m_Cell[nX][nY + 1].bIsInContiguousSea()) || m_pRasterGrid->m_Cell[nX][nY + 1].bIsInThisDownDriftZone(nZone - 1) || m_pRasterGrid->m_Cell[nX][nY + 1].bIsShadowZone() || m_pRasterGrid->m_Cell[nX][nY + 1].bIsCoastline()))
         {
            bSpanBelow = false;
         }

         nX++;
      }
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Locate the down drift boundary and mark it on the grid
//===============================================================================================================================
int CSimulation::nFindDownDriftBoundaryLine(int const nCoast, int const nZone, CGeom2DIPoint const* pPtiStart, int const nStartPoint, int const nShadowBoundaryEndPoint, CGeomILine* pILDownDriftBoundary, vector<CGeom2DIPoint>* pVPtiDownDriftCells, int& nAlongCoastDistanceToShadowEndpoint)
{
   // We will traverse the coastline from the start point of the shadow zone line, going toward the end point. Which direction is this?
   bool bSweepDownCoast = true;

   if (nShadowBoundaryEndPoint < nStartPoint)
      bSweepDownCoast = false;

   // Get the distance (in cells) from the shadow boundary start point to the shadow boundary end point, going along the coast
   nAlongCoastDistanceToShadowEndpoint = tAbs(nShadowBoundaryEndPoint - nStartPoint - 1);

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
   // if ((nDownDriftEndPoint < 0) || (nDownDriftEndPoint >= m_VCoast[nCoast].nGetCoastlineSize()))
   //    LogStream << "nDownDriftEndPoint not within grid" << endl;
   // else
   // {
   //    CGeom2DIPoint PtiDownDriftEndPointTmp = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nDownDriftEndPoint);
   //    LogStream << "nShadowBoundaryEndPoint = " << nStartPoint << " [" << PtiDownDriftEndPointTmp.nGetX() << "][" << PtiDownDriftEndPointTmp.nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiDownDriftEndPointTmp.nGetX()) << ", " << dGridCentroidYToExtCRSY(PtiDownDriftEndPointTmp.nGetY()) << "}" << endl;
   // }
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

   // Now trace the down-drift boundary line: interpolate between cells by a simple DDA line algorithm, see http://en.wikipedia.org/wiki/Digital_differential_analyzer_(graphics_algorithm) Note that Bresenham's algorithm gave occasional gaps
   int const nXStart = pPtiStart->nGetX();
   int const nYStart = pPtiStart->nGetY();
   int const nXEnd = PtiDownDriftEndPoint.nGetX();
   int const nYEnd = PtiDownDriftEndPoint.nGetY();

   // Safety check
   if ((nXStart == nXEnd) && (nYStart == nYEnd))
      return RTN_ERR_DOWNDRIFT_BOUNDARY_NOGOOD;

   double const dXStart = dGridCentroidXToExtCRSX(nXStart);
   double const dYStart = dGridCentroidYToExtCRSY(nYStart);
   double const dXEnd = dGridCentroidXToExtCRSX(nXEnd);
   double const dYEnd = dGridCentroidYToExtCRSY(nYEnd);
   double dXInc = dXEnd - dXStart;
   double dYInc = dYEnd - dYStart;
   double const dLength = tMax(tAbs(dXInc), tAbs(dYInc));

   dXInc /= dLength;
   dYInc /= dLength;

   int const nLineLen = nRound(dLength / m_dCellSide) + 1;

   double dX = dXStart;
   double dY = dYStart;

   // For the downdrift boundary (ext CRS)
   CGeomLine LDownDriftBoundary;

   // Process each interpolated point
   for (int m = 0; m <= nLineLen; m++)
   {
      int nX = nRound(dExtCRSXToGridX(dX));
      int nY = nRound(dExtCRSYToGridY(dY));

      // Sometimes a grid-edge cell is not marked because of rounding error, soi fix this
      if (nX == -1)
         nX = 0;
      if (nX == m_nXGridSize)
         nX = m_nXGridSize-1;
      if (nY == -1)
         nY = 0;
      if (nY == m_nYGridSize)
         nY = m_nYGridSize-1;

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
         pILDownDriftBoundary->Append(nX, nY);

         // Mark the cell (same as the associated shadow zone number)
         m_pRasterGrid->m_Cell[nX][nY].SetDownDriftZoneNumber(nZone);

         // Save cell, for later wave modification
         pVPtiDownDriftCells->push_back(CGeom2DIPoint(nX, nY));

         // LogStream << "Downdrift boundary [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;

         dX += dXInc;
         dY += dYInc;
      }
   }

   // Store the downdrift boundary (external CRS), with the start point first
   m_VCoast[nCoast].AppendDownDriftBoundary(&LDownDriftBoundary);

   return RTN_OK;
}

//===============================================================================================================================
//! Process a single cell which is in the downdrift zone, changing its wave height
//===============================================================================================================================
void CSimulation::ModifyWavesOnDownDriftCell(int const nX, int const nY, int const nSweep, int const nSweepLength)
{
   // Get the pre-existing (i.e. shore-parallel) wave height
   double const dWaveHeight = m_pRasterGrid->m_Cell[nX][nY].dGetWaveHeight();

   // Safety check
   if (bFPIsEqual(dWaveHeight, DBL_NODATA, TOLERANCE))
   {
      // Is not a sea cell
      // LogStream << m_ulIter << ": [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "} ignored, not a sea cell" << endl;

      return;
   }

   if (m_pRasterGrid->m_Cell[nX][nY].bIsShadowZone())
   {
      // This cell is in a wave shadow zone (can happen occasionally), so don't change it
      // LogStream << m_ulIter << ": [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "} ignored, is in a shadow zone" << endl;

      return;
   }

   // Equation 14 from Hurst et al. TODO 056 Check this! Could not get this to work (typo in paper?), so used the equation below instead
   double const dKp = 1 - (0.5 * (1.0 - sin((PI * 90.0 * nSweep) / (180.0 * nSweepLength))));

   // Set the modified wave height
   m_pRasterGrid->m_Cell[nX][nY].SetWaveHeight(dKp * dWaveHeight);

   // LogStream << m_ulIter << ":  modifying waves on downdrift cell [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}, nSweep = " << nSweep << " nSweepLength = " << nSweepLength << " fraction swept = " << nSweep / nSweepLength << " deep water wave height = " << m_pRasterGrid->m_Cell[nX][nY].dGetCellDeepWaterWaveHeight() << " original dWaveHeight = " << dWaveHeight << " dKp = " << dKp << " new wave height = " << dKp * dWaveHeight << endl;
}

//===============================================================================================================================
//! Process a single cell which is in the shadow zone, changing its wave height and orientation
//===============================================================================================================================
void CSimulation::ModifyWavesOnShadowZoneCell(int const nX, int const nY, int const nShadowZoneCoastToCapeSeaHand, CGeom2DIPoint const* pPtiSDStart, CGeom2DIPoint const* pPtiSDEnd)
{
   // OK, we are in the shadow zone and have not already processed this cell. Set it as a +ve number to show that wave values have been changed
   // m_pRasterGrid->m_Cell[nX][nY].SetShadowZoneNumber(nZone);

   // Next calculate wave angle here: first calculate dOmega, the signed angle subtended between this end point and the start point, and this end point and the end of the shadow boundary
   CGeom2DIPoint const PtiThis(nX, nY);

   double const dOmega = 180 * dAngleSubtended(pPtiSDStart, &PtiThis, pPtiSDEnd) / PI;

   // If dOmega is 90 degrees or more in either direction, set both wave angle and wave height to zero
   if (tAbs(dOmega) >= 90)
   {
      m_pRasterGrid->m_Cell[nX][nY].SetWaveAngle(0);
      m_pRasterGrid->m_Cell[nX][nY].SetWaveHeight(0);

      LogStream << m_ulIter << ":  on shadow linking line with coast end [" << pPtiSDStart->nGetX() << "][" << pPtiSDStart->nGetY() << "] = {" << dGridCentroidXToExtCRSX(pPtiSDStart->nGetX()) << ", " << dGridCentroidYToExtCRSY(pPtiSDStart->nGetY()) << "} and shadow boundary end [" << pPtiSDEnd->nGetX() << "][" << pPtiSDEnd->nGetY() << "] = {" << dGridCentroidXToExtCRSX(pPtiSDEnd->nGetX()) << ", " << dGridCentroidYToExtCRSY(pPtiSDEnd->nGetY()) << "}, this point [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl << "angle subtended = " << dOmega << " degrees, m_pRasterGrid->m_Cell[" << nX << "][" << nY << "].dGetCellDeepWaterWaveHeight() = " << m_pRasterGrid->m_Cell[nX][nY].dGetCellDeepWaterWaveHeight() << " degrees, wave orientation = 0 degrees, wave height = 0 m" << endl;
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

      // LogStream << m_ulIter << ":  modifying waves in shadow zone with coast end [" << pPtiSDStart->nGetX() << "][" << pPtiSDStart->nGetY() << "] = {" << dGridCentroidXToExtCRSX(pPtiSDStart->nGetX()) << ", " << dGridCentroidYToExtCRSY(pPtiSDStart->nGetY()) << "} and shadow boundary end [" << pPtiSDEnd->nGetX() << "][" << pPtiSDEnd->nGetY() << "] = {" << dGridCentroidXToExtCRSX(pPtiSDEnd->nGetX()) << ", " << dGridCentroidYToExtCRSY(pPtiSDEnd->nGetY()) << "}, this point [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}, angle subtended = " << dOmega << " degrees, m_pRasterGrid->m_Cell[" << nX << "][" << nY << "].dGetCellDeepWaterWaveHeight() = " << m_pRasterGrid->m_Cell[nX][nY].dGetCellDeepWaterWaveHeight() << " m, dDeltaShadowWaveAngle = " << dDeltaShadowWaveAngle << " degrees, dWaveAngle = " << dWaveAngle << " degrees, dShadowWaveAngle = " << dShadowWaveAngle << " degrees, dWaveHeight = " << dWaveHeight << " m, dKp = " << dKp << ", shadow zone wave height = " << dKp * dWaveHeight << " m" << endl;
   }
}

//===============================================================================================================================
//! Trace a shadow zone boundary by following waves (the boundary will proably be curved)
//===============================================================================================================================
int CSimulation::nFindShadowZoneBoundaryFollowWave(int const nCoast, int const nStartPoint, int& nEndPoint, bool& bHitEdge, bool& bHitCoast, bool& bHitSea, bool& bStillInland, CGeom2DIPoint const* pPtiStart, double dPrevWaveAngle, CGeomILine* pILShadowBoundary)
{
   pILShadowBoundary->Append(pPtiStart);

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

      // Follow the previous cell's wave orientation to find the new boundary cell
      CGeom2DIPoint const PtiNew = PtiFollowWaveAngle(&PtiPrev, dPrevWaveAngle, dCorrection);

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

      LogStream << m_ulIter << ": at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;

      // Having hit sea, have we now hit we hit a coast point? Note that two diagonal(ish) raster lines can cross each other without any intersection, so must also test an adjacent cell for intersection (does not matter which adjacent cell)
      if (bHitSea)
      {
         if (m_pRasterGrid->m_Cell[nX][nY].bIsCoastline() || (bIsWithinValidGrid(nX, nY + 1) && m_pRasterGrid->m_Cell[nX][nY + 1].bIsCoastline()))
         {
            bHitCoast = true;

            if (m_nLogFileDetail >= LOG_FILE_ALL)
               LogStream << m_ulIter << ":\t coast " << nCoast << " possible shadow boundary from start point " << nStartPoint << " [" << m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nStartPoint)->nGetX() << "][" << m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nStartPoint)->nGetY() <<"][ = {" << dGridCentroidXToExtCRSX(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nStartPoint)->nGetX()) << ", " << dGridCentroidYToExtCRSY(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nStartPoint)->nGetY()) << "} hit the coast at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;
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
         LogStream << m_ulIter << ":\t coast " << nCoast << " possible shadow boundary from start point " << nStartPoint << " [" << m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nStartPoint)->nGetX() << "][" << m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nStartPoint)->nGetY() <<"][ = {" << dGridCentroidXToExtCRSX(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nStartPoint)->nGetX()) << ", " << dGridCentroidYToExtCRSY(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nStartPoint)->nGetY()) << "} is still inland after crossing " << MAX_LAND_LENGTH_OF_SHADOW_ZONE_LINE << " cells, abandoning. Starts at [" << pILShadowBoundary->at(0).nGetX() << "][" << pILShadowBoundary->at(0).nGetY() << "] = {" << dGridCentroidXToExtCRSX(pILShadowBoundary->at(0).nGetX()) << ", " << dGridCentroidYToExtCRSY(pILShadowBoundary->at(0).nGetY()) << "} abandoned at [" << pILShadowBoundary->Back().nGetX() << "][" << pILShadowBoundary->Back().nGetY() << "] = {" << dGridCentroidXToExtCRSX(pILShadowBoundary->Back().nGetX()) << ", " << dGridCentroidYToExtCRSY(pILShadowBoundary->Back().nGetY()) << "}" << endl;

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
            LogStream << m_ulIter << ":\t coast " << nCoast << " possible shadow boundary from start point " << nStartPoint << " [" << m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nStartPoint)->nGetX() << "][" << m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nStartPoint)->nGetY() <<"][ = {" << dGridCentroidXToExtCRSX(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nStartPoint)->nGetX()) << ", " << dGridCentroidYToExtCRSY(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nStartPoint)->nGetY()) << "} is too short. Length " << dShadowLen << " m minimum length " << MIN_LENGTH_OF_SHADOW_ZONE_LINE << " m. Starts at [" << pILShadowBoundary->at(0).nGetX() << "][" << pILShadowBoundary->at(0).nGetY() << "] = {" << dGridCentroidXToExtCRSX(pILShadowBoundary->at(0).nGetX()) << ", " << dGridCentroidYToExtCRSY(pILShadowBoundary->at(0).nGetY()) << "} hits coast at [" << pILShadowBoundary->Back().nGetX() << "][" << pILShadowBoundary->Back().nGetY() << "] = {" << dGridCentroidXToExtCRSX(pILShadowBoundary->Back().nGetX()) << ", " << dGridCentroidYToExtCRSY(pILShadowBoundary->Back().nGetY()) << "}" << endl;

         pILShadowBoundary->Clear();

         return RTN_ERR_SHADOW_BOUNDARY_NOGOOD;
      }

      // We've found a valid shadow zone. Now check the last point in the shadow boundary. Note that occasionally this last cell is not 'above' a cell but is above one of its neighbouring cells is: in which case, replace the last point in the shadow boundary with the coordinates of this neighbouring cell
      nEndPoint = m_VCoast[nCoast].nGetCoastPointGivenCell(&PtiEnd);

      if (nEndPoint == INT_NODATA)
      {
         // Could not find a neighbouring cell which is 'under' the coastline
         if (m_nLogFileDetail >= LOG_FILE_ALL)
            LogStream << m_ulIter << ":\t coast " << nCoast << ", no coast point under shadow boundary end point {" << dGridCentroidXToExtCRSX(pILShadowBoundary->Back().nGetX()) << ", " << dGridCentroidYToExtCRSY(pILShadowBoundary->Back().nGetY()) << "}" << endl;

         // TODO 004 Need to fix this, for the moment just abandon this shadow zone and carry on
         pILShadowBoundary->Clear();

         return RTN_ERR_NO_CELL_UNDER_COASTLINE;
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
   int nXEnd;
   int nYEnd;

   do
   {
      nDist++;

      double const dDeltaX = m_dCellSide * sin(dWaveAngle * (PI/180));
      double const dDeltaY = -m_dCellSide * cos(dWaveAngle * (PI/180));

      nX += nRound(dDeltaX);
      nY += nRound(dDeltaY);

      nXEnd = nX;
      nYEnd = nY;

      LogStream << m_ulIter << ":  in shadow zone boundary loop at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;

      // Have we hit the edge of the valid part of the grid?
      if (! bIsWithinValidGrid(nX, nY))
      {
         // Yes we have
         bHitEdge = true;

         if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
            LogStream << m_ulIter << "\t outside valid grid at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "} " << endl;

         if (CREATE_SHADOW_ZONE_IF_HITS_GRID_EDGE)
         {
            // The shadow boundary hits the grid edge but accept it anyway
            LogStream << m_ulIter << "\t hit grid edge at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;

            break;
         }
         else
         {
            LogStream << m_ulIter << "\t abandoning shadow boundary which starts at [" << nXStart << "][" << nYStart << "] = {" << dGridCentroidXToExtCRSX(nXStart) << ", " << dGridCentroidYToExtCRSY(nYStart) << "}" << endl;

            return RTN_ERR_SHADOW_BOUNDARY_NOGOOD;
         }
      }

      // OK so far. Have we hit a sea cell yet?
      if (! bHitSea)
      {
         if (m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea())
            bHitSea = true;
         else
         {
            if (nDist >= MAX_LAND_LENGTH_OF_SHADOW_ZONE_LINE)
            {
               // If we have travelled MAX_LAND_LENGTH_OF_SHADOW_ZONE_LINE cells without hitting sea, then abandon this shadow boundary
               bStillInland = true;

               return RTN_ERR_SHADOW_BOUNDARY_NOGOOD;
            }
         }
      }
      else
      {
         // Having hit sea, have we now hit a coast point? Note that two diagonal(ish) raster lines can cross each other without any intersection, so must also test an adjacent cell for intersection (does not matter which adjacent cell)
         if (m_pRasterGrid->m_Cell[nX][nY].bIsCoastline() || (bIsWithinValidGrid(nX, nY + 1) && m_pRasterGrid->m_Cell[nX][nY + 1].bIsCoastline()))
         {
            bHitCoast = true;

            if (m_nLogFileDetail >= LOG_FILE_ALL)
               LogStream << m_ulIter << ":\t coast " << nCoast << " possible shadow boundary from start point " << nStartPoint << " [" << m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nStartPoint)->nGetX() << "][" << m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nStartPoint)->nGetY() <<"][ = {" << dGridCentroidXToExtCRSX(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nStartPoint)->nGetX()) << ", " << dGridCentroidYToExtCRSY(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nStartPoint)->nGetY()) << "} hit the coast at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;
         }
      }

      if ((nDist > MIN_LENGTH_OF_SHADOW_ZONE_LINE) && m_pRasterGrid->m_Cell[nX][nY].bIsCoastline())
      {
         LogStream << m_ulIter << "\t hit coast at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;

         break;
      }
   } while (true);

   if (bStillInland)
   {
      // Shadow line is still inland after crossing MAX_LAND_LENGTH_OF_SHADOW_ZONE_LINE calls
      if (m_nLogFileDetail >= LOG_FILE_ALL)
         LogStream << m_ulIter << ":\t coast " << nCoast << " possible shadow boundary from start point " << nStartPoint << " [" << m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nStartPoint)->nGetX() << "][" << m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nStartPoint)->nGetY() <<"][ = {" << dGridCentroidXToExtCRSX(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nStartPoint)->nGetX()) << ", " << dGridCentroidYToExtCRSY(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nStartPoint)->nGetY()) << "} is still inland after crossing " << MAX_LAND_LENGTH_OF_SHADOW_ZONE_LINE << " cells, abandoning. Starts at [" << nXStart << "][" << nYStart << "] = {" << dGridCentroidXToExtCRSX(nXStart) << ", " << dGridCentroidYToExtCRSY(nYStart) << "} abandoned at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;

      return RTN_ERR_SHADOW_BOUNDARY_NOGOOD;
   }

   if (bHitCoast)
   {
      // The shadow zone boundary has hit a coast, but is the shadow zone line trivially short?
      CGeom2DIPoint PtiEnd(nXEnd, nYEnd);
      double const dShadowLen = dGetDistanceBetween(pPtiStart, &PtiEnd) * m_dCellSide;

      if (dShadowLen < MIN_LENGTH_OF_SHADOW_ZONE_LINE)
      {
         // Too short, so forget about it
         if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
            LogStream << m_ulIter << ":\t coast " << nCoast << " possible shadow boundary from start point " << nStartPoint << " [" << m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nStartPoint)->nGetX() << "][" << m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nStartPoint)->nGetY() <<"][ = {" << dGridCentroidXToExtCRSX(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nStartPoint)->nGetX()) << ", " << dGridCentroidYToExtCRSY(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nStartPoint)->nGetY()) << "} is too short. Length " << dShadowLen << " m minimum length " << MIN_LENGTH_OF_SHADOW_ZONE_LINE << " m. Starts at [" << nXStart << "][" << nYStart << "] = {" << dGridCentroidXToExtCRSX(nXStart) << ", " << dGridCentroidYToExtCRSY(nYStart) << "} hits coast at [" << nXEnd << "][" << nYEnd << "] = {" << dGridCentroidXToExtCRSX(nXEnd) << ", " << dGridCentroidYToExtCRSY(nYEnd) << "}" << endl;

         return RTN_ERR_SHADOW_BOUNDARY_NOGOOD;
      }

      // We've found a valid shadow zone. Now check the last point in the shadow boundary. Note that occasionally this last cell is not 'above' a cell but is above one of its neighbouring cells is: in which case, replace the last point in the shadow boundary with the coordinates of this neighbouring cell
      nEndPoint = m_VCoast[nCoast].nGetCoastPointGivenCell(&PtiEnd);

      if (nEndPoint == INT_NODATA)
      {
         // Could not find a neighbouring cell which is 'under' the coastline
         if (m_nLogFileDetail >= LOG_FILE_ALL)
            LogStream << m_ulIter << ":\t coast " << nCoast << ", no coast point under shadow boundary end point [" << PtiEnd.nGetX() << "][" << PtiEnd.nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiEnd.nGetX()) << ", " << dGridCentroidYToExtCRSY(PtiEnd.nGetY()) << "}" << endl;

         // Abandon this shadow zone and carry on
         return RTN_ERR_NO_CELL_UNDER_COASTLINE;
      }

      // Safety check
      if ((nXStart == nXEnd) && (nYStart == nYEnd))
         return RTN_ERR_SHADOW_BOUNDARY_NOGOOD;

      double const dXStart = dGridCentroidXToExtCRSX(nXStart);
      double const dYStart = dGridCentroidYToExtCRSY(nYStart);
      double dX = dXStart;
      double dY = dYStart;
      double const dXEnd = dGridCentroidXToExtCRSX(nXEnd);
      double const dYEnd = dGridCentroidYToExtCRSY(nYEnd);

      // Interpolate between start and end cells by a simple DDA line algorithm, see http://en.wikipedia.org/wiki/Digital_differential_analyzer_(graphics_algorithm) Note that Bresenham's algorithm gave occasional gaps
      double dXInc = dXEnd - dXStart;
      double dYInc = dYEnd - dYStart;
      double const dLength = tMax(tAbs(dXInc), tAbs(dYInc));

      dXInc /= dLength;
      dYInc /= dLength;

      // Process each interpolated point
      for (int m = 0; m <= nRound(dLength); m++)
      {
         nX = nRound(dExtCRSXToGridX(dX));
         nY = nRound(dExtCRSYToGridY(dY));

         if (! bIsWithinValidGrid(nX, nY))
         {
            // Safety check
            break;
         }

         // OK, this is part of the shadow boundary so store this coordinate
         CGeom2DIPoint const PtiThis(nX, nY);

         // Make sure we have not already stored this coordinate (can happen, due to rounding)
         if ((pILShadowBoundary->nGetSize() == 0) || (PtiThis != pILShadowBoundary->Back()))
         {
            // Store this coordinate
            pILShadowBoundary->Append(nX, nY);

            LogStream << m_ulIter << ":  append to shadow boundary [" << nX << "][" << nY << "] = {" << dX << ", " << dY << "}" << endl;

            dX += dXInc;
            dY += dYInc;
         }
      }
   }

   return RTN_OK;
}
