/*!
   \file calc_curvature.cpp
   \brief Calculates curvature of 2D vectors
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
#include <iostream>
using std::endl;

#include <cfloat>

#include <cmath>
using std::sqrt;
using std::pow;

#include <numeric>
using std::accumulate;
using std::inner_product;

#include "cme.h"
#include "2d_point.h"
#include "simulation.h"
#include "coast.h"

//===============================================================================================================================
//! Calculates both detailed and smoothed curvature for every point on a coastline
//===============================================================================================================================
void CSimulation::DoCoastCurvature(int const nCoast, int const nHandedness)
{
   if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
      LogStream << m_ulIter << ":\t calculating curvatures for coast " << nCoast << endl;

   int const nCoastSize = m_VCoast[nCoast].nGetCoastlineSize();

   // Start with detailed curvature, do every point on the coastline, apart from the first and last points
   for (int nThisCoastPoint = 1; nThisCoastPoint < (nCoastSize - 1); nThisCoastPoint++)
   {
      // Calculate the signed curvature based on this point, and the points before and after
      double const dCurvature = dCalcCurvature(nHandedness, m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nThisCoastPoint - 1), m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nThisCoastPoint), m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nThisCoastPoint + 1));

      // Set the detailed curvature
      m_VCoast[nCoast].SetDetailedCurvature(nThisCoastPoint, dCurvature);
   }

   // Set the curvature for the first and last coastline points
   double dTemp = m_VCoast[nCoast].dGetDetailedCurvature(1);
   m_VCoast[nCoast].SetDetailedCurvature(0, dTemp);

   dTemp = m_VCoast[nCoast].dGetDetailedCurvature(nCoastSize - 2);
   m_VCoast[nCoast].SetDetailedCurvature(nCoastSize - 1, dTemp);

   // Now create the smoothed curvature
   int const nHalfWindow = m_nCoastSmoothingWindowSize / 2;

   // Apply a running mean smoothing filter, with a variable window size at both ends of the coastline
   for (int i = 0; i < nCoastSize; i++)
   {
      int nTmpWindow = 0;
      double dWindowTot = 0;

      for (int j = -nHalfWindow; j < m_nCoastSmoothingWindowSize - nHalfWindow; j++)
      {
         // For points at both ends of the coastline, use a smaller window
         int const k = i + j;

         if ((k < 0) || (k >= nCoastSize))
            continue;

         dWindowTot += m_VCoast[nCoast].dGetDetailedCurvature(k);
         nTmpWindow++;
      }

      m_VCoast[nCoast].SetSmoothCurvature(i, dWindowTot / static_cast<double>(nTmpWindow));
   }

   // Now calculate the mean and standard deviation of each set of curvature values
   vector<double>* pVDetailed = m_VCoast[nCoast].pVGetDetailedCurvature();

   double dSum = accumulate(pVDetailed->begin(), pVDetailed->end(), 0.0);
   double dMean = dSum / static_cast<double>(pVDetailed->size());

   m_VCoast[nCoast].SetDetailedCurvatureMean(dMean);

   double dSquareSum = inner_product(pVDetailed->begin(), pVDetailed->end(), pVDetailed->begin(), 0.0);
   double dSTD = sqrt(dSquareSum / static_cast<double>(pVDetailed->size()) - dMean * dMean);

   m_VCoast[nCoast].SetDetailedCurvatureSTD(dSTD);

   vector<double>* pVSmooth = m_VCoast[nCoast].pVGetSmoothCurvature();
   dSum = accumulate(pVSmooth->begin(), pVSmooth->end(), 0.0),
   dMean = dSum / static_cast<double>(pVSmooth->size());

   m_VCoast[nCoast].SetSmoothCurvatureMean(dMean);

   dSquareSum = inner_product(pVSmooth->begin(), pVSmooth->end(), pVSmooth->begin(), 0.0), dSTD = sqrt(dSquareSum / static_cast<double>(pVSmooth->size()) - dMean * dMean);
   m_VCoast[nCoast].SetSmoothCurvatureSTD(dSTD);

   double dMaxConvexDetailed = DBL_MAX;
   double dMaxConvexSmoothed = DBL_MAX;

   for (int mm = 0; mm < nCoastSize; mm++)
   {
      if (m_VCoast[nCoast].dGetDetailedCurvature(mm) < dMaxConvexDetailed)
      {
         dMaxConvexDetailed = m_VCoast[nCoast].dGetDetailedCurvature(mm);
      }

      if (m_VCoast[nCoast].dGetSmoothCurvature(mm) < dMaxConvexSmoothed)
      {
         dMaxConvexSmoothed = m_VCoast[nCoast].dGetSmoothCurvature(mm);
         // nMaxConvexSmoothedCoastPoint = mm;
      }
   }

   if (bFPIsEqual(dMaxConvexDetailed, 0.0, TOLERANCE))
   {
      // We have a straight-line coast, so set the point of maximum convexity at the coast mid-point
      int const nMaxConvexCoastPoint = nCoastSize / 2;

      m_VCoast[nCoast].SetDetailedCurvature(nMaxConvexCoastPoint, STRAIGHT_COAST_MAX_DETAILED_CURVATURE);
      m_VCoast[nCoast].SetSmoothCurvature(nMaxConvexCoastPoint, STRAIGHT_COAST_MAX_SMOOTH_CURVATURE);
   }

   // LogStream << "-----------------" << endl;
   // for (int kk = 0; kk < m_VCoast.back().nGetCoastlineSize(); kk++)
   //    LogStream << kk << " [" << m_VCoast.back().pPtiGetCellMarkedAsCoastline(kk)->nGetX() << "][" << m_VCoast.back().pPtiGetCellMarkedAsCoastline(kk)->nGetY() << "] = {" << dGridCentroidXToExtCRSX(m_VCoast.back().pPtiGetCellMarkedAsCoastline(kk)->nGetX()) << ", " << dGridCentroidYToExtCRSY(m_VCoast.back().pPtiGetCellMarkedAsCoastline(kk)->nGetY()) << "} detailed curvature = " << m_VCoast[nCoast].dGetDetailedCurvature(kk) << " smooth curvature = " << m_VCoast[nCoast].dGetSmoothCurvature(kk) << endl;
   // LogStream << "-----------------" << endl;

   // CGeom2DIPoint PtiMax = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nMaxConvexDetailedCoastPoint);

   // if (m_nLogFileDetail >= LOG_FILE_ALL)
   // LogStream << m_ulIter << ": Max detailed convexity (" << m_VCoast[nCoast].dGetDetailedCurvature(nMaxConvexDetailedCoastPoint) << ") at raster coastline point " << nMaxConvexDetailedCoastPoint << " [" << PtiMax.nGetX() << "][" << PtiMax.nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiMax.nGetX()) << ", " << dGridCentroidYToExtCRSY(PtiMax.nGetY()) << "}"  << endl;

   // CGeom2DIPoint PtiMaxSmooth = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nMaxConvexSmoothedCoastPoint);
   // CGeom2DPoint PtMaxSmooth = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nMaxConvexSmoothedCoastPoint);

   // if (m_nLogFileDetail >= LOG_FILE_ALL)
   // LogStream << m_ulIter << ": Max smoothed convexity (" << m_VCoast[nCoast].dGetSmoothCurvature(nMaxConvexSmoothedCoastPoint) << ") near vector coastline point " << nMaxConvexSmoothedCoastPoint << ", at [" << PtiMaxSmooth.nGetX() << "][" << PtiMaxSmooth.nGetY() << "] = {" << PtMaxSmooth.dGetX() << ", " << PtMaxSmooth.dGetY() << "}" << endl;
}

//===============================================================================================================================
//! Returns the 2D curvature of three points: +ve for convex, -ve for concave, or zero if the points are co-linear. Curvature is multiplied by 1000 to give easier-to-read numbers
//===============================================================================================================================
double CSimulation::dCalcCurvature(int const nHandedness, CGeom2DPoint const* pPtBefore, CGeom2DPoint const* pPtThis, CGeom2DPoint const* pPtAfter)
{
   // Calculate the cross product
   double const dCrosProd = (pPtThis->dGetX() - pPtBefore->dGetX()) * (pPtAfter->dGetY() - pPtBefore->dGetY()) - (pPtThis->dGetY() - pPtBefore->dGetY()) * (pPtAfter->dGetX() - pPtBefore->dGetX());

   // Reverse if left-handed
   int const nSide = (nHandedness == RIGHT_HANDED ? 1 : -1);

   return 1000 * nSide * dCrosProd;
}
