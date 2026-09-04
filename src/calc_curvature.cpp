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
//! Calculates both detailed and smoothed curvature (+ve for convex, -ve for concave, zero if the points are co-linear) for every point on a coastline
//===============================================================================================================================
void CSimulation::DoCoastCurvature(int const nCoast, int const nHandedness)
{
   if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
      LogStream << m_ulIter << ":\t calculating curvatures for coast " << nCoast << endl;

   int const nCoastSize = m_VCoast[nCoast].nGetCoastlineSize();

   // Set an arbitrary curvature (zero) for the first and last coastline points
   m_VCoast[nCoast].SetDetailedCurvature(0, 0);
   m_VCoast[nCoast].SetDetailedCurvature(nCoastSize - 1, 0);

   // Start with detailed curvature, do every point on the coastline, apart from the first and last points
   for (int nThisCoastPoint = 1; nThisCoastPoint < (nCoastSize - 1); nThisCoastPoint++)
   {
      // Calculate the signed curvature based on this point, and the points before and after: +ve for convex, -ve for concave, zero if the points are co-linear
      double const dCurvature = dCalcCurvature(nHandedness, m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nThisCoastPoint - 1), m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nThisCoastPoint), m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nThisCoastPoint + 1));

      // Set the detailed curvature
      m_VCoast[nCoast].SetDetailedCurvature(nThisCoastPoint, dCurvature);
   }

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

   double dMaxConvexDetailed = DBL_MIN;
   double dMaxConvexSmoothed = DBL_MIN;
   double dMaxConcaveDetailed = DBL_MAX;
   double dMaxConcaveSmoothed = DBL_MAX;
   vector<int> VnMaxConvexDetailedCoastPoint;
   vector<int> VnMaxConvexSmoothedCoastPoint;
   vector<int> VnMaxConcaveDetailedCoastPoint;
   vector<int> VnMaxConcaveSmoothedCoastPoint;

   // Find maxima of convexity and concavity for detailed and smoothed curvature
   for (int mm = 0; mm < nCoastSize; mm++)
   {
      // Detailed curvature, maximum convexity
      if (m_VCoast[nCoast].dGetDetailedCurvature(mm) > dMaxConvexDetailed)
      {
         dMaxConvexDetailed = m_VCoast[nCoast].dGetDetailedCurvature(mm);

         VnMaxConvexDetailedCoastPoint.clear();
         VnMaxConvexDetailedCoastPoint.push_back(mm);
      }
      else if (bFPIsEqual(m_VCoast[nCoast].dGetDetailedCurvature(mm), dMaxConvexDetailed, TOLERANCE))
         VnMaxConvexDetailedCoastPoint.push_back(mm);

      // Smoothed curvature, maximum convexity
      if (m_VCoast[nCoast].dGetSmoothCurvature(mm) > dMaxConvexSmoothed)
      {
         dMaxConvexSmoothed = m_VCoast[nCoast].dGetSmoothCurvature(mm);

         VnMaxConvexSmoothedCoastPoint.clear();
         VnMaxConvexSmoothedCoastPoint.push_back(mm);
      }
      else if (bFPIsEqual(m_VCoast[nCoast].dGetSmoothCurvature(mm), dMaxConvexSmoothed, TOLERANCE))
         VnMaxConvexSmoothedCoastPoint.push_back(mm);

      // Detailed curvature, maximum concavity
      if (m_VCoast[nCoast].dGetDetailedCurvature(mm) < dMaxConcaveDetailed)
      {
         dMaxConcaveDetailed = m_VCoast[nCoast].dGetDetailedCurvature(mm);

         VnMaxConcaveDetailedCoastPoint.clear();
         VnMaxConcaveDetailedCoastPoint.push_back(mm);
      }
      else if (bFPIsEqual(m_VCoast[nCoast].dGetDetailedCurvature(mm), dMaxConcaveDetailed, TOLERANCE))
         VnMaxConcaveDetailedCoastPoint.push_back(mm);

      // Smoothed curvature, maximum concavity
      if (m_VCoast[nCoast].dGetSmoothCurvature(mm) < dMaxConcaveSmoothed)
      {
         dMaxConcaveSmoothed = m_VCoast[nCoast].dGetSmoothCurvature(mm);

         VnMaxConcaveSmoothedCoastPoint.clear();
         VnMaxConcaveSmoothedCoastPoint.push_back(mm);
      }
      else if (bFPIsEqual(m_VCoast[nCoast].dGetSmoothCurvature(mm), dMaxConcaveSmoothed, TOLERANCE))
         VnMaxConcaveSmoothedCoastPoint.push_back(mm);
   }

   // Is it a straight coastline?
   if (bFPIsEqual(dMaxConvexDetailed, 0.0, TOLERANCE))
   {
      // We have a straight-line coast, so set the point of maximum convexity and maximum concavity at the coast mid-point
      int const nMidCoastPoint = nCoastSize / 2;

      m_VCoast[nCoast].SetDetailedCurvature(nMidCoastPoint, DUMMY_MAX_CONVEX_DETAILED_CURVE);
      m_VCoast[nCoast].SetSmoothCurvature(nMidCoastPoint, DUMMY_MAX_CONVEX_SMOOTH_CURVE);
      dMaxConvexDetailed = DUMMY_MAX_CONVEX_DETAILED_CURVE;
      dMaxConvexSmoothed = DUMMY_MAX_CONVEX_SMOOTH_CURVE;
      VnMaxConvexDetailedCoastPoint.push_back(nMidCoastPoint);
      VnMaxConvexSmoothedCoastPoint.push_back(nMidCoastPoint);
      VnMaxConcaveDetailedCoastPoint.push_back(nMidCoastPoint);
      VnMaxConcaveSmoothedCoastPoint.push_back(nMidCoastPoint);
   }

   // // DEBUG CODE ============================================
   // LogStream << "-----------------" << endl;
   // for (int kk = 0; kk < m_VCoast.back().nGetCoastlineSize(); kk++)
   //    LogStream << kk << " [" << m_VCoast.back().pPtiGetCellMarkedAsCoastline(kk)->nGetX() << "][" << m_VCoast.back().pPtiGetCellMarkedAsCoastline(kk)->nGetY() << "] = {" << dGridCentroidXToExtCRSX(m_VCoast.back().pPtiGetCellMarkedAsCoastline(kk)->nGetX()) << ", " << dGridCentroidYToExtCRSY(m_VCoast.back().pPtiGetCellMarkedAsCoastline(kk)->nGetY()) << "}\t detailed curvature = " << m_VCoast[nCoast].dGetDetailedCurvature(kk) << "\t\t smooth curvature = " << m_VCoast[nCoast].dGetSmoothCurvature(kk) << endl;
   // LogStream << "-----------------" << endl;
   // // DEBUG CODE ============================================

   if (m_nLogFileDetail >= LOG_FILE_ALL)
   {
      // Write out max detailed convexity
      for (int n = 0; n < static_cast<int>(VnMaxConvexDetailedCoastPoint.size()); n++)
         LogStream << m_ulIter << ":\t  max detailed convexity (" << dMaxConvexDetailed << ") at coast point " << VnMaxConvexDetailedCoastPoint[n] << " [" << m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(VnMaxConvexDetailedCoastPoint[n])->nGetX() << "][" << m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(VnMaxConvexDetailedCoastPoint[n])->nGetY() << "] = {" << dGridCentroidXToExtCRSX(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(VnMaxConvexDetailedCoastPoint[n])->nGetX()) << ", " << dGridCentroidYToExtCRSY(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(VnMaxConvexDetailedCoastPoint[n])->nGetY()) << "}"  << endl;

      // Write out max smoothed convexity
      for (int n = 0; n < static_cast<int>(VnMaxConvexSmoothedCoastPoint.size()); n++)
         LogStream << m_ulIter << ":\t  max smoothed convexity (" << dMaxConvexSmoothed << ") at coast point " << VnMaxConvexSmoothedCoastPoint[n] << " [" << m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(VnMaxConvexSmoothedCoastPoint[n])->nGetX() << "][" << m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(VnMaxConvexSmoothedCoastPoint[n])->nGetY() << "] = {" << dGridCentroidXToExtCRSX(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(VnMaxConvexSmoothedCoastPoint[n])->nGetX()) << ", " << dGridCentroidYToExtCRSY(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(VnMaxConvexSmoothedCoastPoint[n])->nGetY()) << "}"  << endl;

      // Write out max detailed concavity
      for (int n = 0; n < static_cast<int>(VnMaxConcaveDetailedCoastPoint.size()); n++)
         LogStream << m_ulIter << ":\t  max detailed concavity (" << dMaxConcaveDetailed << ") at coast point " << VnMaxConcaveDetailedCoastPoint[n] << " [" << m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(VnMaxConcaveDetailedCoastPoint[n])->nGetX() << "][" << m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(VnMaxConcaveDetailedCoastPoint[n])->nGetY() << "] = {" << dGridCentroidXToExtCRSX(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(VnMaxConcaveDetailedCoastPoint[n])->nGetX()) << ", " << dGridCentroidYToExtCRSY(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(VnMaxConcaveDetailedCoastPoint[n])->nGetY()) << "}"  << endl;

      // Write out max smoothed concavity
      for (int n = 0; n < static_cast<int>(VnMaxConcaveSmoothedCoastPoint.size()); n++)
         LogStream << m_ulIter << ":\t  max smoothed concavity (" << dMaxConcaveSmoothed << ") at coast point " << VnMaxConcaveSmoothedCoastPoint[n] << " [" << m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(VnMaxConcaveSmoothedCoastPoint[n])->nGetX() << "][" << m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(VnMaxConcaveSmoothedCoastPoint[n])->nGetY() << "] = {" << dGridCentroidXToExtCRSX(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(VnMaxConcaveSmoothedCoastPoint[n])->nGetX()) << ", " << dGridCentroidYToExtCRSY(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(VnMaxConcaveSmoothedCoastPoint[n])->nGetY()) << "}"  << endl;
   }
}

//===============================================================================================================================
//! Returns the 2D curvature of three points: +ve for convex, -ve for concave, zero if the points are co-linear. Curvature is multiplied by 1000 to give easier-to-read numbers
//===============================================================================================================================
double CSimulation::dCalcCurvature(int const nHandedness, CGeom2DPoint const* pPtBefore, CGeom2DPoint const* pPtThis, CGeom2DPoint const* pPtAfter)
{
   // Calculate the cross product
   double const dCrosProd = (pPtThis->dGetX() - pPtBefore->dGetX()) * (pPtAfter->dGetY() - pPtBefore->dGetY()) - (pPtThis->dGetY() - pPtBefore->dGetY()) * (pPtAfter->dGetX() - pPtBefore->dGetX());

   // Reverse if left-handed
   int const nSide = (nHandedness == RIGHT_HANDED ? 1 : -1);

   // Make it easier to read
   return 1000 * nSide * dCrosProd;
}
