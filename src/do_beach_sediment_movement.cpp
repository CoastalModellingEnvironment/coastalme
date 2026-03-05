/*!
   \file do_beach_sediment_movement.cpp
   \brief Does between-polygon actual (supply-limited) redistribution of transported beach sediment
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

#include <iostream>
using std::endl;

#include <algorithm>
using std::stable_sort;

#include "cme.h"
#include "simulation.h"
#include "coast.h"

namespace
{
//===============================================================================================================================
//! Function used to sort polygons for a single coast before doing the polygon-to-polygon source-target pattern. For both LHS and RHS arguments, the first element is the polygon ID, the second element is the down- or up-coast direction, and third and subsequent elements are adjacent polygon coast IDs in that direction. If the first argument is to be ordered before the second (i.e. the original sequence is to be retained), return true. If the second argument must be ordered before the first (i.e. the arguments must be swapped), return false
//===============================================================================================================================
bool bPolygonAndAdjCompare(const vector<int>& nVLeft, const vector<int>& nVRight)
{
   // Each row (vector) of this vector-of-vectors is:
   // 0: This-polygon ID (in down-coast seq when sorted)
   // 1: This-polygon down-coast (true) or up-coast (false) sediment movement
   // 2 and subsequent: if sediment movement is down-coast, IDs of down-coast adjacent polygons; if sediment movement is up-coast, IDs of up-coast adjacent polygons

   bool const bDownCoastLeft = nVLeft[1];
   bool const bDownCoastRight = nVRight[1];

   if (bDownCoastLeft)
   {
      // LHS polygon is down-coast. First, deal with polygon 0
      if (nVLeft[0] == 0)
         // Do not swap LHS and RHS polygons
         return true;

      if (nVRight[0] == 0)
         // Swap LHS and RHS polygons
         return false;

      if ((nVLeft.size() >= 3) && (nVRight.size() >= 3))
      {
         // Search the adjacent polygons of the LHS argument, is there an "off edge" amongst them?
         bool bLHSOffEdge = false;

         for (unsigned int n = 2; n < nVLeft.size(); n++)
         {
            if (nVLeft[n] == INT_NODATA)
            {
               bLHSOffEdge = true;
               break;
            }
         }

         if (bLHSOffEdge)
            // Yes, there is an "off edge" among the adjacent polygons of the LHS polygon: so swap LHS and RHS polygons
            return false;

         // Search the adjacent polygons of the RHS argument, is there an "off edge" amongst them?
         bool bRHSOffEdge = false;

         for (unsigned int n = 2; n < nVRight.size(); n++)
         {
            if (nVRight[n] == INT_NODATA)
            {
               bRHSOffEdge = true;
               break;
            }
         }

         if (bRHSOffEdge)
            // Yes, there is an "off edge" among the adjacent polygons of the LHS polygon: do not swap LHS and RHS polygons
            return true;

         // Do we have a local change of sediment direction?
         if (! bDownCoastRight)
            // Yes, there is a change of sediment direction between this and the adjacent polygon. Keep the existing sequence
            return true;

         // Now sort out polygon-to-polygon dependencies. We need to put 'target' polygons after 'source' polygons, so that the source is processed before the target. So is the RHS polygon among the adjacent polygons of the LHS polygon?
         bool bLeftFound = false;

         for (unsigned int n = 2; n < nVLeft.size(); n++)
         {
            if (nVLeft[n] == nVRight[0])
               bLeftFound = true;
         }

         if (bLeftFound)
            // Yes, the RHS polygon is among the adjacent polygons of the LHS polygon, so uncons sediment movement is from the LHS polygon to the RHS polygon, and so down-coast (i.e. along the coast in the direction of increasing coastline point numbers). Keep the existing sequence
            return true;

         // Is the LHS polygon among the adjacent polygons of the RHS polygon?
         bool bRightFound = false;

         for (unsigned int n = 2; n < nVRight.size(); n++)
         {
            if (nVRight[n] == nVLeft[0])
               bRightFound = true;
         }

         if (bRightFound)
            // Yes, the LHS polygon is among the adjacent polygons of the RHS polygon, so uncons sediment movement is from the RHS polygon to the LHS polygon, and so up-coast (i.e. along the coast in the direction of decreasing coastline point numbers). Swap them
            return false;
      }

      // Neither polygon has an "off edge", and each polygon does not have the other polygon's coast ID amongst its list of adjacent polygons. So just put the polygon in increasing own-coast sequence
      if (nVLeft[0] < nVRight[0])
         return true;
      else
         return false;
   }
   else
   {
      // LHS polygon is up-coast. First, deal with polygon 0
      if (nVLeft[0] == 0)
         // Swap LHS and RHS polygons
         return false;

      if (nVRight[0] == 0)
         // Do not swap LHS and RHS polygons
         return true;

      if ((nVLeft.size() >= 3) && (nVRight.size() >= 3))
      {
         // Next, search the adjacent polygons of the LHS argument, is there an "off edge" amongst them?
         bool bLHSOffEdge = false;

         for (unsigned int n = 2; n < nVLeft.size(); n++)
         {
            if (nVLeft[n] == INT_NODATA)
            {
               bLHSOffEdge = true;
               break;
            }
         }

         if (bLHSOffEdge)
            // Yes, there is an "off edge" among the adjacent polygons of the LHS polygon: so do not swap LHS and RHS polygons
            return true;

         // Search the adjacent polygons of the RHS argument, is there an "off edge" amongst them?
         bool bRHSOffEdge = false;

         for (unsigned int n = 2; n < nVRight.size(); n++)
         {
            if (nVRight[n] == INT_NODATA)
            {
               bRHSOffEdge = true;
               break;
            }
         }

         if (bRHSOffEdge)
            // Yes, there is an "off edge" among the adjacent polygons of the LHS polygon: swap LHS and RHS polygons
            return false;

         // Do we have a local change of sediment direction?
         if (! bDownCoastRight)
            // Yes, there is a change of sediment direction between this and the adjacent polygon. Keep the existing sequence
            return true;

         // Now sort out polygon-to-polygon dependencies. We need to put 'target' polygons after 'source' polygons, so that the source is processed before the target. So is the RHS polygon among the adjacent polygons of the LHS polygon?
         bool bLeftFound = false;

         for (unsigned int n = 2; n < nVLeft.size(); n++)
         {
            if (nVLeft[n] == nVRight[0])
               bLeftFound = true;
         }

         if (bLeftFound)
            // Yes, the RHS polygon is among the adjacent polygons of the LHS polygon, so uncons sediment movement is from the LHS polygon to the RHS polygon, and so down-coast (i.e. along the coast in the direction of increasing coastline point numbers). Keep the existing sequence of polygons
            return true;

         // Is the LHS polygon among the adjacent polygons of the RHS polygon?
         bool bRightFound = false;

         for (unsigned int n = 2; n < nVRight.size(); n++)
         {
            if (nVRight[n] == nVLeft[0])
               bRightFound = true;
         }

         if (bRightFound)
            // Yes, the LHS polygon is among the adjacent polygons of the RHS polygon, so uncons sediment movement is from the RHS polygon to the LHS polygon, and so up-coast (i.e. along the coast in the direction of decreasing coastline point numbers). Swap the LHS and RHS polygons
            return false;
      }

      // Neither polygon has an "off edge", and each polygon does not have the other polygon's coast ID amongst its list of adjacent polygons. So just put the polygons in decreasing own-coast sequence
      if (nVLeft[0] < nVRight[0])
         return false;
      else
         return true;
   }

   // Default return value, should never get here
   return true;
}
} // namespace

//===============================================================================================================================
//! Does between-polygon and within-polygon actual (supply-limited) redistribution of transported beach sediment
//===============================================================================================================================
int CSimulation::nDoAllActualBeachErosionAndDeposition(void)
{
   if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
      LogStream << m_ulIter << ": Calculating unconsolidated sediment transport" << endl;

   int const nCoastSize = static_cast<int>(m_VCoast.size());

   for (int nCoast = 0; nCoast < nCoastSize; nCoast++)
   {
      // Update the values of pre-existing unconsolidated sediment, for all three size classes, to include unconsolidated sediment derived from platform erosion, cliff collapse, and sediment input events
      AllPolygonsUpdateStoredUncons(nCoast);
   }

   if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
   {
      WritePolygonSedimentBeforeMovement();
      WritePolygonPotentialErosion();
   }

   // Here we set up a 3D vector (i.e. a vector-of-vectors-of-vectors) to sort the polygons. The innermost (1D) vector (nVPolyAndAdj) comprises these elements, for a single polygon:
   // 0 This-polygon coast ID (in down-coast seq when sorted)
   // 1 This-polygon down-coast (true) or up-coast (false) sediment movement
   // 2 and subsequent: if sediment movement is down-coast, coast IDs of down-coast adjacent polygons; if sediment movement is up-coast, coast IDs of up-coast adjacent polygons
   vector<int> nVPolyAndAdj;

   // The middle-most (2D) vector comprises nVPolyAndAdj vectors, for a single coast
   vector<vector<int>> nVVCoastPolyAndAdjacent;

   // The outermost (3D) vector comprises nVVCoastPolyAndAdjacent for all coasts
   vector<vector<vector<int>>> nVVVAllCoastPolyAndAdjacent;

   // Now preparte to route actually-eroded sand/coarse sediment to adjacent polygons, or off-grid. Sort polygons first
   for (int nCoast = 0; nCoast < nCoastSize; nCoast++)
   {
      // For each coast, update the values of pre-existing unconsolidated sediment, for all three size classes, to include unconsolidated sediment derived from platform erosion, cliff collapse, and sediment input events
      AllPolygonsUpdateStoredUncons(nCoast);

      nVVCoastPolyAndAdjacent.clear();
      for (int nn = 0; nn < m_VCoast[nCoast].nGetNumPolygons(); nn++)
      {
         CGeomCoastPolygon const* pPolygon = m_VCoast[nCoast].pGetPolygon(nn);
         nVPolyAndAdj.clear();

         // The first [0] nVPolyAndAdj array item is the polygon's down-coast ID
         nVPolyAndAdj.push_back(pPolygon->nGetPolygonCoastID());

         if (pPolygon->bDownCoastThisIter())
         {
            // Sediment is leaving this polygon in a down-coast direction. Set this as the second [1] nVPolyAndAdj array item
            nVPolyAndAdj.push_back(true);

            // Third [2] and subsequent nVPolyAndAdj array items are the down-coast seqs of adjacent down-coast polygons
            for (int nAdj = 0; nAdj < pPolygon->nGetNumDownCoastAdjacentPolygons(); nAdj++)
            {
               // Save the coast ID of each down-coast adjacent polygon
               int const nAdjPolyID = pPolygon->nGetDownCoastAdjacentPolygon(nAdj);
               nVPolyAndAdj.push_back(nAdjPolyID);
            }
         }
         else
         {
            // Sediment is leaving this polygon in an up-coast direction. Set this as the second [1] nVPolyAndAdj array item
            nVPolyAndAdj.push_back(false);

            // Third [2] and susequent nVPolyAndAdj array items are the down-coast IDs of adjacent up-coast polygons
            for (int nAdj = 0; nAdj < pPolygon->nGetNumUpCoastAdjacentPolygons(); nAdj++)
            {
               // Save the coast ID of each up-coast adjacent polygon
               int const nAdjPolyID = pPolygon->nGetUpCoastAdjacentPolygon(nAdj);
               nVPolyAndAdj.push_back(nAdjPolyID);
            }
         }
         // Save info for this polygon
         nVVCoastPolyAndAdjacent.push_back(nVPolyAndAdj);
      }
      // Save info for this coast
      nVVVAllCoastPolyAndAdjacent.push_back(nVVCoastPolyAndAdjacent);
   }

   if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
   {
      // Write out the unsorted polygon sequence for all coasts
      WritePolygonUnsortedSequence(nVVVAllCoastPolyAndAdjacent);
   }

   for (int nCoast = 0; nCoast < nCoastSize; nCoast++)
   {
      // OK, now sort the array using bPolygonAndAdjCompare(), so that 'target' polygons are processed after 'source' polygons NOTE: crashes if just use "sort", related to this? https://stackoverflow.com/questions/18291620/why-will-stdsort-crash-if-the-comparison-function-is-not-as-operator
      stable_sort(nVVVAllCoastPolyAndAdjacent[nCoast].begin(), nVVVAllCoastPolyAndAdjacent[nCoast].end(), bPolygonAndAdjCompare);

      // And check for circularities i.e. where poly X -> poly Y -> poly X. Note that we only look for two-way circularities. i.e. we ignore poly A -> poly B -> Poly C -> poly A patterns. These are probably pretty rare, however
      vector<int> VnSourcePolygons;

      for (int n = 0; n < static_cast<int>(nVVVAllCoastPolyAndAdjacent[nCoast].size()); n++)
      {
         int const nThisPoly = nVVVAllCoastPolyAndAdjacent[nCoast][n][0];
         VnSourcePolygons.push_back(nThisPoly);

         for (int m = 2; m < static_cast<int>(nVVVAllCoastPolyAndAdjacent[nCoast][n].size()); m++)
         {
            // Check the adjacent polygon(s) for circularities
            int const nToFind = nVVVAllCoastPolyAndAdjacent[nCoast][n][m];
            vector<int>::iterator const it = find(VnSourcePolygons.begin(), VnSourcePolygons.end(), nToFind);

            if (it != VnSourcePolygons.end())
            {
               // Uh-oh: this adjacent polygon is in the list of previously-processed source polygons. So store the coast ID numbers of the polygons with circularity in both polygons
               CGeomCoastPolygon* pPoly = m_VCoast[nCoast].pGetPolygon(nThisPoly);
               pPoly->AddCircularity(nToFind);

               pPoly = m_VCoast[nCoast].pGetPolygon(nToFind);
               pPoly->AddCircularity(nThisPoly);
            }
         }
      }
   }

   if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
   {
      // Write out the sorted polygon-to-polygon sequence for all coasts
      WritePolygonSortedSequence(nVVVAllCoastPolyAndAdjacent);
   }

   // Now process, for each coast, all polygons and do the actual (supply-limited) unconsolidated sediment movement
   for (int nCoast = 0; nCoast < nCoastSize; nCoast++)
   {
      int nRet;
      int const nNumPolygons = m_VCoast[nCoast].nGetNumPolygons();

      for (int nPoly = 0; nPoly < nNumPolygons; nPoly++)
      {
         int const nPolygon = nVVVAllCoastPolyAndAdjacent[nCoast][nPoly][0];
         CGeomCoastPolygon* pPolygon = m_VCoast[nCoast].pGetPolygon(nPolygon);

         // // DEBUG CODE =====================
         // int nUpCoastProfile = pPolygon->nGetUpCoastProfile();
         // CGeomProfile* pUpCoastProfile = m_VCoast[nCoast].pGetProfile(nUpCoastProfile);
         // int nDownCoastProfile = pPolygon->nGetDownCoastProfile();
         // CGeomProfile* pDownCoastProfile = m_VCoast[nCoast].pGetProfile(nDownCoastProfile);
         // int nNumUpCoastCell = pUpCoastProfile->nGetNumCellsInProfile();
         // int nNumDownCoastCell = pDownCoastProfile->nGetNumCellsInProfile();
         // LogStream << "pUpCoastProfile->nGetNumCellsInProfile() = " << nNumUpCoastCell << " pDownCoastProfile->nGetNumCellsInProfile() = " << nNumDownCoastCell << endl;
         // // DEBUG CODE =====================

         // // DEBUG CODE ============================================================================================================================================
         // // Get total depths of sand consolidated and unconsolidated for every cell
         // if (m_ulIter == 5)
         // {
         // double dTmpSandUncons = 0;
         // for (int nX1 = 0; nX1 < m_nXGridSize; nX1++)
         // {
         // for (int nY1 = 0; nY1 < m_nYGridSize; nY1++)
         // {
         // dTmpSandUncons += m_pRasterGrid->m_Cell[nX1][nY1].dGetUnconsSandDepthAllLayers();
         // }
         // }
         //
         // LogStream << endl;
         // LogStream << "*****************************" << endl;
         // LogStream << m_ulIter << ": before beach movement on nPoly = " << nPoly << " TOTAL UNCONSOLIDATED SAND ON ALL CELLS = " << dTmpSandUncons * m_dCellArea << endl;
         // }
         // // DEBUG CODE ============================================================================================================================================

         // Do deposition first: does this polygon have coarse deposition?
         double dCoarseDepositionTarget = pPolygon->dGetToDoBeachDepositionUnconsCoarse();

         if (dCoarseDepositionTarget > 0)
         {
            // It does, first tho', if we have some coarse sediment which we were unable to deposit on the previously-processed polygon (which could be the last-processed polygon of the previous timestep), then add this in
            if (m_dDepositionCoarseDiff > MASS_BALANCE_TOLERANCE)
            {
               // We had some coarse unconsolidated sediment which we were unable to deoosit on the last polygon (which could have been the last polygon of the previous iteration). So add it to the total to be deposited here
               if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
                  LogStream << m_ulIter << ": nPoly = " << nPolygon << " dCoarseDepositionTarget was = " << dCoarseDepositionTarget * m_dCellArea << " adding m_dDepositionCoarseDiff = " << m_dDepositionCoarseDiff * m_dCellArea;

               dCoarseDepositionTarget += m_dDepositionCoarseDiff;
               m_dDepositionCoarseDiff = 0;

               if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
                  LogStream << " dCoarseDepositionTarget now = " << dCoarseDepositionTarget << endl;
            }

            // OK, do deposition of coarse sediment: calculate a net increase in depth of coarse-sized unconsolidated sediment on the cells within the polygon. Note that some cells may decrease in elevation (i.e. have some coarse-sized sediment erosion) however
            double dCoarseDeposited = 0;
            nRet = nDoUnconsDepositionOnPolygon(nCoast, pPolygon, TEXTURE_COARSE, dCoarseDepositionTarget, dCoarseDeposited);

            if (nRet != RTN_OK)
               return nRet;

            double const dCoarseNotDeposited = dCoarseDepositionTarget - dCoarseDeposited;

            if (dCoarseNotDeposited > 0)
            {
               m_dDepositionCoarseDiff += dCoarseNotDeposited;
            }
         }

         // Does this polygon have sand deposition?
         double dSandDepositionTarget = pPolygon->dGetToDoBeachDepositionUnconsSand();

         if (dSandDepositionTarget > 0)
         {
            // It does, first tho', if we have some sand sediment which we were unable to deposit on the previously-processed polygon (which could be the last-processed polygon of the previous timestep), then add this in
            if (m_dDepositionSandDiff > MASS_BALANCE_TOLERANCE)
            {
               // We had some sand unconsolidated sediment which we were unable to deoosit on the last polygon (which could have been the last polygon of the previous iteration). So add it to the total to be deposited here
               if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
                  LogStream << m_ulIter << ":\t nPolygon = " << nPolygon << " dSandDepositionTarget was = " << dSandDepositionTarget * m_dCellArea << " adding m_dDepositionSandDiff = " << m_dDepositionSandDiff * m_dCellArea;

               dSandDepositionTarget += m_dDepositionSandDiff;
               m_dDepositionSandDiff = 0;

               if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
                  LogStream << " dSandDepositionTarget now = " << dSandDepositionTarget << endl;
            }

            // Now do deposition of sand sediment: calculate a net increase in depth of sand-sized unconsolidated sediment on the cells within the polygon. Note that some cells may decrease in elevation (i.e. have some sand-sized sediment erosion) however
            double dSandDeposited = 0;
            nRet = nDoUnconsDepositionOnPolygon(nCoast, pPolygon, TEXTURE_SAND, dSandDepositionTarget, dSandDeposited);

            if (nRet != RTN_OK)
               return nRet;

            double const dSandNotDeposited = dSandDepositionTarget - dSandDeposited;

            if (dSandNotDeposited > 0)
            {
               m_dDepositionSandDiff += dSandNotDeposited;
            }

            // // DEBUG CODE #####################
            // if (m_ulIter == 5)
            // {
            // LogStream << m_ulIter << ": after sand deposition on nPoly = " << nPoly << " dSandDepositionTarget = " << dSandDepositionTarget * m_dCellArea << " dSandDeposited = " << dSandDeposited * m_dCellArea << " dSandNotDeposited = " << dSandNotDeposited * m_dCellArea << " m_dDepositionSandDiff = " << m_dDepositionSandDiff * m_dCellArea << endl;
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
            // LogStream << m_ulIter << ": after sand deposition on nPoly = " << nPoly << " TOTAL UNCONSOLIDATED SAND ON ALL CELLS = " << dTmpSandUncons * m_dCellArea << endl;
            // }
            // // DEBUG CODE #####################
         }

         // Now do erosion
         double const dPotentialErosion = -pPolygon->dGetPotentialErosion();

         if (dPotentialErosion > 0)
         {
            // There is some erosion on this polygon: process this in the sequence fine, sand, coarse. Is there any fine sediment on this polygon?
            double const dExistingUnconsFine = pPolygon->dGetPreExistingUnconsFine();

            if (dExistingUnconsFine > 0)
            {
               // Yes there is, so crudely partition this potential value for this size class by erodibility, the result will almost always be much greater than actual (supply limited) erosion
               double const dFinePotentialErosion = dPotentialErosion * m_dFineErodibilityNormalized;

               // Now reduce this further, by considering the total depth of fine sediment on the polygon
               double const dFineErosionTarget = tMin(dFinePotentialErosion, dExistingUnconsFine);

               // OK, do the supply-limited erosion of fine sediment
               double dFineEroded = 0;
               nRet = nDoUnconsErosionOnPolygon(nCoast, pPolygon, TEXTURE_FINE, dFineErosionTarget, dFineEroded);

               if (nRet != RTN_OK)
                  return nRet;

               if (dFineEroded > 0)
               {
                  // We eroded some fine sediment, so add to the this-iteration total. Note that total this gets added in to the suspended load elsewhere, so no need to do it here
                  m_dThisIterBeachErosionFine += dFineEroded;

                  // Also add to the suspended load
                  m_dThisIterFineSedimentToSuspension += dFineEroded;

                  // Store the amount of unconsolidated fine beach sediment eroded for this polygon
                  pPolygon->SetBeachErosionUnconsFine(-dFineEroded);
               }
            }

            // Is there any sand-sized sediment on this polygon?
            double const dExistingUnconsSand = pPolygon->dGetPreExistingUnconsSand();
            double dSandEroded = 0;

            if (dExistingUnconsSand > 0)
            {
               // There is: so crudely partition this potential value for this size class by erodibility, the result will almost always be much greater than actual (supply limited) erosion
               double const dSandPotentialErosion = dPotentialErosion * m_dSandErodibilityNormalized;

               // Now reduce this further, by considering the total depth of sand sediment on the polygon
               double const dSandErosionTarget = tMin(dSandPotentialErosion, dExistingUnconsSand);

               // OK, do the supply-limited erosion of sand sediment
               nRet = nDoUnconsErosionOnPolygon(nCoast, pPolygon, TEXTURE_SAND, dSandErosionTarget, dSandEroded);

               if (nRet != RTN_OK)
                  return nRet;

               if (dSandEroded > 0)
               {
                  // We eroded some sand sediment, so add to the this-iteration total
                  m_dThisIterBeachErosionSand += dSandEroded;

                  // Store the amount eroded for this polygon
                  pPolygon->SetBeachErosionUnconsSand(-dSandEroded);
               }

               // // DEBUG CODE #####################
               // if (m_ulIter == 5)
               // {
               // LogStream << m_ulIter << ": nPoly = " << nPoly << " dSandErosionTarget = " << dSandErosionTarget * m_dCellArea << " m_dDepositionSandDiff = " << m_dDepositionSandDiff * m_dCellArea << " dSandEroded = " << dSandEroded * m_dCellArea << " m_dThisIterBeachErosionSand = " << m_dThisIterBeachErosionSand * m_dCellArea << endl;
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
               // LogStream << m_ulIter << ": after sand erosion on nPoly = " << nPoly << " TOTAL UNCONSOLIDATED SAND ON ALL CELLS = " << dTmpSandUncons * m_dCellArea << endl;
               // }
               // // DEBUG CODE #####################
            }

            // Is there any coarse sediment on this polygon?
            double const dExistingUnconsCoarse = pPolygon->dGetPreExistingUnconsCoarse();
            double dCoarseEroded = 0;

            if (dExistingUnconsCoarse > 0)
            {
               // There is: so crudely partition this potential value for this size class by erodibility, the result will almost always be much greater than actual (supply limited) erosion
               double const dCoarsePotentialErosion = dPotentialErosion * m_dCoarseErodibilityNormalized;

               // Now reduce this further, by considering the total depth of coarse sediment on the polygon
               double const dCoarseErosionTarget = tMin(dCoarsePotentialErosion, dExistingUnconsCoarse);

               // OK, do the supply-limited erosion of coarse sediment
               nRet = nDoUnconsErosionOnPolygon(nCoast, pPolygon, TEXTURE_COARSE, dCoarseErosionTarget, dCoarseEroded);

               if (nRet != RTN_OK)
                  return nRet;

               if (dCoarseEroded > 0)
               {
                  // We eroded some coarse sediment, so add to the this-iteration toal
                  m_dThisIterBeachErosionCoarse += dCoarseEroded;

                  // Store the amount eroded for this polygon
                  pPolygon->SetBeachErosionUnconsCoarse(-dCoarseEroded);
               }
            }

            // OK we now have the actual values of sediment eroded from this polygon, so next determine where this eroded sand and coarse sediment goes (have to consider fine sediment too, because this goes off-grid on grid-edge polygons). Only do this if some sand or coarse was eroded on this polygon
            if ((dSandEroded + dCoarseEroded) > 0)
            {
               if (pPolygon->bDownCoastThisIter())
               {
                  // Moving eroded sediment down-coast
                  int const nNumAdjPoly = pPolygon->nGetNumDownCoastAdjacentPolygons();

                  for (int nn = 0; nn < nNumAdjPoly; nn++)
                  {
                     int const nAdjPoly = pPolygon->nGetDownCoastAdjacentPolygon(nn);

                     // if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
                     // LogStream << m_ulIter << ": polygon " << nPoly << " moves sediment down-coast to polygon " << nAdjPoly << endl;

                     if (nAdjPoly == INT_NODATA)
                     {
                        // Sediment is leaving the grid
                        if (pPolygon->bIsCoastStartPolygon())
                        {
                           // Error: uncons sediment movement is down-coast and off-edge, but this polygon is at the up-coast end of the coastline
                           if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
                              LogStream << m_ulIter << ": " << ERR << "in sediment export. Unconsolidated sediment movement is DOWN-COAST, and sediment is leaving the grid, but polygon " << nPolygon << " is at the up-coast end of the coastline. This will result in mass balance problems." << endl;
                        }
                        else if (pPolygon->bIsCoastEndPolygon())
                        {
                           // This is the polygon at the down-coast end of the coastline, and uncons sediment movement is down-coast. Decide what to do based on the user setting m_nUnconsSedimentHandlingAtGridEdges
                           if (m_nUnconsSedimentHandlingAtGridEdges == GRID_EDGE_CLOSED)
                           {
                              // Closed grid edges: no uncons sediment moves off-grid, nothing is removed from this polygon, so cannot adjust sediment export
                              if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
                                 LogStream << m_ulIter << ": when adjusting sediment export, polygon " << nPolygon << " is at the down-coast end of the coastline, and actual sediment movement is DOWN-COAST. Since grid edges are closed, no sand or coarse unconsolidated sediment goes off-grid so cannot adjust sediment export. This will result in mass balance problems." << endl;
                           }
                           else if (m_nUnconsSedimentHandlingAtGridEdges == GRID_EDGE_OPEN)
                           {
                              // Open grid edges, so this sediment goes off-grid
                              m_dThisIterLeftGridUnconsSand += dSandEroded;
                              m_dThisIterLeftGridUnconsCoarse += dCoarseEroded;
                           }
                           else if (m_nUnconsSedimentHandlingAtGridEdges == GRID_EDGE_RECIRCULATE)
                           {
                              // Re-circulating grid edges, so adjust the sediment exported to the polygon at the up-coast end of this coastline
                              int const nOtherEndPoly = 0;
                              CGeomCoastPolygon* pOtherEndPoly = m_VCoast[nCoast].pGetPolygon(nOtherEndPoly);

                              if (dSandEroded > 0)
                              {
                                 // Add to the still-to-do total of unconsolidated sand to be deposited on the polygon at the up-coast end of this coastline
                                 pOtherEndPoly->AddToDoBeachDepositionUnconsSand(dSandEroded);
                              }

                              if (dCoarseEroded > 0)
                              {
                                 // Add to the still-to-do total of unconsolidated coarse sediment to be deposited on the polygon at the up-coast end of this coastline
                                 pOtherEndPoly->AddToDoBeachDepositionUnconsCoarse(dCoarseEroded);
                              }
                           }
                        }
                     }
                     else
                     {
                        // This polygon is not at the grid edge
                        CGeomCoastPolygon* pAdjPolygon = m_VCoast[nCoast].pGetPolygon(nAdjPoly);
                        double const dBoundaryShare = pPolygon->dGetDownCoastAdjacentPolygonBoundaryShare(nn);

                        if (dSandEroded > 0)
                        {
                           // if (m_ulIter == 5)
                           // LogStream << m_ulIter << ": B polygon not at grid edge nPoly = " << nPoly << " nAdjPoly = " << nAdjPoly << ", dGetToDoBeachDepositionUnconsSand() on nAdjPoly is = " << pAdjPolygon->dGetToDoBeachDepositionUnconsSand() * m_dCellArea << endl;

                           // Add to the still-to-do total of unconsolidated sand to be deposited on the adjacent polygon
                           pAdjPolygon->AddToDoBeachDepositionUnconsSand(dSandEroded * dBoundaryShare);

                           // if (m_ulIter == 5)
                           // LogStream << m_ulIter << ": B after nAdjPoly = " << nAdjPoly << " AddToDoBeachDepositionUnconsSand(" << dSandEroded * dBoundaryShare * m_dCellArea << ") dGetToDoBeachDepositionUnconsSand() now = " << pAdjPolygon->dGetToDoBeachDepositionUnconsSand() * m_dCellArea << endl;
                        }

                        if (dCoarseEroded > 0)
                        {
                           // if (m_nLogFileDetail >= LOG_FILE_ALL)
                           // LogStream << m_ulIter << ": polygon = " << nPoly << " adjacent polygon = " << nAdjPoly << ", beach deposition 1 of uncons coarse was = " << pAdjPolygon->dGetToDoBeachDepositionUnconsCoarse() * m_dCellArea;

                           // Add to the still-to-do total of unconsolidated coarse sediment to be deposited on the adjacent polygon
                           pAdjPolygon->AddToDoBeachDepositionUnconsCoarse(dCoarseEroded * dBoundaryShare);

                           // if (m_nLogFileDetail >= LOG_FILE_ALL)
                           // LogStream << " after AddToDoBeachDepositionUnconsCoarse(" << dCoarseEroded * dBoundaryShare << ") beach deposition 1 of uncons coarse now = " << pAdjPolygon->dGetToDoBeachDepositionUnconsCoarse() * m_dCellArea << endl;
                        }
                     }
                  }

                  // if (m_nLogFileDetail >= LOG_FILE_ALL)
                  // LogStream << m_ulIter << ": 1 uncons sand eroded = " << dSandEroded * m_dCellArea << " 1 uncons coarse eroded = " << dCoarseEroded * m_dCellArea << endl;
               }
               else
               {
                  // Moving eroded sediment up-coast
                  int const nNumAdjPoly = pPolygon->nGetNumUpCoastAdjacentPolygons();

                  for (int nn = 0; nn < nNumAdjPoly; nn++)
                  {
                     int const nAdjPoly = pPolygon->nGetUpCoastAdjacentPolygon(nn);
                     // if (m_nLogFileDetail >= LOG_FILE_ALL)
                     // LogStream << m_ulIter << ": polygon " << nPoly << " moves sediment up-coast to polygon " << nAdjPoly << endl;

                     if (nAdjPoly == INT_NODATA)
                     {
                        // Sediment is leaving the grid
                        if (pPolygon->bIsCoastEndPolygon())
                        {
                           // Error: uncons sediment movement is down-coast and off-edge, but this polygon is at the up-coast end of the coastline
                           if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
                              LogStream << m_ulIter << ": " << ERR << "in sediment export. Unconsolidated sediment movement is UP-COAST, and sediment is leaving the grid, but polygon " << nPolygon << " is at the down-coast end of the coastline. This will result in mass balance problems." << endl;
                        }
                        else if (pPolygon->bIsCoastStartPolygon())
                        {
                           // This is the polygon at the up-coast end of the coastline, and uncons sediment movement is up-coast. Decide what to do based on the user setting m_nUnconsSedimentHandlingAtGridEdges
                           if (m_nUnconsSedimentHandlingAtGridEdges == GRID_EDGE_CLOSED)
                           {
                              // Closed grid edges: no uncons sediment moves off-grid, nothing is removed from this polygon, so cannot adjust sediment export
                              if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
                                 LogStream << m_ulIter << ": when adjusting sediment export, polygon " << nPolygon << " is at the up-coast end of the coastline, and actual sediment movement is UP-COAST. Since grid edges are closed, no sand or coarse unconsolidated sediment goes off-grid so cannot adjust sediment export" << endl;
                           }
                           else if (m_nUnconsSedimentHandlingAtGridEdges == GRID_EDGE_OPEN)
                           {
                              // Open grid edges, so this sediment goes off-grid
                              m_dThisIterLeftGridUnconsSand += dSandEroded;
                              m_dThisIterLeftGridUnconsCoarse += dCoarseEroded;
                           }
                           else if (m_nUnconsSedimentHandlingAtGridEdges == GRID_EDGE_RECIRCULATE)
                           {
                              // Re-circulating grid edges, so adjust the sediment exported to the polygon at the up-coast end of this coastline TODO 016 Check whether this causes mass balance problems, depending on the sequence of polygon processing
                              int const nOtherEndPoly = 0;
                              CGeomCoastPolygon* pOtherEndPoly = m_VCoast[nCoast].pGetPolygon(nOtherEndPoly);

                              if (dSandEroded > 0)
                              {
                                 // Add to the still-to-do total of unconsolidated sand to be deposited on the polygon at the up-coast end of this coastline
                                 pOtherEndPoly->AddToDoBeachDepositionUnconsSand(dSandEroded);
                              }

                              if (dCoarseEroded > 0)
                              {
                                 // Add to the still-to-do total of unconsolidated coarse sediment to be deposited on the polygon at the up-coast end of this coastline
                                 pOtherEndPoly->AddToDoBeachDepositionUnconsCoarse(dCoarseEroded);
                              }
                           }
                        }
                     }
                     else
                     {
                        // This polygon is not at the grid edge
                        CGeomCoastPolygon* pAdjPolygon = m_VCoast[nCoast].pGetPolygon(nAdjPoly);
                        double const dBoundaryShare = pPolygon->dGetUpCoastAdjacentPolygonBoundaryShare(nn);

                        if (dSandEroded > 0)
                        {
                           // if (m_ulIter == 5)
                           // LogStream << m_ulIter << ": A polygon not at grid edge nPoly = " << nPoly << " nAdjPoly = " << nAdjPoly << ", dGetToDoBeachDepositionUnconsSand() on nAdjPoly is = " << pAdjPolygon->dGetToDoBeachDepositionUnconsSand() * m_dCellArea << endl;

                           // Add to the still-to-do total of unconsolidated sand to be deposited on the the adjacent polygon
                           pAdjPolygon->AddToDoBeachDepositionUnconsSand(dSandEroded * dBoundaryShare);

                           // if (m_ulIter == 5)
                           // LogStream << m_ulIter << ": A after nAdjPoly = " << nAdjPoly << " AddToDoBeachDepositionUnconsSand(" << dSandEroded * dBoundaryShare * m_dCellArea << ") dGetToDoBeachDepositionUnconsSand() now = " << pAdjPolygon->dGetToDoBeachDepositionUnconsSand() * m_dCellArea << endl;
                        }

                        if (dCoarseEroded > 0)
                        {
                           // if (m_ulIter == 5)
                           // LogStream << m_ulIter << ": polygon not at grid edge nPoly = " << nPoly << " nAdjPoly = " << nAdjPoly << ", beach deposition of uncons coarse was = " << pAdjPolygon->dGetToDoBeachDepositionUnconsCoarse() * m_dCellArea << endl;

                           // Add to the still-to-do total of unconsolidated coarse sediment to be deposited on the adjacent polygon
                           pAdjPolygon->AddToDoBeachDepositionUnconsCoarse(+dCoarseEroded * dBoundaryShare);

                           // if (m_ulIter == 5)
                           // LogStream << " After AddToDoBeachDepositionUnconsCoarse(" << dCoarseEroded * dBoundaryShare << ") uncons coarse now = " << pAdjPolygon->dGetToDoBeachDepositionUnconsCoarse() * m_dCellArea << endl;
                        }
                     }
                  }
               }
            }

            // if (m_nLogFileDetail >= LOG_FILE_ALL)
            // LogStream << m_ulIter << ": sand eroded on poly = " << dSandEroded * m_dCellArea << " coarse eroded on poly = " << dCoarseEroded * m_dCellArea << endl;

         } // if (dPotentialErosion > 0)

         // // DEBUG CODE ============================================================================================================================================
         // // Get total depths of unconsolidated sand for every cell
         // if (m_ulIter == 5)
         // {
         // double dTmpSandUncons = 0;
         // for (int nX1 = 0; nX1 < m_nXGridSize; nX1++)
         // {
         // for (int nY1 = 0; nY1 < m_nYGridSize; nY1++)
         // {
         // dTmpSandUncons += m_pRasterGrid->m_Cell[nX1][nY1].dGetUnconsSandDepthAllLayers();
         // }
         // }
         //
         // LogStream << endl;
         // LogStream << m_ulIter << ": TOTAL UNCONSOLIDATED SAND ON ALL CELLS after beach movement on nPoly = " << nPoly << " dTmpSandUncons = " << dTmpSandUncons * m_dCellArea << " (dTmpSandUncons - m_dStartIterUnconsSandAllCells) =  " << (dTmpSandUncons - m_dStartIterUnconsSandAllCells) * m_dCellArea << endl;
         //
         //    // Now get the total in all still-to-do fields of all polygons
         // double dToDoTot = 0;
         //
         // for (int nn = 0; nn < nPolygons; nn++)
         // {
         // CGeomCoastPolygon* pThisPolygon = m_VCoast[nCoast].pGetPolygon(nn);
         // dToDoTot += pThisPolygon->dGetToDoBeachDepositionUnconsSand();
         // }
         //
         // LogStream << endl  << m_ulIter << ": dToDoTot = " << dToDoTot * m_dCellArea << endl;
         // LogStream  << m_ulIter << ": dTmpSandUncons + dToDoTot = " << (dTmpSandUncons + dToDoTot) * m_dCellArea << endl << endl;
         //
         // LogStream << m_ulIter << ": m_dThisIterLeftGridUnconsSand = " << m_dThisIterLeftGridUnconsSand * m_dCellArea << endl;
         //
         // LogStream << "*****************************" << endl;
         // }
         // // DEBUG CODE ============================================================================================================================================

      } // for (int n = 0; n < nNumPolygons; n++)

      // OK we have processed all polygons, But if there are adjacent-polygon circularities (i.e. Polygon A -> Polygon B -> Polygon A) then we may have some still-to-do deposition on at least one polygon. So look through all polygons and check their still-to-do lists
      int const nPolygons = m_VCoast[nCoast].nGetNumPolygons();

      // for (int nn = 0; nn < nPolygons; nn++)
      for (int nn = nPolygons - 1; nn >= 0; nn--)
      {
         int const nThisPoly = nVVVAllCoastPolyAndAdjacent[nCoast][nn][0];
         CGeomCoastPolygon* pThisPolygon = m_VCoast[nCoast].pGetPolygon(nThisPoly);

         double const dSandToDepositOnPoly = pThisPolygon->dGetToDoBeachDepositionUnconsSand();

         if (dSandToDepositOnPoly > 0)
         {
            // There is some still-to-do deposition of sand sediment on this polygon: calculate a net increase in depth of sand-sized unconsolidated sediment on the cells within the polygon. Note that some cells may decrease in elevation (i.e. have some sand-sized sediment erosion) however
            double dSandDeposited = 0;
            nRet = nDoUnconsDepositionOnPolygon(nCoast, pThisPolygon, TEXTURE_SAND, dSandToDepositOnPoly, dSandDeposited);

            if (nRet != RTN_OK)
               return nRet;

            double const dSandNotDeposited = dSandToDepositOnPoly - dSandDeposited;

            if (dSandNotDeposited > 0)
               m_dDepositionSandDiff += dSandNotDeposited;

            if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
               LogStream << m_ulIter << ":\t coast " << nCoast << " re-processing nThisPoly = " << nThisPoly << " dSandDeposited = " << dSandDeposited * m_dCellArea << " dSandNotDeposited = " << dSandNotDeposited * m_dCellArea << " m_dDepositionSandDiff = " << m_dDepositionSandDiff * m_dCellArea << endl;
         }

         double const dCoarseToDepositOnPoly = pThisPolygon->dGetToDoBeachDepositionUnconsCoarse();

         if (dCoarseToDepositOnPoly > 0)
         {
            // There is some still-to-do deposition of coarse sediment on this polygon: calculate a net increase in depth of coarse-sized unconsolidated sediment on the cells within the polygon. Note that some cells may decrease in elevation (i.e. have some coarse-sized sediment erosion) however
            double dCoarseDeposited = 0;
            nRet = nDoUnconsDepositionOnPolygon(nCoast, pThisPolygon, TEXTURE_COARSE, dCoarseToDepositOnPoly, dCoarseDeposited);

            if (nRet != RTN_OK)
               return nRet;

            double const dCoarseNotDeposited = dCoarseToDepositOnPoly - dCoarseDeposited;

            if (dCoarseNotDeposited > 0)
               m_dDepositionCoarseDiff += dCoarseNotDeposited;

            if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
               LogStream << m_ulIter << ":\t coast " << nCoast << " re-processing nThisPoly = " << nThisPoly << " dCoarseDeposited = " << dCoarseDeposited * m_dCellArea << " dCoarseNotDeposited = " << dCoarseNotDeposited * m_dCellArea << " m_dDepositionCoarseDiff = " << m_dDepositionCoarseDiff * m_dCellArea << endl;
         }
      }
   }

   if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
   {
      // Write out actual sediment movement for all coasts and all polygons
      WritePolygonActualMovement(nVVVAllCoastPolyAndAdjacent);
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Before simulating beach erosion, update the per-polygon values of pre-existing unconsolidated sediment for all three size classes, to include unconsolidated sediment derived from platform erosion, cliff collapse, and sediment input events
//===============================================================================================================================
void CSimulation::AllPolygonsUpdateStoredUncons(int const nCoast)
{
   int const nNumPolygons = m_VCoast[nCoast].nGetNumPolygons();

   for (int nPoly = 0; nPoly < nNumPolygons; nPoly++)
   {
      CGeomCoastPolygon* pPolygon = m_VCoast[nCoast].pGetPolygon(nPoly);

      // Only include unconsolidated fine sediment from sediment input events, don't include any unconsolidated fine sediment from platform erosion and cliff collapse since these have gone to suspension
      double const dFine = pPolygon->dGetSedimentInputUnconsFine();
      pPolygon->SetPreExistingUnconsFine(dFine);

      // Include unconsolidated sand sediment derived from platform erosion, cliff collapse, and sediment input events
      double const dSand = pPolygon->dGetPreExistingUnconsSand() + pPolygon->dGetPlatformErosionUnconsSand() + pPolygon->dGetCliffCollapseUnconsSandDeposition() + pPolygon->dGetSedimentInputUnconsSand();
      pPolygon->SetPreExistingUnconsSand(dSand);

      // Include unconsolidated coarse sediment derived from platform erosion, cliff collapse, and sediment input events
      double const dCoarse = pPolygon->dGetPreExistingUnconsCoarse() + pPolygon->dGetPlatformErosionUnconsCoarse() + pPolygon->dGetCliffCollapseUnconsCoarseDeposition() + pPolygon->dGetSedimentInputUnconsCoarse();
      pPolygon->SetPreExistingUnconsCoarse(dCoarse);
   }
}
