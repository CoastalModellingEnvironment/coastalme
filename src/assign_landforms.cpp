/*!
   \file assign_landforms.cpp
   \brief Assigns landform categories to coastlines and coastal cells, and to all other dryland cells
   \details TODO 001 A more detailed description of these routines.
   \author David Favis-Mortlock
   \author Andres Payo
   \author Wilf Chun
   \date 2026
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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "cme.h"
#include "simulation.h"
#include "coast.h"
#include "cliff.h"
#include "drift.h"
#include "intervention.h"
#include "cell_layer.h"

//===============================================================================================================================
//! Each timestep, classify landforms for all cells (cells on the coastline will be changed later)
//===============================================================================================================================
int CSimulation::nAssignLandformsForAllCells(void)
{
   int const NO_ACTION = -1;

   // First pass: collect information about cells that need to be changed. This avoids race conditions from reading neighbour cells while writing to current cells
   vector<vector<int>> VCellToUpdate(m_nXGridSize, vector<int>(m_nYGridSize, NO_ACTION));

   // Read-only phase: determine what changes need to be made
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif

   for (int nX = 0; nX < m_nXGridSize; nX++)
   {
      for (int nY = 0; nY < m_nYGridSize; nY++)
      {
         // Get this cell's existing landform category
         CRWCellLandform const* pLandform = m_pRasterGrid->m_Cell[nX][nY].pGetLandform();
         int const nCat = pLandform->nGetLandformCategory();

         // OK, use these rules to set landform categories
         if (nCat == LF_INTERVENTION_STRUCT)
            // Already set to intervention, so don't change
            continue;

         if (nCat == LF_INTERVENTION_NON_STRUCT)
            // Already set to intervention, so don't change
            continue;

         if (nCat == LF_SEDIMENT_INPUT_UNCONSOLIDATED)
            // Already set to sediment input, so don't change
            continue;

         if (nCat == LF_SEDIMENT_INPUT_CONSOLIDATED)
            // Already set to sediment input, so don't change
            continue;

         // Now maybe change landform category
         if (m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea())
         {
            // This is a sea cell, is it surrounded by drift cells?
            if (bSurroundedByDriftCells(nX, nY))
            {
               // Set to beach
               VCellToUpdate[nX][nY] = LF_DRIFT_BEACH;
               continue;
            }

            if (m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevIncTalus() > m_dThisIterSWL)
            {
               // This sea cell has elevation is above this-iteration SWL
               if (nCat != LF_ISLAND)
                  // Set to island
                  VCellToUpdate[nX][nY] = LF_ISLAND;
               continue;
            }

            // Is this cell already marked as sea?
            if (nCat != LF_SEA)
               // Mark it as sea
               VCellToUpdate[nX][nY] = LF_SEA;
            continue;
         }

         // OK, this is not sea. Are we down to basement?
         if (m_pRasterGrid->m_Cell[nX][nY].bBasementElevIsMissingValue())
         {
            // Down to basement
            if (nCat != LF_UNKNOWN)
               // Set to unknown landform
               VCellToUpdate[nX][nY] = LF_UNKNOWN;
            continue;
         }

         int const nTopLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopNonZeroLayerAboveBasement();
         CRWCellLayer* pTopLayer = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer);

         if (pTopLayer->bHasTalus())
         {
            // There is talus here
            VCellToUpdate[nX][nY] = LF_DRIFT_TALUS;
            continue;
         }

         if (pTopLayer->bHasUncons())
         {
            // This is unconsolidated sediment here, so set to beach
            VCellToUpdate[nX][nY] = LF_DRIFT_BEACH;
            continue;
         }

         // Default
         VCellToUpdate[nX][nY] = LF_HINTERLAND;
      }
   }

   // Write phase: apply the changes
   for (int nX = 0; nX < m_nXGridSize; nX++)
   {
      for (int nY = 0; nY < m_nYGridSize; nY++)
      {
         int const nAction = VCellToUpdate[nX][nY];

         if (nAction == NO_ACTION)
            // No change
            continue;

         CRWCellLandform* pLandform = m_pRasterGrid->m_Cell[nX][nY].pGetLandform();

         switch (nAction)
         {
         case LF_UNKNOWN:
            // Set to unknown landform
            pLandform->SetLandformCategory(LF_UNKNOWN);
            break;

         case LF_HINTERLAND:
            // Set to hinterland
            pLandform->SetLandformCategory(LF_HINTERLAND);
            break;

         case LF_CLIFF:
            // Set to cliff
            pLandform->SetLandformCategory(LF_CLIFF);
            break;

         case LF_DRIFT_TALUS:
            // Set to talus
            pLandform->SetLandformCategory(LF_DRIFT_TALUS);
            break;

         case LF_DRIFT_BEACH:
            // Set to beach
            pLandform->SetLandformCategory(LF_DRIFT_BEACH);
            break;

         case LF_DRIFT_DUNES:
            // TODO not yet implemented
            // Set to dunes
            pLandform->SetLandformCategory(LF_DRIFT_DUNES);
            break;

         case LF_INTERVENTION_STRUCT:
            // Set to structural intervention
            pLandform->SetLandformCategory(LF_INTERVENTION_STRUCT);
            break;

         case LF_INTERVENTION_NON_STRUCT:
            // Set to non-structural intervention
            pLandform->SetLandformCategory(LF_INTERVENTION_NON_STRUCT);
            break;

         case LF_ISLAND:
            // Set to island
            pLandform->SetLandformCategory(LF_ISLAND);
            break;

         case LF_SEDIMENT_INPUT_UNCONSOLIDATED:
            // Set to unconsolidated sediment input
            pLandform->SetLandformCategory(LF_SEDIMENT_INPUT_UNCONSOLIDATED);
            break;

         case LF_SEDIMENT_INPUT_CONSOLIDATED:
            // Set to consolidated sediment input
            pLandform->SetLandformCategory(LF_SEDIMENT_INPUT_CONSOLIDATED);
            break;
         }
      }
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Each timestep, classify coastal landforms and assign a coastal landform object to every point on every coastline. If, for a given cell, the coastal landform class has not changed then it inherits values from the previous timestep
//===============================================================================================================================
int CSimulation::nAssignLandformsForAllCoasts(void)
{
   // For each coastline, put a coastal landform at every point along the coastline
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      for (int nCoastPoint = 0; nCoastPoint < m_VCoast[nCoast].nGetCoastlineSize(); nCoastPoint++)
      {
         // Get the coords of the grid cell marked as coastline for the coastal landform object
         int const nX = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint)->nGetX();
         int const nY = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint)->nGetY();

         // Store the coastline number and the number of the coastline point in the cell so we can get these quickly later
         m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetPointOnCoast(nCoastPoint);

         // OK, start assigning coastal landforms. First, is there an intervention on this cell?
         int nCat = m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLandformCategory();
         if ((nCat == LF_INTERVENTION_STRUCT) || (nCat == LF_INTERVENTION_NON_STRUCT))
         {
            // There is, so create an intervention object on the vector coastline with these attributes
            CACoastLandform* pIntervention = new CRWIntervention(&m_VCoast[nCoast], nCoast, nCoastPoint, nCat);
            m_VCoast[nCoast].AppendCoastLandform(pIntervention);

            // LogStream << nCoastPoint << " [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "} " << m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLandformCategory() << " " << m_pRasterGrid->m_Cell[nX][nY].dGetInterventionHeight() << endl;

            continue;
         }

         // OK the landform on this coast cell is something other than an intervention. First check for talus
         int const nTopLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopNonZeroLayerAboveBasement();
         CRWCellLayer* pTopLayer = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer);

         if (pTopLayer->bHasTalus())
         {
            // There is talus on this cell
            m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetLandformCategory(LF_DRIFT_TALUS);

            CACoastLandform* pDrift = new CRWDrift(&m_VCoast[nCoast], nCoast, nCoastPoint, LF_DRIFT_TALUS);
            m_VCoast[nCoast].AppendCoastLandform(pDrift);

            continue;
         }

         // Next, do some safety checks. Note that layer 0 is the first layer above basement
         int const nLayer = m_pRasterGrid->m_Cell[nX][nY].nGetLayerAtElev(m_dThisIterSWL);
         if (nLayer == ELEV_IN_BASEMENT)
         {
            // Should never happen
            LogStream << m_ulIter << ": SWL (" << m_dThisIterSWL << ") is in basement on cell [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}, cannot assign coastal landform for coastline " << nCoast << endl;

            return RTN_ERR_CANNOT_ASSIGN_COASTAL_LANDFORM;
         }

         if (nLayer == ELEV_ABOVE_SEDIMENT_TOP)
         {
            // Again, should never happen
            LogStream << m_ulIter << ":\t SWL (" << m_dThisIterSWL << ") is above sediment-top elevation inc. any talus (" << m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevIncTalus() << ") on cell [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}, cannot assign coastal landform for coastline " << nCoast << endl;

            // TODO DFM bodge ========================
            // We have unconsolidated sediment at SWL, so this is a drift cell: create a drift object on the vector coastline with these attributes
            CACoastLandform* pDrift = new CRWDrift(&m_VCoast[nCoast], nCoast, nCoastPoint, LF_DRIFT_BEACH);
            m_VCoast[nCoast].AppendCoastLandform(pDrift);

            m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetLandformCategory(LF_DRIFT_BEACH);
            continue;
            // TODO DFM bodge ========================

            // return RTN_ERR_CANNOT_ASSIGN_COASTAL_LANDFORM;
         }

         // OK, now check what we have at SWL on this cell: is it unconsolidated or consolidated sediment?
         double const dConsSedTop = m_pRasterGrid->m_Cell[nX][nY].dGetConsSedTopElevForLayerAboveBasement(nLayer);
         if (dConsSedTop >= m_dThisIterSWL)
         {
            // We have consolidated sediment at or above SWL on this cell. Are we considering cliff collapse?
            if (m_bDoCliffCollapse)
            {
               // OK we are considering cliff collapse, and we have consolidated sediment at SWL, so this is a cliff cell. Get the existing landform category for this cell
               nCat = m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLandformCategory();
               if (nCat == LF_CLIFF)
               {
                  // This cell was a cliff in some previous timestep. Is the pre-existing notch still below the top of the consolidated sediment?
                  double dNotchApexElev = m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->dGetCliffNotchApexElev();
                  double const dSedTopElevNoTalus = m_pRasterGrid->m_Cell[nX][nY].dGetConsSedTopElevOmitTalus();
                  if (dNotchApexElev < dSedTopElevNoTalus)
                  {
                     // Yes, the notch is still below the top of the consolidated sediment, so get the pre-existing data stored in the cell
                     double const dAccumWaveEnergy = m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->dGetAccumWaveEnergy();
                     double const dNotchIncision = m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->dGetCliffNotchIncisionDepth();

                     // Set this as a cliff cell on the coastline
                     m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetLandformCategory(LF_CLIFF);

                     // Create a cliff object on the vector coastline with these attributes
                     CACoastLandform* pCliff = new CRWCliff(&m_VCoast[nCoast], nCoast, nCoastPoint, m_dCellSide, dNotchIncision, dNotchApexElev, dAccumWaveEnergy);
                     m_VCoast[nCoast].AppendCoastLandform(pCliff);

                     if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
                        LogStream << m_ulIter << ":\t continues to be a cliff at [" << nX << "][" << nY << "] dAccumWaveEnergy = " << dAccumWaveEnergy << " dSedTopElevNoTalus = " << dSedTopElevNoTalus << " dSedTopElevIncTalus = " << m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevIncTalus() << " dNotchApexElev = " << dNotchApexElev << " dNotchIncision = " << dNotchIncision << endl;
                  }
                  else
                  {
                     // This was a cliff in the previous timestep, but the notch is no longer below the top of the consolidated sediment. Create a cliff object on the vector coastline without a notch
                     double const dAccumWaveEnergy = m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->dGetAccumWaveEnergy();
                     double const dNotchIncision = DBL_NODATA;
                     dNotchApexElev = DBL_NODATA;

                     // This is a cliff cell on the coastline without a notch
                     m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetLandformCategory(LF_CLIFF);
                     m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetCliffNotchApexElev(dNotchApexElev);
                     m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetCliffNotchIncisionDepth(dNotchIncision);

                     // Create a cliff object on the vector coastline with these attributes
                     CACoastLandform* pCliff = new CRWCliff(&m_VCoast[nCoast], nCoast, nCoastPoint, m_dCellSide, dNotchIncision, dNotchApexElev, dAccumWaveEnergy);
                     m_VCoast[nCoast].AppendCoastLandform(pCliff);

                     double const dSedTopElevIncTalus = m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevIncTalus();

                     if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
                        LogStream << m_ulIter << ":\t PROBLEM cliff with notch above sediment top (inc any talus) at [" << nX << "][" << nY << "] dAccumWaveEnergy = " << dAccumWaveEnergy << " dNotchApexElev = " << dNotchApexElev << " dSedTopElevNoTalus = " << dSedTopElevNoTalus << " dSedTopElevIncTalus = " << dSedTopElevIncTalus << " dNotchIncision = " << dNotchIncision << endl;
                  }
               }
               else
               {
                  // This was not a cliff in the previous timestep, but it is now
                  m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetLandformCategory(LF_CLIFF);

                  // Get the pre-existing wave energy stored in the cell
                  double const dAccumWaveEnergy = m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->dGetAccumWaveEnergy();

                  // The DBL_NODATA Values indicate that the cliff object does not have an incised notch
                  double dNotchIncision = DBL_NODATA;
                  double dNotchApexElev = DBL_NODATA;

                  // Would the new notch apex elevation be below the top of the cell's consolidated sediment?
                  double const dSedTopElevNoTalus = m_pRasterGrid->m_Cell[nX][nY].dGetConsSedTopElevOmitTalus();
                  if (m_dThisIterNewNotchApexElev < dSedTopElevNoTalus)
                  {
                     // Yes it would, so this new cliff object has a notch
                     dNotchIncision = 0;
                     dNotchApexElev = m_dThisIterNewNotchApexElev - SED_ELEV_TOLERANCE;
                  }

                  // Create a cliff object on the vector coastline with these attributes
                  CACoastLandform* pCliff = new CRWCliff(&m_VCoast[nCoast], nCoast, nCoastPoint, m_dCellSide, dNotchIncision, dNotchApexElev, dAccumWaveEnergy);
                  m_VCoast[nCoast].AppendCoastLandform(pCliff);

                  // if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
                  // {
                  //    if (bFPIsEqual(dNotchIncision, 0.0, TOLERANCE))
                  //       LogStream << m_ulIter << ":\t coastline cliff created at [" << nX << "][" << nY << "] dAccumWaveEnergy = " << dAccumWaveEnergy << " dSedTopElevNoTalus = " << dSedTopElevNoTalus << " dSedTopElevIncTalus = " << m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevIncTalus() << " dNotchApexElev = " << dNotchApexElev << " dNotchIncision = " << dNotchIncision << endl;
                  //    else
                  //       LogStream << m_ulIter << ":\t coastline no-notch cliff created at [" << nX << "][" << nY << "] dAccumWaveEnergy = " << dAccumWaveEnergy << " dSedTopElevNoTalus = " << dSedTopElevNoTalus << " dSedTopElevIncTalus = " << m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevIncTalus() << " dNotchApexElev = " << dNotchApexElev << " dNotchIncision = " << dNotchIncision << endl;
                  // }
               }
            }
            else
            {
               // We have consolidated sediment at SWL but we are not considering cliff collapse. Get the pre-existing wave energy stored in the cell
               double const dAccumWaveEnergy = m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->dGetAccumWaveEnergy();

               // Create a cliff object on the vector coastline
               CACoastLandform* pCliff = new CRWCliff(&m_VCoast[nCoast], nCoast, nCoastPoint, m_dCellSide, DBL_NODATA, DBL_NODATA, dAccumWaveEnergy);
               m_VCoast[nCoast].AppendCoastLandform(pCliff);

               if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
                  LogStream << m_ulIter << ":\t coastline cliff created (cliff collapse not considered) at " << nX << "][" << nY << "] dAccumWaveEnergy = " << dAccumWaveEnergy << " dAllSedTopElevNoTalus = " << m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevOmitTalus() << " dAllSedTopElevIncTalus = " << m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevIncTalus() << " dConsSedTopElevNoTalus = " << m_pRasterGrid->m_Cell[nX][nY].dGetConsSedTopElevOmitTalus() << " dConsSedTopElevIncTalus = " << m_pRasterGrid->m_Cell[nX][nY].dGetConsSedTopElevIncTalus() << " dNotchApexElev = " << DBL_NODATA << " dNotchIncision = " << DBL_NODATA << endl;
            }
         }
         else
         {
            // We have unconsolidated sediment at SWL, so this is a drift cell: create a drift object on the vector coastline with these attributes
            CACoastLandform* pDrift = new CRWDrift(&m_VCoast[nCoast], nCoast, nCoastPoint, LF_DRIFT_BEACH);
            m_VCoast[nCoast].AppendCoastLandform(pDrift);

            m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetLandformCategory(LF_DRIFT_BEACH);

// #ifdef _DEBUG
//             LogStream << m_ulIter << ":\t drift created at [" << nX << "][" << nY << "]" << endl;
// #endif
         }
      }
   }

   // // DEBUG CODE ============================================================================================================================================
   // for (int i = 0; i < static_cast<int>(m_VCoast.size()); i++)
   // {
   //    for (int j = 0; j < m_VCoast[i].nGetCoastlineSize(); j++)
   //    {
   //       int nX = m_VCoast[i].pPtiGetCellMarkedAsCoastline(j)->nGetX();
   //       int nY = m_VCoast[i].pPtiGetCellMarkedAsCoastline(j)->nGetY();
   //
   //       LogStream << m_ulIter << ": coast cell " << j << " at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "} has landform category = ";
   //
   //       int nCat = m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLandformCategory();
   //       int nSubCat = m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFSubCategory();
   //
   //       switch (nCat)
   //       {
   //          case LF_HINTERLAND:
   //             LogStream << "hinterland";
   //             break;
   //
   //          case LF_SEA:
   //             LogStream << "sea";
   //             break;
   //
   //          case LF_CLIFF:
   //             LogStream << "cliff";
   //             break;
   //
   //          case LF_DRIFT:
   //             LogStream << "drift";
   //             break;
   //
   //          case LF_INTERVENTION:
   //             LogStream << "intervention";
   //             break;
   //
   //          case LF_UNKNOWN:
   //             LogStream << "none";
   //             break;
   //
   //          default:
   //             LogStream << "NONE";
   //             break;
   //       }
   //
   //       LogStream << " and landform subcategory = ";
   //
   //       switch (nSubCat)
   //       {
   //          case LF_CLIFF:
   //             LogStream << "cliff";
   //             break;
   //
   //          case LF_DRIFT_TALUS:
   //             LogStream << "talus";
   //             break;
   //
   //          case LF_DRIFT_BEACH:
   //             LogStream << "beach";
   //             break;
   //
   //          case LF_DRIFT_DUNES:
   //             LogStream << "dunes";
   //             break;
   //
   //          case LF_INTERVENTION_STRUCT:
   //             LogStream << "structural intervention";
   //             break;
   //
   //          case LF_INTERVENTION_NON_STRUCT:
   //             LogStream << "non-structural intervention";
   //             break;
   //
   //          case LF_UNKNOWN:
   //             LogStream << "none";
   //             break;
   //
   //          default:
   //             LogStream << "NONE";
   //             break;
   //       }
   //       LogStream << endl;
   //    }
   // }
   // // DEBUG CODE ============================================================================================================================================

   return RTN_OK;
}

//===============================================================================================================================
//! At the end of each timestep, this routine stores the attributes from a single coastal landform object in the grid cell 'under' the object, ready for the next timestep
//===============================================================================================================================
int CSimulation::nLandformToGrid(int const nCoast, int const nPoint)
{
   // What is the coastal landform here?
   CACoastLandform* pCoastLandform = m_VCoast[nCoast].pGetCoastLandform(nPoint);
   int const nCat = pCoastLandform->nGetLandFormCategory();
   if (nCat == LF_CLIFF)
   {
      // It's a cliff
      CRWCliff const* pCliff = reinterpret_cast<CRWCliff*>(pCoastLandform);

      int const nX = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nPoint)->nGetX();
      int const nY = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nPoint)->nGetY();

      if (! pCliff->bHasCollapsed())
      {
         // The cliff has not collapsed. Get attribute values from the cliff object
         double const dNotchBaseElev = pCliff->dGetNotchApexElev();
         double const dNotchIncision = pCliff->dGetNotchIncision();

         // And store some attribute values in the cliff cell
         m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetLandformCategory(LF_CLIFF);
         m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetCliffNotchApexElev(dNotchBaseElev);
         m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetCliffNotchIncisionDepth(dNotchIncision);
      }
      else
      {
         // The cliff has collapsed: all sediment above the base of the erosional notch is gone from this cliff object via cliff collapse, so this cell is no longer a cliff
//          int const nTopLayer = m_pRasterGrid->m_Cell[nX][nY].nGetNumOfTopLayerAboveBasement();
//
//          // Safety check
//          if (nTopLayer == INT_NODATA)
//             return RTN_ERR_NO_TOP_LAYER;

         // // Update the cell's layer elevations
         // m_pRasterGrid->m_Cell[nX][nY].CalcAllLayerElevsAndD50();
         //
         // // And update the cell's sea depth
         // m_pRasterGrid->m_Cell[nX][nY].SetSeaDepth();
      }

      // Always accumulate wave energy
      m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetAccumWaveEnergy(pCliff->dGetTotAccumWaveEnergy());
   }
   // else if (nCat == LF_DRIFT)
   // {
   //    // It's drift, so calculate D50 TODO 002 Why might we need this?
   // }

   return RTN_OK;
}

//===============================================================================================================================
//! Returns true if this cell has eight drift cells surrounding it
//===============================================================================================================================
bool CSimulation::bSurroundedByDriftCells(int const nX, int const nY)
{
   int nXTmp;
   int nYTmp;
   int nAdjacent = 0;

   // North
   nXTmp = nX;
   nYTmp = nY - 1;

   if (bIsWithinValidGrid(nXTmp, nYTmp))
   {
      CRWCellLandform const* pLandform = m_pRasterGrid->m_Cell[nXTmp][nYTmp].pGetLandform();
      int const nCat = pLandform->nGetLandformCategory();

      if ((nCat == LF_DRIFT_BEACH) || (nCat == LF_DRIFT_TALUS) || (nCat == LF_DRIFT_DUNES) || (nCat == LF_CLIFF))
         nAdjacent++;
   }

   // North-east
   nXTmp = nX + 1;
   nYTmp = nY - 1;

   if (bIsWithinValidGrid(nXTmp, nYTmp))
   {
      CRWCellLandform const* pLandform = m_pRasterGrid->m_Cell[nXTmp][nYTmp].pGetLandform();
      int const nCat = pLandform->nGetLandformCategory();

      if ((nCat == LF_DRIFT_BEACH) || (nCat == LF_DRIFT_TALUS) || (nCat == LF_DRIFT_DUNES) || (nCat == LF_CLIFF))
         nAdjacent++;
   }

   // East
   nXTmp = nX + 1;
   nYTmp = nY;

   if (bIsWithinValidGrid(nXTmp, nYTmp))
   {
      CRWCellLandform const* pLandform = m_pRasterGrid->m_Cell[nXTmp][nYTmp].pGetLandform();
      int const nCat = pLandform->nGetLandformCategory();

      if ((nCat == LF_DRIFT_BEACH) || (nCat == LF_DRIFT_TALUS) || (nCat == LF_DRIFT_DUNES) || (nCat == LF_CLIFF))
         nAdjacent++;
   }

   // South-east
   nXTmp = nX + 1;
   nYTmp = nY + 1;

   if (bIsWithinValidGrid(nXTmp, nYTmp))
   {
      CRWCellLandform const* pLandform = m_pRasterGrid->m_Cell[nXTmp][nYTmp].pGetLandform();
      int const nCat = pLandform->nGetLandformCategory();

      if ((nCat == LF_DRIFT_BEACH) || (nCat == LF_DRIFT_TALUS) || (nCat == LF_DRIFT_DUNES) || (nCat == LF_CLIFF))
         nAdjacent++;
   }

   // South
   nXTmp = nX;
   nYTmp = nY + 1;

   if (bIsWithinValidGrid(nXTmp, nYTmp))
   {
      CRWCellLandform const* pLandform = m_pRasterGrid->m_Cell[nXTmp][nYTmp].pGetLandform();
      int const nCat = pLandform->nGetLandformCategory();

      if ((nCat == LF_DRIFT_BEACH) || (nCat == LF_DRIFT_TALUS) || (nCat == LF_DRIFT_DUNES) || (nCat == LF_CLIFF))
         nAdjacent++;
   }

   // South-west
   nXTmp = nX - 1;
   nYTmp = nY + 1;

   if (bIsWithinValidGrid(nXTmp, nYTmp))
   {
      CRWCellLandform const* pLandform = m_pRasterGrid->m_Cell[nXTmp][nYTmp].pGetLandform();
      int const nCat = pLandform->nGetLandformCategory();

      if ((nCat == LF_DRIFT_BEACH) || (nCat == LF_DRIFT_TALUS) || (nCat == LF_DRIFT_DUNES) || (nCat == LF_CLIFF))
         nAdjacent++;
   }

   // West
   nXTmp = nX - 1;
   nYTmp = nY;

   if (bIsWithinValidGrid(nXTmp, nYTmp))
   {
      CRWCellLandform const* pLandform = m_pRasterGrid->m_Cell[nXTmp][nYTmp].pGetLandform();
      int const nCat = pLandform->nGetLandformCategory();

      if ((nCat == LF_DRIFT_BEACH) || (nCat == LF_DRIFT_TALUS) || (nCat == LF_DRIFT_DUNES) || (nCat == LF_CLIFF))
         nAdjacent++;
   }

   // North-west
   nXTmp = nX - 1;
   nYTmp = nY - 1;

   if (bIsWithinValidGrid(nXTmp, nYTmp))
   {
      CRWCellLandform const* pLandform = m_pRasterGrid->m_Cell[nXTmp][nYTmp].pGetLandform();
      int const nCat = pLandform->nGetLandformCategory();

      if ((nCat == LF_DRIFT_BEACH) || (nCat == LF_DRIFT_TALUS) || (nCat == LF_DRIFT_DUNES) || (nCat == LF_CLIFF))
         nAdjacent++;
   }

   if (nAdjacent == 8)
   {
      // This cell has eight LF_DRIFT neighbours
      return true;
   }

   return false;
}
