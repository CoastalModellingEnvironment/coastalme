/*!
   \file do_surge_flood.cpp
   \brief Does flood/surge stuff (not yet implemented)
   \details TODO 001 A more detailed description of these routines.
   \author Manuel Cobos Budia
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
using std::isnan;

#include <iostream>
using std::endl;
using std::ios;

#include <stack>
using std::stack;

#include "cme.h"
#include "2d_point.h"
#include "i_line.h"
#include "line.h"
#include "simulation.h"
#include "raster_grid.h"
#include "coast.h"
#include "2di_point.h"

// //===============================================================================================================================
// //! Use the sealevel, wave set-up and run-up to evaluate flood hydraulically connected TODO 007 Finish surge and runup stuff
// //===============================================================================================================================
// void CSimulation::FloodFillLand(int const nXStart, int const nYStart)
// {
//    // The flood is at a user-specified location. So get the location from values read from the shapefile
//    long const unsigned int nLocIDs = m_VdFloodLocationX.size();
//    double dDiffTotWaterLevel = 0;
//
//    double dAuxWaterLevelDiff = 0;
//
//    if (nLocIDs == 0)
//    {
//       int pointCounter = 0;
//
//       for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
//       {
//          int const nCoastSize = m_VCoast[nCoast].pLGetCoastlineExtCRS()->nGetSize();
//
//          for (int nCoastPoint = 0; nCoastPoint < nCoastSize; nCoastPoint++)
//          {
//             dAuxWaterLevelDiff = m_VCoast[nCoast].dGetLevel(nCoastPoint, m_nLevel);
//
//             if (!isnan(dAuxWaterLevelDiff))
//             {
//                if (tAbs(dAuxWaterLevelDiff) < 1) // Limiting the maximum value that can be found (dAuxWaterLevelDiff != DBL_NODATA)
//                {
//                   pointCounter++;
//                   dDiffTotWaterLevel += dAuxWaterLevelDiff;
//                }
//             }
//          }
//       }
//
//       if (pointCounter > 0)
//          dDiffTotWaterLevel /= pointCounter;
//       else
//          dDiffTotWaterLevel = 0;
//    }
//    else
//    {
//       for (long unsigned int n = 0; n < nLocIDs; n++)
//       {
//          double const dPointGridXExtCRS = m_VdFloodLocationX[n];
//          double const dPointGridYExtCRS = m_VdFloodLocationY[n];
//          double dMinDiffTotWaterLevelAtCoast = 1e10;
//
//          for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
//          {
//             int const nCoastSize = m_VCoast[nCoast].pLGetCoastlineExtCRS()->nGetSize();
//             double dMinDistSquare = 1e10;
//
//             for (int nCoastPoint = 0; nCoastPoint < nCoastSize; nCoastPoint++)
//             {
//                double const dCoastPointXExtCRS = m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nCoastPoint)->dGetX();
//                double const dCoastPointYExtCRS = m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nCoastPoint)->dGetY();
//
//                double const dDistSquare = (dCoastPointXExtCRS - dPointGridXExtCRS) * (dCoastPointXExtCRS - dPointGridXExtCRS) + (dCoastPointYExtCRS - dPointGridYExtCRS) * (dCoastPointYExtCRS - dPointGridYExtCRS);
//
//                if (dDistSquare < dMinDistSquare)
//                {
//                   dAuxWaterLevelDiff = m_VCoast[nCoast].dGetLevel(nCoastPoint, m_nLevel);
//
//                   if (!isnan(dAuxWaterLevelDiff))
//                   {
//                      dMinDistSquare = dDistSquare;
//                      dMinDiffTotWaterLevelAtCoast = dAuxWaterLevelDiff;
//                   }
//                }
//             }
//          }
//
//          dDiffTotWaterLevel += dMinDiffTotWaterLevelAtCoast;
//       }
//
//       dDiffTotWaterLevel /= static_cast<double>(nLocIDs);
//    }
//
//    m_dThisIterDiffTotWaterLevel = dDiffTotWaterLevel;
//
//    switch (m_nLevel)
//    {
//    case 0: // WAVESETUP + STORMSURGE:
//       m_dThisIterDiffWaveSetupSurgeWaterLevel = m_dThisIterDiffTotWaterLevel;
//       break;
//
//    case 1: // WAVESETUP + STORMSURGE + RUNUP:
//       m_dThisIterDiffWaveSetupSurgeRunupWaterLevel = m_dThisIterDiffTotWaterLevel;
//       break;
//    }
//
//    // Create an empty stack
//    stack<CGeom2DIPoint> PtiStackFlood;
//
//    // Start at the given edge cell, push this onto the stack
//    PtiStackFlood.push(CGeom2DIPoint(nXStart, nYStart));
//
//    // Then do the cell-by-cell fill loop until there are no more cell coordinates on the stack
//    while (! PtiStackFlood.empty())
//    {
//       CGeom2DIPoint const Pti = PtiStackFlood.top();
//       PtiStackFlood.pop();
//
//       int nX = Pti.nGetX();
//       int const nY = Pti.nGetY();
//
//       while (nX >= 0)
//       {
//          if (m_pRasterGrid->m_Cell[nX][nY].bIsCellFloodCheck())
//             break;
//
//          if (! m_pRasterGrid->m_Cell[nX][nY].bElevLessThanSWL())
//             break;
//
//          nX--;
//       }
//
//       nX++;
//
//       bool bSpanAbove = false;
//       bool bSpanBelow = false;
//
//       while (nX < m_nXGridSize)
//       {
//          if (m_pRasterGrid->m_Cell[nX][nY].bIsCellFloodCheck())
//             break;
//
//          if (! m_pRasterGrid->m_Cell[nX][nY].bElevLessThanSWL())
//             break;
//
//          // Flood this cell
//          m_pRasterGrid->m_Cell[nX][nY].SetCheckFloodCell();
//          m_pRasterGrid->m_Cell[nX][nY].SetInContiguousFlood();
//
//          switch (m_nLevel)
//          {
//          case 0: // WAVESETUP + STORMSURGE:
//             m_pRasterGrid->m_Cell[nX][nY].SetFloodBySetupSurge();
//             break;
//
//          case 1: // WAVESETUP + STORMSURGE + RUNUP:
//             m_pRasterGrid->m_Cell[nX][nY].SetFloodBySetupSurgeRunup();
//             break;
//          }
//
//          if ((! bSpanAbove) && (nY > 0) && (m_pRasterGrid->m_Cell[nX][nY - 1].bElevLessThanSWL()) && (!m_pRasterGrid->m_Cell[nX][nY - 1].bIsCellFloodCheck()))
//          {
//             PtiStackFlood.push(CGeom2DIPoint(nX, nY - 1));
//             bSpanAbove = true;
//          }
//
//          else if (bSpanAbove && (nY > 0) && (!m_pRasterGrid->m_Cell[nX][nY - 1].bElevLessThanSWL()))
//          {
//             bSpanAbove = false;
//          }
//
//          if ((! bSpanBelow) && (nY < m_nYGridSize - 1) && (m_pRasterGrid->m_Cell[nX][nY + 1].bElevLessThanSWL()) && (!m_pRasterGrid->m_Cell[nX][nY + 1].bIsCellFloodCheck()))
//          {
//             PtiStackFlood.push(CGeom2DIPoint(nX, nY + 1));
//             bSpanBelow = true;
//          }
//
//          else if (bSpanBelow && (nY < m_nYGridSize - 1) && (!m_pRasterGrid->m_Cell[nX][nY + 1].bElevLessThanSWL()))
//          {
//             bSpanBelow = false;
//          }
//
//          nX++;
//       }
//    }
// }

// //===============================================================================================================================
// //! Locates all the potential coastline start points on the edges of the raster grid, then traces vector coastline(s) from these start points
// //===============================================================================================================================
// int CSimulation::nTraceAllFloodCoasts(void)
// {
//    vector<bool> VbPossibleStartCellLHEdge;
//    vector<bool> VbTraced;
//    vector<int> VnSearchDirection;
//    vector<CGeom2DIPoint> V2DIPossibleStartCell;
//
//    // Go along the list of edge cells and look for possible coastline start cells
//    for (unsigned int n = 0; n < m_VPtiAllEdgeCell.size() - 1; n++)
//    {
//       if (m_bOmitSearchNorthEdge && (m_VPtiAllEdgeCellEdge[n] == NORTH || m_VPtiAllEdgeCellEdge[n + 1] == NORTH))
//          continue;
//
//       if (m_bOmitSearchSouthEdge && (m_VPtiAllEdgeCellEdge[n] == SOUTH || m_VPtiAllEdgeCellEdge[n + 1] == SOUTH))
//          continue;
//
//       if (m_bOmitSearchWestEdge && (m_VPtiAllEdgeCellEdge[n] == WEST || m_VPtiAllEdgeCellEdge[n + 1] == WEST))
//          continue;
//
//       if (m_bOmitSearchEastEdge && (m_VPtiAllEdgeCellEdge[n] == EAST || m_VPtiAllEdgeCellEdge[n + 1] == EAST))
//          continue;
//
//       int const nXThis = m_VPtiAllEdgeCell[n].nGetX();
//       int const nYThis = m_VPtiAllEdgeCell[n].nGetY();
//       int const nXNext = m_VPtiAllEdgeCell[n + 1].nGetX();
//       int const nYNext = m_VPtiAllEdgeCell[n + 1].nGetY();
//
//       // Get "Is it sea?" information for 'this' and 'next' cells
//       bool const bThisCellIsSea = m_pRasterGrid->m_Cell[nXThis][nYThis].bIsInContiguousSeaFlood();
//       bool const bNextCellIsSea = m_pRasterGrid->m_Cell[nXNext][nYNext].bIsInContiguousSeaFlood();
//
//       // Are we at a coast?
//       if ((! bThisCellIsSea) && bNextCellIsSea)
//       {
//          // 'This' cell is just inland, has it already been flagged as a possible start for a coastline (even if this subsequently 'failed' as a coastline)?
//          // if (! m_pRasterGrid->m_Cell[nXThis][nYThis].bIsPossibleCoastStartCell())
//          {
//             // It has not, so flag it
//             m_pRasterGrid->m_Cell[nXThis][nYThis].SetPossibleFloodStartCell();
//
//             // And save it
//             V2DIPossibleStartCell.push_back(CGeom2DIPoint(nXThis, nYThis));
//             VbPossibleStartCellLHEdge.push_back(true);
//             VnSearchDirection.push_back(nGetOppositeDirection(m_VPtiAllEdgeCellEdge[n]));
//             VbTraced.push_back(false);
//          }
//       }
//
//       else if (bThisCellIsSea && (! bNextCellIsSea))
//       {
//          // The 'next' cell is just inland, has it already been flagged as a possible start for a coastline (even if this subsequently 'failed' as a coastline)?
//          // if (! m_pRasterGrid->m_Cell[nXNext][nYNext].bIsPossibleCoastStartCell())
//          {
//             // It has not, so flag it
//             m_pRasterGrid->m_Cell[nXNext][nYNext].SetPossibleFloodStartCell();
//
//             // And save it
//             V2DIPossibleStartCell.push_back(CGeom2DIPoint(nXNext, nYNext));
//             VbPossibleStartCellLHEdge.push_back(false);
//             VnSearchDirection.push_back(nGetOppositeDirection(m_VPtiAllEdgeCellEdge[n + 1]));
//             VbTraced.push_back(false);
//          }
//       }
//    }
//
//    bool bAtLeastOneCoastTraced = false;
//
//    for (unsigned int n = 0; n < V2DIPossibleStartCell.size(); n++)
//    {
//       if (!VbTraced[n])
//       {
//          int nRet = 0;
//
//          if (VbPossibleStartCellLHEdge[n])
//          {
//             nRet = nTraceFloodCoastLine(n, VnSearchDirection[n], LEFT_HANDED, &VbTraced, &V2DIPossibleStartCell);
//          }
//
//          else
//          {
//             nRet = nTraceFloodCoastLine(n, VnSearchDirection[n], RIGHT_HANDED, &VbTraced, &V2DIPossibleStartCell);
//          }
//
//          if (nRet == RTN_OK)
//          {
//             // We have a valid coastline starting from this possible start cell
//             VbTraced[n] = true;
//             bAtLeastOneCoastTraced = true;
//          }
//       }
//    }
//
//    if (bAtLeastOneCoastTraced)
//       return RTN_OK;
//    else
//       return RTN_ERR_TRACING_FLOOD;
// }

// //===============================================================================================================================
// //! Traces a coastline (which is defined to be just above still water level) on the grid using the 'wall follower' rule for maze traversal (http://en.wikipedia.org/wiki/Maze_solving_algorithm#Wall_follower). The vector coastlines are then smoothed
// //===============================================================================================================================
// int CSimulation::nTraceFloodCoastLine(unsigned int const nTraceFromStartCellIndex, int const nStartSearchDirection, int const nHandedness, vector<bool>* pVbTraced, vector<CGeom2DIPoint> const* pV2DIPossibleStartCell)
// {
//    bool bHitStartCell = false;
//    bool bAtCoast = false;
//    bool bHasLeftStartEdge = false;
//    bool bTooLong = false;
//    bool bOffEdge = false;
//    bool bRepeating = false;
//
//    int const nStartX = pV2DIPossibleStartCell->at(nTraceFromStartCellIndex).nGetX();
//    int const nStartY = pV2DIPossibleStartCell->at(nTraceFromStartCellIndex).nGetY();
//    int nX = nStartX;
//    int nY = nStartY;
//    int nSearchDirection = nStartSearchDirection;
//    int nRoundLoop = -1;
//    // nThisLen = 0;
//    // nLastLen = 0,
//    // nPreLastLen = 0;
//
//    // Temporary coastline as integer points (grid CRS)
//    CGeomILine ILTempGridCRS;
//
//    // Mark the start cell as coast and add it to the vector object
//    m_pRasterGrid->m_Cell[nStartX][nStartY].SetAsFloodline(true);
//    CGeom2DIPoint const PtiStart(nStartX, nStartY);
//    ILTempGridCRS.Append(&PtiStart);
//
//    // Start at this grid-edge point and trace the rest of the coastline using the 'wall follower' rule for maze traversal, trying to keep next to cells flagged as sea
//    do
//    {
//       //       // DEBUG CODE ==============================================================================================================
//       // LogStream << "Now at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;
//       // LogStream << "ILTempGridCRS is now:" << endl;
//       // for (int n = 0; n < ILTempGridCRS.nGetSize(); n++)
//       // LogStream << "[" << ILTempGridCRS[n].nGetX() << "][" << ILTempGridCRS[n].nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILTempGridCRS[n].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILTempGridCRS[n].nGetY()) << "}" << endl;
//       // LogStream <<  "=================" << endl;
//       //       // DEBUG CODE ==============================================================================================================
//
//       // Safety check
//       if (++nRoundLoop > m_nCoastMax)
//       {
//          bTooLong = true;
//
//          // LogStream << m_ulIter << ":\t coast " << nCoast << " abandoning coastline tracing from [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "}, exceeded maximum search length (" << m_nCoastMax << ")" << endl;
//
//          // for (int n = 0; n < ILTempGridCRS.nGetSize(); n++)
//          // LogStream << "[" << ILTempGridCRS[n].nGetX() << "][" << ILTempGridCRS[n].nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILTempGridCRS[n].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILTempGridCRS[n].nGetY()) << "}" << endl;
//          // LogStream << endl;
//
//          break;
//       }
//
//       // Another safety check
//       if ((nRoundLoop > 10) && (ILTempGridCRS.nGetSize() < 2))
//       {
//          // We've been 10 times round the loop but the coast is still less than 2 coastline points in length, so we must be repeating
//          bRepeating = true;
//
//          break;
//       }
//
//       // OK so far: so have we left the start edge?
//       if (! bHasLeftStartEdge)
//       {
//          // We have not yet left the start edge
//          if (((nStartSearchDirection == SOUTH) && (nY > nStartY)) || ((nStartSearchDirection == NORTH) && (nY < nStartY)) ||
//              ((nStartSearchDirection == EAST) && (nX > nStartX)) || ((nStartSearchDirection == WEST) && (nX < nStartX)))
//             bHasLeftStartEdge = true;
//
//          // Flag this cell to ensure that it is not chosen as a coastline start cell later
//          m_pRasterGrid->m_Cell[nX][nY].SetPossibleFloodStartCell();
//          // LogStream << "Flagging [" << nX << "][" << nY << "] as possible coast start cell NOT YET LEFT EDGE" << endl;
//       }
//
//       // Leave the loop if the vector coastline has left the start edge, then we find a coast cell which is a possible start cell from which a coastline has not yet been traced
//       // if (bHasLeftStartEdge && bAtCoast)
//       {
//          for (unsigned int nn = 0; nn < pVbTraced->size(); nn++)
//          {
//             if ((nn != nTraceFromStartCellIndex) && (! pVbTraced->at(nn)))
//             {
//                // LogStream << "[" << pV2DIPossibleStartCell->at(nn).nGetX() << "][" << pV2DIPossibleStartCell->at(nn).nGetY() << "]" << endl;
//
//                if (bAtCoast && (nX == pV2DIPossibleStartCell->at(nn).nGetX()) && (nY == pV2DIPossibleStartCell->at(nn).nGetY()))
//                {
//                   if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
//                      LogStream << m_ulIter << ": Possible flood coastline found, traced from [" << nStartX << "][" << nStartY << "] and hit another possible flood coast start cell at [" << nX << "][" << nY << "]" << endl;
//
//                   pVbTraced->at(nn) = true;
//                   bHitStartCell = true;
//                   break;
//                }
//             }
//          }
//
//          // LogStream << endl;
//       }
//
//       if (bHitStartCell)
//          break;
//
//       // OK now sort out the next iteration of the search
//       int nXSeaward = 0;
//       int nYSeaward = 0;
//       int nSeawardNewDirection = 0;
//       int nXStraightOn = 0;
//       int nYStraightOn = 0;
//       int nXAntiSeaward = 0;
//       int nYAntiSeaward = 0;
//       int nAntiSeawardNewDirection = 0;
//       int nXGoBack = 0;
//       int nYGoBack = 0;
//       int nGoBackNewDirection = 0;
//
//       CGeom2DIPoint const Pti(nX, nY);
//
//       // Set up the variables
//       switch (nHandedness)
//       {
//       case RIGHT_HANDED:
//
//          // The sea is to the right-hand side of the coast as we traverse it. We are just inland, so we need to keep heading right to find the sea
//          switch (nSearchDirection)
//          {
//          case NORTH:
//             // The sea is towards the RHS (E) of the coast, so first try to go right (to the E)
//             nXSeaward = nX + 1;
//             nYSeaward = nY;
//             nSeawardNewDirection = EAST;
//
//             // If can't do this, try to go straight on (to the N)
//             nXStraightOn = nX;
//             nYStraightOn = nY - 1;
//
//             // If can't do either of these, try to go anti-seaward i.e. towards the LHS (W)
//             nXAntiSeaward = nX - 1;
//             nYAntiSeaward = nY;
//             nAntiSeawardNewDirection = WEST;
//
//             // As a last resort, go back (to the S)
//             nXGoBack = nX;
//             nYGoBack = nY + 1;
//             nGoBackNewDirection = SOUTH;
//
//             break;
//
//          case EAST:
//             // The sea is towards the RHS (S) of the coast, so first try to go right (to the S)
//             nXSeaward = nX;
//             nYSeaward = nY + 1;
//             nSeawardNewDirection = SOUTH;
//
//             // If can't do this, try to go straight on (to the E)
//             nXStraightOn = nX + 1;
//             nYStraightOn = nY;
//
//             // If can't do either of these, try to go anti-seaward i.e. towards the LHS (N)
//             nXAntiSeaward = nX;
//             nYAntiSeaward = nY - 1;
//             nAntiSeawardNewDirection = NORTH;
//
//             // As a last resort, go back (to the W)
//             nXGoBack = nX - 1;
//             nYGoBack = nY;
//             nGoBackNewDirection = WEST;
//
//             break;
//
//          case SOUTH:
//             // The sea is towards the RHS (W) of the coast, so first try to go right (to the W)
//             nXSeaward = nX - 1;
//             nYSeaward = nY;
//             nSeawardNewDirection = WEST;
//
//             // If can't do this, try to go straight on (to the S)
//             nXStraightOn = nX;
//             nYStraightOn = nY + 1;
//
//             // If can't do either of these, try to go anti-seaward i.e. towards the LHS (E)
//             nXAntiSeaward = nX + 1;
//             nYAntiSeaward = nY;
//             nAntiSeawardNewDirection = EAST;
//
//             // As a last resort, go back (to the N)
//             nXGoBack = nX;
//             nYGoBack = nY - 1;
//             nGoBackNewDirection = NORTH;
//
//             break;
//
//          case WEST:
//             // The sea is towards the RHS (N) of the coast, so first try to go right (to the N)
//             nXSeaward = nX;
//             nYSeaward = nY - 1;
//             nSeawardNewDirection = NORTH;
//
//             // If can't do this, try to go straight on (to the W)
//             nXStraightOn = nX - 1;
//             nYStraightOn = nY;
//
//             // If can't do either of these, try to go anti-seaward i.e. towards the LHS (S)
//             nXAntiSeaward = nX;
//             nYAntiSeaward = nY + 1;
//             nAntiSeawardNewDirection = SOUTH;
//
//             // As a last resort, go back (to the E)
//             nXGoBack = nX + 1;
//             nYGoBack = nY;
//             nGoBackNewDirection = EAST;
//
//             break;
//          }
//
//          break;
//
//       case LEFT_HANDED:
//
//          // The sea is to the left-hand side of the coast as we traverse it. We are just inland, so we need to keep heading left to find the sea
//          switch (nSearchDirection)
//          {
//          case NORTH:
//             // The sea is towards the LHS (W) of the coast, so first try to go left (to the W)
//             nXSeaward = nX - 1;
//             nYSeaward = nY;
//             nSeawardNewDirection = WEST;
//
//             // If can't do this, try to go straight on (to the N)
//             nXStraightOn = nX;
//             nYStraightOn = nY - 1;
//
//             // If can't do either of these, try to go anti-seaward i.e. towards the RHS (E)
//             nXAntiSeaward = nX + 1;
//             nYAntiSeaward = nY;
//             nAntiSeawardNewDirection = EAST;
//
//             // As a last resort, go back (to the S)
//             nXGoBack = nX;
//             nYGoBack = nY + 1;
//             nGoBackNewDirection = SOUTH;
//
//             break;
//
//          case EAST:
//             // The sea is towards the LHS (N) of the coast, so first try to go left (to the N)
//             nXSeaward = nX;
//             nYSeaward = nY - 1;
//             nSeawardNewDirection = NORTH;
//
//             // If can't do this, try to go straight on (to the E)
//             nXStraightOn = nX + 1;
//             nYStraightOn = nY;
//
//             // If can't do either of these, try to go anti-seaward i.e. towards the RHS (S)
//             nXAntiSeaward = nX;
//             nYAntiSeaward = nY + 1;
//             nAntiSeawardNewDirection = SOUTH;
//
//             // As a last resort, go back (to the W)
//             nXGoBack = nX - 1;
//             nYGoBack = nY;
//             nGoBackNewDirection = WEST;
//
//             break;
//
//          case SOUTH:
//             // The sea is towards the LHS (E) of the coast, so first try to go left (to the E)
//             nXSeaward = nX + 1;
//             nYSeaward = nY;
//             nSeawardNewDirection = EAST;
//
//             // If can't do this, try to go straight on (to the S)
//             nXStraightOn = nX;
//             nYStraightOn = nY + 1;
//
//             // If can't do either of these, try to go anti-seaward i.e. towards the RHS (W)
//             nXAntiSeaward = nX - 1;
//             nYAntiSeaward = nY;
//             nAntiSeawardNewDirection = WEST;
//
//             // As a last resort, go back (to the N)
//             nXGoBack = nX;
//             nYGoBack = nY - 1;
//             nGoBackNewDirection = NORTH;
//
//             break;
//
//          case WEST:
//             // The sea is towards the LHS (S) of the coast, so first try to go left (to the S)
//             nXSeaward = nX;
//             nYSeaward = nY + 1;
//             nSeawardNewDirection = SOUTH;
//
//             // If can't do this, try to go straight on (to the W)
//             nXStraightOn = nX - 1;
//             nYStraightOn = nY;
//
//             // If can't do either of these, try to go anti-seaward i.e. towards the RHS (N)
//             nXAntiSeaward = nX;
//             nYAntiSeaward = nY - 1;
//             nAntiSeawardNewDirection = NORTH;
//
//             // As a last resort, go back (to the E)
//             nXGoBack = nX + 1;
//             nYGoBack = nY;
//             nGoBackNewDirection = EAST;
//
//             break;
//          }
//
//          break;
//       }
//
//       // Now do the actual search for this timestep: first try going in the direction of the sea. Is this seaward cell still within the grid?
//       if (bIsWithinValidGrid(nXSeaward, nYSeaward))
//       {
//          // It is, so check if the cell in the seaward direction is a sea cell
//          if (m_pRasterGrid->m_Cell[nXSeaward][nYSeaward].bIsInContiguousSeaFlood())
//          {
//             // There is sea in this seaward direction, so we are on the coast
//             bAtCoast = true;
//
//             // Has the current cell already marked been marked as a coast cell?
//             if (! m_pRasterGrid->m_Cell[nX][nY].bIsFloodline())
//             {
//                // Not already marked, is this an intervention cell with the top above SWL?
//                if ((bIsInterventionCell(nX, nY)) && (!m_pRasterGrid->m_Cell[nX][nY].bElevLessThanSWL()))
//                {
//                   // It is, so mark as coast and add it to the vector object
//                   m_pRasterGrid->m_Cell[nX][nY].SetAsFloodline(true);
//                   ILTempGridCRS.Append(&Pti);
//                }
//
//                else if (! m_pRasterGrid->m_Cell[nX][nY].bElevLessThanSWL())
//                {
//                   // The sediment top is above SWL so mark as coast and add it to the vector object
//                   m_pRasterGrid->m_Cell[nX][nY].SetAsFloodline(true);
//                   ILTempGridCRS.Append(&Pti);
//                }
//             }
//          }
//
//          else
//          {
//             // The seaward cell is not a sea cell, so we will move to it next time
//             nX = nXSeaward;
//             nY = nYSeaward;
//
//             // And set a new search direction, to keep turning seaward
//             nSearchDirection = nSeawardNewDirection;
//             continue;
//          }
//       }
//
//       // OK, we couldn't move seaward (but we may have marked the current cell as coast) so next try to move straight on. Is this straight-ahead cell still within the grid?
//       if (bIsWithinValidGrid(nXStraightOn, nYStraightOn))
//       {
//          // It is, so check if there is sea immediately in front
//          if (m_pRasterGrid->m_Cell[nXStraightOn][nYStraightOn].bIsInContiguousSeaFlood())
//          {
//             // Sea is in front, so we are on the coast
//             bAtCoast = true;
//
//             // Has the current cell already marked been marked as a floodline cell?
//             if (! m_pRasterGrid->m_Cell[nX][nY].bIsFloodline())
//             {
//                // Not already marked, is this an intervention cell with the top above SWL?
//                if ((bIsInterventionCell(nX, nY)) && (!m_pRasterGrid->m_Cell[nX][nY].bElevLessThanSWL()))
//                {
//                   // It is, so mark as coast and add it to the vector object
//                   m_pRasterGrid->m_Cell[nX][nY].SetAsFloodline(true);
//                   ILTempGridCRS.Append(&Pti);
//                }
//
//                else if (! m_pRasterGrid->m_Cell[nX][nY].bElevLessThanSWL())
//                {
//                   // The sediment top is above SWL so mark as coast and add it to the vector object
//                   m_pRasterGrid->m_Cell[nX][nY].SetAsFloodline(true);
//                   ILTempGridCRS.Append(&Pti);
//                }
//             }
//          }
//
//          else
//          {
//             // The straight-ahead cell is not a sea cell, so we will move to it next time
//             nX = nXStraightOn;
//             nY = nYStraightOn;
//
//             // The search direction remains unchanged
//             continue;
//          }
//       }
//
//       // Couldn't move either seaward or straight on (but we may have marked the current cell as coast) so next try to move in the anti-seaward direction. Is this anti-seaward cell still within the grid?
//       if (bIsWithinValidGrid(nXAntiSeaward, nYAntiSeaward))
//       {
//          // It is, so check if there is sea in this anti-seaward cell
//          if (m_pRasterGrid->m_Cell[nXAntiSeaward][nYAntiSeaward].bIsInContiguousSeaFlood())
//          {
//             // There is sea on the anti-seaward side, so we are on the coast
//             bAtCoast = true;
//
//             // Has the current cell already marked been marked as a floodline cell?
//             if (! m_pRasterGrid->m_Cell[nX][nY].bIsFloodline())
//             {
//                // Not already marked, is this an intervention cell with the top above SWL?
//                if ((bIsInterventionCell(nX, nY)) && (!m_pRasterGrid->m_Cell[nX][nY].bElevLessThanSWL()))
//                {
//                   // It is, so mark as coast and add it to the vector object
//                   m_pRasterGrid->m_Cell[nX][nY].SetAsFloodline(true);
//                   ILTempGridCRS.Append(&Pti);
//                }
//
//                else if (! m_pRasterGrid->m_Cell[nX][nY].bElevLessThanSWL())
//                {
//                   // The sediment top is above SWL so mark as coast and add it to the vector object
//                   m_pRasterGrid->m_Cell[nX][nY].SetAsFloodline(true);
//                   ILTempGridCRS.Append(&Pti);
//                }
//             }
//          }
//
//          else
//          {
//             // The anti-seaward cell is not a sea cell, so we will move to it next time
//             nX = nXAntiSeaward;
//             nY = nYAntiSeaward;
//
//             // And set a new search direction, to keep turning seaward
//             nSearchDirection = nAntiSeawardNewDirection;
//             continue;
//          }
//       }
//
//       // Could not move to the seaward side, move straight ahead, or move to the anti-seaward side, so we must be in a single-cell dead end! As a last resort, turn round and move back to where we just came from, but first check that this is a valid cell
//       if (bIsWithinValidGrid(nXGoBack, nYGoBack))
//       {
//          nX = nXGoBack;
//          nY = nYGoBack;
//
//          // And change the search direction
//          nSearchDirection = nGoBackNewDirection;
//       }
//
//       else
//       {
//          // Our final choice is not a valid cell, so give up
//          bOffEdge = true;
//          break;
//       }
//    } while (true);
//
//    // OK, we have finished tracing this coastline on the grid. But is the coastline too long or too short?
//    int nCoastSize = ILTempGridCRS.nGetSize();
//
//    if (bOffEdge)
//    {
//       if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
//          LogStream << m_ulIter << ": Ignoring possible flood coastline from [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "} since hit off-edge cell at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}, coastline size is " << nCoastSize << endl;
//
//       // Unmark these cells as coast cells
//       for (int n = 0; n < nCoastSize; n++)
//          m_pRasterGrid->m_Cell[ILTempGridCRS[n].nGetX()][ILTempGridCRS[n].nGetY()].SetAsFloodline(false);
//
//       return RTN_ERR_TRACING_FLOOD;
//    }
//
//    if (bTooLong)
//    {
//       // Around loop too many times, so abandon this coastline
//       if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
//       {
//          LogStream << m_ulIter << ":\t abandoning possible flood coastline from [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "} since round loop " << nRoundLoop << " times (m_nCoastMax = " << m_nCoastMax << "), coastline size is " << nCoastSize;
//
//          if (nCoastSize > 0)
//             LogStream << ", it ended at [" << ILTempGridCRS[nCoastSize - 1].nGetX() << "][" << ILTempGridCRS[nCoastSize - 1].nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILTempGridCRS[nCoastSize - 1].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILTempGridCRS[nCoastSize - 1].nGetY()) << "}";
//
//          LogStream << endl;
//       }
//
//       // Unmark these cells as coast cells
//       for (int n = 0; n < nCoastSize; n++)
//          m_pRasterGrid->m_Cell[ILTempGridCRS[n].nGetX()][ILTempGridCRS[n].nGetY()].SetAsFloodline(false);
//
//       return RTN_ERR_TRACING_FLOOD;
//    }
//
//    if (bRepeating)
//    {
//       if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
//       {
//          LogStream << m_ulIter << ": Ignoring possible flood coastline from [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "} since repeating, coastline size is " << nCoastSize;
//
//          if (nCoastSize > 0)
//             LogStream << ", it ended at [" << ILTempGridCRS[nCoastSize - 1].nGetX() << "][" << ILTempGridCRS[nCoastSize - 1].nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILTempGridCRS[nCoastSize - 1].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILTempGridCRS[nCoastSize - 1].nGetY()) << "}";
//
//          LogStream << endl;
//       }
//
//       // Unmark these cells as coast cells
//       for (int n = 0; n < nCoastSize; n++)
//          m_pRasterGrid->m_Cell[ILTempGridCRS[n].nGetX()][ILTempGridCRS[n].nGetY()].SetAsFloodline(false);
//
//       return RTN_ERR_TRACING_FLOOD;
//    }
//
//    if (nCoastSize == 0)
//    {
//       // Zero-length coastline, so abandon it
//       if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
//          LogStream << m_ulIter << ": abandoning zero-length flood coastline from [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "}" << endl;
//
//       return RTN_ERR_TRACING_FLOOD;
//    }
//
//    if (nCoastSize < m_nCoastMin)
//    {
//       // The vector coastline is too small, so abandon it
//       if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
//          LogStream << m_ulIter << ":\t Ignoring possible flood coastline from [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "} to [" << ILTempGridCRS[nCoastSize - 1].nGetX() << "][" << ILTempGridCRS[nCoastSize - 1].nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILTempGridCRS[nCoastSize - 1].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILTempGridCRS[nCoastSize - 1].nGetY()) << "} since size (" << nCoastSize << ") is less than minimum (" << m_nCoastMin << ")" << endl;
//
//       // Unmark these cells as coast cells
//       for (int n = 0; n < nCoastSize; n++)
//          m_pRasterGrid->m_Cell[ILTempGridCRS[n].nGetX()][ILTempGridCRS[n].nGetY()].SetAsFloodline(false);
//
//       return RTN_ERR_TRACING_FLOOD;
//    }
//
//    // OK this new coastline is fine
//    int const nEndX = nX;
//    int const nEndY = nY;
//    int const nCoastEndX = ILTempGridCRS[nCoastSize - 1].nGetX();
//    int const nCoastEndY = ILTempGridCRS[nCoastSize - 1].nGetY();
//
//    if ((nCoastEndX != nEndX) || (nCoastEndY != nEndY))
//    {
//       // The grid-edge cell at nEndX, nEndY is not already at end of ILTempGridCRS. But is the final cell in ILTempGridCRS already at the edge of the grid?
//       if (! m_pRasterGrid->m_Cell[nCoastEndX][nCoastEndY].bIsBoundingBoxEdge())
//       {
//          // The final cell in ILTempGridCRS is not a grid-edge cell, so add the grid-edge cell and mark the cell as coastline
//          ILTempGridCRS.Append(nEndX, nEndY);
//          nCoastSize++;
//
//          m_pRasterGrid->m_Cell[nEndX][nEndY].SetAsFloodline(true);
//       }
//    }
//
//    // Need to specify start edge and end edge for smoothing routines
//    // int
//    // nStartEdge = m_pRasterGrid->m_Cell[nStartX][nStartY].nGetBoundingBoxEdge(),
//    // nEndEdge = m_pRasterGrid->m_Cell[nEndX][nEndY].nGetBoundingBoxEdge();
//
//    // Next, convert the grid coordinates in ILTempGridCRS (integer values stored as doubles) to external CRS coordinates (which will probably be non-integer, again stored as doubles). This is done now, so that smoothing is more effective
//    CGeomLine LTempExtCRS;
//
//    for (int j = 0; j < nCoastSize; j++)
//    {
//       LTempExtCRS.Append(dGridCentroidXToExtCRSX(ILTempGridCRS[j].nGetX()), dGridCentroidYToExtCRSY(ILTempGridCRS[j].nGetY()));
//    }
//
//    // Now do some smoothing of the vector output, if desired
//    // if (m_nCoastSmooth == SMOOTH_RUNNING_MEAN)
//    // LTempExtCRS = LSmoothCoastRunningMean(&LTempExtCRS);
//    // else if (m_nCoastSmooth == SMOOTH_SAVITZKY_GOLAY)
//    // LTempExtCRS = LSmoothCoastSavitzkyGolay(&LTempExtCRS, nStartEdge, nEndEdge);
//
//    //    // DEBUG CODE ==================================================================================================
//    // LogStream << "==================================" << endl;
//    // for (int j = 0; j < nCoastSize; j++)
//    // {
//    // LogStream << "{" << dGridCentroidXToExtCRSX(ILTempGridCRS[j].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILTempGridCRS[j].nGetY()) << "}" << "\t{" << LTempExtCRS.dGetXAt(j) << ", " << LTempExtCRS.dGetYAt(j) << "}" << endl;
//    // }
//    // LogStream << "==================================" << endl;
//    //    // DEBUG CODE ==================================================================================================
//
//    // Create a new coastline object and append to it the vector of coastline objects
//    CRWCoast const CoastTmp(this);
//    int nCoast;
//
//    switch (m_nLevel)
//    {
//    case 0:
//       m_VFloodWaveSetupSurge.push_back(CoastTmp);
//       nCoast = static_cast<int>(m_VFloodWaveSetupSurge.size()) - 1;
//       m_VFloodWaveSetupSurge[nCoast].SetCoastlineExtCRS(&LTempExtCRS);
//       break;
//
//    case 1:
//       m_VFloodWaveSetupSurgeRunup.push_back(CoastTmp);
//       nCoast = static_cast<int>(m_VFloodWaveSetupSurgeRunup.size()) - 1;
//       m_VFloodWaveSetupSurgeRunup[nCoast].SetCoastlineExtCRS(&LTempExtCRS);
//    }
//
//    return RTN_OK;
// }
