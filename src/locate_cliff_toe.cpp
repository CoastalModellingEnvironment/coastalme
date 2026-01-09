/*!
  \file locate_cliff_toe.cpp
  \brief Locates and traces the cliff toe on the raster grid
  \details TODO 001 A more detailed description of these routines.
  \author Wilf Chun
  \author David Favis-Mortlock
  \author Andres Payo
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
using std::sqrt;

#include <cstdio>
using std::size_t;

#include <stdint.h>

#include <sstream>
using std::stringstream;

// #include <algorithm>
// using std::max;
// using std::min;

#include <utility>
// using std::make_pair;
using std::pair;

// #include "cme.h"
// #include "cell.h"
// #include "raster_grid.h"
// #include "simulation.h"
// #include "2di_point.h"
// #include "line.h"

// /*===============================================================================================================================
// //! Locates and traces the cliff toe
// ===============================================================================================================================*/
// int CSimulation::nLocateCliffToe(void)
// {
//    nCalcSlopeAtAllCells();
//    nLocateCliffCell();
//    nRemoveSmallCliffIslands(50);
//    nTraceSeawardCliffEdge();
//    nValidateCliffToeEdges();
//
//    return RTN_OK;
// }

// /*===============================================================================================================================
// //! Calculates slope of the top surface of all sediment for all cells, using finite difference method, for cliff toe detection
// ===============================================================================================================================*/
// void CSimulation::nCalcSlopeAtAllCells(void)
// {
//    for (int nX = 0; nX < m_nXGridSize; nX++)
//    {
//       for (int nY = 0; nY < m_nYGridSize; nY++)
//       {
//          // Look at the four surounding cells
//          if ((nX > 0) && (nX < m_nXGridSize - 1) && (nY > 0) && (nY < m_nYGridSize - 1))
//          {
//             double const dElevLeft = m_pRasterGrid->m_Cell[nX - 1][nY].dGetAllSedTopElevOmitTalus();
//             double const dElevRight = m_pRasterGrid->m_Cell[nX + 1][nY].dGetAllSedTopElevOmitTalus();
//             double const dElevUp = m_pRasterGrid->m_Cell[nX][nY - 1].dGetAllSedTopElevOmitTalus();
//             double const dElevDown = m_pRasterGrid->m_Cell[nX][nY + 1].dGetAllSedTopElevOmitTalus();
//
//             // Calculate slope using finite difference method
//             double const dSlopeX = (dElevRight - dElevLeft) / (2.0 * m_dCellSide);
//             double const dSlopeY = (dElevDown - dElevUp) / (2.0 * m_dCellSide);
//             double const dSlope = sqrt(dSlopeX * dSlopeX + dSlopeY * dSlopeY);
//             m_pRasterGrid->m_Cell[nX][nY].SetSlopeForCliffToe(dSlope);
//          }
//       }
//    }
// }

// /*===============================================================================================================================
// //! Finds cells with slope greater than threshold, for cliff toe detection
// ===============================================================================================================================*/
// void CSimulation::nLocateCliffCell(void)
// {
//    for (int nX = 0; nX < m_nXGridSize; nX++)
//    {
//       for (int nY = 0; nY < m_nYGridSize; nY++)
//       {
//          double const dSlope = m_pRasterGrid->m_Cell[nX][nY].dGetSlopeForCliffToe();
//          if (dSlope >= m_dSlopeThresholdForCliffToe)
//          {
//             m_pRasterGrid->m_Cell[nX][nY].SetAsCliffToe(true);
//          }
//       }
//    }
// }

// /*===============================================================================================================================
// //! Removes small cliff toe islands if they are below a given number of cells, using flood fill algorithm
// ===============================================================================================================================*/
// void CSimulation::nRemoveSmallCliffIslands(int const dMinCliffCellThreshold)
// {
//    // Create a 2D array to track which cells have been visited during flood fill
//    vector<vector<bool>> bVisited(static_cast<unsigned int>(m_nXGridSize), vector<bool>(static_cast<unsigned int>(m_nYGridSize), false));
//
//    // Vector to store cells that belong to small cliff islands to be removed
//    vector<pair<int, int>> VSmallIslandCells;
//
//    // Loop through all cells to find unvisited cliff cells
//    for (unsigned int nX = 0; nX < static_cast<unsigned int>(m_nXGridSize); nX++)
//    {
//       for (unsigned int nY = 0; nY < static_cast<unsigned int>(m_nYGridSize); nY++)
//       {
//          // Check if this is an unvisited cliff cell
//          if ((! bVisited[nX][nY]) && m_pRasterGrid->m_Cell[nX][nY].bIsCliffToe())
//          {
//             // Found the start of a new cliff region - use flood fill to find all connected cliff cells
//             vector<pair<int, int>> VCurrentCliffRegion;
//
//             // Stack for iterative flood fill algorithm
//             vector<pair<int, int>> VStack;
//             VStack.push_back(make_pair(nX, nY));
//
//             // Flood fill to find all connected cliff cells
//             while (! VStack.empty())
//             {
//                pair<int, int> const currentCell = VStack.back();
//                VStack.pop_back();
//
//                size_t const nCurX = static_cast<size_t>(currentCell.first);
//                size_t const nCurY = static_cast<size_t>(currentCell.second);
//
//                // Skip if already visited or out of bounds
//                if (nCurX >= static_cast<unsigned int>(m_nXGridSize) ||
//                    nCurY >= static_cast<unsigned int>(m_nYGridSize) || bVisited[nCurX][nCurY] ||
//                    (! m_pRasterGrid->m_Cell[nCurX][nCurY].bIsCliffToe()))
//                {
//                   continue;
//                }
//
//                // Mark as visited and add to current cliff region
//                bVisited[nCurX][nCurY] = true;
//                VCurrentCliffRegion.push_back(make_pair(nCurX, nCurY));
//
//                // Add neighboring cells to stack (8-connectivity: N, NE, E, SE, S, SW, W, NW)
//                VStack.push_back(make_pair(nCurX - 1, nCurY));     // North
//                VStack.push_back(make_pair(nCurX - 1, nCurY + 1)); // Northeast
//                VStack.push_back(make_pair(nCurX, nCurY + 1));     // East
//                VStack.push_back(make_pair(nCurX + 1, nCurY + 1)); // Southeast
//                VStack.push_back(make_pair(nCurX + 1, nCurY));     // South
//                VStack.push_back(make_pair(nCurX + 1, nCurY - 1)); // Southwest
//                VStack.push_back(make_pair(nCurX, nCurY - 1));     // West
//                VStack.push_back(make_pair(nCurX - 1, nCurY - 1)); // Northwest
//             }
//
//             // Calculate area of this cliff region (number of cells * cell area)
//             int const dCliffRegionArea = static_cast<int>(VCurrentCliffRegion.size());
//
//             // If area is below threshold, mark all cells in this region for removal
//             if (dCliffRegionArea < dMinCliffCellThreshold)
//             {
//                VSmallIslandCells.insert(VSmallIslandCells.end(), VCurrentCliffRegion.begin(), VCurrentCliffRegion.end());
//             }
//          }
//       }
//    }
//
//    // Remove cliff designation from all small island cells
//    for (const auto &cell : VSmallIslandCells)
//    {
//       m_pRasterGrid->m_Cell[cell.first][cell.second].SetAsCliffToe(false);
//    }
// }

// /*===============================================================================================================================
// //! Traces the seaward extent of cliff cells using wall-following algorithm
// ===============================================================================================================================*/
// void CSimulation::nTraceSeawardCliffEdge(void)
// {
//    // Clear previous cliff edges
//    m_VCliffToe.clear();
//
//    // Find all possible cliff edge start points by scanning for cliff/non-cliff transitions
//    vector<CGeom2DIPoint> V2DIPossibleStartCell;
//    vector<bool> VbPossibleStartCellHandedness;
//    vector<int> VnSearchDirection;
//
//    // Scan all cells to find cliff toe edges by checking elevation patterns
//    for (int nX = 2; nX < m_nXGridSize - 2; nX++)
//    {
//       for (int nY = 2; nY < m_nYGridSize - 2; nY++)
//       {
//          if (m_pRasterGrid->m_Cell[nX][nY].bIsCliffToe())
//          {
//             // East direction (check if this is a seaward-facing cliff toe)
//             if (! m_pRasterGrid->m_Cell[nX][nY + 1].bIsCliffToe())
//             {
//                V2DIPossibleStartCell.push_back(CGeom2DIPoint(nX, nY));
//                VbPossibleStartCellHandedness.push_back(true);
//                VnSearchDirection.push_back(EAST);
//             }
//
//             // South direction
//             if (! m_pRasterGrid->m_Cell[nX + 1][nY].bIsCliffToe())
//             {
//                V2DIPossibleStartCell.push_back(CGeom2DIPoint(nX, nY));
//                VbPossibleStartCellHandedness.push_back(true);
//                VnSearchDirection.push_back(SOUTH);
//             }
//
//             // West direction
//             if (! m_pRasterGrid->m_Cell[nX][nY - 1].bIsCliffToe())
//             {
//                V2DIPossibleStartCell.push_back(CGeom2DIPoint(nX, nY));
//                VbPossibleStartCellHandedness.push_back(true);
//                VnSearchDirection.push_back(WEST);
//             }
//
//             // North direction
//             if (! m_pRasterGrid->m_Cell[nX - 1][nY].bIsCliffToe())
//             {
//                V2DIPossibleStartCell.push_back(CGeom2DIPoint(nX, nY));
//                VbPossibleStartCellHandedness.push_back(true);
//                VnSearchDirection.push_back(NORTH);
//             }
//          }
//       }
//    }
//
//    // Create a 2D array to track which cells have been used in cliff edge tracing
//    vector<vector<bool>> bUsedInCliffTrace(m_nXGridSize, vector<bool>(m_nYGridSize, false));
//
//    int nCliffEdgesTraced = 0;
//
//    // For each possible start point, attempt to trace a cliff edge
//    for (size_t nStartPoint = 0; nStartPoint < V2DIPossibleStartCell.size(); nStartPoint++)
//    {
//       int const nXStart = V2DIPossibleStartCell[nStartPoint].nGetX();
//       int const nYStart = V2DIPossibleStartCell[nStartPoint].nGetY();
//
//       // Skip if this cell has already been used in another cliff trace
//       if (bUsedInCliffTrace[nXStart][nYStart])
//       {
//          continue;
//       }
//
//       // Begin cliff edge tracing using wall following algorithm
//       vector<CGeom2DIPoint> VCliffEdge;
//       int nSearchDirection = VnSearchDirection[nStartPoint];
//
//       int nX = nXStart;
//       int nY = nYStart;
//       int nLength = 0;
//       const int nMaxLength = m_nXGridSize * m_nYGridSize; // Safety limit
//
//       do
//       {
//          // Add current point to cliff edge
//          VCliffEdge.push_back(CGeom2DIPoint(nX, nY));
//          bUsedInCliffTrace[nX][nY] = true;
//          nLength++;
//
//          // Find next cell using wall following algorithm (right-hand rule)
//          int nXSeaward;
//          int nYSeaward;
//          int nXStraightOn;
//          int nYStraightOn;
//          int nXAntiSeaward, nYAntiSeaward, nXGoBack, nYGoBack;
//
//          // Calculate candidate positions based on search direction
//          switch (nSearchDirection)
//          {
//          case NORTH:
//             nXSeaward = nX + 1;
//             nYSeaward = nY;            // Right (East)
//             nXStraightOn = nX;
//             nYStraightOn = nY - 1;     // Straight (North)
//             nXAntiSeaward = nX - 1;
//             nYAntiSeaward = nY;        // Left (West)
//             nXGoBack = nX;
//             nYGoBack = nY + 1;         // Back (South)
//
//             break;
//
//          case EAST:
//             nXSeaward = nX;
//             nYSeaward = nY + 1;        // Right (South)
//             nXStraightOn = nX + 1;
//             nYStraightOn = nY;         // Straight (East)
//             nXAntiSeaward = nX;
//             nYAntiSeaward = nY - 1;    // Left (North)
//             nXGoBack = nX - 1;
//             nYGoBack = nY;             // Back (West)
//
//             break;
//
//          case SOUTH:
//             nXSeaward = nX - 1;
//             nYSeaward = nY;            // Right (West)
//             nXStraightOn = nX;
//             nYStraightOn = nY + 1;     // Straight (South)
//             nXAntiSeaward = nX + 1;
//             nYAntiSeaward = nY;        // Left (East)
//             nXGoBack = nX;
//             nYGoBack = nY - 1;         // Back (North)
//
//             break;
//
//          case WEST:
//             nXSeaward = nX;
//             nYSeaward = nY - 1;        // Right (North)
//             nXStraightOn = nX - 1;
//             nYStraightOn = nY;         // Straight (West)
//             nXAntiSeaward = nX;
//             nYAntiSeaward = nY + 1;    // Left (South)
//             nXGoBack = nX + 1;
//             nYGoBack = nY;             // Back (East)
//
//             break;
//
//          default:
//             nXSeaward = nXStraightOn = nXAntiSeaward = nXGoBack = nX;
//             nYSeaward = nYStraightOn = nYAntiSeaward = nYGoBack = nY;
//
//             break;
//          }
//
//          // Try to move using wall following priority: seaward, straight, anti-seaward, back
//          bool bFoundNextCell = false;
//
//          // 1. Try seaward (right turn)
//          if (bIsWithinValidGrid(nXSeaward, nYSeaward) && m_pRasterGrid->m_Cell[nXSeaward][nYSeaward].bIsCliffToe())
//          {
//             nX = nXSeaward;
//             nY = nYSeaward;
//
//             // Update search direction (turn right)
//             switch (nSearchDirection)
//             {
//             case NORTH:
//                nSearchDirection = EAST;
//                break;
//
//             case EAST:
//                nSearchDirection = SOUTH;
//                break;
//
//             case SOUTH:
//                nSearchDirection = WEST;
//                break;
//
//             case WEST:
//                nSearchDirection = NORTH;
//                break;
//             }
//
//             bFoundNextCell = true;
//          }
//
//          // 2. Try straight ahead
//          else if (bIsWithinValidGrid(nXStraightOn, nYStraightOn) && m_pRasterGrid->m_Cell[nXStraightOn][nYStraightOn].bIsCliffToe())
//          {
//             nX = nXStraightOn;
//             nY = nYStraightOn;
//
//             // Direction stays the same
//             bFoundNextCell = true;
//          }
//
//          // 3. Try anti-seaward (left turn)
//          else if (bIsWithinValidGrid(nXAntiSeaward, nYAntiSeaward) && m_pRasterGrid->m_Cell[nXAntiSeaward][nYAntiSeaward].bIsCliffToe())
//          {
//             nX = nXAntiSeaward;
//             nY = nYAntiSeaward;
//
//             // Update search direction (turn left)
//             switch (nSearchDirection)
//             {
//             case NORTH:
//                nSearchDirection = WEST;
//                break;
//
//             case EAST:
//                nSearchDirection = NORTH;
//                break;
//
//             case SOUTH:
//                nSearchDirection = EAST;
//                break;
//
//             case WEST:
//                nSearchDirection = SOUTH;
//                break;
//             }
//
//             bFoundNextCell = true;
//          }
//
//          // 4. Try going back (U-turn)
//          else if (bIsWithinValidGrid(nXGoBack, nYGoBack) && m_pRasterGrid->m_Cell[nXGoBack][nYGoBack].bIsCliffToe())
//          {
//             nX = nXGoBack;
//             nY = nYGoBack;
//
//             // Update search direction (U-turn)
//             nSearchDirection = nGetOppositeDirection(nSearchDirection);
//             bFoundNextCell = true;
//          }
//
//          if (! bFoundNextCell)
//          {
//             break;   // No valid next cell found, end this cliff edge trace
//          }
//
//       } while ((nX != nXStart || nY != nYStart) && nLength < nMaxLength);
//
//       // Only keep cliff edges that have a reasonable length
//       if (VCliffEdge.size() > 2)
//       {
//          nCliffEdgesTraced++;
//
//          // Convert grid coordinates to external CRS for smoothing
//          CGeomLine CliffEdgeExtCRS;
//          for (const auto &point : VCliffEdge)
//          {
//             CliffEdgeExtCRS.Append(dGridCentroidXToExtCRSX(point.nGetX()), dGridCentroidYToExtCRSY(point.nGetY()));
//             bUsedInCliffTrace[point.nGetX()][point.nGetY()] = true;
//          }
//
//          // Apply cliff edge specific smoothing
//          if (m_nCliffEdgeSmooth == SMOOTH_RUNNING_MEAN)
//             CliffEdgeExtCRS = LSmoothCoastRunningMean(&CliffEdgeExtCRS);
//
//          else if (m_nCliffEdgeSmooth == SMOOTH_SAVITZKY_GOLAY)
//             CliffEdgeExtCRS = LSmoothCoastSavitzkyGolay(&CliffEdgeExtCRS, 0, 0);
//
//          // Store the smoothed cliff edge in external CRS
//          m_VCliffToe.push_back(CliffEdgeExtCRS);
//       }
//    }
// }

// /*===============================================================================================================================
// //! Validates cliff edges to ensure they represent cliff toes, truncating edges that trace cliff tops
// ===============================================================================================================================*/
// void CSimulation::nValidateCliffToeEdges(void)
// {
//    vector<CGeomLine> ValidatedCliffEdges;
//
//    // Process each traced cliff edge
//    for (size_t nEdge = 0; nEdge < m_VCliffToe.size(); nEdge++)
//    {
//       CGeomLine &CliffEdge = m_VCliffToe[nEdge];
//
//       // Try validating in forward direction first
//       CGeomLine const ForwardValidated = nValidateCliffToeDirection(CliffEdge, false);
//
//       // If we got a very short result (broke early), try reverse direction
//       CGeomLine ReverseValidated;
//       if (ForwardValidated.nGetSize() < CliffEdge.nGetSize() * 0.2)
//       {
//          ReverseValidated = nValidateCliffToeDirection(CliffEdge, true);
//       }
//
//       // Use whichever result is longer
//       CGeomLine BestValidated;
//       if (ReverseValidated.nGetSize() > ForwardValidated.nGetSize())
//       {
//          BestValidated = ReverseValidated;
//       }
//       else
//       {
//          BestValidated = ForwardValidated;
//       }
//
//       // Only keep edges that have a reasonable length after validation
//       if (BestValidated.nGetSize() > 2)
//       {
//          ValidatedCliffEdges.push_back(BestValidated);
//       }
//    }
//
//    // Replace the original cliff edges with the validated ones
//    m_VCliffToe = ValidatedCliffEdges;
// }

// /*===============================================================================================================================
// //! Validates cliff toe direction
// ===============================================================================================================================*/
// CGeomLine CSimulation::nValidateCliffToeDirection(CGeomLine& CliffEdge, bool bReverse)
// {
//    CGeomLine ValidatedEdge;
//
//    int nConsecutiveFailures = 0;
//    int const nMaxConsecutiveFailures = 2;
//    int const nSize = CliffEdge.nGetSize();
//
//    // Check each point along the cliff edge
//    for (int i = 0; i < nSize - 1; i++)
//    {
//       // Determine which point to process based on direction
//       int const nPoint = bReverse ? (nSize - 1 - i) : i;
//       int const nNextPoint = bReverse ? (nSize - 2 - i) : (i + 1);
//
//       // Skip if we've gone beyond bounds
//       if (nNextPoint < 0 || nNextPoint >= nSize)
//          continue;
//
//       // Get current and next point in grid coordinates
//       int nX = static_cast<int>((CliffEdge.dGetXAt(nPoint) - m_dGeoTransform[0]) / m_dGeoTransform[1]);
//       int nY = static_cast<int>((CliffEdge.dGetYAt(nPoint) - m_dGeoTransform[3]) / m_dGeoTransform[5]);
//
//       int nNextX = static_cast<int>((CliffEdge.dGetXAt(nNextPoint) - m_dGeoTransform[0]) /  m_dGeoTransform[1]);
//       int nNextY = static_cast<int>((CliffEdge.dGetYAt(nNextPoint) - m_dGeoTransform[3]) /  m_dGeoTransform[5]);
//
//       // Ensure coordinates are within grid bounds
//       nX = max(0, min(m_nXGridSize - 1, nX));
//       nY = max(0, min(m_nYGridSize - 1, nY));
//       nNextX = max(0, min(m_nXGridSize - 1, nNextX));
//       nNextY = max(0, min(m_nYGridSize - 1, nNextY));
//
//       // Calculate direction of travel
//       int const nDirX = nNextX - nX;
//       int const nDirY = nNextY - nY;
//
//       // Get perpendicular directions (90 degrees to travel direction)
//       int const nLeftX = nX - nDirY;         // Rotate direction 90 degrees counter-clockwise
//       int const nLeftY = nY + nDirX;
//       int const nRightX = nX + nDirY;        // Rotate direction 90 degrees clockwise
//       int const nRightY = nY - nDirX;
//
//       bool bIsValidToe = true;
//
//       // Check if perpendicular cells are within bounds
//       if (bIsWithinValidGrid(nLeftX, nLeftY) && bIsWithinValidGrid(nRightX, nRightY))
//       {
//          bool const bLeftIsCliff = m_pRasterGrid->m_Cell[nLeftX][nLeftY].bIsCliffToe();
//          bool const bRightIsCliff = m_pRasterGrid->m_Cell[nRightX][nRightY].bIsCliffToe();
//
//          // One should be cliff and one should be not cliff for a valid cliff edge
//          if (bLeftIsCliff != bRightIsCliff)
//          {
//             // Get the elevation of these two adjacent cells
//             double const dLeftElev = m_pRasterGrid->m_Cell[nLeftX][nLeftY].dGetAllSedTopElevOmitTalus();
//             double const dRightElev = m_pRasterGrid->m_Cell[nRightX][nRightY].dGetAllSedTopElevOmitTalus();
//
//             // Determine which is the cliff side and which is the non-cliff side
//             double dCliffElev;
//             double dNonCliffElev;
//             if (bLeftIsCliff)
//             {
//                dCliffElev = dLeftElev;
//                dNonCliffElev = dRightElev;
//             }
//             else
//             {
//                dCliffElev = dRightElev;
//                dNonCliffElev = dLeftElev;
//             }
//
//             // If the non-cliff cell is higher than the cliff cell, this is not the toe
//             if (dNonCliffElev > dCliffElev)
//             {
//                bIsValidToe = false;
//             }
//          }
//       }
//
//       if (bIsValidToe)
//       {
//          // Reset consecutive failure counter
//          nConsecutiveFailures = 0;
//
//          // Add this point to the validated edge
//          ValidatedEdge.Append(CliffEdge.dGetXAt(nPoint), CliffEdge.dGetYAt(nPoint));
//       }
//       else
//       {
//          // Increment consecutive failure counter
//          nConsecutiveFailures++;
//
//          // If we haven't hit the threshold yet, still add the point
//          if (nConsecutiveFailures < nMaxConsecutiveFailures)
//          {
//             ValidatedEdge.Append(CliffEdge.dGetXAt(nPoint), CliffEdge.dGetYAt(nPoint));
//          }
//          else
//          {
//             // Too many consecutive failures - truncate here
//             break;
//          }
//       }
//    }
//
//    return ValidatedEdge;
// }
