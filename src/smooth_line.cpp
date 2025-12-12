/*!
   \file smooth_line.cpp
   \brief Smoothing routines for CGeomLine objects
   \details The Savitzky-Golay routines are modified from C originals by Jean-Pierre Moreau (jpmoreau@wanadoo.fr, http://jean-pierre.moreau.pagesperso-orange.fr/index.html), to whom we are much indebted
   \author David Favis-Mortlock
   \author Andres Payo
   \author Wilf Chun
   \date 2025
   \copyright GNU General Public License
*/

/* ==============================================================================================================================
   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public  License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
==============================================================================================================================*/
#include <assert.h>

#include <cmath>
using std::pow;

#include <iostream>
using std::cerr;
using std::endl;
using std::ios;

#include "cme.h"
#include "simulation.h"
#include "2d_point.h"
#include "2di_point.h"
#include "line.h"

// Visible throughout this file. Note that arrays here are used from index 1
typedef double Matrix[SAVGOL_POLYNOMIAL_MAX_ORDER + 2][SAVGOL_POLYNOMIAL_MAX_ORDER + 2];

namespace
{
void LUDecomp(Matrix, int const, int const, int[], int*, int*);
void LULinearSolve(Matrix const, int const, int const[], double[]);
} // namespace

//===============================================================================================================================
//! Calculates the Savitzky-Golay smoothing coefficients for a given size of smoothing window. Derived from a C original by Jean-Pierre Moreau (jpmoreau@wanadoo.fr, http://jean-pierre.moreau.pagesperso-orange.fr/index.html), to whom we are much indebted
//===============================================================================================================================
void CSimulation::CalcSavitzkyGolayCoeffs(void)
{
   m_VnSavGolIndexCoast.resize(m_nCoastSmoothingWindowSize + 1, 0);

   // Note that m_nCoastSmoothingWindowSize must be odd (have already checked this)
   int const nHalfWindow = m_nCoastSmoothingWindowSize / 2;

   // Calculate the shift index for this value of nHalfWindow
   int j = 3;

   for (int i = 2; i <= nHalfWindow + 1; i++)
   {
      m_VnSavGolIndexCoast[i] = i - j;
      j += 2;
   }

   j = 2;

   for (int i = nHalfWindow + 2; i <= m_nCoastSmoothingWindowSize; i++)
   {
      m_VnSavGolIndexCoast[i] = i - j;
      j += 2;
   }

   // Now calculate the Savitzky-Golay filter coefficients
   m_VdSavGolFCRWCoast.resize(m_nCoastSmoothingWindowSize + 1, 0);

   CalcSavitzkyGolay(&(m_VdSavGolFCRWCoast.at(0)), m_nCoastSmoothingWindowSize, nHalfWindow, nHalfWindow, 0, m_nSavGolCoastPoly);
}

//===============================================================================================================================
//! Does smoothing of a CGeomLine coastline vector (is in external CRS coordinates) using a Savitzky-Golay filter. Derived from a C original by Jean-Pierre Moreau (jpmoreau@wanadoo.fr, http://jean-pierre.moreau.pagesperso-orange.fr/index.html), to whom we are much indebted
//===============================================================================================================================
CGeomLine CSimulation::LSmoothCoastSavitzkyGolay(CGeomLine* pLineIn, int const nStartEdge, int const nEndEdge) const
{
   // Note that m_nCoastSmoothingWindowSize must be odd (have already checked this)
   int const nHalfWindow = m_nCoastSmoothingWindowSize / 2;

   // Make a copy of the unsmoothed CGeomLine (must be blank)
   int const nSize = pLineIn->nGetSize();
   CGeomLine LTemp;
   LTemp.Resize(nSize);

   // Apply the Savitzky-Golay smoothing filter
   for (int i = 0; i < nSize; i++)
   {
      CGeom2DPoint const PtThis(pLineIn->dGetXAt(i), pLineIn->dGetYAt(i));
      CGeom2DIPoint const PtiThis = PtiExtCRSToGridRound(&PtThis);
      int nXThis = PtiThis.nGetX();
      int nYThis = PtiThis.nGetY();

      // Safety check
      if (! bIsWithinValidGrid(nXThis, nYThis))
         KeepWithinValidGrid(nXThis, nYThis);

      // Don't smooth intervention cells
      if (bIsInterventionCell(nXThis, nYThis))
      {
         LTemp[i] = pLineIn->pPtGetAt(i);
         continue;
      }

      if (i < nHalfWindow)
      {
         // For the first few values of LTemp, just apply a running mean with a variable-sized window
         int nTmpWindow = 0;
         double dWindowTotX = 0, dWindowTotY = 0;

         for (int j = -nHalfWindow; j < m_nCoastSmoothingWindowSize - nHalfWindow; j++)
         {
            int const k = i + j;

            if ((k > 0) && (k < nSize))
            {
               dWindowTotX += pLineIn->dGetXAt(k);
               dWindowTotY += pLineIn->dGetYAt(k);
               nTmpWindow++;
            }
         }

         switch (nStartEdge)
         {
         case NORTH:
         case SOUTH:
            // Don't apply the filter in the Y direction
            LTemp[i] = CGeom2DPoint(dWindowTotX / nTmpWindow, pLineIn->dGetYAt(i));
            // LTemp.SetXAt(i, dWindowTotX / static_cast<double>(nTmpWindow));
            // LTemp.SetYAt(i, pLineIn->dGetYAt(i));
            break;

         case EAST:
         case WEST:
            // Don't apply the filter in the X direction
            LTemp[i] = CGeom2DPoint(pLineIn->dGetXAt(i), dWindowTotY / nTmpWindow);
            // LTemp.SetXAt(i, pLineIn->dGetXAt(i));
            // LTemp.SetYAt(i, dWindowTotY / static_cast<double>(nTmpWindow));
            break;
         }
      }

      else if (i >= (nSize - nHalfWindow))
      {
         // For the last few values of PtVTemp, just apply a running mean with a variable-sized window
         int nTmpWindow = 0;
         double dWindowTotX = 0, dWindowTotY = 0;

         for (int j = -nHalfWindow; j < m_nCoastSmoothingWindowSize - nHalfWindow; j++)
         {
            int const k = i + j;

            if ((k > 0) && (k < nSize))
            {
               dWindowTotX += pLineIn->dGetXAt(k);
               dWindowTotY += pLineIn->dGetYAt(k);
               nTmpWindow++;
            }
         }

         switch (nEndEdge)
         {
         case NORTH:
         case SOUTH:
            // Don't apply the filter in the Y direction
            LTemp[i] = CGeom2DPoint(dWindowTotX / nTmpWindow, pLineIn->dGetYAt(i));
            // LTemp.SetXAt(i, dWindowTotX / static_cast<double>(nTmpWindow));
            // LTemp.SetYAt(i, pLineIn->dGetYAt(i));
            break;

         case EAST:
         case WEST:
            // Don't apply the filter in the X direction
            LTemp[i] = CGeom2DPoint(pLineIn->dGetXAt(i), dWindowTotY / nTmpWindow);
            // LTemp.SetXAt(i, pLineIn->dGetXAt(i));
            // LTemp.SetYAt(i, dWindowTotY / static_cast<double>(nTmpWindow));
            break;
         }
      }

      else
      {
         // For all other PtVTemp values, calc Savitzky-Golay weighted values for both X and Y
         for (int j = 0; j < m_nCoastSmoothingWindowSize; j++)
         {
            int const k = i + m_VnSavGolIndexCoast[j + 1];

            if ((k >= 0) && (k < nSize)) // Skip points that do not exist, note starts from 1
            {
               double dX = LTemp.dGetXAt(i);
               dX += m_VdSavGolFCRWCoast[j + 1] * pLineIn->dGetXAt(k);
               // LTemp.SetXAt(i, dX);

               double dY = LTemp.dGetYAt(i);
               dY += m_VdSavGolFCRWCoast[j + 1] * pLineIn->dGetYAt(k);

               LTemp[i] = CGeom2DPoint(dX, dY);
               // LTemp.SetYAt(i, dY);
            }
         }
      }
   }

   // Return the smoothed CGeomLine
   return LTemp;
}

//===============================================================================================================================
//! Does running-mean smoothing of a CGeomLine coastline vector (is in external CRS coordinates)
//===============================================================================================================================
CGeomLine CSimulation::LSmoothCoastRunningMean(CGeomLine* pLineIn) const
{
   // Note that m_nCoastSmoothingWindowSize must be odd (have already checked this)
   int const nHalfWindow = m_nCoastSmoothingWindowSize / 2;
   double const dHalfWindow = nHalfWindow;

   // Make a copy of the unsmoothed CGeomLine
   int const nSize = pLineIn->nGetSize();
   CGeomLine LTemp;
   LTemp = *pLineIn;

   // Apply the running mean smoothing filter, with a variable window size at both ends of the line
   for (int i = 0; i < nSize; i++)
   {
      // Don't smooth intervention cells
      CGeom2DPoint const PtThis(pLineIn->dGetXAt(i), pLineIn->dGetYAt(i));
      CGeom2DIPoint const PtiThis = PtiExtCRSToGridRound(&PtThis);
      int nXThis = tMin(PtiThis.nGetX(), m_nXGridSize - 1);
      int nYThis = tMin(PtiThis.nGetY(), m_nYGridSize - 1);

      // Safety checks
      nXThis = tMax(nXThis, 0);
      nYThis = tMax(nYThis, 0);

      if (bIsInterventionCell(nXThis, nYThis))
      {
         LTemp[i] = pLineIn->pPtGetAt(i);
         continue;
      }

      // bool bNearStartEdge = false, bNearEndEdge = false;
      // int consTant = 0;
      double nTmpWindow = 0;
      double dWindowTotX = 0, dWindowTotY = 0;

      if (i < nHalfWindow)
      {
         for (int j = 0; j <= i; j++)
         {
            // // For points at both ends of the coastline, use a smaller window
            double const weight = (dHalfWindow - tAbs(i - j)) / dHalfWindow;
            dWindowTotX += pLineIn->dGetXAt(j) * weight;
            dWindowTotY += pLineIn->dGetYAt(j) * weight;
            nTmpWindow += weight;
         }
      }

      else if (i >= nSize - nHalfWindow)
      {
         for (int j = nSize - 1; j >= i; j--)
         {
            double const weight = (dHalfWindow - tAbs(i - j)) / dHalfWindow;
            dWindowTotX += pLineIn->dGetXAt(j) * weight;
            dWindowTotY += pLineIn->dGetYAt(j) * weight;
            nTmpWindow += weight;
         }
      } // namespace name

      else
      {
         for (int j = i - nHalfWindow; j < i + nHalfWindow; j++)
         {
            double const weight = (dHalfWindow - tAbs(i - j)) / dHalfWindow;
            dWindowTotX += pLineIn->dGetXAt(j) * weight;
            dWindowTotY += pLineIn->dGetYAt(j) * weight;
            nTmpWindow += weight;
         }
      }

      // if (bNearStartEdge)
      // {
      //    // We are near the start edge
      // switch (nStartEdge)
      // {
      // case NORTH:
      //       // Don't apply the filter in the y direction
      // LTemp.SetXAt(i, dWindowTotX / static_cast<double>(nTmpWindow));
      // break;
      // case SOUTH:
      //       // Don't apply the filter in the y direction
      // LTemp.SetXAt(i, dWindowTotX / static_cast<double>(nTmpWindow));
      // break;

      // case EAST:
      //       // Don't apply the filter in the x direction
      // LTemp.SetYAt(i, dWindowTotY / static_cast<double>(nTmpWindow));
      // break;
      // case WEST:
      //       // Don't apply the filter in the x direction
      // LTemp.SetYAt(i, dWindowTotY / static_cast<double>(nTmpWindow));
      // break;
      // }
      // }
      // else if (bNearEndEdge)
      // {
      //    // We are near the end edge
      // switch (nEndEdge)
      // {
      // case NORTH:
      //       // Don't apply the filter in the y direction
      // LTemp.SetXAt(i, dWindowTotX / static_cast<double>(nTmpWindow));
      // break;
      // case SOUTH:
      //       // Don't apply the filter in the y direction
      // LTemp.SetXAt(i, dWindowTotX / static_cast<double>(nTmpWindow));
      // break;

      // case EAST:
      //       // Don't apply the filter in the x direction
      // LTemp.SetYAt(i, dWindowTotY / static_cast<double>(nTmpWindow));
      // break;
      // case WEST:
      //       // Don't apply the filter in the x direction
      // LTemp.SetYAt(i, dWindowTotY / static_cast<double>(nTmpWindow));
      // break;
      // }
      // }
      // else
      // {
      //    // Not near any edge, apply both x and y filters
      LTemp[i] = CGeom2DPoint(dWindowTotX / nTmpWindow, dWindowTotY / nTmpWindow);
      // LTemp.SetXAt(i, dWindowTotX / static_cast<double>(nTmpWindow));
      // LTemp.SetYAt(i, dWindowTotY / static_cast<double>(nTmpWindow));
      // }
   }

   // Return the smoothed CGeomLine
   return LTemp;
}

//===============================================================================================================================
//! Does running-mean smoothing of the slope of a coastline-normal profile
//===============================================================================================================================
vector<double> CSimulation::dVSmoothProfileSlope(vector<double>* pdVSlope) const
{
   // Make a copy of the unsmoothed profile slope vector
   int const nSize = static_cast<int>(pdVSlope->size());
   vector<double> dVSmoothed = *pdVSlope;

   // Note that m_nProfileSmoothWindow must be odd (have already checked this)
   int const nHalfWindow = m_nProfileSmoothWindow / 2;

   // Apply the running mean smoothing filter, with a variable window size at both ends of the line
   for (int i = 0; i < nSize; i++)
   {
      int nTmpWindow = 0;
      double dWindowTot = 0;

      for (int j = -nHalfWindow; j < m_nProfileSmoothWindow - nHalfWindow; j++)
      {
         // For points at both ends of the profile, use a smaller window
         int const k = i + j;

         if ((k < 0) || (k >= nSize))
            continue;

         dWindowTot += pdVSlope->at(k);
         nTmpWindow++;
      }

      dVSmoothed[i] = dWindowTot / static_cast<double>(nTmpWindow);

      // If necessary, constrain the slope as in SCAPE
      if (dVSmoothed[i] >= 0)
         dVSmoothed[i] = tMin(dVSmoothed[i], m_dProfileMaxSlope);

      else
         dVSmoothed[i] = tMax(dVSmoothed[i], -m_dProfileMaxSlope);
   }

   // Return the smoothed vector
   return dVSmoothed;
}

//===============================================================================================================================
//! CalcSavitzkyGolay uses LULinearSolve and LUDecomp. It returns dFilterCoeffsArray[nWindowSize], a set of Savitzky-Golay filter coefficients in wrap-around order. nLeft is the number of leftward data points used, nRight is the number of rightward data points, making the total number of data points used = nLeft + nRight + 1. nDerivOrder is the order of the derivative desired (e.g. nDerivOrder = 0 for smoothed function), nSmoothPolyOrder is the order of the smoothing polynomial also equal to the highest conserved moment; usual values are nSmoothPolyOrder = 2 or nSmoothPolyOrder = 4. Derived from a C original by Jean-Pierre Moreau (jpmoreau@wanadoo.fr, http://jean-pierre.moreau.pagesperso-orange.fr/index.html), to whom we are much indebted
//===============================================================================================================================
void CSimulation::CalcSavitzkyGolay(double dFilterCoeffsArray[], int const nWindowSize, int const nLeft, int const nRight, int const nDerivOrder, int const nSmoothPolyOrder)
{
   // Some sanity checks
   if ((nWindowSize < nLeft + nRight + 1) || (nLeft < 0) || (nRight < 0) || (nDerivOrder > nSmoothPolyOrder) || (nSmoothPolyOrder > static_cast<int>(SAVGOL_POLYNOMIAL_MAX_ORDER)) || (nLeft + nRight < nSmoothPolyOrder))
   {
      cerr << ERR << "Error in arguments to CalcSavitzkyGolay" << endl;
      return;
   }

   // Initialise arrays (including index 0 for tidiness)
   int nIndexArray[SAVGOL_POLYNOMIAL_MAX_ORDER + 2];
   Matrix dMatrix;
   double dSolutionArray[SAVGOL_POLYNOMIAL_MAX_ORDER + 2];

   for (int i = 0; i <= SAVGOL_POLYNOMIAL_MAX_ORDER + 1; i++)
   {
      for (int j = 0; j <= SAVGOL_POLYNOMIAL_MAX_ORDER + 1; j++)
         dMatrix[i][j] = 0;

      dSolutionArray[i] = 0;
      nIndexArray[i] = 0;
   }

   for (int ipj = 0; ipj <= 2 * nSmoothPolyOrder; ipj++)
   {
      // Set up the equations for the desired least squares fit
      double dSum = 0;

      if (ipj == 0)
         dSum = 1;

      for (int k = 1; k <= nRight; k++)
         dSum += pow(k, ipj);

      for (int k = 1; k <= nLeft; k++)
         dSum += pow(-k, ipj);

      int const mm = tMin(ipj, 2 * nSmoothPolyOrder - ipj);
      int imj = -mm;

      do
      {
         dMatrix[1 + (ipj + imj) / 2][1 + (ipj - imj) / 2] = dSum;
         imj += 2;
      } while (imj <= mm);
   }

   // Solve them using LU decomposition
   int nDCode = 0, nICode = 0;
   LUDecomp(dMatrix, nSmoothPolyOrder + 1, SAVGOL_POLYNOMIAL_MAX_ORDER + 1, nIndexArray, &nDCode, &nICode);

   // Right-hand side vector is unit vector, depending on which derivative we want
   dSolutionArray[nDerivOrder + 1] = 1;

   // Backsubstitute, giving one row of the inverse matrix
   LULinearSolve(dMatrix, nSmoothPolyOrder + 1, nIndexArray, dSolutionArray);

   // Zero the output array (it may be bigger than the number of coefficients). Again include index 0 for tidiness
   // for (int n = 0; n < SAVGOL_POLYNOMIAL_MAX_ORDER+1; n++)
   // for (int n = 0; n <= nWindowSize; n++)
   // dFilterCoeffsArray[n] = 0;

   for (int n = -nLeft; n <= nRight; n++)
   {
      // Each Savitzky-Golay coefficient is the dot product of powers of an integer with the inverse matrix row
      double dSum = dSolutionArray[1];
      double dFac = 1;

      for (int m = 1; m <= nSmoothPolyOrder; m++)
      {
         dFac *= n;
         dSum += dSolutionArray[m + 1] * dFac;
      }

      // Store in wrap-around order
      int const nInd = ((nWindowSize - n) % nWindowSize) + 1;
      dFilterCoeffsArray[nInd] = dSum;
   }
}

namespace
{
//===============================================================================================================================
//! Given an N x N matrix A, this routine replaces it by the LU decomposition of a rowwise permutation of itself. A and N are input. nIndexArray is an output vector which records the row permutation effected by the partial pivoting; D is output as -1 or 1, depending on whether the number of row interchanges was even or odd, respectively. This routine is used in combination with LULinearSolve to solve linear equations or to invert a matrix. Returns with nICode = 1 if matrix is singular. Derived from a C original by Jean-Pierre Moreau (jpmoreau@wanadoo.fr, http://jean-pierre.moreau.pagesperso-orange.fr/index.html), to whom we are much indebted
//===============================================================================================================================
void LUDecomp(Matrix A, int const N, int const np, int nIndexArray[], int* nDCode, int* nICode)
{
   if (N >= np)
   {
      // cerr << ERR << "in LUDecomp" << endl;
      return;
   }

   double const TINY = 1e-12;
   double AMAX, DUM, SUM;
   double* VV = new double[np];
   int I, IMAX = 0, J, K;

   *nDCode = 1;
   *nICode = 0;

   for (I = 1; I <= N; I++)
   {
      AMAX = 0.0;

      for (J = 1; J <= N; J++)
         if (tAbs(A[I][J]) > AMAX)
            AMAX = tAbs(A[I][J]);

      if (AMAX < TINY)
      {
         *nICode = 1;
         delete[] VV;
         return;
      }

      VV[I] = 1.0 / AMAX;
   }

   for (J = 1; J <= N; J++)
   {
      for (I = 1; I < J; I++)
      {
         SUM = A[I][J];

         for (K = 1; K < I; K++)
            SUM -= A[I][K] * A[K][J];

         A[I][J] = SUM;
      }

      AMAX = 0.0;

      for (I = J; I <= N; I++)
      {
         SUM = A[I][J];

         for (K = 1; K < J; K++)
            SUM -= A[I][K] * A[K][J];

         A[I][J] = SUM;
         DUM = VV[I] * tAbs(SUM);

         if (DUM >= AMAX)
         {
            IMAX = I;
            AMAX = DUM;
         }
      }

      if (J != IMAX)
      {
         for (K = 1; K <= N; K++)
         {
            DUM = A[IMAX][K];
            A[IMAX][K] = A[J][K];
            A[J][K] = DUM;
         }

         *nDCode = -(*nDCode);
         VV[IMAX] = VV[J];
      }

      nIndexArray[J] = IMAX;

      if (tAbs(A[J][J]) < TINY)
         A[J][J] = TINY;

      if (J != N)
      {
         DUM = 1.0 / A[J][J];

         for (I = J + 1; I <= N; I++)
            A[I][J] *= DUM;
      }
   }

   delete[] VV;
}
} // namespace

namespace
{
//===============================================================================================================================
//! Solves the set of N linear equations A . X = B.  Here A is input, not as the matrix A but rather as its LU decomposition, determined by the routine LUDecomp. nIndexArray is input as the permutation vector returned by LUDecomp. B is input as the right-hand side vector B, and returns with the solution vector X. A, N and nIndexArray are not modified by this routine and can be used for successive calls with different right-hand sides. This routine is also efficient for plain matrix inversion. Derived from a C original by Jean-Pierre Moreau (jpmoreau@wanadoo.fr, http://jean-pierre.moreau.pagesperso-orange.fr/index.html), to whom we are much indebted
//===============================================================================================================================
void LULinearSolve(Matrix const A, int const N, int const nIndexArray[], double B[])
{
   int II = 0;
   double SUM;

   for (int I = 1; I <= N; I++)
   {
      int const LL = nIndexArray[I];
      SUM = B[LL];
      B[LL] = B[I];

      if (II != 0)
         for (int J = II; J < I; J++)
            SUM -= A[I][J] * B[J];

      else if (! bFPIsEqual(SUM, 0.0, TOLERANCE))
         II = I;

      B[I] = SUM;
   }

   for (int I = N; I > 0; I--)
   {
      SUM = B[I];

      if (I < N)
         for (int J = I + 1; J <= N; J++)
            SUM -= A[I][J] * B[J];

      B[I] = SUM / A[I][I];
   }
}
} // namespace
