/*!
   \file interpolate.cpp
   \brief Returns interpolated value at x from parallel arrays
   \details TODO 001 A more detailed description of these routines.
   \author Modified by David Favis-Mortlock and Andres Payo
   \date 2026
   \copyright GNU Lesser General Public License
*/

/* ===============================================================================================================================
   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
===============================================================================================================================*/
#include <assert.h>

#include <cfloat>

#include <vector>
using std::vector;

#include "cme.h"
#include "simulation.h"

//===============================================================================================================================
//! From https://cplusplus.com/forum/general/216928/
//! Returns interpolated value at x from parallel arrays (VdXdata, VdYdata). Assumes that VdXdata has at least two elements, is sorted and is strictly monotonically increasing. The boolean argument extrapolate determines behaviour beyond ends of array (if needed). For this version, both lots of data are doubles
//===============================================================================================================================
double CSimulation::dGetInterpolatedValue(vector<double> const* pVdXdata, vector<double> const* pVdYdata, double dX, bool bExtrapolate)
{
   int const size = static_cast<int>(pVdXdata->size());

   int i = 0; // Find left end of interval for interpolation

   if (dX >= pVdXdata->at(size - 2)) // Special case: beyond right end
   {
      i = size - 2;
   }

   else
   {
      while (dX > pVdXdata->at(i + 1))
         i++;
   }

   double const dXL = pVdXdata->at(i);
   double dYL = pVdYdata->at(i);
   double const dXR = pVdXdata->at(i + 1);
   double dYR = pVdYdata->at(i + 1); // Points on either side (unless beyond ends)

   if (! bExtrapolate) // If beyond ends of array and not extrapolating
   {
      if (dX < dXL)
         dYR = dYL;

      if (dX > dXR)
         dYL = dYR;
   }

   double const ddYdX = (dYR - dYL) / (dXR - dXL); // Gradient

   return (dYL + ddYdX * (dX - dXL)); // Linear interpolation
}

//===============================================================================================================================
//! From https://cplusplus.com/forum/general/216928/
//! Returns interpolated value at x from parallel arrays (VdXdata, VdYdata). Assumes that VdXdata has at least two elements, is sorted and is strictly monotonically increasing. The boolean argument extrapolate determines behaviour beyond ends of array (if needed). For this version, one lot of data is integer and the other is double
//===============================================================================================================================
double CSimulation::dGetInterpolatedValue(vector<int> const* pVnXdata, vector<double> const* pVdYdata, int nX, bool bExtrapolate)
{
   unsigned int const nSize = static_cast<unsigned int>(pVnXdata->size());

   int i = 0; // Find left end of interval for interpolation

   if (nX >= pVnXdata->at(nSize - 2)) // Special case: beyond right end
   {
      i = nSize - 2;
   }

   else
   {
      while (nX > pVnXdata->at(i + 1))
         i++;
   }

   int const nXL = pVnXdata->at(i);
   int const nXR = pVnXdata->at(i + 1);

   double dYL = pVdYdata->at(i);
   double dYR = pVdYdata->at(i + 1); // Points on either side (unless beyond ends)

   if (! bExtrapolate) // If beyond ends of array and not extrapolating
   {
      if (nX < nXL)
         dYR = dYL;

      if (nX > nXR)
         dYL = dYR;
   }

   double const ddYdX = (dYR - dYL) / static_cast<double>(nXR - nXL); // Gradient

   return dYL + ddYdX * static_cast<double>(nX - nXL); // Linear interpolation
}

//===============================================================================================================================
//! This is used by VdInterpolateCShoreProfileOutput, it returns the index of the value in pVdX which is less than or equal to the absolute difference between dValueIn and the pVdX value
//===============================================================================================================================
int CSimulation::nFindIndex(vector<double> const* pVdX, double const dValueIn)
{
   double dLastValue = DBL_MAX;
   int nIndexFound = 0;

   for (unsigned int i = 0; i < pVdX->size(); ++i)
   {
      double const dThisValue = tAbs(dValueIn - pVdX->at(i));

      if (dThisValue <= dLastValue)
      {
         dLastValue = dThisValue;
         nIndexFound = i;
      }
   }

   return nIndexFound;
}

//===============================================================================================================================
//! Returns a linearly interpolated vector of doubles, to make CShore profile output compatible with CME. The array pVdY has been output by CShore and so always has length CSHORE_ARRAY_OUT_SIZE, whereas all other arrays have sizes which depend on CME at runtime
//===============================================================================================================================
vector<double> CSimulation::VdInterpolateCShoreProfileOutput(vector<double> const* pVdX, vector<double> const* pVdY, vector<double> const* pVdXNew)
{
   int const nXSize = static_cast<int>(pVdX->size());
   int const nXNewSize = static_cast<int>(pVdXNew->size());

   // assert(nXSize > 0);
   // assert(nXNewSize > 0);

   double dX;
   double dY;
   vector<double> VdYNew(nXNewSize, 0.0);

   for (int i = 0; i < nXNewSize; ++i)
   {
      int const idx = nFindIndex(pVdX, pVdXNew->at(i));

      if (pVdX->at(idx) > pVdXNew->at(i))
      {
         if (idx > 0)
         {
            dX = pVdX->at(idx) - pVdX->at(idx - 1);
            dY = pVdY->at(idx) - pVdY->at(idx - 1);
         }

         else
         {
            dX = pVdX->at(idx + 1) - pVdX->at(idx);
            dY = pVdY->at(idx + 1) - pVdY->at(idx);
         }
      }

      else
      {
         if (idx < nXSize - 1)
         {
            dX = pVdX->at(idx + 1) - pVdX->at(idx);
            dY = pVdY->at(idx + 1) - pVdY->at(idx);
         }

         else
         {
            dX = pVdX->at(idx) - pVdX->at(idx - 1);
            dY = pVdY->at(idx) - pVdY->at(idx - 1);
         }
      }

      // Safety check: this crashes (divide by zero) if there are identical consecutive values in pVdX, and thus if dX becomes 0. To prevent this, if dX is near zero, set to a small non-zero number
      if (bFPIsEqual(dX, 0.0, TOLERANCE))
         dX = 1e-10;

      double const dM = dY / dX;
      double const dB = pVdY->at(idx) - pVdX->at(idx) * dM;

      // VdYNew[i] = (pVdXNew->at(i) * dM) + dB;
      VdYNew[nXNewSize - 1 - i] = (pVdXNew->at(i) * dM) + dB;
   }

   return VdYNew;
}
