/*!
   \file calc_external_forcing.cpp
   \brief Calculates external forcings
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

#include "cme.h"
#include "simulation.h"

//===============================================================================================================================
//! Calculate external forcings: change in still water level, tide level and deep water waves height, orientation and period. If we are considering cliff collapse, then also calculate mean high water elevation
//===============================================================================================================================
int CSimulation::nCalcExternalForcing(void)
{
   // Increment long-term SWL change (increment may be zero)
   m_dAccumulatedSeaLevelChange += m_dDeltaSWLPerTimestep;

   // Increment this-iteration mean SWL with long-term SWL change
   m_dThisIterMeanSWL = m_dInitialMeanSWL + m_dAccumulatedSeaLevelChange;

   int const nSize = static_cast<int>(m_VdTideData.size());
   if (nSize == 0)
   {
      // No tide data
      m_dThisIterSWL = m_dThisIterMeanSWL;
      m_dThisIterMHWElev = m_dThisIterMeanSWL;

      if (m_bDoCliffCollapse)
         CalcMHWAndNewCliffNotchElevation(-1);
   }
   else
   {
      // We have tide data
      static int snTideDataCount = 0;

      // Wrap the tide data, i.e. start again with the first line of the tide data if we do not have enough
      if (snTideDataCount > nSize - 1)
         snTideDataCount = 0;

      // This-iteration SWL includes both tidal change and long-term SWL change
      m_dThisIterSWL = m_dThisIterMeanSWL + m_VdTideData[snTideDataCount];

      if (m_bDoCliffCollapse)
         CalcMHWAndNewCliffNotchElevation(snTideDataCount);

      snTideDataCount++;
   }

   LogStream << m_ulIter << ": SWL = " << m_dThisIterSWL << endl;

   m_bHighestSWLSoFar = false;
   m_bLowestSWLSoFar = false;

   // Maybe update min and max still water levels so far during this simulation
   if (m_dThisIterSWL > m_dMaxSWLSoFar)
   {
      m_bHighestSWLSoFar = true;
      m_dMaxSWLSoFar = m_dThisIterSWL;
   }

   if (m_dThisIterSWL < m_dMinSWLSoFar)
   {
      m_bLowestSWLSoFar = true;
      m_dMinSWLSoFar = m_dThisIterSWL;
   }

   // Update the wave height, orientation and period for this time step and start again with the first record if we do not have enough
   if (m_bHaveWaveStationData)
   {
      static int snWaveStationDataCount = 0;

      if (snWaveStationDataCount > m_nDeepWaterWaveDataNumTimeSteps - 1)
      {
         // Wrap the tide data, i.e. start again with the first record if we do not have enough
         snWaveStationDataCount = 0;
      }

      // Do we have just a single wave station?
      if (m_bSingleDeepWaterWaveValues)
      {
         // Yes, just a single wave station
         m_dAllCellsDeepWaterWaveHeight = m_VdTSDeepWaterWaveStationHeight[snWaveStationDataCount];
         m_dAllCellsDeepWaterWaveAngle = m_VdTSDeepWaterWaveStationAngle[snWaveStationDataCount];
         m_dAllCellsDeepWaterWavePeriod = m_VdTSDeepWaterWaveStationPeriod[snWaveStationDataCount];
      }
      else
      {
         // More than one wave station, so update this time step's deep water wave values for use in the nInterpolateAllDeepWaterWaveValues() routine. Note that the order on the vector is determined by the points ID i.e. to ensure that stations match with time series
         int const nNumberDeepWaterWaveStations = static_cast<int>(m_VnDeepWaterWaveStationID.size());
         int const nTot = nNumberDeepWaterWaveStations * snWaveStationDataCount;

         for (int j = 0; j < nNumberDeepWaterWaveStations; j++)
         {
            m_VdThisIterDeepWaterWaveStationHeight[j] = m_VdTSDeepWaterWaveStationHeight[(m_VnDeepWaterWaveStationID[j] - 1) + nTot];
            m_VdThisIterDeepWaterWaveStationAngle[j] = m_VdTSDeepWaterWaveStationAngle[(m_VnDeepWaterWaveStationID[j] - 1) + nTot];
            m_VdThisIterDeepWaterWaveStationPeriod[j] = m_VdTSDeepWaterWaveStationPeriod[(m_VnDeepWaterWaveStationID[j] - 1) + nTot];
         }
      }

      snWaveStationDataCount++;
   }

   return RTN_OK;
}
