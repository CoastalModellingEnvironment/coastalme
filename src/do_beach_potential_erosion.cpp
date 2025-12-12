/*!
   \file do_beach_potential_erosion.cpp
   \brief Calculates potential (i.e. not constrained by the availability of unconsolidated sediment) beach erosion of unconsolidated sediment on coastal polygons
   \details TODO 001 A more detailed description of these routines.
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

#include <algorithm>
using std::sort;

#include <utility>
using std::make_pair;
using std::pair;

#include "cme.h"
#include "simulation.h"
#include "coast.h"

namespace
{
//===============================================================================================================================
//! Function used to sort polygon length values. If the first argument must be ordered before the second, return true
//===============================================================================================================================
bool bPolygonLengthPairCompare(const pair<int, double>& prLeft, const pair<int, double>& prRight)
{
   // Sort in ascending order (i.e. most concave first)
   return prLeft.second < prRight.second;
}
} // namespace

//===============================================================================================================================
//! Uses either the CERC equation or the Kamphuis (1990) equation to calculate potential (unconstrained) sediment movement between polygons
//===============================================================================================================================
void CSimulation::DoAllPotentialBeachErosion(void)
{
   // Do this for each coast
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      int const nNumPolygons = m_VCoast[nCoast].nGetNumPolygons();

      // Create a vector of pairs: the first value of the pair is the profile index, the second is the seaward length of that profile
      vector<pair<int, double>> prVPolygonLength;

      for (int nPoly = 0; nPoly < nNumPolygons; nPoly++)
      {
         CGeomCoastPolygon const* pPolygon = m_VCoast[nCoast].pGetPolygon(nPoly);
         double const dSeawardLength = pPolygon->dGetLength();
         prVPolygonLength.push_back(make_pair(nPoly, dSeawardLength));
      }

      // Sort this pair vector in ascending order, so that the polygons with the shortest length (i.e. the most concave polygons) are first
      sort(prVPolygonLength.begin(), prVPolygonLength.end(), bPolygonLengthPairCompare);

      // Do this for every coastal polygon in sequence of coastline concavity
      for (int n = 0; n < nNumPolygons; n++)
      {
         int const nThisPoly = prVPolygonLength[n].first;

         CGeomCoastPolygon* pPolygon = m_VCoast[nCoast].pGetPolygon(nThisPoly);

         // Calculate the average breaking wave height and angle along this polygon's segment of coastline
         int const nStartNormal = pPolygon->nGetUpCoastProfile();
         int const nEndNormal = pPolygon->nGetDownCoastProfile();
         int const nCoastStartPoint = m_VCoast[nCoast].pGetProfile(nStartNormal)->nGetCoastPoint();
         int const nCoastEndPoint = m_VCoast[nCoast].pGetProfile(nEndNormal)->nGetCoastPoint();
         int nCoastPoints = 0;
         int nActiveZonePoints = 0;

         double dAvgBreakingWaveHeight = 0;
         double dAvgBreakingWaveAngle = 0;
         double dAvgDeepWaterWavePeriod = 0;
         double dAvgFluxOrientation = 0;
         double dAvgBreakingDepth = 0;
         double dAvgBreakingDist = 0;

         // Calculate the average tangent to the polygon's coast segment, the average wave breaking height, the average depth of breaking, and the average distance of breaking, for this coast segment
         for (int nCoastPoint = nCoastStartPoint; nCoastPoint < nCoastEndPoint - 1; nCoastPoint++)
         {
            nCoastPoints++;
            dAvgFluxOrientation += m_VCoast[nCoast].dGetFluxOrientation(nCoastPoint);

            double const dThisBreakingWaveHeight = m_VCoast[nCoast].dGetBreakingWaveHeight(nCoastPoint);

            if (! bFPIsEqual(dThisBreakingWaveHeight, DBL_NODATA, TOLERANCE))
            {
               // We are in the active zone
               nActiveZonePoints++;
               dAvgBreakingWaveHeight += dThisBreakingWaveHeight;

               double const dThisBreakingWaveAngle = m_VCoast[nCoast].dGetBreakingWaveAngle(nCoastPoint);
               double const dThisDeepWaterWavePeriod = m_VCoast[nCoast].dGetCoastDeepWaterWavePeriod(nCoastPoint);

               dAvgBreakingWaveAngle += dThisBreakingWaveAngle;
               dAvgDeepWaterWavePeriod += dThisDeepWaterWavePeriod;

               dAvgBreakingDepth += m_VCoast[nCoast].dGetDepthOfBreaking(nCoastPoint);

               dAvgBreakingDist += (m_VCoast[nCoast].nGetBreakingDistance(nCoastPoint) * m_dCellSide);
            }
         }

         // Safety check
         if (nCoastPoints == 0)
            nCoastPoints = 1;

         // Calc the averages
         dAvgFluxOrientation /= nCoastPoints;

         if (nActiveZonePoints > 0)
         {
            // Only calculate sediment movement if the polygon has at least one coast point in its coastline segment which is in the active zone
            dAvgBreakingWaveHeight /= nActiveZonePoints;
            dAvgBreakingWaveAngle /= nActiveZonePoints;
            dAvgDeepWaterWavePeriod /= nActiveZonePoints;
            dAvgBreakingDepth /= nActiveZonePoints;
            dAvgBreakingDist /= nActiveZonePoints;

            // Get the coast handedness, and (based on the average tangent) calculate the direction towards which a coastline-normal profile points
            int const nSeaHand = m_VCoast[nCoast].nGetSeaHandedness();
            double dNormalOrientation;

            if (nSeaHand == RIGHT_HANDED)
               dNormalOrientation = dKeepWithin360(dAvgFluxOrientation - 90);
            else
               dNormalOrientation = dKeepWithin360(dAvgFluxOrientation + 90);

            // Determine dThetaBr, the angle between the coastline-normal orientation and the breaking wave orientation (the direction FROM which the waves move). This tells us whether the sediment movement is up-coast (-ve) or down-coast (+ve)
            double dThetaBr = dNormalOrientation - dAvgBreakingWaveAngle;

            if (dThetaBr > 270)
               dThetaBr = dAvgBreakingWaveAngle + 360.0 - dNormalOrientation;
            else if (dThetaBr < -270)
               dThetaBr = dNormalOrientation + 360.0 - dAvgBreakingWaveAngle;

            bool bDownCoast = true;

            if (dThetaBr < 0)
               bDownCoast = false;

            // And save the direction of sediment movement in the polygon object
            pPolygon->SetDownCoastThisIter(bDownCoast);

            // Now that we have the direction of sediment movement, normalize dThetaBr to be always +ve so that subsequent calculations are clearer
            dThetaBr = tAbs(dThetaBr);

            // Safety check: not sure why we need this, but get occasional big values
            dThetaBr = tMin(dThetaBr, 90.0);

            // Calculate the immersed weight of sediment transport
            double dImmersedWeightTransport = 0;

            if (m_nBeachErosionDepositionEquation == UNCONS_SEDIMENT_EQUATION_CERC)
            {
               // Use the CERC equation (Komar and Inman, 1970; USACE, 1984), this describes the immersive weight transport of sand (i.e. sand transport in suspension). Depth-integrated alongshore volumetric sediment transport is a function of breaking wave height Hb and angle αb:
               //
               //    Qls = Kls * Hb^(5/2) * sin(2 * αb)
               //
               // where Kls is a transport coefficient which varies between 0.4 to 0.79

               dImmersedWeightTransport = m_dKLS / (16 * pow(m_dBreakingWaveHeightDepthRatio, 0.5)) * m_dSeaWaterDensity * pow(m_dG, 1.5) * pow(dAvgBreakingWaveHeight, 2.5) * sin((PI / 180) * 2 * dThetaBr);
            }
            else if (m_nBeachErosionDepositionEquation == UNCONS_SEDIMENT_EQUATION_KAMPHUIS)
            {
               // Use the Kamphuis (1990) equation to estimate the immersive weight transport of sand in kg/s:
               //
               //    Qls = 2.33 * (Tp^(1.5)) * (tanBeta^(0.75)) * (d50^(-0.25)) * (Hb^2) * (sin(2 * αb)^(0.6))
               //
               // where:
               // Tp = peak wave period
               // tanBeta = beach slope, defined as the ratio of the water depth at the breaker line and the distance from the still water beach line to the breaker line
               // d50 = median particle size in surf zone (m)

               if (dAvgBreakingDist > 0)
               {
                  double const dD50 = pPolygon->dGetAvgUnconsD50();

                  if (dD50 > 0)
                  {
                     double const dBeachSlope = dAvgBreakingDepth / dAvgBreakingDist;

                     // Note that we use a calibration constant here (m_dKamphuis)
                     dImmersedWeightTransport = m_dKamphuis * 2.33 * pow(dAvgDeepWaterWavePeriod, 1.5) * pow(dBeachSlope, 0.75) * pow(dD50, -0.25) * pow(dAvgBreakingWaveHeight, 2) * pow(sin((PI / 180) * 2 * dThetaBr), 0.6);
                  }
               }
            }

            // Convert from immersed weight rate to bulk volumetric (sand and voids) transport rate (m3/s)
            double dSedimentVol = m_dInmersedToBulkVolumetric * dImmersedWeightTransport;

            // Convert to a volume during this timestep
            dSedimentVol *= (m_dTimeStep * 3600);

            // Convert to a depth in m
            double dSedimentDepth = dSedimentVol / m_dCellArea;

            // If this is just a tiny depth, do nothing
            if (dSedimentDepth < SED_ELEV_TOLERANCE)
               dSedimentDepth = 0;

            // LogStream << m_ulIter << ": polygon = " << nThisPoly << " nActiveZonePoints = " << nActiveZonePoints << " dAvgBreakingWaveHeight = " << dAvgBreakingWaveHeight << " dAvgFluxOrientation = " << dAvgFluxOrientation << " dNormalOrientation = " << dNormalOrientation << " dAvgBreakingWaveAngle = " << dAvgBreakingWaveAngle <<  " potential sediment transport this timestep = " << dSedimentDepth << " m " << (bDownCoast ? "DOWN" : "UP") << " coast" << endl;

            // Store the potential erosion value for this polygon
            pPolygon->AddPotentialErosion(-dSedimentDepth);
            // LogStream << "\tPotential erosion on polygon " << nThisPoly << " -dSedimentDepth = " << -dSedimentDepth << endl;
         }

         // else
         // LogStream << m_ulIter << ": polygon = " << nThisPoly << " NOT IN ACTIVE ZONE dAvgFluxOrientation = " << dAvgFluxOrientation << endl;
      }
   }
}
