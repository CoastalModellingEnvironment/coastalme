/*!
   \file create_profiles.cpp
   \brief Creates profiles which are approximately normal to the coastline, these will become inter-polygon boundaries
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

#include <cstdio>
#include <cmath>

#include <iostream>
using std::cerr;
using std::endl;
using std::ios;

#include <algorithm>
using std::find;
using std::sort;

#include <utility>
using std::make_pair;
using std::pair;

#include <random>
using std::normal_distribution;

// #include <string>
// using std::to_string;

#include "cme.h"
#include "simulation.h"
#include "coast.h"
#include "2d_point.h"
#include "2di_point.h"

namespace
{
//===============================================================================================================================
//! Helper function used when sorting unsigned coastline curvature values, to locate start points of normal profiles. If the first argument must be ordered before the second, return true
//===============================================================================================================================
bool bCurvaturePairCompareDescending(const pair<int, double>& prLeft, const pair<int, double>& prRight)
{
   // Sort curvature (low values are straight, high values are curved)
   return prLeft.second < prRight.second;
}
} // namespace

//===============================================================================================================================
//! Create coastline-normal profiles for all coastlines. The first profiles are created 'around' the most concave bits of coast. Also create 'special' profiles at the start and end of the coast, and put these onto the raster grid
//===============================================================================================================================
int CSimulation::nCreateAllProfiles(void)
{
   if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
      LogStream << endl << m_ulIter << ": Creating profiles" << endl;

   for (unsigned int nCoast = 0; nCoast < m_VCoast.size(); nCoast++)
   {
      int nProfile = 0;
      int const nCoastSize = m_VCoast[nCoast].nGetCoastlineSize();

      // Create a bool vector to mark coast points which have been searched
      vector<bool> bVCoastPointDone(nCoastSize, false);

      // Now create a vector of pairs: the first value of the pair is the coastline point, the second is the coastline's curvature at that point
      vector<pair<int, double>> prVCurvature;

      for (int nCoastPoint = 0; nCoastPoint < nCoastSize; nCoastPoint++)
      {
         double dCurvature;

         int const nCat = m_VCoast[nCoast].pGetCoastLandform(nCoastPoint)->nGetLandFormCategory();
         if ((nCat == LF_INTERVENTION_STRUCT) || (nCat == LF_INTERVENTION_NON_STRUCT))
         {
            // This is an intervention coast point, which is likely to have some sharp angles. So store the detailed curvature
            dCurvature = m_VCoast[nCoast].dGetDetailedCurvature(nCoastPoint);
         }
         else
         {
            // Not an intervention coast point, so store the smoothed curvature
            dCurvature = m_VCoast[nCoast].dGetSmoothCurvature(nCoastPoint);
         }

         prVCurvature.push_back(make_pair(nCoastPoint, dCurvature));
      }

      // Sort this pair vector in descending order, so that the most concave and convex curvature points are first
      sort(prVCurvature.begin(), prVCurvature.end(), bCurvaturePairCompareDescending);

      // // DEBUG CODE =======================================================================================================================
      // for (int n = 0; n < prVCurvature.size(); n++)
      //    LogStream << prVCurvature[n].first << "\t" << prVCurvature[n].second << endl;
      // LogStream << endl << endl;
      // // DEBUG CODE =======================================================================================================================

      // And mark points at and near the start and end of the coastline so that they don't get searched (will be creating 'special' start- and end-of-coast profiles at these end points later)
      for (int n = 0; n < m_nCoastProfileSpacing; n++)
      {
         if (n < nCoastSize)
            bVCoastPointDone[n] = true;

         int const m = nCoastSize - n - 1;

         if (m >= 0)
            bVCoastPointDone[m] = true;
      }

      // Now locate the start points for all coastline-normal profiles (except the grid-edge ones), at points of maximum convexity. Then create the profiles
      LocateAndCreateProfiles(nCoast, nProfile, &bVCoastPointDone, &prVCurvature);

      // Did we fail to create any coast-normal profiles? If so, quit
      if (nProfile == 0)
      {
         // If this is the only coastline, then we have a problem
         if (nCoast == 0)
         {
            string strErr = ERR + "timestep " + strDblToStr(m_ulIter) + ": could not create coast-normal profiles for coastline " + strDblToStr(nCoast);

            // This is the only coastline
            if (m_ulIter == 1)
               strErr += ". Check the SWL";

            strErr += "\n";
            cerr << strErr;
            LogStream << strErr;

            return RTN_ERR_NO_PROFILES_1;
         }
         else
         {
            // This is not the only coastline, so attempt to continue
            string const strErr = WARN + "timestep " + strDblToStr(m_ulIter) + ": could not create coast-normal profiles for coastline " + strDblToStr(nCoast) + " continuing however\n";
            cerr << strErr;
            LogStream << strErr;

            continue;
         }
      }

      // Locate and create a 'special' profile at the grid edge, first at the beginning of the coastline. Then put this onto the raster grid
      int nRet = nLocateAndCreateGridEdgeProfile(true, nCoast, nProfile);
      if (nRet != RTN_OK)
         return nRet;

      // Locate a second 'special' profile at the grid edge, this time at end of the coastline. Then put this onto the raster grid
      nRet = nLocateAndCreateGridEdgeProfile(false, nCoast, ++nProfile);
      if (nRet != RTN_OK)
         return nRet;

      // Insert pointers to profiles at coastline points in the profile-all-coastpoint index
      m_VCoast[nCoast].InsertProfilesInProfileCoastPointIndex();

      // // DEBUG CODE ===================================================================================================
      // LogStream << endl << "===========================================================================================" << endl;
      // LogStream << "PROFILES BEFORE ADDING BEFORE- AND AFTER-PROFILE NUMBERS" << endl;
      // int nNumProfiles = m_VCoast[nCoast].nGetNumProfiles();
      // for (int nn = 0; nn < nNumProfiles; nn++)
      // {
      // CGeomProfile* pProfile = m_VCoast[nCoast].pGetProfile(nn);
      //
      // LogStream << nn << " nCoastID = " << pProfile->nGetProfileID() << " nGlobalID = " << pProfile->nGetProfileID() << " nGetCoastPoint = " << pProfile->nGetCoastPoint() << " pGetUpCoastAdjacentProfile = " << pProfile->pGetUpCoastAdjacentProfile() << " pGetDownCoastAdjacentProfile = " << pProfile->pGetDownCoastAdjacentProfile() << endl;
      // }
      // LogStream << "===================================================================================================" << endl << endl;
      // // DEBUG CODE ===================================================================================================

      CGeomProfile* pLastProfile;
      CGeomProfile* pThisProfile;

      // Go along the coastline and give each profile the number of the adjacent up-coast profile and the adjacent down-coast profile
      for (int nCoastPoint = 0; nCoastPoint < nCoastSize; nCoastPoint++)
      {
         if (m_VCoast[nCoast].bIsProfileAtCoastPoint(nCoastPoint))
         {
            pThisProfile = m_VCoast[nCoast].pGetProfileAtCoastPoint(nCoastPoint);

            if (nCoastPoint == 0)
            {
               // This is the profile at the upcoast end of the coast, so set this profile's up-coast profile to NULL
               pThisProfile->SetUpCoastAdjacentProfile(NULL);

               pLastProfile = pThisProfile;

               continue;
            }

            pLastProfile->SetDownCoastAdjacentProfile(pThisProfile);
            pThisProfile->SetUpCoastAdjacentProfile(pLastProfile);

            if (nCoastPoint == nCoastSize - 1)
            {
               // This is the profile at the downcoast end of the coast, so set this profile's down-coast profile to NULL
               pThisProfile->SetDownCoastAdjacentProfile(NULL);
            }

            pLastProfile = pThisProfile;
         }
      }

      // // DEBUG CODE ============================
      // for (int nCoastPoint = 0; nCoastPoint < nCoastSize; nCoastPoint++)
      // {
      //    if (m_VCoast[nCoast].bIsProfileAtCoastPoint(nCoastPoint))
      //    {
      //       pThisProfile = m_VCoast[nCoast].pGetProfileAtCoastPoint(nCoastPoint);
      //
      //       LogStream << m_ulIter << ":\t NULL PROFILE CHECK nCoastPoint = " << nCoastPoint << " nCoastSize = " << nCoastSize << " this profile = " << pThisProfile->nGetProfileID() << " upcoast profile = ";
      //       if (pThisProfile->pGetUpCoastAdjacentProfile() == NULL)
      //          LogStream << "NULL";
      //       else
      //          LogStream << pThisProfile->pGetUpCoastAdjacentProfile()->nGetProfileID();
      //       LogStream << " downcoast profile = ";
      //       if (pThisProfile->pGetDownCoastAdjacentProfile() == NULL)
      //          LogStream << "NULL";
      //       else
      //          LogStream << pThisProfile->pGetDownCoastAdjacentProfile()->nGetProfileID();
      //       LogStream << endl;
      //    }
      // }
      // // DEBUG CODE ============================

      // And create an index to this coast's profiles in along-coastline sequence
      m_VCoast[nCoast].CreateProfileDownCoastIndex();

      // // DEBUG CODE =======================================================================================================================
      // for (int n = 0; n < m_VCoast[nCoast].nGetNumProfiles(); n++)
      // {
      // CGeomProfile* pProfile = m_VCoast[nCoast].pGetProfileWithDownCoastSeq(n);
      // CGeomProfile* pUpCoastProfile = pProfile->pGetUpCoastAdjacentProfile();
      // CGeomProfile* pDownCoastProfile = pProfile->pGetDownCoastAdjacentProfile();
      // int nUpCoastProfile = INT_NODATA;
      // int nDownCoastProfile = INT_NODATA;
      // if (pUpCoastProfile != 0)
      // nUpCoastProfile = pUpCoastProfile->nGetProfileID();
      // if (pDownCoastProfile != 0)
      // nDownCoastProfile = pDownCoastProfile->nGetProfileID();
      // LogStream << "nCoastID = " << pProfile->nGetProfileID() << "\t up-coast profile = " << nUpCoastProfile << "\t down-coast profile = " << nDownCoastProfile << endl;
      // }
      // LogStream << endl;
      // // DEBUG CODE =======================================================================================================================

      //       // DEBUG CODE =======================================================================================================================
      // int nProf = 0;
      // for (int n = 0; n < nCoastSize; n++)
      // {
      //          // LogStream << n << "\t";
      //
      // //          LogStream << m_VCoast[nCoast].dGetDetailedCurvature(n) << "\t";
      // //
      // //          LogStream << m_VCoast[nCoast].dGetSmoothCurvature(n) << "\t";
      // //
      // //          if (m_VCoast[nCoast].pGetCoastLandform(n)->nGetLandFormCategory() == LF_INTERVENTION)
      // //             LogStream << "I\t";
      //
      // if (m_VCoast[nCoast].bIsProfileAtCoastPoint(n))
      // {
      // CGeomProfile* pProfile = m_VCoast[nCoast].pGetProfileAtCoastPoint(n);
      //
      // LogStream << "profile " << pProfile->nGetProfileID() << " at coast point " << n << " adjacent up-coast profile = " << pProfile->pGetUpCoastAdjacentProfile() << " adjacent down-coast profile = " << pProfile->pGetDownCoastAdjacentProfile() << endl;
      //
      // nProf++;
      // }
      // }
      // LogStream << endl;
      // LogStream << "nProf = " << nProf << endl;
      //       // DEBUG CODE =======================================================================================================================

      // // DEBUG CODE =======================================================================================================================
      // LogStream << "=====================" << endl;
      // for (int n = 0; n < m_VCoast[nCoast].nGetNumProfiles(); n++)
      // {
      // CGeomProfile* pProfile = m_VCoast[nCoast].pGetProfileWithDownCoastSeq(n);
      // int nSecondProfileStartCoastPoint = pProfile->nGetCoastPoint();
      //
      // LogStream << n << "\t nCoastID = " << pProfile->nGetProfileID() << "\tnSecondProfileStartCoastPoint = " << nSecondProfileStartCoastPoint << endl;
      // }
      // LogStream << endl;
      // LogStream << "=====================" << endl;
      // // DEBUG CODE =======================================================================================================================
   }

   return RTN_OK;
}

//===============================================================================================================================
//! For a single coastline, locate the start points for all coastline-normal profiles (except the grid-edge profiles). Then create the profiles
//===============================================================================================================================
void CSimulation::LocateAndCreateProfiles(int const nCoast, int& nProfile, vector<bool>* pbVCoastPointDone, vector<pair<int, double>> const* prVCurvature)
{
   // LogStream << m_ulIter << "\t in LocateAndCreateProfiles() nCoast = " << nCoast << " nProfile = " << nProfile << endl;
   int const nCoastSize = m_VCoast[nCoast].nGetCoastlineSize();

   // Work along the vector of curvature pairs starting at the most concave and most convex end
   for (int n = nCoastSize - 1; n >= 0; n--)
   {
      // Have we searched all the coastline points?
      int nStillToSearch = 0;

      for (int m = 0; m < nCoastSize; m++)
         if (! pbVCoastPointDone->at(m))
            nStillToSearch++;

      if (nStillToSearch == 0)
      {
         // OK we are done here
         // LogStream << "nStillToSearch = " << nStillToSearch << endl;
         return;
      }

      // This convex-or-concave point on the coastline is a potential location for a normal
      int const nNormalPoint = prVCurvature->at(n).first;

      // Ignore each end of the coastline
      if ((nNormalPoint == 0) || (nNormalPoint == nCoastSize - 1))
      {
         // LogStream <<"Ignoring start or end of coastline, nNormalPoint = " << nNormalPoint << endl;
         continue;
      }

      if (! pbVCoastPointDone->at(nNormalPoint))
      {
         // We have not already searched this coast point. Is it an intervention coast point?
         bool bIntervention = false;

         int const nCat = m_VCoast[nCoast].pGetCoastLandform(nNormalPoint)->nGetLandFormCategory();
         if ((nCat == LF_INTERVENTION_STRUCT) || (nCat == LF_INTERVENTION_NON_STRUCT))
         {
            // It is an intervention
            bIntervention = true;

            // We don't want profiles on very smooth sections of structural interventions
            if ((nCat == LF_INTERVENTION_STRUCT) && (prVCurvature->at(n).second < 0.1))
               continue;
         }

         CGeom2DIPoint const PtiThis = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nNormalPoint);

         // Create a profile here
         // LogStream << m_ulIter << ":\t creating profile at nNormalPoint = " << nNormalPoint << endl;
         int const nRet = nCreateProfile(nCoast, nCoastSize, nNormalPoint, nProfile, bIntervention, &PtiThis);

         // Mark this coast point as searched
         pbVCoastPointDone->at(nNormalPoint) = true;

         if (nRet != RTN_OK)
         {
            // This potential profile is no good (has hit coast, or hit dry land, etc.) so forget about it
            // LogStream << m_ulIter << ":\t profile at coastpoint " << nNormalPoint << " is no good" << endl;

            continue;
         }

         // // DEBUG CODE =================
         // LogStream << "After nCreateProfile() ===========" << endl;
         // CGeomProfile* pProfile = m_VCoast[nCoast].pGetProfile(nProfile);
         // LogStream << pProfile->nGetProfileID() << "\t";
         //
         // int nPointsInProfile = pProfile->nGetProfileSize();
         //
         // for (int nPoint = 0; nPoint < nPointsInProfile; nPoint++)
         // {
         //    CGeom2DPoint Pt = *pProfile->pPtGetPointInProfile(nPoint);
         //    LogStream << " {" << Pt.dGetX() << ", " << Pt.dGetY() << "}";
         // }
         // LogStream << endl << "===========" << endl;
         // // DEBUG CODE =================

//          // DEBUG CODE ===================================================================================================
//          LogStream << endl << "===========================================================================================" << endl;
//          LogStream << "PROFILES JUST AFTER CREATION" << endl;
//          int nNumProfiles = m_VCoast[nCoast].nGetNumProfiles();
//          for (int nn = 0; nn < nNumProfiles; nn++)
//          {
//             CGeomProfile* pProfile = m_VCoast[nCoast].pGetProfile(nn);
//
//             LogStream << nn << " nCoastID = " << pProfile->nGetProfileID() << " nGlobalID = " << pProfile->nGetProfileID() << " nGetCoastPoint = " << pProfile->nGetCoastPoint() << " pGetUpCoastAdjacentProfile = " << pProfile->pGetUpCoastAdjacentProfile() << " pGetDownCoastAdjacentProfile = " << pProfile->pGetDownCoastAdjacentProfile() << endl;
//          }
//          LogStream << "===================================================================================================" << endl << endl;
//          // DEBUG CODE ===================================================================================================
//
//          // DEBUG CODE ===================================================================================================
//          LogStream << "++++++++++++++++++++++" << endl;
//          LogStream << endl << "Just created profile " << nProfile << endl;
//          int nProf = 0;
//          for (int nnn = 0; nnn < nCoastSize; nnn++)
//          {
//             if (m_VCoast[nCoast].bIsProfileAtCoastPoint(nnn))
//             {
//                CGeomProfile* pProfile = m_VCoast[nCoast].pGetProfileAtCoastPoint(nnn);
//
//                LogStream << "profile " << pProfile->nGetProfileID() << " at coast point " << nnn << " adjacent up-coast profile = " << pProfile->pGetUpCoastAdjacentProfile() << " adjacent down-coast profile = " << pProfile->pGetDownCoastAdjacentProfile() << endl;
//
//                nProf++;
//             }
//          }
//          LogStream << endl;
//          LogStream << "nProf = " << nProf << endl;
//          LogStream << "++++++++++++++++++++++" << endl;
//          // DEBUG CODE ===================================================================================================

         // CGeom2DPoint PtThis = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nNormalPoint);
         // if (m_nLogFileDetail >= LOG_FILE_ALL)
         //    LogStream << m_ulIter << ":\t  coast " << nCoast << " profile " << nProfile << " created at coast point " << nNormalPoint << " [" << PtiThis.nGetX() << "][" << PtiThis.nGetY() << "] = {" << PtThis.dGetX() << ", " << PtThis.dGetY() << "} (smoothed curvature = " << m_VCoast[nCoast].dGetSmoothCurvature(nNormalPoint) << ", detailed curvature = " << m_VCoast[nCoast].dGetDetailedCurvature(nNormalPoint) << ")" << endl;

         // // DEBUG CODE =================================================================================
         // if (m_pRasterGrid->m_Cell[PtiThis.nGetX()][PtiThis.nGetY()].bIsCoastline())
         //    LogStream << m_ulIter << ":\t cell[" << PtiThis.nGetX() << "][" << PtiThis.nGetY() << "] IS coastline, coast number = " << m_pRasterGrid->m_Cell[PtiThis.nGetX()][PtiThis.nGetY()].nGetCoastline() << endl;
         // else
         //    LogStream << m_ulIter << ":\t ******* cell[" << PtiThis.nGetX() << "][" << PtiThis.nGetY() << "] IS NOT coastline" << endl;
         // // DEBUG CODE =================================================================================

         // This profile is fine
         nProfile++;

         // We need to mark points on either side of this profile so that we don't get profiles which are too close together. However best-placed profiles on narrow intervention structures may need to be quite close
         double dNumToMark = m_nCoastProfileSpacing;

         if (bIntervention)
            dNumToMark = m_nCoastProfileInterventionSpacing;

         // If we have a random factor for profile spacing, then modify the profile spacing
         if (m_dCoastNormalRandSpacingFactor > 0)
         {
            // Draw a sample from the unit normal distribution using random number generator 0
            double const dRand = m_dGetFromUnitNormalDist(m_Rand[0]);

            double const dTmp = dRand * m_dCoastNormalRandSpacingFactor * dNumToMark;
            dNumToMark += dTmp;

            // // Make sure number to mark is not too small or too big TODO 011
            // if (bIntervention)
            // {
            //    dNumToMark = tMin(dNumToMark, m_nCoastProfileInterventionSpacing * 0.75);
            //    dNumToMark = tMax(dNumToMark, m_nCoastProfileInterventionSpacing * 1.25);
            // }
            // else
            // {
            //    dNumToMark = tMin(dNumToMark, m_nCoastProfileSpacing * 0.75);
            //    dNumToMark = tMax(dNumToMark, m_nCoastProfileSpacing * 1.25);
            // }

            // TODO 014 Assume that the above is the profile spacing on straight bits of coast. Try gradually increasing the profile spacing with increasing concavity, and decreasing the profile spacing with increasing convexity. Could use a Michaelis-Menten S-curve relationship for this i.e.
            // double fReN = pow(NowCell[nX][nY].dGetReynolds(m_dNu), m_dDepN);
            // double fC1 = m_dC1Laminar - ((m_dC1Diff * fReN) / (fReN + m_dReMidN));
         }

         // Mark points on either side of the profile
         for (int m = 1; m < dNumToMark; m++)
         {
            int nTmpPoint = nNormalPoint + m;

            if (nTmpPoint < nCoastSize)
               pbVCoastPointDone->at(nTmpPoint) = true;

            nTmpPoint = nNormalPoint - m;

            if (nTmpPoint >= 0)
               pbVCoastPointDone->at(nTmpPoint) = true;
         }
      }
   }
}

//===============================================================================================================================
//! Creates a single coastline-normal profile (which may be an intervention profile)
//===============================================================================================================================
int CSimulation::nCreateProfile(int const nCoast, int const nCoastSize, int const nProfileStartPoint, int const nProfile, bool const bIntervention, CGeom2DIPoint const* pPtiStart)
{
   // LogStream << m_ulIter << ":\t in nCreateProfile() nProfileStartPoint = " << nProfileStartPoint << " nProfile = " << nProfile << " start point [" << pPtiStart->nGetX() << "][" << pPtiStart->nGetY() << "]" << endl;

   // OK, we have flagged the start point of this new coastline-normal profile, so create it. Make the start of the profile the centroid of the actual cell that is marked as coast (not the cell under the smoothed vector coast, they may well be different)
   CGeom2DPoint PtStart;      // In external CRS
   PtStart.SetX(dGridCentroidXToExtCRSX(pPtiStart->nGetX()));
   PtStart.SetY(dGridCentroidYToExtCRSY(pPtiStart->nGetY()));

   CGeom2DPoint PtEnd;        // In external CRS
   CGeom2DIPoint PtiEnd;      // In grid CRS
   int const nRet = nGetCoastNormalEndPoint(nCoast, nProfileStartPoint, nCoastSize, &PtStart, m_dCoastNormalLength, &PtEnd, &PtiEnd, bIntervention);
   if (nRet == RTN_ERR_NO_SOLUTION_FOR_ENDPOINT)
   {
      // Could not solve end-point equation, so forget about this profile
      LogStream << m_ulIter << ":\t  could not solve end-point equation for profile" << endl;

      return nRet;
   }

   int const nXEnd = PtiEnd.nGetX();
   int const nYEnd = PtiEnd.nGetY();

   // Safety check: is the end point in the contiguous sea?
   if (! m_pRasterGrid->m_Cell[nXEnd][nYEnd].bIsInContiguousSea())
   {
      if (m_nLogFileDetail >= LOG_FILE_ALL)
         LogStream << m_ulIter << ":\t possible profile from coast " << nCoast << " point " << nProfileStartPoint << " [" << pPtiStart->nGetX() << "][" << pPtiStart->nGetY() << "] = {" << PtStart.dGetX() << ", " << PtStart.dGetY() << "} since has inland end point [" << nXEnd << "][" << nYEnd << "] = {" << dGridCentroidXToExtCRSX(nXEnd) << ", " << dGridCentroidYToExtCRSY(nYEnd) << "}, abandoning" << endl;

      return RTN_ERR_PROFILE_ENDPOINT_IS_INLAND;
   }

   // Warning: is the water depth at the end point less than the depth of closure?
   if (m_pRasterGrid->m_Cell[nXEnd][nYEnd].dGetSeaDepth() < m_dDepthOfClosure)
   {
      if (m_nLogFileDetail >= LOG_FILE_ALL)
      {
         // LogStream << m_ulIter << ":\t possible profile from coast " << nCoast << " point " << nProfileStartPoint  << " [" << pPtiStart->nGetX() << "][" << pPtiStart->nGetY() << "] = {" << PtStart.dGetX() << ", " << PtStart.dGetY() << "} is too short for depth of closure " << m_dDepthOfClosure << " at end point [" << nXEnd << "][" << nYEnd << "] = {" << dGridCentroidXToExtCRSX(nXEnd) << ", " << dGridCentroidYToExtCRSY(nYEnd) << "}, continuing however" << endl;
      }

      // return RTN_ERR_PROFILE_END_INSUFFICIENT_DEPTH;
   }

   // No problems, so create the new profile
   CGeomProfile* pProfile = new CGeomProfile(nCoast, nProfileStartPoint, nProfile, bIntervention);

   // And create the profile's coastline-normal vector. Only two points (start and end points, both external CRS) are stored
   vector<CGeom2DPoint> VNormal;
   VNormal.push_back(PtStart);
   VNormal.push_back(PtEnd);

   // Set the start and end points (external CRS) of the profile
   pProfile->SetPointsInProfile(&VNormal);

   // Append an empty line segment to the new profile's CGeomMultiLine then set nProfile as the only co-incident profile of the only line segment
   pProfile->AppendLineSegment();
   pProfile->AppendPairToFinalLineSegment(make_pair(nProfile, 0));

   assert(pProfile->nGetProfileSize() == (1 + pProfile->nGetNumLineSegments()));

   // Save the profile, note that several fields in the profile are still blank
   m_VCoast[nCoast].AppendProfile(pProfile);

   // // DEBUG CODE =================
   // LogStream << "in nCreateProfile() ===========" << endl;
   // // CGeomProfile* pProfile = m_VCoast[nCoast].pGetProfile(nProfile);
   // LogStream << pProfile->nGetProfileID() << "\t";
   //
   // int nPointsInProfile = pProfile->nGetProfileSize();
   //
   // for (int nPoint = 0; nPoint < nPointsInProfile; nPoint++)
   // {
   // CGeom2DPoint Pt = *pProfile->pPtGetPointInProfile(nPoint);
   // LogStream << " {" << Pt.dGetX() << ", " << Pt.dGetY() << "}";
   // }
   // LogStream << endl << "===========" << endl;
   // // DEBUG CODE =================

   assert(pProfile->nGetProfileSize() > 0);

   // LogStream << m_ulIter << ":\t  coast " << nCoast << " profile " << nProfile << " created at coast point " << nProfileStartPoint << " from [" << pPtiStart->nGetX() << "][" << pPtiStart->nGetY() << "] = {" << PtStart.dGetX() << ", " << PtStart.dGetY() << "} to [" << PtiEnd.nGetX() << "][" << PtiEnd.nGetY() << "] = {" << PtEnd.dGetX() << ", " << PtEnd.dGetY() << "}" << (pProfile->bIsIntervention() ? ", from intervention" : "") << endl;

   return RTN_OK;
}

//===============================================================================================================================
//! Creates a 'special' profile at each end of a coastline, at the edge of the raster grid. This profile is not necessarily normal to the coastline since it goes along the grid's edge
//===============================================================================================================================
int CSimulation::nLocateAndCreateGridEdgeProfile(bool const bCoastStart, int const nCoast, int& nProfile)
{
   int const nCoastSize = m_VCoast[nCoast].nGetCoastlineSize();
   int const nHandedness = m_VCoast[nCoast].nGetSeaHandedness();
   int const nProfileLen = nRound(m_dCoastNormalLength / m_dCellSide);                 // Profile length in grid CRS
   // int nProfileStartEdge;

   CGeom2DIPoint PtiProfileStart;                                                      // In grid CRS
   vector<CGeom2DIPoint> VPtiNormalPoints;                                             // In grid CRS

   if (bCoastStart)
   {
      // At start of coast
      PtiProfileStart = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(0);             // Grid CRS
      // nProfileStartEdge = m_VCoast[nCoast].nGetStartEdge();
   }
   else
   {
      // At end of coast
      PtiProfileStart = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastSize - 1); // Grid CRS
      // nProfileStartEdge = m_VCoast[nCoast].nGetEndEdge();
   }

   VPtiNormalPoints.push_back(PtiProfileStart);

   // Find the start cell in the list of edge cells
   auto it = find(m_VPtiAllEdgeCell.begin(), m_VPtiAllEdgeCell.end(), PtiProfileStart);
   if (it == m_VPtiAllEdgeCell.end())
   {
      // Not found. This can happen because of rounding problems, i.e. the cell which was stored as the first cell of the raster coastline
      if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
         LogStream << m_ulIter << ": " << ERR << " when constructing start-of-coast profile, [" << PtiProfileStart.nGetX() << "][" << PtiProfileStart.nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiProfileStart.nGetX()) << ", " << dGridCentroidYToExtCRSY(PtiProfileStart.nGetY()) << "} not found in list of edge cells" << endl;

      return RTN_ERR_COAST_CANT_FIND_EDGE_CELL;
   }

   // Found
   int nPos = static_cast<int>(it - m_VPtiAllEdgeCell.begin());

   // Now construct the edge profile, searching for edge cells
   for (int n = 0; n < nProfileLen; n++)
   {
      if (bCoastStart)
      {
         if (nHandedness == LEFT_HANDED)
         {
            // At the start of a coast which is left-handed. The list of all edge cells is in clockwise sequence, so go in this direction
            nPos++;

            if (nPos >= static_cast<int>(m_VPtiAllEdgeCell.size()))
            {
               // We've reached the end of the list of edge cells, so carry on from the beginning of the list
               nPos = 0;
            }
         }
         else
         {
            // At the start of a coast which is right-handed. The list of all edge cells is in clockwise sequence, so go in the opposite direction
            nPos--;

            if (nPos < 0)
            {
               // We've reached the beginning of the list of edge cells, so carry on from the end of the list
               nPos = m_VPtiAllEdgeCell.size()-1;
            }
         }
      }
      else
      {
         if (nHandedness == LEFT_HANDED)
         {
            // At the end of a coast which is left-handed. The list of all edge cells is in clockwise sequence, so go in the opposite direction
            nPos--;

            if (nPos < 0)
            {
               // We've reached the beginning of the list of edge cells, so carry on from the end of the list
               nPos = m_VPtiAllEdgeCell.size()-1;
            }
         }
         else
         {
            // At the end of a coast which is right-handed. The list of all edge cells is in clockwise sequence, so go in this direction
            nPos++;

            if (nPos >= static_cast<int>(m_VPtiAllEdgeCell.size()))
            {
               // We've reached the end of the list of edge cells, so carry on from the beginning of the list
               nPos = 0;
            }
         }
      }

      CGeom2DIPoint const Pti = m_VPtiAllEdgeCell[nPos];

      int const nX = Pti.nGetX();
      int const nY = Pti.nGetY();

      // Have we hit a profile belonging to another coast?
      if (m_pRasterGrid->m_Cell[nX][nY].bIsProfile())
      {
         // We have hit a profile, so assume it belongs to another coast
         break;
      }

      // Have we hit another coast?
      if (m_pRasterGrid->m_Cell[nX][nY].bIsCoastline())
      {
         // We have hit a coast cell, so assume it belongs to another coast
         break;
      }

      // All OK so far, so append this grid-edge cell, making sure that there is no gap between this and the previously-appended cell (if there is, will get problems with cell-by-cell fill)
      AppendEnsureNoGap(&VPtiNormalPoints, &m_VPtiAllEdgeCell[nPos]);

      // Have we hit a corner point?
      it = find(m_VPtiBoundingBoxCorner.begin(), m_VPtiBoundingBoxCorner.end(), Pti);
      if (it != m_VPtiBoundingBoxCorner.end())
      {
         // We've reached the end of a grid side before the profile is long enough. OK, we can live with this
         break;
      }
   }

   int nProfileStartPoint;
   CGeomProfile* pProfile;
   CGeom2DIPoint const PtiDummy(INT_NODATA, INT_NODATA);

   if (bCoastStart)
   {
      nProfileStartPoint = 0;

      // Create the new start-of-coast profile
      pProfile = new CGeomProfile(nCoast, nProfileStartPoint, nProfile, false);

      // Mark this as a start-of-coast profile
      pProfile->SetStartOfCoast(true);
   }
   else
   {
      nProfileStartPoint = nCoastSize - 1;

      // Create the new end-of-coast profile
      pProfile = new CGeomProfile(nCoast, nProfileStartPoint, nProfile, false);

      // Mark this as an end-of-coast profile
      pProfile->SetEndOfCoast(true);
   }

   // Create the list of cells 'under' this grid-edge profile. Note that more than two cells are stored
   int const nPointsSize = static_cast<int>(VPtiNormalPoints.size());
   for (int n = 0; n < nPointsSize; n++)
   {
      int const nX = VPtiNormalPoints[n].nGetX();
      int const nY = VPtiNormalPoints[n].nGetY();

      // Mark each cell in the raster grid
      m_pRasterGrid->m_Cell[nX][nY].SetCoastAndProfileID(nCoast, nProfile);

      // Store the raster grid coordinates in the profile object
      pProfile->AppendCellInProfile(nX, nY);

      if ((n == 0) || (n == nPointsSize-1))
      {
         // Store the external CRS coordinates of the first and last points in the profile object
         CGeom2DPoint const Pt(dGridCentroidXToExtCRSX(nX), dGridCentroidYToExtCRSY(nY));       // In external CRS
         pProfile->AppendPointInProfile(&Pt);
      }
   }

   int const nEndX = VPtiNormalPoints.back().nGetX();
   int const nEndY = VPtiNormalPoints.back().nGetY();

   // Get the deep water wave height and orientation values at the end of the profile
   double const dDeepWaterWaveHeight = m_pRasterGrid->m_Cell[nEndX][nEndY].dGetCellDeepWaterWaveHeight();
   double const dDeepWaterWaveAngle = m_pRasterGrid->m_Cell[nEndX][nEndY].dGetCellDeepWaterWaveAngle();
   double const dDeepWaterWavePeriod = m_pRasterGrid->m_Cell[nEndX][nEndY].dGetCellDeepWaterWavePeriod();

   // And store them in this profile
   pProfile->SetProfileDeepWaterWaveHeight(dDeepWaterWaveHeight);
   pProfile->SetProfileDeepWaterWaveAngle(dDeepWaterWaveAngle);
   pProfile->SetProfileDeepWaterWavePeriod(dDeepWaterWavePeriod);

   // Append an empty line segment to the new profile's CGeomMultiLine then set nProfile as the only co-incident profile of the only line segment
   pProfile->AppendLineSegment();
   pProfile->AppendPairToFinalLineSegment(make_pair(nProfile, 0));

   assert(pProfile->nGetProfileSize() == (1 + pProfile->nGetNumLineSegments()));

   // Store the grid-edge profile
   m_VCoast[nCoast].AppendProfile(pProfile);
   m_VCoast[nCoast].SetProfileAtCoastPoint(nProfileStartPoint, pProfile);

   if (m_nLogFileDetail >= LOG_FILE_ALL)
      LogStream << m_ulIter << ":\t  coast " << nCoast << " grid-edge profile " << nProfile << " created at coast " << (bCoastStart ? "start" : "end") << " point " << (bCoastStart ? 0 : nCoastSize - 1) << ", from [" << PtiProfileStart.nGetX() << "][" << PtiProfileStart.nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiProfileStart.nGetX()) << ", " << dGridCentroidYToExtCRSY(PtiProfileStart.nGetY()) << "} to [" << VPtiNormalPoints.back().nGetX() << "][" << VPtiNormalPoints.back().nGetY() << "] = {" << dGridCentroidXToExtCRSX(VPtiNormalPoints.back().nGetX()) << ", " << dGridCentroidYToExtCRSY(VPtiNormalPoints.back().nGetY()) << "}" << endl;

   assert(pProfile->nGetProfileSize() > 0);

   return RTN_OK;
}

//===============================================================================================================================
//! Finds the end point of a coastline-normal line, given the start point on the vector coastline. If however the start point is on the grid edge, then the end point is also on the grid edge, and the line joining the start and end points is not usually normal to the vector coast. All input coordinates are in the external CRS
//===============================================================================================================================
int CSimulation::nGetCoastNormalEndPoint(int const nCoast, int const nFirstProfileStartCoastPoint, int const nCoastSize, CGeom2DPoint const* pPtStart, double const dLineLength, CGeom2DPoint* pPtEnd, CGeom2DIPoint* pPtiEnd, bool const bIntervention)
{
   int nAvgSize = 21; // TODO 011 This should be a user input

   if (bIntervention)
      nAvgSize = 3;

   double dXEnd1 = 0;
   double dXEnd2 = 0;
   double dYEnd1 = 0;
   double dYEnd2 = 0;

   CGeom2DPoint PtBefore;
   CGeom2DPoint PtAfter;

   // This is not an intervention profile. It could be a cliff collapse profile, which could be a grid-edge profile
   double const dXStart = pPtStart->dGetX();
   double const dYStart = pPtStart->dGetY();

   int const nXStart = nRound(dExtCRSXToGridX(dXStart));
   int const nYStart = nRound(dExtCRSYToGridY(dYStart));
   int const nLineLength = nConvertMetresToNumCells(dLineLength);

   // LogStream << nXStart << ", " << nYStart << endl;
   if ((nXStart == 0) || (nXStart == m_nXGridSize-1))
   {
      // Yes it is a grid-edge profile
      dXEnd1 = dGridXToExtCRSX(nXStart);
      dXEnd2 = dXEnd1;

      dYEnd1 = dGridYToExtCRSY(nYStart + nLineLength);
      dYEnd2 = dGridYToExtCRSY(nYStart - nLineLength);
   }
   else if ((nYStart == 0) || (nYStart == m_nYGridSize-1))
   {
      // Yes it is a grid-edge profile
      dYEnd1 = dGridYToExtCRSY(nYStart);
      dYEnd2 = dYEnd1;

      dXEnd1 = dGridXToExtCRSX(nXStart + nLineLength);
      dXEnd2 = dGridXToExtCRSX(nXStart - nLineLength);
   }
   else
   {
      // This is not a grid-edge profile, so put a maximum of  points before the start point into a vector
      vector<CGeom2DPoint> VPtBeforeToAverage;

      for (int n = 1; n <= nAvgSize; n++)
      {
         int const nPoint = nFirstProfileStartCoastPoint - n;
         if (nPoint < 0)
            break;

         VPtBeforeToAverage.push_back(*m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nPoint));
      }

      // Put a maximum of nAvgSize points after the start point into a vector
      vector<CGeom2DPoint> VPtAfterToAverage;

      for (int n = 1; n <= nAvgSize; n++)
      {
         int const nPoint = nFirstProfileStartCoastPoint + n;
         if (nPoint > nCoastSize - 1)
            break;

         VPtAfterToAverage.push_back(*m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nPoint));
      }

      // Now average each of these vectors of points: results are in PtBefore and PtAfter (coordinates in external CRS)
      PtBefore = PtAverage(&VPtBeforeToAverage);
      PtAfter = PtAverage(&VPtAfterToAverage);

      // Get the y = a * x + b equation of the straight line linking the coastline points before and after 'this' coastline point. For this linking line, slope a = (y2 - y1) / (x2 - x1)
      double const dYDiff = PtAfter.dGetY() - PtBefore.dGetY();
      double const dXDiff = PtAfter.dGetX() - PtBefore.dGetX();

      if (bFPIsEqual(dYDiff, 0.0, TOLERANCE))
      {
         // The linking line runs W-E or E-W, so a straight line at right angles to this runs N-S or S-N. Calculate the two possible end points for this coastline-normal profile
         dXEnd1 = dXEnd2 = pPtStart->dGetX();
         dYEnd1 = pPtStart->dGetY() + dLineLength;
         dYEnd2 = pPtStart->dGetY() - dLineLength;
      }
      else if (bFPIsEqual(dXDiff, 0.0, TOLERANCE))
      {
         // The linking line runs N-S or S-N, so a straight line at right angles to this runs W-E or E-W. Calculate the two possible end points for this coastline-normal profile
         dYEnd1 = dYEnd2 = pPtStart->dGetY();
         dXEnd1 = pPtStart->dGetX() + dLineLength;
         dXEnd2 = pPtStart->dGetX() - dLineLength;
      }
      else
      {
         // The linking line runs neither W-E nor N-S so we have to work a bit harder to find the end-point of the coastline-normal profile
         double const dA = dYDiff / dXDiff;

         // Now calculate the equation of the straight line which is perpendicular to this linking line
         double const dAPerp = -1 / dA;
         double const dBPerp = pPtStart->dGetY() - (dAPerp * pPtStart->dGetX());

         // Calculate the end point of the profile: first do some substitution then rearrange as a quadratic equation i.e. in the form Ax^2 + Bx + C = 0 (see http://math.stackexchange.com/questions/228841/how-do-i-calculate-the-intersections-of-a-straight-line-and-a-circle)
         double const dQuadA = 1 + (dAPerp * dAPerp);
         double const dQuadB = 2 * ((dBPerp * dAPerp) - (dAPerp * pPtStart->dGetY()) - pPtStart->dGetX());
         double const dQuadC = ((pPtStart->dGetX() * pPtStart->dGetX()) + (pPtStart->dGetY() * pPtStart->dGetY()) + (dBPerp * dBPerp) - (2 * pPtStart->dGetY() * dBPerp) - (dLineLength * dLineLength));

         // Solve for x and y using the quadratic formula x = (−B ± sqrt(B^2 − 4AC)) / 2A
         double const dDiscriminant = (dQuadB * dQuadB) - (4 * dQuadA * dQuadC);

         if (dDiscriminant < 0)
         {
            LogStream << ERR << "timestep " << m_ulIter << ": discriminant < 0 when finding profile end point on coastline " << nCoast << ", from coastline point " << nFirstProfileStartCoastPoint << "), ignored" << endl;
            return RTN_ERR_NO_SOLUTION_FOR_ENDPOINT;
         }

         dXEnd1 = (-dQuadB + sqrt(dDiscriminant)) / (2 * dQuadA);
         dYEnd1 = (dAPerp * dXEnd1) + dBPerp;
         dXEnd2 = (-dQuadB - sqrt(dDiscriminant)) / (2 * dQuadA);
         dYEnd2 = (dAPerp * dXEnd2) + dBPerp;
      }
   }

   // We have two possible solutions, so decide which of the two endpoints to use then create the profile end-point (coordinates in external CRS)
   int const nSeaHand = m_VCoast[nCoast].nGetSeaHandedness(); // Assumes handedness is either 0 or 1 (i.e. not -1)
   *pPtEnd = PtChooseEndPoint(nSeaHand, &PtBefore, &PtAfter, dXEnd1, dYEnd1, dXEnd2, dYEnd2);

   // Check that pPtiEnd is not off the grid. Note that pPtiEnd is not necessarily a cell centroid
   pPtiEnd->SetXY(nRound(dExtCRSXToGridX(pPtEnd->dGetX())), nRound(dExtCRSYToGridY(pPtEnd->dGetY())));

   if (! bIsWithinValidGrid(pPtiEnd))
   {
      // LogStream << m_ulIter << ": profile endpoint is outside grid [" << pPtiEnd->nGetX() << "][" << pPtiEnd->nGetY() << "] = {" << pPtEnd->dGetX() << ", " << pPtEnd->dGetY() << "}. The profile starts at coastline point " << nFirstProfileStartCoastPoint << " = {" << pPtStart->dGetX() << ", " << pPtStart->dGetY() << "}" << endl;

      // The end point is off the grid, so constrain it to be within the valid grid
      CGeom2DIPoint const PtiStart(nRound(dExtCRSXToGridX(pPtStart->dGetX())), nRound(dExtCRSYToGridY(pPtStart->dGetY())));
      KeepWithinValidGrid(&PtiStart, pPtiEnd);

      pPtEnd->SetX(dGridCentroidXToExtCRSX(pPtiEnd->nGetX()));
      pPtEnd->SetY(dGridCentroidYToExtCRSY(pPtiEnd->nGetY()));

      // LogStream << m_ulIter << ":\t coast " << nCoast << " profile endpoint constrained to be within grid, is now [" << pPtiEnd->nGetX() << "][" << pPtiEnd->nGetY() << "] = {" << pPtEnd->dGetX() << ", " << pPtEnd->dGetY() << "}. The profile starts at coastline point " << nFirstProfileStartCoastPoint << " = {" << pPtStart->dGetX() << ", " << pPtStart->dGetY() << "}" << endl;
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Choose which end point to use for the coastline-normal profile
//===============================================================================================================================
CGeom2DPoint CSimulation::PtChooseEndPoint(int const nHand, CGeom2DPoint const* PtBefore, CGeom2DPoint const* PtAfter, double const dXEnd1, double const dYEnd1, double const dXEnd2, double const dYEnd2)
{
   CGeom2DPoint PtChosen;

   // All coordinates here are in the external CRS, so the origin of the grid is the bottom left
   if (nHand == RIGHT_HANDED)
   {
      // The sea is to the right of the linking line. So which way is the linking line oriented? First check the N-S component
      if (PtAfter->dGetY() > PtBefore->dGetY())
      {
         // We are going S to N and the sea is to the right: the normal endpoint is to the E. We want the larger of the two x values
         if (dXEnd1 > dXEnd2)
         {
            PtChosen.SetX(dXEnd1);
            PtChosen.SetY(dYEnd1);
         }
         else
         {
            PtChosen.SetX(dXEnd2);
            PtChosen.SetY(dYEnd2);
         }
      }
      else if (PtAfter->dGetY() < PtBefore->dGetY())
      {
         // We are going N to S and the sea is to the right: the normal endpoint is to the W. We want the smaller of the two x values
         if (dXEnd1 < dXEnd2)
         {
            PtChosen.SetX(dXEnd1);
            PtChosen.SetY(dYEnd1);
         }
         else
         {
            PtChosen.SetX(dXEnd2);
            PtChosen.SetY(dYEnd2);
         }
      }
      else
      {
         // No N-S component i.e. the linking line is exactly W-E. So check the W-E component
         if (PtAfter->dGetX() > PtBefore->dGetX())
         {
            // We are going W to E and the sea is to the right: the normal endpoint is to the s. We want the smaller of the two y values
            if (dYEnd1 < dYEnd2)
            {
               PtChosen.SetX(dXEnd1);
               PtChosen.SetY(dYEnd1);
            }
            else
            {
               PtChosen.SetX(dXEnd2);
               PtChosen.SetY(dYEnd2);
            }
         }
         else // Do not check for (PtAfter->dGetX() == PtBefore->dGetX()), since this would mean the two points are co-incident
         {
            // We are going E to W and the sea is to the right: the normal endpoint is to the N. We want the larger of the two y values
            if (dYEnd1 > dYEnd2)
            {
               PtChosen.SetX(dXEnd1);
               PtChosen.SetY(dYEnd1);
            }
            else
            {
               PtChosen.SetX(dXEnd2);
               PtChosen.SetY(dYEnd2);
            }
         }
      }
   }
   else // nHand == LEFT_HANDED
   {
      // The sea is to the left of the linking line. So which way is the linking line oriented? First check the N-S component
      if (PtAfter->dGetY() > PtBefore->dGetY())
      {
         // We are going S to N and the sea is to the left: the normal endpoint is to the W. We want the smaller of the two x values
         if (dXEnd1 < dXEnd2)
         {
            PtChosen.SetX(dXEnd1);
            PtChosen.SetY(dYEnd1);
         }
         else
         {
            PtChosen.SetX(dXEnd2);
            PtChosen.SetY(dYEnd2);
         }
      }
      else if (PtAfter->dGetY() < PtBefore->dGetY())
      {
         // We are going N to S and the sea is to the left: the normal endpoint is to the E. We want the larger of the two x values
         if (dXEnd1 > dXEnd2)
         {
            PtChosen.SetX(dXEnd1);
            PtChosen.SetY(dYEnd1);
         }
         else
         {
            PtChosen.SetX(dXEnd2);
            PtChosen.SetY(dYEnd2);
         }
      }
      else
      {
         // No N-S component i.e. the linking line is exactly W-E. So check the W-E component
         if (PtAfter->dGetX() > PtBefore->dGetX())
         {
            // We are going W to E and the sea is to the left: the normal endpoint is to the N. We want the larger of the two y values
            if (dYEnd1 > dYEnd2)
            {
               PtChosen.SetX(dXEnd1);
               PtChosen.SetY(dYEnd1);
            }
            else
            {
               PtChosen.SetX(dXEnd2);
               PtChosen.SetY(dYEnd2);
            }
         }
         else // Do not check for (PtAfter->dGetX() == PtBefore->dGetX()), since this would mean the two points are co-incident
         {
            // We are going E to W and the sea is to the left: the normal endpoint is to the S. We want the smaller of the two y values
            if (dYEnd1 < dYEnd2)
            {
               PtChosen.SetX(dXEnd1);
               PtChosen.SetY(dYEnd1);
            }
            else
            {
               PtChosen.SetX(dXEnd2);
               PtChosen.SetY(dYEnd2);
            }
         }
      }
   }

   return PtChosen;
}

