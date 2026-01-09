/*!
   \file create_profiles.cpp
   \brief Creates profiles which are approximately normal to the coastline, these will become inter-polygon boundaries
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

#include "cme.h"
#include "simulation.h"
#include "coast.h"
#include "2d_point.h"
#include "2di_point.h"

namespace
{
//===============================================================================================================================
//! Function used to sort coastline curvature values when locating start points of normal profiles
//===============================================================================================================================
bool bCurvaturePairCompareDescending(const pair<int, double>& prLeft, const pair<int, double>& prRight)
{
   // Sort in descending order (i.e. most concave first)
   return prLeft.second > prRight.second;
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
      // DEBUG CODE ==============================================================================================================
      if ((m_ulIter == 851) && (nCoast == 1))
      {
         LogStream <<  "=================" << endl;
      }
      // DEBUG CODE ==============================================================================================================

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
         if ((nCat != LF_INTERVENTION_STRUCT) && (nCat != LF_INTERVENTION_NON_STRUCT))
         {
            // Not an intervention coast point, so store the smoothed curvature
            dCurvature = m_VCoast[nCoast].dGetSmoothCurvature(nCoastPoint);
         }
         else
         {
            // This is an intervention coast point, which is likely to have some sharp angles. So store the detailed curvature
            dCurvature = m_VCoast[nCoast].dGetDetailedCurvature(nCoastPoint);
         }

         prVCurvature.push_back(make_pair(nCoastPoint, dCurvature));
      }

      // Sort this pair vector in descending order, so that the most convex curvature points are first
      sort(prVCurvature.begin(), prVCurvature.end(), bCurvaturePairCompareDescending);

      // // DEBUG CODE =======================================================================================================================
      // for (int n = 0; n < prVCurvature.size(); n++)
      // {
      // LogStream << prVCurvature[n].first << "\t" << prVCurvature[n].second << endl;
      // }
      // LogStream << endl << endl;
      // // DEBUG CODE =======================================================================================================================

      // And mark points at and near the start and end of the coastline so that they don't get searched (will be creating 'special' start- and end-of-coast profiles at these end points later)
      for (int n = 0; n < m_nCoastNormalSpacing; n++)
      {
         if (n < nCoastSize)
            bVCoastPointDone[n] = true;

         int const m = nCoastSize - n - 1;

         if (m >= 0)
            bVCoastPointDone[m] = true;
      }

      // Now locate the start points for all coastline-normal profiles (except the grid-edge ones), at points of maximum convexity. Then create the profiles
      LocateAndCreateProfiles(nCoast, nProfile, &bVCoastPointDone, &prVCurvature);

      // Did we fail to create any normal profiles? If so, quit
      if (nProfile < 0)
      {
         string strErr = ERR + "timestep " + strDblToStr(m_ulIter) + ": could not create profiles for coastline " + strDblToStr(nCoast);

         if (m_ulIter == 1)
            strErr += ". Check the SWL";

         strErr += "\n";

         cerr << strErr;
         LogStream << strErr;

         return RTN_ERR_NO_PROFILES_1;
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
            // nThisProfile = pThisProfile->nGetProfileID();

            if (nCoastPoint == 0)
            {
               pThisProfile->SetUpCoastAdjacentProfile(NULL);

               // LogStream << "nCoastPoint = " << nCoastPoint << " ThisProfile = " << nThisProfile << " ThisProfile UpCoast = " << pThisProfile->pGetUpCoastAdjacentProfile() << " ThisProfile DownCoast = " << pThisProfile->pGetDownCoastAdjacentProfile() << endl;

               pLastProfile = pThisProfile;
               // nLastProfile = nThisProfile;
               continue;
            }

            pLastProfile->SetDownCoastAdjacentProfile(pThisProfile);
            pThisProfile->SetUpCoastAdjacentProfile(pLastProfile);

            if (nCoastPoint == nCoastSize - 1)
               pThisProfile->SetDownCoastAdjacentProfile(NULL);

            pLastProfile = pThisProfile;
         }
      }

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
      // int nStartPoint = pProfile->nGetCoastPoint();
      //
      // LogStream << n << "\t nCoastID = " << pProfile->nGetProfileID() << "\tnStartPoint = " << nStartPoint << endl;
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

   // Work along the vector of curvature pairs starting at the convex end
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

      // This convex point on the coastline is a potential location for a normal
      int const nNormalPoint = prVCurvature->at(n).first;

      // Ignore each end of the coastline
      if ((nNormalPoint == 0) || (nNormalPoint == nCoastSize - 1))
      {
         // LogStream <<"Ignoring start or end of coastline, nNormalPoint = " << nNormalPoint << endl;
         continue;
      }

      // TODO 089 When choosing locations for profiles, do coast first then interventions

      if (! pbVCoastPointDone->at(nNormalPoint))
      {
         // We have not already searched this coast point. Is it an intervention coast point?
         bool bIntervention = false;

         int const nCat = m_VCoast[nCoast].pGetCoastLandform(nNormalPoint)->nGetLandFormCategory();
         if ((nCat == LF_INTERVENTION_STRUCT) || (nCat == LF_INTERVENTION_NON_STRUCT))
         {
            // It is an intervention
            bIntervention = true;
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
            LogStream << m_ulIter << ":\t profile at coastpoint " << nNormalPoint << " is no good" << endl;
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
         //    LogStream << m_ulIter << ":\t coast " << nCoast << " profile " << nProfile << " created at coast point " << nNormalPoint << " [" << PtiThis.nGetX() << "][" << PtiThis.nGetY() << "] = {" << PtThis.dGetX() << ", " << PtThis.dGetY() << "} (smoothed curvature = " << m_VCoast[nCoast].dGetSmoothCurvature(nNormalPoint) << ", detailed curvature = " << m_VCoast[nCoast].dGetDetailedCurvature(nNormalPoint) << ")" << endl;

         // // DEBUG CODE =================================================================================
         // if (m_pRasterGrid->m_Cell[PtiThis.nGetX()][PtiThis.nGetY()].bIsCoastline())
         //    LogStream << m_ulIter << ":\t cell[" << PtiThis.nGetX() << "][" << PtiThis.nGetY() << "] IS coastline, coast number = " << m_pRasterGrid->m_Cell[PtiThis.nGetX()][PtiThis.nGetY()].nGetCoastline() << endl;
         // else
         //    LogStream << m_ulIter << ":\t ******* cell[" << PtiThis.nGetX() << "][" << PtiThis.nGetY() << "] IS NOT coastline" << endl;
         // // DEBUG CODE =================================================================================

         // This profile is fine
         nProfile++;

         // We need to mark points on either side of this profile so that we don't get profiles which are too close together. However best-placed profiles on narrow intervention structures may need to be quite closes
         double dNumToMark = m_nCoastNormalSpacing;

         if (bIntervention)
            dNumToMark = m_nCoastNormalInterventionSpacing;

         // If we have a random factor for profile spacing, then modify the profile spacing
         if (m_dCoastNormalRandSpacingFactor > 0)
         {
            // Draw a sample from the unit normal distribution using random number generator 0
            double const dRand = m_dGetFromUnitNormalDist(m_Rand[0]);

            double const dTmp = dRand * m_dCoastNormalRandSpacingFactor * dNumToMark;
            dNumToMark += dTmp;

            // Make sure number to mark is not too small or too big TODO 011
            if (bIntervention)
            {
               dNumToMark = tMin(dNumToMark, m_nCoastNormalInterventionSpacing * 0.75);
               dNumToMark = tMax(dNumToMark, m_nCoastNormalInterventionSpacing * 1.25);
            }

            else
            {
               dNumToMark = tMin(dNumToMark, m_nCoastNormalSpacing * 0.75);
               dNumToMark = tMax(dNumToMark, m_nCoastNormalSpacing * 1.25);
            }

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
         LogStream << m_ulIter << ": coast " << nCoast << ", possible profile with start point " << nProfileStartPoint << " has inland end point at [" << nXEnd << "][" << nYEnd << "] = {" << dGridCentroidXToExtCRSX(nXEnd) << ", " << dGridCentroidYToExtCRSY(nYEnd) << "}, ignoring" << endl;

      return RTN_ERR_PROFILE_ENDPOINT_IS_INLAND;
   }

   // Warning: is the water depth at the end point less than the depth of closure?
   if (m_pRasterGrid->m_Cell[nXEnd][nYEnd].dGetSeaDepth() < m_dDepthOfClosure)
   {
      if (m_nLogFileDetail >= LOG_FILE_ALL)
      {
         LogStream << m_ulIter << ": coast " << nCoast << ", possible profile with start point " << nProfileStartPoint << " is too short for depth of closure " << m_dDepthOfClosure << " at end point [" << nXEnd << "][" << nYEnd << "] = {" << dGridCentroidXToExtCRSX(nXEnd) << ", " << dGridCentroidYToExtCRSY(nYEnd) << "}, continuing however. Suggest increasing the length of coastline normals" << endl;
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

   // Create the profile's CGeomMultiLine then set nProfile as the only co-incident profile of the only line segment
   pProfile->AppendLineSegment();
   pProfile->AppendCoincidentProfileToLineSegments(make_pair(nProfile, 0));

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

   // assert(pProfile->nGetProfileSize() > 0);

   LogStream << m_ulIter << ":\t coast " << nCoast << " profile " << nProfile << " created at coast point " << nProfileStartPoint << " from [" << pPtiStart->nGetX() << "][" << pPtiStart->nGetY() << "] = {" << PtStart.dGetX() << ", " << PtStart.dGetY() << "} to [" << PtiEnd.nGetX() << "][" << PtiEnd.nGetY() << "] = {" << PtEnd.dGetX() << ", " << PtEnd.dGetY() << "}" << (pProfile->bIsIntervention() ? ", from intervention" : "") << endl;

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
               // We've reached the end of the list of edge cells before the profile is long enough. OK, we can live with this
               break;
            }
         }
         else
         {
            // At the start of a coast which is right-handed. The list of all edge cells is in clockwise sequence, so go in the opposite direction
            nPos--;

            if (nPos < 0)
            {
               // We've reached the beginning of the list of edge cells before the profile is long enough. OK, we can live with this
               break;
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
               // We've reached the beginning of the list of edge cells before the profile is long enough. OK, we can live with this
               break;
            }
         }
         else
         {
            // At the end of a coast which is right-handed. The list of all edge cells is in clockwise sequence, so go in this direction
            nPos++;

            if (nPos >= static_cast<int>(m_VPtiAllEdgeCell.size()))
            {
               // We've reached the end of the list of edge cells before the profile is long enough. OK, we can live with this
               break;
            }
         }
      }

      // Have we hit a corner point?
      CGeom2DIPoint const Pti = m_VPtiAllEdgeCell[nPos];
      it = find(m_VPtiBoundingBoxCorner.begin(), m_VPtiBoundingBoxCorner.end(), Pti);
      if (it != m_VPtiBoundingBoxCorner.end())
      {
         // We've reached the end of a grid side before the profile is long enough. OK, we can live with this
         break;
      }

      // All OK, so append this grid-edge cell, making sure that there is no gap between this and the previously-appended cell (if there is, will get problems with cell-by-cell fill)
      AppendEnsureNoGap(&VPtiNormalPoints, &m_VPtiAllEdgeCell[nPos]);
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
   for (unsigned int n = 0; n < VPtiNormalPoints.size(); n++)
   {
      int const nX = VPtiNormalPoints[n].nGetX();
      int const nY = VPtiNormalPoints[n].nGetY();

      // Mark each cell in the raster grid
      m_pRasterGrid->m_Cell[nX][nY].SetCoastAndProfileID(nCoast, nProfile);

      // Store the raster grid coordinates in the profile object
      pProfile->AppendCellInProfile(nX, nY);

      CGeom2DPoint const Pt(dGridCentroidXToExtCRSX(nX), dGridCentroidYToExtCRSY(nY)); // In external CRS

      // Store the external coordinates in the profile object. Note that for this grid-edge profile, the coordinates of the cells and the coordinates of points on the profile itself are identical, this is not the case for ordinary profiles
      pProfile->AppendPointInProfile(&Pt);
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

   // Create the profile's CGeomMultiLine then set nProfile as the only co-incident profile of the only line segment
   pProfile->AppendLineSegment();
   pProfile->AppendCoincidentProfileToLineSegments(make_pair(nProfile, 0));

   // Store the grid-edge profile
   m_VCoast[nCoast].AppendProfile(pProfile);
   m_VCoast[nCoast].SetProfileAtCoastPoint(nProfileStartPoint, pProfile);

   if (m_nLogFileDetail >= LOG_FILE_ALL)
      LogStream << m_ulIter << ":\t coast " << nCoast << " grid-edge profile " << nProfile << " created at coast " << (bCoastStart ? "start" : "end") << " point " << (bCoastStart ? 0 : nCoastSize - 1) << ", from [" << PtiProfileStart.nGetX() << "][" << PtiProfileStart.nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiProfileStart.nGetX()) << ", " << dGridCentroidYToExtCRSY(PtiProfileStart.nGetY()) << "} to [" << VPtiNormalPoints.back().nGetX() << "][" << VPtiNormalPoints.back().nGetY() << "] = {" << dGridCentroidXToExtCRSX(VPtiNormalPoints.back().nGetX()) << ", " << dGridCentroidYToExtCRSY(VPtiNormalPoints.back().nGetY()) << "}" << endl;

   // assert(pProfile->nGetProfileSize() > 0);

   return RTN_OK;
}

//===============================================================================================================================
//! Finds the end point of a coastline-normal line, given the start point on the vector coastline. If however the start point is on the grid edge (only applicable to cliff collapse [rofiles), then the end point is also on the grid edge, and the line joining the start and end points is not necessarily normal to the vector coast. All input coordinates are in the external CRS
//===============================================================================================================================
int CSimulation::nGetCoastNormalEndPoint(int const nCoast, int const nStartCoastPoint, int const nCoastSize, CGeom2DPoint const* pPtStart, double const dLineLength, CGeom2DPoint* pPtEnd, CGeom2DIPoint* pPtiEnd, bool const bIntervention)
{
   int const AVGSIZE = 21; // TODO 011 This should be a user input

   double dXEnd1 = 0;
   double dXEnd2 = 0;
   double dYEnd1 = 0;
   double dYEnd2 = 0;

   CGeom2DPoint PtBefore;
   CGeom2DPoint PtAfter;

   if (bIntervention)
   {
      // This is an intervention profile, so just use one point on either side (coordinates in external CRS). TODO Note this this assumes that this intervention profile is not at the start or end of the coastline
      PtBefore = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nStartCoastPoint - 1);
      PtAfter = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nStartCoastPoint + 1);
   }
   else
   {
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
         // This is not a grid-edge profile, so put a maximum of AVGSIZE points before the start point into a vector
         vector<CGeom2DPoint> VPtBeforeToAverage;

         for (int n = 1; n <= AVGSIZE; n++)
         {
            int const nPoint = nStartCoastPoint - n;
            if (nPoint < 0)
               break;

            VPtBeforeToAverage.push_back(*m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nPoint));
         }

         // Put a maximum of AVGSIZE points after the start point into a vector
         vector<CGeom2DPoint> VPtAfterToAverage;

         for (int n = 1; n <= AVGSIZE; n++)
         {
            int const nPoint = nStartCoastPoint + n;
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
               LogStream << ERR << "timestep " << m_ulIter << ": discriminant < 0 when finding profile end point on coastline " << nCoast << ", from coastline point " << nStartCoastPoint << "), ignored" << endl;
               return RTN_ERR_NO_SOLUTION_FOR_ENDPOINT;
            }

            dXEnd1 = (-dQuadB + sqrt(dDiscriminant)) / (2 * dQuadA);
            dYEnd1 = (dAPerp * dXEnd1) + dBPerp;
            dXEnd2 = (-dQuadB - sqrt(dDiscriminant)) / (2 * dQuadA);
            dYEnd2 = (dAPerp * dXEnd2) + dBPerp;
         }
      }
   }

   // We have two possible solutions, so decide which of the two endpoints to use then create the profile end-point (coordinates in external CRS)
   int const nSeaHand = m_VCoast[nCoast].nGetSeaHandedness(); // Assumes handedness is either 0 or 1 (i.e. not -1)
   *pPtEnd = PtChooseEndPoint(nSeaHand, &PtBefore, &PtAfter, dXEnd1, dYEnd1, dXEnd2, dYEnd2);

   // Check that pPtiEnd is not off the grid. Note that pPtiEnd is not necessarily a cell centroid
   pPtiEnd->SetXY(nRound(dExtCRSXToGridX(pPtEnd->dGetX())), nRound(dExtCRSYToGridY(pPtEnd->dGetY())));

   if (! bIsWithinValidGrid(pPtiEnd))
   {
      // LogStream << m_ulIter << ": profile endpoint is outside grid [" << pPtiEnd->nGetX() << "][" << pPtiEnd->nGetY() << "] = {" << pPtEnd->dGetX() << ", " << pPtEnd->dGetY() << "}. The profile starts at coastline point " << nStartCoastPoint << " = {" << pPtStart->dGetX() << ", " << pPtStart->dGetY() << "}" << endl;

      // The end point is off the grid, so constrain it to be within the valid grid
      CGeom2DIPoint const PtiStart(nRound(dExtCRSXToGridX(pPtStart->dGetX())), nRound(dExtCRSYToGridY(pPtStart->dGetY())));
      KeepWithinValidGrid(&PtiStart, pPtiEnd);

      pPtEnd->SetX(dGridCentroidXToExtCRSX(pPtiEnd->nGetX()));
      pPtEnd->SetY(dGridCentroidYToExtCRSY(pPtiEnd->nGetY()));

      // LogStream << m_ulIter << ":\t coast " << nCoast << " profile endpoint constrained to be within grid, is now [" << pPtiEnd->nGetX() << "][" << pPtiEnd->nGetY() << "] = {" << pPtEnd->dGetX() << ", " << pPtEnd->dGetY() << "}. The profile starts at coastline point " << nStartCoastPoint << " = {" << pPtStart->dGetX() << ", " << pPtStart->dGetY() << "}" << endl;
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

//===============================================================================================================================
//! Checks all coastline-normal profiles for intersection, and modifies those that intersect
//===============================================================================================================================
void CSimulation::CheckForIntersectingProfiles(void)
{
   LogStream << endl << m_ulIter << ": Checking for profile intersection" << endl;

   // Do once for every coastline object
   int const nCoastLines = static_cast<int>(m_VCoast.size());

   for (int nCoast = 0; nCoast < nCoastLines; nCoast++)
   {
      int const nCoastSize = m_VCoast[nCoast].nGetCoastlineSize();

      // Do once for every profile, in along-coast sequence
      for (int nCoastPoint = 0; nCoastPoint < nCoastSize; nCoastPoint++)
      {
         if (! m_VCoast[nCoast].bIsProfileAtCoastPoint(nCoastPoint))
            continue;

         // There is a profile at this coast point
         CGeomProfile* pFirstProfile = m_VCoast[nCoast].pGetProfileAtCoastPoint(nCoastPoint);
         int const nFirstProfile = pFirstProfile->nGetProfileID();

         // Only check this profile if it is problem free, and is not a start- or end-of-coast profile. Continue checking if it has been truncated, however
         if (! pFirstProfile->bProfileOKIncTruncated())
         {
            // LogStream << m_ulIter << ": nCoastPoint = " << nCoastPoint << " pFirstProfile = " << pFirstProfile->nGetProfileID() << " is not OK (could be a start- or end-of-coast profile), abandoning" << endl;
            continue;
         }

         // OK we have found a first profile. Now go along the coast in alternate directions: first down-coast (in the direction of increasing coast point numbers) then up-coast
         for (int nDirection = DIRECTION_DOWNCOAST; nDirection <= DIRECTION_UPCOAST; nDirection++)
         {
            int nStartPoint;

            if (nDirection == DIRECTION_DOWNCOAST)
               nStartPoint = nCoastPoint + 1;
            else
               nStartPoint = nCoastPoint - 1;

            for (int nSecondCoastPoint = nStartPoint; (nDirection == DIRECTION_DOWNCOAST) ? nSecondCoastPoint < nCoastSize : nSecondCoastPoint >= 0; (nDirection == DIRECTION_DOWNCOAST) ? nSecondCoastPoint++ : nSecondCoastPoint--)
            // // In this direction, look at profiles which are increasingly close to the first profile
            // int nStartPoint;
            // if (nDirection == DIRECTION_DOWNCOAST)
            // nStartPoint = 0;
            // else
            // nStartPoint = nCoastSize - 1;
            //
            // for (int nSecondCoastPoint = nStartPoint; (nDirection == DIRECTION_DOWNCOAST) ? nSecondCoastPoint < nCoastPoint : nSecondCoastPoint > nCoastPoint; (nDirection == DIRECTION_DOWNCOAST) ? nSecondCoastPoint++ : nSecondCoastPoint--)
            {
               if (m_VCoast[nCoast].bIsProfileAtCoastPoint(nSecondCoastPoint))
               {
                  // There is a profile at the second coast point, so get a pointer to it
                  CGeomProfile* pSecondProfile = m_VCoast[nCoast].pGetProfileAtCoastPoint(nSecondCoastPoint);
                  int const nSecondProfile = pSecondProfile->nGetProfileID();

                  // LogStream << m_ulIter << ": " << (nDirection == DIRECTION_DOWNCOAST ? "down" : "up") << "-coast search, nCoastPoint = " << nCoastPoint << " nSecondCoastPoint = " << nSecondCoastPoint << " (profiles " << pFirstProfile->nGetProfileID() << " and " << pSecondProfile->nGetProfileID() << ")" << endl;

                  // Only check this profile if it is problem free, and is not a start- or end-of-coast profile. Continue checking if it has been truncated, however
                  if (! pSecondProfile->bProfileOKIncTruncated())
                  {
                     // LogStream << m_ulIter << ": second profile = " << pSecondProfile->nGetProfileID() << " is not OK (could be a start- or end-of-coast profile), abandoning" << endl;
                     continue;
                  }

                  // Only check these two profiles for intersection if they are are not co-incident in the final line segment of both profiles (i.e. the profiles have not already intersected)
                  if ((pFirstProfile->bFindProfileInCoincidentProfilesOfLastLineSegment(nSecondProfile)) || (pSecondProfile->bFindProfileInCoincidentProfilesOfLastLineSegment(nFirstProfile)))
                  {
                     // LogStream << m_ulIter << ": profiles " << pFirstProfile->nGetProfileID() << " and " << pSecondProfile->nGetProfileID() << " are are not co-incident in the final line segment of both profiles (i.e. the profiles have not already intersected), abandoning" << endl;
                     continue;
                  }

                  // OK go for it
                  int nProf1LineSeg = 0;
                  int nProf2LineSeg = 0;
                  double dIntersectX = 0;
                  double dIntersectY = 0;
                  double dAvgEndX = 0;
                  double dAvgEndY = 0;

                  if (bCheckForIntersection(pFirstProfile, pSecondProfile, nProf1LineSeg, nProf2LineSeg, dIntersectX, dIntersectY, dAvgEndX, dAvgEndY))
                  {
                     // The profiles intersect. Decide which profile to truncate, and which to retain
                     int nPoint = -1;

                     if (pFirstProfile->bIsIntervention())
                     {
                        LogStream << m_ulIter << ": profiles " << nFirstProfile << " and " << nSecondProfile << " intersect, truncate " << nFirstProfile << " since it is an intervention profile" << endl;

                        // Truncate the first profile, since it is an intervention profile
                        TruncateOneProfileRetainOtherProfile(nCoast, pFirstProfile, pSecondProfile, dIntersectX, dIntersectY, nProf1LineSeg, nProf2LineSeg, false);
                     }
                     else if (pSecondProfile->bIsIntervention())
                     {
                        LogStream << m_ulIter << ": profiles " << nFirstProfile << " and " << nSecondProfile << " intersect, truncate " << nSecondProfile << " since it is an intervention profile" << endl;

                        // Truncate the second profile, since it is an intervention profile
                        TruncateOneProfileRetainOtherProfile(nCoast, pSecondProfile, pFirstProfile, dIntersectX, dIntersectY, nProf2LineSeg, nProf1LineSeg, false);
                     }
                     // Is the point of intersection already present in the first profile (i.e. because there has already been an intersection at this point between the first profile and some other profile)?
                     else if (pFirstProfile->bIsPointInProfile(dIntersectX, dIntersectY, nPoint))
                     {
                        LogStream << m_ulIter << ":\t  coast " << nCoast << " profiles " << nFirstProfile << " and " << nSecondProfile << " intersect, but point {" << dIntersectX << ", " << dIntersectY << "} is already present in profile " << nFirstProfile << " as point " << nPoint << endl;

                        // Truncate the second profile and merge it with the first profile
                        TruncateOneProfileRetainOtherProfile(nCoast, pSecondProfile, pFirstProfile, dIntersectX, dIntersectY, nProf2LineSeg, nProf1LineSeg, true);
                     }
                     // Is the point of intersection already present in the second profile?
                     else if (pSecondProfile->bIsPointInProfile(dIntersectX, dIntersectY, nPoint))
                     {
                        LogStream << m_ulIter << ":\t coast " << nCoast << " profiles " << nFirstProfile << " and " << nSecondProfile << " intersect, but point {" << dIntersectX << ", " << dIntersectY << "} is already present in profile " << nSecondProfile << " as point " << nPoint << endl;

                        // Truncate the first profile and merge it with the second profile
                        TruncateOneProfileRetainOtherProfile(nCoast, pFirstProfile, pSecondProfile, dIntersectX, dIntersectY, nProf1LineSeg, nProf2LineSeg, true);
                     }
                     else
                     {
                        // The point of intersection is not already present in either profile, so get the number of line segments of each profile
                        int const nFirstProfileLineSegments = pFirstProfile->nGetNumLineSegments();
                        int const nSecondProfileLineSegments = pSecondProfile->nGetNumLineSegments();

                        // assert(nProf1LineSeg < nFirstProfileLineSegments);
                        // assert(nProf2LineSeg < nSecondProfileLineSegments);

                        // Next check whether the point of intersection is on the final line segment of both profiles
                        if ((nProf1LineSeg == (nFirstProfileLineSegments - 1)) && (nProf2LineSeg == (nSecondProfileLineSegments - 1)))
                        {
                           // Yes, the point of intersection is on the final line segment of both profiles, so merge the profiles seaward of the point of intersection
                           MergeProfilesAtFinalLineSegments(nCoast, pFirstProfile, pSecondProfile, nFirstProfileLineSegments, nSecondProfileLineSegments, dIntersectX, dIntersectY, dAvgEndX, dAvgEndY);

                           // LogStream << m_ulIter << ": " << ((nDirection == DIRECTION_DOWNCOAST) ? "down" : "up") << "-coast search, end-segment intersection between profiles " << nFirstProfile << " and " << nSecondProfile << " at [" << dIntersectX << ", " << dIntersectY << "] in line segment [" << nProf1LineSeg << "] of " << nFirstProfileLineSegments << " segments, and line segment [" << nProf2LineSeg << "] of " << nSecondProfileLineSegments << " segments, respectively" << endl;

                           // // DEBUG CODE =============================================================================================
                           // int nSizeTmp = pFirstProfile->nGetProfileSize();
                           // CGeom2DPoint PtEndTmp = *pFirstProfile->pPtGetPointInProfile(nSizeTmp-1);
                           //
                           // LogStream << m_ulIter << ": end of first profile (" << nFirstProfile << ") is point " << nSizeTmp-1 << " at [" << dExtCRSXToGridX(PtEndTmp.dGetX()) << "][" << dExtCRSYToGridY(PtEndTmp.dGetY()) << "} = {" << PtEndTmp.dGetX() << ", " << PtEndTmp.dGetY() << "}" << endl;
                           //
                           // nSizeTmp = pSecondProfile->nGetProfileSize();
                           // PtEndTmp = *pSecondProfile->pPtGetPointInProfile(nSizeTmp-1);
                           //
                           // LogStream << m_ulIter << ": end of second profile (" << nSecondProfile << ") is point " << nSizeTmp-1 << " at [" << dExtCRSXToGridX(PtEndTmp.dGetX()) << "][" << dExtCRSYToGridY(PtEndTmp.dGetY()) << "} = {" << PtEndTmp.dGetX() << ", " << PtEndTmp.dGetY() << "}" << endl;
                           // // DEBUG CODE =============================================================================================
                        }
                        else
                        {
                           // The profiles intersect, but the point of intersection is not on the final line segment of both profiles. One of the profiles will be truncated, the other profile will be retained
                           // LogStream << m_ulIter << ": " << ((nDirection == DIRECTION_DOWNCOAST) ? "down" : "up") << "-coast search, intersection (NOT both end segments) between profiles " << nFirstProfile << " and " << nSecondProfile << " at [" << dIntersectX << ", " << dIntersectY << "] in line segment [" << nProf1LineSeg << "] of " << nFirstProfileLineSegments << ", and line segment [" << nProf2LineSeg << "] of " << nSecondProfileLineSegments << ", respectively" << endl;

                           // Decide which profile to truncate, and which to retain
                           if (pFirstProfile->bIsIntervention())
                           {
                              // Truncate the first profile, since it is an intervention profile
                              // LogStream << m_ulIter << ": pFirstProfile = " << pFirstProfile->nGetProfileID() << " pSecondProfile = " << pSecondProfile->nGetProfileID() << ", pFirstProfile is an intervention profile, so truncate pFirstProfile (" << pFirstProfile->nGetProfileID() << ")" << endl;

                              TruncateOneProfileRetainOtherProfile(nCoast, pFirstProfile, pSecondProfile, dIntersectX, dIntersectY, nProf1LineSeg, nProf2LineSeg, false);
                           }
                           else if (pSecondProfile->bIsIntervention())
                           {
                              // Truncate the second profile, since it is an intervention profile
                              // LogStream << m_ulIter << ": pFirstProfile = " << pFirstProfile->nGetProfileID() << " pSecondProfile = " << pSecondProfile->nGetProfileID() << ", pSecondProfile is an intervention profile, so truncate pSecondProfile (" << pSecondProfile->nGetProfileID() << ")" << endl;

                              TruncateOneProfileRetainOtherProfile(nCoast, pSecondProfile, pFirstProfile, dIntersectX, dIntersectY, nProf2LineSeg, nProf1LineSeg, false);
                           }
                           else if (nFirstProfileLineSegments < nSecondProfileLineSegments)
                           {
                              // Truncate the first profile, since it has a smaller number of line segments
                              // LogStream << m_ulIter << ": pFirstProfile = " << pFirstProfile->nGetProfileID() << " pSecondProfile = " << pSecondProfile->nGetProfileID() << ", pFirstProfile has a smaller number of line segments, so truncate pFirstProfile (" << pFirstProfile->nGetProfileID() << ")" << endl;

                              TruncateOneProfileRetainOtherProfile(nCoast, pFirstProfile, pSecondProfile, dIntersectX, dIntersectY, nProf1LineSeg, nProf2LineSeg, false);
                           }
                           else if (nFirstProfileLineSegments > nSecondProfileLineSegments)
                           {
                              // Truncate the second profile, since it has a smaller number of line segments
                              // LogStream << m_ulIter << ": pFirstProfile = " << pFirstProfile->nGetProfileID() << " pSecondProfile = " << pSecondProfile->nGetProfileID() << ", pSecondProfile has a smaller number of line segments, so truncate pSecondProfile (" << pSecondProfile->nGetProfileID() << ")" << endl;

                              TruncateOneProfileRetainOtherProfile(nCoast, pSecondProfile, pFirstProfile, dIntersectX, dIntersectY, nProf2LineSeg, nProf1LineSeg, false);
                           }
                           else
                           {
                              // Both profiles have the same number of line segments, so choose randomly. Draw a sample from the unit normal distribution using random number generator 1
                              double const dRand = m_dGetFromUnitNormalDist(m_Rand[0]);

                              if (dRand >= 0.0)
                              {
                                 // LogStream << m_ulIter << ": pFirstProfile = " << pFirstProfile->nGetProfileID() << " pSecondProfile = " << pSecondProfile->nGetProfileID() << ", same number of line segment, randomly truncate pFirstProfile" << endl;

                                 TruncateOneProfileRetainOtherProfile(nCoast, pFirstProfile, pSecondProfile, dIntersectX, dIntersectY, nProf1LineSeg, nProf2LineSeg, false);
                              }
                              else
                              {
                                 // LogStream << m_ulIter << ": pFirstProfile = " << pFirstProfile->nGetProfileID() << " pSecondProfile = " << pSecondProfile->nGetProfileID() << ", same number of line segment, randomly truncate pSecondProfile" << endl;

                                 TruncateOneProfileRetainOtherProfile(nCoast, pSecondProfile, pFirstProfile, dIntersectX, dIntersectY, nProf2LineSeg, nProf1LineSeg, false);
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
}

//===============================================================================================================================
//! Check all coastline-normal profiles and modify the profiles if they intersect, then mark valid profiles on the raster grid
//===============================================================================================================================
int CSimulation::nCheckAndMarkAllProfiles(void)
{
   // Check to see which coastline-normal profiles intersect. Then modify intersecting profiles so that the sections of each profile seaward of the point of intersection are 'shared' i.e. are multi-lines. This creates the boundaries of the triangular polygons
   CheckForIntersectingProfiles();

   // Again check the normal profiles for insufficient length: is the water depth at the end point less than the depth of closure? We do this again because some profiles may have been shortened as a result of intersection. Do once for every coastline object
   for (unsigned int nCoast = 0; nCoast < m_VCoast.size(); nCoast++)
   {
      for (int n = 0; n < m_VCoast[nCoast].nGetNumProfiles(); n++)
      {
         CGeomProfile* pProfile = m_VCoast[nCoast].pGetProfileWithDownCoastSeq(n);
         int const nProfile = pProfile->nGetProfileID();

         if (pProfile->bProfileOK())
         {
            int const nSize = pProfile->nGetProfileSize();

            // Safety check
            if (nSize == 0)
            {
               // pProfile->SetTooShort(true);
               m_VCoast[nCoast].pGetProfile(nProfile)->SetTooShort(true);
               LogStream << "Profile " << nProfile << " is too short, size = " << nSize << endl;
               continue;
            }

            CGeom2DPoint const* pPtEnd = pProfile->pPtGetPointInProfile(nSize - 1);
            CGeom2DIPoint const PtiEnd = PtiExtCRSToGridRound(pPtEnd);
            int nXEnd = PtiEnd.nGetX();
            int nYEnd = PtiEnd.nGetY();

            // Safety checks: the point may be outside the grid, so keep it within the grid
            nXEnd = tMin(nXEnd, m_nXGridSize - 1);
            nYEnd = tMin(nYEnd, m_nYGridSize - 1);
            nXEnd = tMax(nXEnd, 0);
            nYEnd = tMax(nYEnd, 0);

            if (m_pRasterGrid->m_Cell[nXEnd][nYEnd].dGetSeaDepth() < m_dDepthOfClosure)
            {
               if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
                  LogStream << m_ulIter << ": coast " << nCoast << ", profile " << nProfile << " is invalid, is too short for depth of closure " << m_dDepthOfClosure << " at end point [" << nXEnd << "][" << nYEnd << "] = {" << pPtEnd->dGetX() << ", " << pPtEnd->dGetY() << "}, flagging as too short" << endl;

               // pProfile->SetTooShort(true);
               m_VCoast[nCoast].pGetProfile(nProfile)->SetTooShort(true);
            }
         }
      }

      // For this coast, put all valid coastline-normal profiles (apart from the profiles at the start and end of the coast, since they have already been done) onto the raster grid. But if the profile is not long enough, crosses a coastline, hits dry land, or hits another profile, then mark the profile as invalid
      int nValidProfiles = 0;
      MarkProfilesOnGrid(nCoast, nValidProfiles);

      if (nValidProfiles == 0)
      {
         // Problem! No valid profiles, so quit
         cerr << m_ulIter << ": " << ERR << "no coastline-normal profiles created" << endl;
         return RTN_ERR_NO_PROFILES_2;
      }

      // // DEBUG CODE ===========================================================================================================
      // if (m_ulIter == 109)
      // {
      // string strOutFile = m_strOutPath;
      // strOutFile += "00_profile_raster_";
      // strOutFile += to_string(m_ulIter);
      // strOutFile += ".tif";
      //
      // GDALDriver* pDriver = GetGDALDriverManager()->GetDriverByName("gtiff");
      // GDALDataset* pDataSet = pDriver->Create(strOutFile.c_str(), m_nXGridSize, m_nYGridSize, 1, GDT_Float64, m_papszGDALRasterOptions);
      // pDataSet->SetProjection(m_strGDALBasementDEMProjection.c_str());
      // pDataSet->SetGeoTransform(m_dGeoTransform);
      //
      // int nn = 0;
      // double* pdRaster = new double[m_nXGridSize * m_nYGridSize];
      // for (int nY = 0; nY < m_nYGridSize; nY++)
      // {
      // for (int nX = 0; nX < m_nXGridSize; nX++)
      // {
      // if (m_pRasterGrid->m_Cell[nX][nY].bIsCoastline())
      // pdRaster[nn] = -1;
      // else
      // {
      //
      // //         pdRaster[nn] = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonID();
      // pdRaster[nn] = m_pRasterGrid->m_Cell[nX][nY].nGetProfileID();
      // }
      //
      // nn++;
      // }
      // }
      //
      // GDALRasterBand* pBand = pDataSet->GetRasterBand(1);
      // pBand->SetNoDataValue(m_dMissingValue);
      // int nRet1 = pBand->RasterIO(GF_Write, 0, 0, m_nXGridSize, m_nYGridSize, pdRaster, m_nXGridSize, m_nYGridSize, GDT_Float64, 0, 0, NULL);
      //
      // if (nRet1 == CE_Failure)
      // return RTN_ERR_GRIDCREATE;
      //
      // GDALClose(pDataSet);
      // delete[] pdRaster;
      // }
      // // DEBUG CODE ===========================================================================================================
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Checks all line segments of a pair of coastline-normal profiles for intersection. If the lines intersect, returns true with the numbers of the line segments at which intersection occurs in nProfile1LineSegment and nProfile1LineSegment, the intersection point in dXIntersect and dYIntersect, and the 'average' seaward endpoint of the two intersecting profiles at dXAvgEnd and dYAvgEnd
//===============================================================================================================================
bool CSimulation::bCheckForIntersection(CGeomProfile* const pVProfile1, CGeomProfile* const pVProfile2, int& nProfile1LineSegment, int& nProfile2LineSegment, double& dXIntersect, double& dYIntersect, double& dXAvgEnd, double& dYAvgEnd)
{
   // For both profiles, look at all line segments
   int const nProfile1NumSegments = pVProfile1->nGetNumLineSegments();
   int const nProfile2NumSegments = pVProfile2->nGetNumLineSegments();
   // nProfile1Size = pVProfile1->nGetProfileSize(),
   // nProfile2Size = pVProfile2->nGetProfileSize();

   // assert(nProfile1Size == nProfile1NumSegments+1);
   // assert(nProfile2Size == nProfile2NumSegments+1);

   for (int i = 0; i < nProfile1NumSegments; i++)
   {
      for (int j = 0; j < nProfile2NumSegments; j++)
      {
         // In external coordinates
         double const dX1 = pVProfile1->pPtVGetPoints()->at(i).dGetX();
         double const dY1 = pVProfile1->pPtVGetPoints()->at(i).dGetY();
         double const dX2 = pVProfile1->pPtVGetPoints()->at(i + 1).dGetX();
         double const dY2 = pVProfile1->pPtVGetPoints()->at(i + 1).dGetY();

         double const dX3 = pVProfile2->pPtVGetPoints()->at(j).dGetX();
         double const dY3 = pVProfile2->pPtVGetPoints()->at(j).dGetY();
         double const dX4 = pVProfile2->pPtVGetPoints()->at(j + 1).dGetX();
         double const dY4 = pVProfile2->pPtVGetPoints()->at(j + 1).dGetY();

         // Uses Cramer's Rule to solve the equations. Modified from code at http://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect (in turn based on Andre LeMothe's "Tricks of the Windows Game Programming Gurus")
         double const dDiffX1 = dX2 - dX1;
         double const dDiffY1 = dY2 - dY1;
         double const dDiffX2 = dX4 - dX3;
         double const dDiffY2 = dY4 - dY3;

         double dS = -999;
         double dT = -999;
         double dTmp = 0;

         dTmp = -dDiffX2 * dDiffY1 + dDiffX1 * dDiffY2;

         if (! bFPIsEqual(dTmp, 0.0, TOLERANCE))
            dS = (-dDiffY1 * (dX1 - dX3) + dDiffX1 * (dY1 - dY3)) / dTmp;

         dTmp = -dDiffX2 * dDiffY1 + dDiffX1 * dDiffY2;

         if (! bFPIsEqual(dTmp, 0.0, TOLERANCE))
            dT = (dDiffX2 * (dY1 - dY3) - dDiffY2 * (dX1 - dX3)) / dTmp;

         if (dS >= 0 && dS <= 1 && dT >= 0 && dT <= 1)
         {
            // Collision detected, calculate intersection coordinates
            dXIntersect = dX1 + (dT * dDiffX1);
            dYIntersect = dY1 + (dT * dDiffY1);

            // And calc the average end-point coordinates
            dXAvgEnd = (dX2 + dX4) / 2;
            dYAvgEnd = (dY2 + dY4) / 2;

            // Get the line segments at which intersection occurred
            nProfile1LineSegment = i;
            nProfile2LineSegment = j;

            // LogStream << "\t" << "INTERSECTION dX2 = " << dX2 << " dX4 = " << dX4 << " dY2 = " << dY2 << " dY4 = " << dY4 << endl;
            return true;
         }
      }
   }

   // No intersection
   return false;
}

//===============================================================================================================================
//! For this coastline, marks all coastline-normal profiles (apart from the two 'special' ones at the start and end of the coast) onto the raster grid, i.e. rasterizes multi-line vector objects onto the raster grid. Note that this doesn't work if the vector has already been interpolated to fit on the grid i.e. if distances between vector points are just one cell apart
//===============================================================================================================================
void CSimulation::MarkProfilesOnGrid(int const nCoast, int& nValidProfiles)
{
   // How many profiles on this coast?
   int const nProfiles = m_VCoast[nCoast].nGetNumProfiles();

   if (nProfiles == 0)
   {
      // This can happen if the coastline is very short, so just give a warning and carry on with the next coastline
      if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
         LogStream << WARN << m_ulIter << ": coast " << nCoast << " has no profiles" << endl;

      return;
   }

   static bool bDownCoast = true;

   // Now do this for every profile, alternate between up-coast and down-coast directions
   for (int n = 0; n < nProfiles; n++)
   {
      CGeomProfile* pProfile;

      if (bDownCoast)
         pProfile = m_VCoast[nCoast].pGetProfileWithDownCoastSeq(n);
      else
         pProfile = m_VCoast[nCoast].pGetProfileWithUpCoastSeq(n);

      // Don't do this for the first and last profiles (i.e. the profiles at the start and end of the coast) since these are put onto the grid elsewhere
      if (pProfile->bIsGridEdge())
         continue;

      int const nProfile = pProfile->nGetProfileID();

      // If this profile has a problem, then forget about it
      // if (! pProfile->bProfileOK())
      // {
      // LogStream << m_ulIter << ": in MarkProfilesOnGrid() profile " << nProfile << " is not OK" << endl;
      // continue;
      // }

      int const nPoints = pProfile->nGetProfileSize();

      if (nPoints < 2)
      {
         // Need at least two points in the profile, so this profile is invalid: mark it
         m_VCoast[nCoast].pGetProfile(nProfile)->SetTooShort(true);

         if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
            LogStream << m_ulIter << ": coast " << nCoast << ", profile " << nProfile << " is invalid, has only " << nPoints << " points" << endl;

         continue;
      }

      // OK, go for it: set up temporary vectors to hold the x-y coords (in grid CRS) of the cells which we will mark
      vector<CGeom2DIPoint> VCellsToMark;
      vector<bool> bVShared;
      bool bTooShort = false;
      bool bTruncatedSameCoast = false;
      bool bHitCoast = false;
      bool bHitLand = false;
      bool bHitIntervention = false;
      bool bHitAnotherProfile = false;

      CreateRasterizedProfile(nCoast, pProfile, &VCellsToMark, &bVShared, bTooShort, bTruncatedSameCoast, bHitCoast, bHitLand, bHitIntervention, bHitAnotherProfile);

      if ((bTruncatedSameCoast && (! ACCEPT_TRUNCATED_PROFILES)) || bTooShort || bHitCoast || bHitLand || bHitIntervention || bHitAnotherProfile || VCellsToMark.size() == 0)
         continue;

      // This profile is fine
      nValidProfiles++;

      for (unsigned int k = 0; k < VCellsToMark.size(); k++)
      {
         // Ignore duplicate points
         if ((k > 0) && (VCellsToMark[k] == m_VCoast[nCoast].pGetProfile(nProfile)->pPtiGetLastCellInProfile()))
            continue;

         // Mark each cell in the raster grid
         int const nXTmp = VCellsToMark[k].nGetX();
         int const nYTmp = VCellsToMark[k].nGetY();
         m_pRasterGrid->m_Cell[nXTmp][nYTmp].SetCoastAndProfileID(nCoast, nProfile);

         // Store the raster grid coordinates in the profile object
         m_VCoast[nCoast].pGetProfile(nProfile)->AppendCellInProfile(nXTmp, nYTmp);

         // Mark the shared (i.e. multi-line) parts of the profile (if any)
         // if (bVShared[k])
         // {
         // m_VCoast[nCoast].pGetProfile(nProfile)->AppendPointShared(true);
         // //                  LogStream << m_ulIter << ": profile " << j << " point " << k << " marked as shared" << endl;
         // }
         // else
         // {
         // m_VCoast[nCoast].pGetProfile(nProfile)->AppendPointShared(false);
         // //                  LogStream << m_ulIter << ": profile " << nProfile << " point " << k << " marked as NOT shared" << endl;
         // }
      }

      // Get the deep water wave height and orientation values at the end of the profile
      double const dDeepWaterWaveHeight = m_pRasterGrid->m_Cell[VCellsToMark.back().nGetX()][VCellsToMark.back().nGetY()].dGetCellDeepWaterWaveHeight();
      double const dDeepWaterWaveAngle = m_pRasterGrid->m_Cell[VCellsToMark.back().nGetX()][VCellsToMark.back().nGetY()].dGetCellDeepWaterWaveAngle();
      double const dDeepWaterWavePeriod = m_pRasterGrid->m_Cell[VCellsToMark.back().nGetX()][VCellsToMark.back().nGetY()].dGetCellDeepWaterWavePeriod();

      // And store them for this profile
      m_VCoast[nCoast].pGetProfile(nProfile)->SetProfileDeepWaterWaveHeight(dDeepWaterWaveHeight);
      m_VCoast[nCoast].pGetProfile(nProfile)->SetProfileDeepWaterWaveAngle(dDeepWaterWaveAngle);
      m_VCoast[nCoast].pGetProfile(nProfile)->SetProfileDeepWaterWavePeriod(dDeepWaterWavePeriod);
   }

   bDownCoast = ! bDownCoast;
}

//===============================================================================================================================
//! Given a pointer to a coastline-normal profile, returns an output vector of cells which are 'under' every line segment of the profile. If there is a problem with the profile (e.g. a rasterized cell is dry land or coast, or the profile has to be truncated) then we pass this back as an error code
//===============================================================================================================================
void CSimulation::CreateRasterizedProfile(int const nCoast, CGeomProfile* pProfile, vector<CGeom2DIPoint>* pVIPointsOut, vector<bool>* pbVShared, bool& bTooShort, bool& bTruncatedSameCoast, bool& bHitCoast, bool& bHitLand, bool& bHitIntervention, bool& bHitAnotherProfile)
{
   int const nProfile = pProfile->nGetProfileID();
   int nSeg = 0;
   int const nNumSegments = pProfile->nGetNumLineSegments();

   pVIPointsOut->clear();

   // LogStream << m_ulIter << ": in CreateRasterizedProfile() *pPtiStart for profile " << nProfile << " is [" << pPtiStart->nGetX() << "][" << pPtiStart->nGetY() << "]" << endl;
   int nXStartLast = INT_NODATA;
   int nYStartLast = INT_NODATA;
   int nXEndLast = INT_NODATA;
   int nYEndLast = INT_NODATA;

   // Do for every segment of this profile
   for (nSeg = 0; nSeg < nNumSegments; nSeg++)
   {
      // Do once for every line segment
      CGeom2DIPoint PtiSegStart;

      if (nSeg == 0)
      {
         // If this is the first segment, use the coastline start point to prevent external CRS to grid CRS rounding errors
         int const nCoastPoint = pProfile->nGetCoastPoint();
         PtiSegStart = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint);
      }
      else
      {
         CGeom2DPoint const* pPtSegStart = pProfile->pPtGetPointInProfile(nSeg);

         // Convert from the external CRS to grid CRS
         PtiSegStart = PtiExtCRSToGridRound(pPtSegStart);
      }

      CGeom2DPoint const* pPtSegEnd = pProfile->pPtGetPointInProfile(nSeg + 1); // This is OK

      // Convert from the external CRS to grid CRS
      CGeom2DIPoint const PtiSegEnd = PtiExtCRSToGridRound(pPtSegEnd);

      // Safety check
      if (PtiSegStart == PtiSegEnd)
         continue;

      int const nXStart = PtiSegStart.nGetX();
      int const nYStart = PtiSegStart.nGetY();
      int const nXEnd = PtiSegEnd.nGetX();
      int const nYEnd = PtiSegEnd.nGetY();

      bool bShared = false;

      if (pProfile->nGetNumCoincidentProfilesInLineSegment(nSeg) > 1)
      {
         bShared = true;

         // If this is the second or more of several coincident line segments (i.e. it has the same start and end points as the previous line segment) then ignore it
         if ((nXStart == nXStartLast) && (nYStart == nYStartLast) && (nXEnd == nXEndLast) && (nYEnd == nYEndLast))
            continue;
      }

      // Interpolate between cells by a simple DDA line algorithm, see http://en.wikipedia.org/wiki/Digital_differential_analyzer_(graphics_algorithm) Note that Bresenham's algorithm gave occasional gaps
      double dXInc = nXEnd - nXStart;
      double dYInc = nYEnd - nYStart;
      double const dLength = tMax(tAbs(dXInc), tAbs(dYInc));

      dXInc /= dLength;
      dYInc /= dLength;

      double dX = nXStart;
      double dY = nYStart;

      // Process each interpolated point
      for (int m = 0; m <= nRound(dLength); m++)
      {
         int const nX = nRound(dX);
         int const nY = nRound(dY);

         // Do some checking of this interpolated point, but only if this is not a grid-edge profile (these profiles are always valid)
         if (! pProfile->bIsGridEdge())
         {
            // Is the interpolated point within the valid raster grid?
            if (! bIsWithinValidGrid(nX, nY))
            {
               // It is outside the valid grid, so mark this profile and quit the loop
               bTruncatedSameCoast = true;

               if (! ACCEPT_TRUNCATED_PROFILES)
                  pProfile->SetTruncatedSameCoast(true);

               if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
                  LogStream << m_ulIter << ": profile " << nProfile << " is invalid, truncated at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;

               break;
            }

            // Check again: is this cell (or an adjacent cell: does not matter which) already marked as 'under' a profile?
            int nYTmp = nY+1;
            if (nY+1 >= m_nYGridSize)
               nYTmp = nY-1;

            if (m_pRasterGrid->m_Cell[nX][nY].bIsProfile() || m_pRasterGrid->m_Cell[nX][nYTmp].bIsProfile())
            {
               // This cell or an adjacent cell, is 'under' a profile, so now check if the profile belongs to another coast
               int const nHitProfileCoast1 = m_pRasterGrid->m_Cell[nX][nY].nGetProfileCoastID();
               int const nHitProfileCoast2 = m_pRasterGrid->m_Cell[nX][nYTmp].nGetProfileCoastID();

               if ((nHitProfileCoast1 == nCoast) || (nHitProfileCoast2 == nCoast))
               {
                  // The profile belongs to the same coast, mark this profile as invalid
                  bHitAnotherProfile = true;
                  pProfile->SetHitAnotherProfile(true);
                  return;
               }
            }

            // If this is the first line segment of the profile, then once we are clear of the coastline (when m > PROFILE_CHECK_DIST_FROM_COAST), check if this profile hits land at this interpolated point. NOTE Get problems here since if the coastline vector has been heavily smoothed, this can result is 'false positives' profiles marked as invalid which are not actually invalid, because the profile hits land when m = 0 or m = 1. This results in some cells being flagged as profile cells which are actually inland
            if (m > PROFILE_CHECK_DIST_FROM_COAST)
            {
               // Check this cell. Two diagonal(ish) raster lines can cross each other without any intersection, so must also test an adjacent cell for intersection (does not matter which adjacent cell)
               if ((m_pRasterGrid->m_Cell[nX][nY].bIsCoastline()) || (bIsWithinValidGrid(nX, nY + 1) && m_pRasterGrid->m_Cell[nX][nY + 1].bIsCoastline()))
               {
                  // We've hit a coastline so set a switch and mark the profile, then quit
                  bHitCoast = true;
                  pProfile->SetHitCoast(true);
                  int const nHitCoast = m_pRasterGrid->m_Cell[nX][nY].nGetCoastline();

                  if (m_nLogFileDetail >= LOG_FILE_ALL)
                     LogStream << m_ulIter << ":\t coast " << nCoast << " profile " << nProfile << " is invalid, hit coast " << nHitCoast << " at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;

                  return;
               }

               if (! m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea())
               {
                  // We've hit dry land, so set a switch and mark the profile
                  bHitLand = true;
                  pProfile->SetHitLand(true);

                  LogStream << m_ulIter << ":\t coast " << nCoast << " profile " << nProfile << " HIT LAND at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}, elevation = " << m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevIncTalus() << ", SWL = " << m_dThisIterSWL << endl;

                  return;
               }

               if (m_pRasterGrid->m_Cell[nX][nY].nGetInterventionClass() != INT_NODATA)
               {
                  // We've hit an intervention, so set a switch and mark the profile
                  bHitIntervention = true;
                  pProfile->SetHitIntervention(true);

                  LogStream << m_ulIter << ":\t coast " << nCoast << " profile " << nProfile << " HIT INTERVENTION at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}, elevation = " << m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevIncTalus() << ", SWL = " << m_dThisIterSWL << endl;

                  return;
               }
            }

            // Now check to see if we hit another profile which is not a coincident normal to this normal
            if (m_pRasterGrid->m_Cell[nX][nY].bIsProfile())
            {
               // We've hit a raster cell which is already marked as 'under' a normal profile. Get the number of the profile which marked this cell, and the coast to hich this profile belongs
               int const nHitProfile = m_pRasterGrid->m_Cell[nX][nY].nGetProfileID();
               int const nHitProfileCoast = m_pRasterGrid->m_Cell[nX][nY].nGetProfileCoastID();

               // Do both profiles belong to the same coast?
               if (nCoast == nHitProfileCoast)
               {
                  // Both profiles belong to the same coast. Is this the number of a coincident profile of this profile?
                  if (! pProfile->bFindProfileInCoincidentProfilesOfLastLineSegment(nHitProfile))
                  {
                     // It isn't a coincident profile, so we have just hit an unrelated profile. Mark this profile as invalid and move on
                     pProfile->SetHitAnotherProfile(true);
                     bHitAnotherProfile = true;

                     return;
                  }
               }
            }
         }

         // Append this point to the output vector
         pVIPointsOut->push_back(CGeom2DIPoint(nX, nY)); // Is in raster grid coordinates
         pbVShared->push_back(bShared);

         // And increment for next time
         dX += dXInc;
         dY += dYInc;
      }

      nXStartLast = nXStart;
      nYStartLast = nYStart;
      nXEndLast = nXEnd;
      nYEndLast = nYEnd;

      if (bTruncatedSameCoast)
         break;
   }

   if (bTruncatedSameCoast)
   {
      if (nSeg < (nNumSegments - 1))
         // We are truncating the profile, so remove any line segments after this one
         pProfile->TruncateLineSegments(nSeg);

      // Shorten the vector input. Ignore CPPCheck errors here, since we know that pVIPointsOut is not empty
      int const nLastX = pVIPointsOut->at(pVIPointsOut->size() - 1).nGetX();
      int const nLastY = pVIPointsOut->at(pVIPointsOut->size() - 1).nGetY();

      pProfile->pPtGetPointInProfile(nSeg + 1)->SetX(dGridCentroidXToExtCRSX(nLastX));
      pProfile->pPtGetPointInProfile(nSeg + 1)->SetY(dGridCentroidYToExtCRSY(nLastY));
   }

   // // DEBUG CODE =====================================================================================
   // LogStream << "====================" << endl;
   // LogStream << m_ulIter << ": for profile " << nProfile << " pPtiStart = [" << pPtiStart->nGetX() << "][" << pPtiStart->nGetY() << "] pPtiEnd = [" << pPtiEnd->nGetX() << "][" << pPtiEnd->nGetY() << "] pVIPointsOut->size() = " << pVIPointsOut->size() << endl;
   // // for (int n = 0; n < static_cast<int>(pVIPointsOut->size()); n++)
   // //    LogStream << "\t[" << pVIPointsOut->at(n).nGetX() << "][" << pVIPointsOut->at(n).nGetY() << "]" << endl;
   // LogStream << "====================" << endl;
   // // DEBUG CODE =====================================================================================

   if (pVIPointsOut->size() < 3)
   {
      // Coastline-normal profiles cannot be very short (e.g. with less than 3 cells), since we cannot calculate along-profile slope properly for such short profiles
      bTooShort = true;
      pProfile->SetTooShort(true);

      if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
      {
         // Ignore CPPCheck errors here, since we know that pVIPointsOut is not empty
         LogStream << m_ulIter << ": profile " << nProfile << " is invalid, is too short, only " << pVIPointsOut->size() << " points, HitLand?" << bHitLand << ". From [" << pVIPointsOut->at(0).nGetX() << "][" << pVIPointsOut->at(0).nGetY() << "] = {" << dGridCentroidXToExtCRSX(pVIPointsOut->at(0).nGetX()) << ", " << dGridCentroidYToExtCRSY(pVIPointsOut->at(0).nGetY()) << "} to [" << pVIPointsOut->at(pVIPointsOut->size() - 1).nGetX() << "][" << pVIPointsOut->at(pVIPointsOut->size() - 1).nGetY() << "] = {" << dGridCentroidXToExtCRSX(pVIPointsOut->at(pVIPointsOut->size() - 1).nGetX()) << ", " << dGridCentroidYToExtCRSY(pVIPointsOut->at(pVIPointsOut->size() - 1).nGetY()) << "}" << endl;
      }
   }
}

//===============================================================================================================================
//! Merges two profiles which intersect at their final (most seaward) line segments, seaward of their point of intersection
//===============================================================================================================================
void CSimulation::MergeProfilesAtFinalLineSegments(int const nCoast, CGeomProfile* pFirstProfile, CGeomProfile* pSecondProfile, int const nFirstProfileLineSegments, int const nSecondProfileLineSegments, double const dIntersectX, double const dIntersectY, double const dAvgEndX, double const dAvgEndY)
{
   // The point of intersection is on the final (most seaward) line segment of both profiles. Put together a vector of coincident profile numbers (with no duplicates) for both profiles
   int nCombinedLastSeg = 0;
   vector<pair<int, int>> prVCombinedProfilesCoincidentProfilesLastSeg;

   for (unsigned int n = 0; n < pFirstProfile->pprVGetPairedCoincidentProfilesForLineSegment(nFirstProfileLineSegments - 1)->size(); n++)
   {
      pair<int, int> prTmp;
      prTmp.first = pFirstProfile->pprVGetPairedCoincidentProfilesForLineSegment(nFirstProfileLineSegments - 1)->at(n).first;
      prTmp.second = pFirstProfile->pprVGetPairedCoincidentProfilesForLineSegment(nFirstProfileLineSegments - 1)->at(n).second;

      bool bFound = false;

      for (unsigned int m = 0; m < prVCombinedProfilesCoincidentProfilesLastSeg.size(); m++)
      {
         if (prVCombinedProfilesCoincidentProfilesLastSeg[m].first == prTmp.first)
         {
            bFound = true;
            break;
         }
      }

      if (! bFound)
      {
         prVCombinedProfilesCoincidentProfilesLastSeg.push_back(prTmp);
         nCombinedLastSeg++;
      }
   }

   for (unsigned int n = 0; n < pSecondProfile->pprVGetPairedCoincidentProfilesForLineSegment(nSecondProfileLineSegments - 1)->size(); n++)
   {
      pair<int, int> prTmp;
      prTmp.first = pSecondProfile->pprVGetPairedCoincidentProfilesForLineSegment(nSecondProfileLineSegments - 1)->at(n).first;
      prTmp.second = pSecondProfile->pprVGetPairedCoincidentProfilesForLineSegment(nSecondProfileLineSegments - 1)->at(n).second;

      bool bFound = false;

      for (unsigned int m = 0; m < prVCombinedProfilesCoincidentProfilesLastSeg.size(); m++)
      {
         if (prVCombinedProfilesCoincidentProfilesLastSeg[m].first == prTmp.first)
         {
            bFound = true;
            break;
         }
      }

      if (! bFound)
      {
         prVCombinedProfilesCoincidentProfilesLastSeg.push_back(prTmp);
         nCombinedLastSeg++;
      }
   }

   // Increment the number of each line segment
   for (int m = 0; m < nCombinedLastSeg; m++)
      prVCombinedProfilesCoincidentProfilesLastSeg[m].second++;

   vector<pair<int, int>> prVFirstProfileCoincidentProfilesLastSeg = *pFirstProfile->pprVGetPairedCoincidentProfilesForLineSegment(nFirstProfileLineSegments - 1);
   vector<pair<int, int>> prVSecondProfileCoincidentProfilesLastSeg = *pSecondProfile->pprVGetPairedCoincidentProfilesForLineSegment(nSecondProfileLineSegments - 1);
   int const nNumFirstProfileCoincidentProfilesLastSeg = static_cast<int>(prVFirstProfileCoincidentProfilesLastSeg.size());
   int const nNumSecondProfileCoincidentProfilesLastSeg = static_cast<int>(prVSecondProfileCoincidentProfilesLastSeg.size());

   // LogStream << m_ulIter << ": END-SEGMENT INTERSECTION between profiles " << nFirstProfile << " and " << nSecondProfile << " at line segment " << nFirstProfileLineSegments-1 << "/" << nFirstProfileLineSegments-1 << ", and line segment " << nSecondProfileLineSegments-1 << "/" << nSecondProfileLineSegments-1 << ", respectively. Both truncated at [" << dIntersectX << ", " << dIntersectY << "] then profiles {" << nFirstProfile << "} and {" << nSecondProfile << "} extended to [" << dAvgEndX << ", " << dAvgEndY << "]" << endl;

   // Truncate the first profile, and all co-incident profiles, at the point of intersection
   for (int n = 0; n < nNumFirstProfileCoincidentProfilesLastSeg; n++)
   {
      int const nThisProfile = prVFirstProfileCoincidentProfilesLastSeg[n].first;
      CGeomProfile* pThisProfile = m_VCoast[nCoast].pGetProfile(nThisProfile);
      int const nProfileLength = pThisProfile->nGetProfileSize();

      // This is the final line segment of the first 'main' profile. We are assuming that it is also the final line segment of all co-incident profiles. This is fine, altho' each profile may well have a different number of line segments landwards i.e. the number of the line segment may be different for each co-incident profile
      pThisProfile->SetPointInProfile(nProfileLength - 1, dIntersectX, dIntersectY);
   }

   // Truncate the second profile, and all co-incident profiles, at the point of intersection
   for (int n = 0; n < nNumSecondProfileCoincidentProfilesLastSeg; n++)
   {
      int const nThisProfile = prVSecondProfileCoincidentProfilesLastSeg[n].first;
      CGeomProfile* pThisProfile = m_VCoast[nCoast].pGetProfile(nThisProfile);
      int const nProfileLength = pThisProfile->nGetProfileSize();

      // This is the final line segment of the second 'main' profile. We are assuming that it is also the final line segment of all co-incident profiles. This is fine, altho' each profile may well have a different number of line segments landwards i.e. the number of the line segment may be different for each co-incident profile
      pThisProfile->SetPointInProfile(nProfileLength - 1, dIntersectX, dIntersectY);
   }

   // Append a new straight line segment to the existing line segment(s) of the first profile, and to all co-incident profiles
   for (int nThisLineSeg = 0; nThisLineSeg < nNumFirstProfileCoincidentProfilesLastSeg; nThisLineSeg++)
   {
      int const nThisProfile = prVFirstProfileCoincidentProfilesLastSeg[nThisLineSeg].first;
      CGeomProfile* pThisProfile = m_VCoast[nCoast].pGetProfile(nThisProfile);

      // Update this profile
      pThisProfile->AppendPointInProfile(dAvgEndX, dAvgEndY);

      // Append details of the combined profiles
      pThisProfile->AppendLineSegment();

      for (int m = 0; m < nCombinedLastSeg; m++)
         pThisProfile->AppendCoincidentProfileToLineSegments(prVCombinedProfilesCoincidentProfilesLastSeg[m]);
   }

   // Append a new straight line segment to the existing line segment(s) of the second profile, and to all co-incident profiles
   for (int nThisLineSeg = 0; nThisLineSeg < nNumSecondProfileCoincidentProfilesLastSeg; nThisLineSeg++)
   {
      int const nThisProfile = prVSecondProfileCoincidentProfilesLastSeg[nThisLineSeg].first;
      CGeomProfile* pThisProfile = m_VCoast[nCoast].pGetProfile(nThisProfile);

      // Update this profile
      pThisProfile->AppendPointInProfile(dAvgEndX, dAvgEndY);

      // Append details of the combined profiles
      pThisProfile->AppendLineSegment();

      for (int m = 0; m < nCombinedLastSeg; m++)
         pThisProfile->AppendCoincidentProfileToLineSegments(prVCombinedProfilesCoincidentProfilesLastSeg[m]);
   }

   // // DEBUG CODE ****************************************************************
   // int nFirstProfileLineSeg= pFirstProfile->nGetNumLineSegments();
   // int nSecondProfileLineSeg = pSecondProfile->nGetNumLineSegments();
   //
   // LogStream << "\tProfile {" << nFirstProfile << "} now has " << nFirstProfileLineSeg << " line segments" << endl;
   // for (int m = 0; m < nFirstProfileLineSeg; m++)
   // {
   // vector<pair<int, int> > prVCoincidentProfiles = *pFirstProfile->pprVGetPairedCoincidentProfilesForLineSegment(m);
   // LogStream << "\tCo-incident profiles and line segments for line segment " << m << " of profile {" << nFirstProfile << "} are {";
   // for (int nn = 0; nn < prVCoincidentProfiles.size(); nn++)
   // LogStream << " " << prVCoincidentProfiles[nn].first << "[" << prVCoincidentProfiles[nn].second << "] ";
   // LogStream << " }" << endl;
   // }
   // LogStream << "\tProfile {" << nSecondProfile << "} now has " << nSecondProfileLineSeg << " line segments" << endl;
   // for (int m = 0; m < nSecondProfileLineSeg; m++)
   // {
   // vector<pair<int, int> > prVCoincidentProfiles = *pSecondProfile->pprVGetPairedCoincidentProfilesForLineSegment(m);
   // LogStream << "\tCo-incident profiles and line segments for line segment " << m << " of profile {" << nSecondProfile << "} are {";
   // for (int nn = 0; nn < prVCoincidentProfiles.size(); nn++)
   // LogStream << " " << prVCoincidentProfiles[nn].first << "[" << prVCoincidentProfiles[nn].second << "] ";
   // LogStream << " }" << endl;
   // }
   // // DEBUG CODE ******************************************************************
}

//===============================================================================================================================
//! Truncates one intersecting profile at the point of intersection, and retains the other profile
//===============================================================================================================================
void CSimulation::TruncateOneProfileRetainOtherProfile(int const nCoast, CGeomProfile* pProfileToTruncate, CGeomProfile* pProfileToRetain, double dIntersectX, double dIntersectY, int nProfileToTruncateIntersectLineSeg, int nProfileToRetainIntersectLineSeg, bool const bAlreadyPresent)
{
   // // Occasionally, profiles cross each other, with the crossing not detected. So check for intersection between pProfileToTruncate and all profiles (starting from the last i.e. in reverse order) in all segments (starting from the last i.e. in reverse order) of pProfileToRetain's CGeomMultiLine
   // bool bFound = false;
   // int nNumSegProfRetain = pProfileToRetain->nGetNumLineSegments();
   // for (int nSeg = nNumSegProfRetain-1; nSeg >= 0; nSeg--)
   // {
   // if (bFound)
   // break;
   //
   // int nNumProfInSeg = pProfileToRetain->nGetNumCoincidentProfilesInLineSegment(nSeg);
   // for (int nProf = nNumProfInSeg-1; nProf >= 0; nProf--)
   // {
   // int nThisProf =  pProfileToRetain->nGetCoincidentProfileForLineSegment(nSeg, nProf);
   // CGeomProfile* pThisProf = m_VCoast[nCoast].pGetProfile(nThisProf);
   //
   // int nProfToTruncLineSeg = 0;
   // int nThisProfLineSeg = 0;
   // double dTmpIntersectX = 0;
   // double dTmpIntersectY = 0;
   // double dAvgEndX = 0;
   // double dAvgEndY = 0;
   //
   // if (bCheckForIntersection(pProfileToTruncate, pThisProf, nProfToTruncLineSeg, nThisProfLineSeg, dTmpIntersectX, dTmpIntersectY, dAvgEndX, dAvgEndY))
   // {
   //          // An intersection was found: so the profile with which pProfileToTruncate intersects becomes the new pProfileToRetain, and dIntersectX, dIntersectY, nProfileToTruncateIntersectLineSeg, and nProfileToRetainIntersectLineSeg are also changed
   // pProfileToRetain = pThisProf;
   // dIntersectX = dTmpIntersectX;
   // dIntersectY = dTmpIntersectY;
   // nProfileToRetainIntersectLineSeg = nThisProfLineSeg;
   // nProfileToTruncateIntersectLineSeg = nProfToTruncLineSeg;
   //
   // bFound = true;
   // break;
   // }
   // }
   // }

   // Insert the intersection point into the main retain-profile if it is not already in the profile, and do the same for all co-incident profiles of the main retain-profile. Also add details of the to-truncate profile (and all its coincident profiles) to every line segment of the main to-retain profile which is seaward of the point of intersection
   int const nRet = nInsertPointIntoProfilesIfNeededThenUpdate(nCoast, pProfileToRetain, dIntersectX, dIntersectY, nProfileToRetainIntersectLineSeg, pProfileToTruncate, nProfileToTruncateIntersectLineSeg, bAlreadyPresent);

   if (nRet != RTN_OK)
   {
      // LogStream << m_ulIter << ": error in nInsertPointIntoProfilesIfNeededThenUpdate()" << endl;
      return;
   }

   // Get all profile points of the main retain-profile seawards from the intersection point, and do the same for the corresponding line segments (including coincident profiles). This also includes details of the main to-truncate profile (and all its coincident profiles)
   vector<CGeom2DPoint> PtVProfileLastPart;
   vector<vector<pair<int, int>>> prVLineSegLastPart;

   if (bAlreadyPresent)
   {
      PtVProfileLastPart = pProfileToRetain->PtVGetThisPointAndAllAfter(nProfileToRetainIntersectLineSeg);
      prVLineSegLastPart = pProfileToRetain->prVVGetAllLineSegAfter(nProfileToRetainIntersectLineSeg);
   }

   else
   {
      PtVProfileLastPart = pProfileToRetain->PtVGetThisPointAndAllAfter(nProfileToRetainIntersectLineSeg + 1);
      prVLineSegLastPart = pProfileToRetain->prVVGetAllLineSegAfter(nProfileToRetainIntersectLineSeg + 1);
   }

   // assert(PtVProfileLastPart.size() > 1);
   // assert(prVLineSegLastPart.size() > 0);

   // Truncate the truncate-profile at the point of intersection, and do the same for all its co-incident profiles. Then append the profile points of the main to-retain profile seaward from the intersection point, and do the same for the corresponding line segments (including coincident profiles)
   TruncateProfileAndAppendNew(nCoast, pProfileToTruncate, nProfileToTruncateIntersectLineSeg, &PtVProfileLastPart, &prVLineSegLastPart);

   // assert(m_VCoast[nCoast].pGetProfile(nProfileToTruncate)->nGetProfileSize() > 1);
   // assert(pProfileToRetain->nGetNumLineSegments() > 0);
   // assert(m_VCoast[nCoast].pGetProfile(nProfileToTruncate)->nGetNumLineSegments() > 0);
}

//===============================================================================================================================
//! Inserts an intersection point into the profile that is to be retained, if that point is not already present in the profile, then does the same for all co-incident profiles. Finally adds the numbers of the to-truncate profile (and all its coincident profiles) to the seaward line segments of the to-retain profile and all its coincident profiles
//===============================================================================================================================
int CSimulation::nInsertPointIntoProfilesIfNeededThenUpdate(int const nCoast, CGeomProfile* pProfileToRetain, double const dIntersectX, double const dIntersectY, int const nProfileToRetainIntersectLineSeg, CGeomProfile* pProfileToTruncate, int const nProfileToTruncateIntersectLineSeg, bool const bAlreadyPresent)
{
   // // DEBUG CODE ****************************************************************
   // // Get the index numbers of all coincident profiles for the 'main' to-retain profile for the line segment in which intersection occurred
   // vector<pair<int, int> > prVRetainCoincidentProfilesCHECK1 = *m_VCoast[nCoast].pGetProfile(nMainProfile)->pprVGetPairedCoincidentProfilesForLineSegment(nProfileToRetainIntersectLineSeg);
   // int nNumRetainCoincidentCHECK1 = prVRetainCoincidentProfilesCHECK1.size();
   // for (int nn = 0; nn < nNumRetainCoincidentCHECK1; nn++)
   // {
   // int nThisProfile = prVRetainCoincidentProfilesCHECK1[nn].first;
   // LogStream << "\tBEFORE nInsertPointIntoProfilesIfNeededThenUpdate(): " << (nThisProfile == nMainProfile ? "MAIN" : "COINCIDENT") << " to-retain profile {" << nThisProfile << "} has " << m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetProfileSize() << " points ";
   // for (int nn = 0; nn < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetProfileSize(); nn++)
   // LogStream << "[" << m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nn)->dGetX() << ", " << m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nn)->dGetY() << "] ";
   // LogStream << "), and " << m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetNumLineSegments() << " line segments, co-incident profiles and their line segments are ";
   // for (int mm = 0; mm < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetNumLineSegments(); mm++)
   // {
   // vector<pair<int, int> > prVTmp = *m_VCoast[nCoast].pGetProfile(nThisProfile)->pprVGetPairedCoincidentProfilesForLineSegment(mm);
   // LogStream << "{ ";
   // for (int nn = 0; nn < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetNumCoincidentProfilesInLineSegment(mm); nn++)
   // LogStream << prVTmp[nn].first << "[" << prVTmp[nn].second << "] ";
   // LogStream << "} ";
   // }
   // LogStream << endl;
   //
   // for (int nPoint = 0; nPoint < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetProfileSize()-1; nPoint++)
   // {
   // CGeom2DPoint
   // Pt1 = *m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nPoint),
   // Pt2 = *m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nPoint+1);
   //
   // if (Pt1 == Pt2)
   // LogStream << m_ulIter << ": IDENTICAL POINTS before changes, in profile {" << nThisProfile << "} points = " << nPoint << " and " << nPoint+1 << endl;
   // }
   // }
   //
   // // Get the index numbers of all coincident profiles for the 'main' to-truncate profile for the line segment in which intersection occurred
   // vector<pair<int, int> > prVTruncateCoincidentProfilesCHECK1 = *m_VCoast[nCoast].pGetProfile(nProfileToTruncate)->pprVGetPairedCoincidentProfilesForLineSegment(nProfileToTruncateIntersectLineSeg);
   // int nNumTruncateCoincidentCHECK1 = prVTruncateCoincidentProfilesCHECK1.size();
   // for (int nn = 0; nn < nNumTruncateCoincidentCHECK1; nn++)
   // {
   // int nThisProfile = prVTruncateCoincidentProfilesCHECK1[nn].first;
   // LogStream << "\tBEFORE nInsertPointIntoProfilesIfNeededThenUpdate(): " << (nThisProfile == nProfileToTruncate ? "MAIN" : "COINCIDENT") << " to-truncate profile {" << nThisProfile << "} has " << m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetProfileSize() << " points ";
   // for (int nn = 0; nn < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetProfileSize(); nn++)
   // LogStream << "[" << m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nn)->dGetX() << ", " << m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nn)->dGetY() << "] ";
   // LogStream << "), and " << m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetNumLineSegments() << " line segments, co-incident profiles and their line segments are ";
   // for (int mm = 0; mm < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetNumLineSegments(); mm++)
   // {
   // vector<pair<int, int> > prVTmp = *m_VCoast[nCoast].pGetProfile(nThisProfile)->pprVGetPairedCoincidentProfilesForLineSegment(mm);
   // LogStream << "{ ";
   // for (int nn = 0; nn < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetNumCoincidentProfilesInLineSegment(mm); nn++)
   // LogStream << prVTmp[nn].first << "[" << prVTmp[nn].second << "] ";
   // LogStream << "} ";
   // }
   // LogStream << endl;
   //
   // for (int nPoint = 0; nPoint < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetProfileSize()-1; nPoint++)
   // {
   // CGeom2DPoint
   // Pt1 = *m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nPoint),
   // Pt2 = *m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nPoint+1);
   //
   // if (Pt1 == Pt2)
   // LogStream << m_ulIter << ": IDENTICAL POINTS before changes, in profile {" << nThisProfile << "} points = " << nPoint << " and " << nPoint+1 << endl;
   // }
   // }
   // // DEBUG CODE ******************************************************************

   int const nProfileToRetain = pProfileToRetain->nGetProfileID();

   // Get the index numbers of all coincident profiles for the 'main' to-retain profile for the line segment in which intersection occurs
   vector<pair<int, int>> prVCoincidentProfiles = *pProfileToRetain->pprVGetPairedCoincidentProfilesForLineSegment(nProfileToRetainIntersectLineSeg);
   int const nNumCoincident = static_cast<int>(prVCoincidentProfiles.size());
   vector<int> nLineSegAfterIntersect(nNumCoincident, -1); // The line segment after the point of intersection, for each co-incident profile

   // Do this for the main profile and all profiles which are co-incident for this line segment
   for (int nn = 0; nn < nNumCoincident; nn++)
   {
      int const nThisProfile = prVCoincidentProfiles[nn].first;  // The number of this profile
      int const nThisLineSeg = prVCoincidentProfiles[nn].second; // The line segment of this profile
      CGeomProfile* pThisProfile = m_VCoast[nCoast].pGetProfile(nThisProfile);

      // Is the intersection point already present in the to-retain profile?
      if (! bAlreadyPresent)
      {
         // It is not already present, so insert it and also update the associated multi-line
         if (! pThisProfile->bInsertIntersection(dIntersectX, dIntersectY, nThisLineSeg))
         {
            // Error
            LogStream << WARN << m_ulIter << ": cannot insert a line segment after the final line segment (" << nThisLineSeg << ") for " << (nThisProfile == nProfileToRetain ? "main" : "co-incident") << " profile (" << nThisProfile << "), abandoning" << endl;

            return RTN_ERR_CANNOT_INSERT_POINT;
         }

         // LogStream << "\tIntersection point NOT already in " << (nThisProfile == nProfileToRetain ? "main" : "co-incident") << " profile {" << nThisProfile << "}, inserted it as point " << nThisLineSeg+1 << endl;
      }

      // Get the line segment after intersection
      nLineSegAfterIntersect[nn] = nThisLineSeg + 1;
   }

   // for (int nn = 0; nn < nNumCoincident; nn++)
   // LogStream << "\tFor profile " << prVCoincidentProfiles[nn].first << " line segment [" << nLineSegAfterIntersect[nn] << "] is immediately after the intersection point" << endl;

   // Get the coincident profiles for the to-truncate profile, at the line segment where intersection occurs
   vector<pair<int, int>> prVToTruncateCoincidentProfiles = *pProfileToTruncate->pprVGetPairedCoincidentProfilesForLineSegment(nProfileToTruncateIntersectLineSeg);
   int const nNumToTruncateCoincident = static_cast<int>(prVToTruncateCoincidentProfiles.size());

   // Now add the number of the to-truncate profile, and all its coincident profiles, to all line segments which are seaward of the point of intersection. Do this for the main profile and all profiles which are co-incident for this line segment
   for (int nn = 0; nn < nNumCoincident; nn++)
   {
      int const nThisProfile = prVCoincidentProfiles[nn].first; // The number of this profile
      CGeomProfile* pThisProfile = m_VCoast[nCoast].pGetProfile(nThisProfile);

      // Get the number of line segments for this to-retain profile (will have just increased, if we just inserted a point)
      int const nNumLineSegs = pThisProfile->nGetNumLineSegments();

      // Do for all line segments seaward of the point of intersection
      for (int nLineSeg = nLineSegAfterIntersect[nn], nIncr = 0; nLineSeg < nNumLineSegs; nLineSeg++, nIncr++)
      {
         // // This can happen occasionally
         // if (nThisProfile == nProfileToTruncateIntersectLineSeg)
         // {
         // LogStream << "\t*** ERROR nThisProfile = " << nThisProfile << " nProfileToTruncateIntersectLineSeg = " << nProfileToTruncateIntersectLineSeg << ", ignoring" << endl;
         // pThisProfile->SetHitAnotherProfile(true);
         // continue;
         // }

         // Add the number of the to-truncate profile, and all its coincident profiles, to this line segment
         for (int m = 0; m < nNumToTruncateCoincident; m++)
         {
            int const nProfileToAdd = prVToTruncateCoincidentProfiles[m].first;
            int const nProfileToAddLineSeg = prVToTruncateCoincidentProfiles[m].second;

            // LogStream << "\tAdding " << (nProfileToAdd == nProfileToTruncateIntersectLineSeg ? "main" : "co-incident") << " truncate-profile " << nProfileToAdd << ", line segment [" << nProfileToAddLineSeg + nIncr << "] to line segment " << nLineSeg << " of " << (nThisProfile == nProfileToRetain ? "main" : "co-incident") << " to-retain profile " << nThisProfile << endl;

            pThisProfile->AddCoincidentProfileToExistingLineSegment(nLineSeg, nProfileToAdd, nProfileToAddLineSeg + nIncr);
         }
      }
   }

   // // DEBUG CODE ****************************************************************
   // Get the index numbers of all coincident profiles for the 'main' profile for the line segment in which intersection occurred
   // vector<pair<int, int> > prVCoincidentProfilesCHECK2 = *m_VCoast[nCoast].pGetProfile(nProfileToRetain)->pprVGetPairedCoincidentProfilesForLineSegment(nProfileToRetainIntersectLineSeg);
   // int nNumCoincidentCHECK2 = prVCoincidentProfilesCHECK2.size();
   // for (int nn = 0; nn < nNumCoincidentCHECK2; nn++)
   // {
   // int nThisProfile = prVCoincidentProfilesCHECK2[nn].first;
   // LogStream << "\tAFTER nInsertPointIntoProfilesIfNeededThenUpdate(): " << (nThisProfile == nProfileToRetain ? "MAIN" : "COINCIDENT") << " to-retain profile {" << nThisProfile << "} has " << m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetProfileSize() << " points ";
   // for (int nn = 0; nn < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetProfileSize(); nn++)
   // LogStream << "[" << m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nn)->dGetX() << ", " << m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nn)->dGetY() << "] ";
   // LogStream << "), and " << m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetNumLineSegments() << " line segments, co-incident profiles and their line segments are ";
   // for (int nLineSeg = 0; nLineSeg < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetNumLineSegments(); nLineSeg++)
   // {
   // vector<pair<int, int> > prVTmp = *m_VCoast[nCoast].pGetProfile(nThisProfile)->pprVGetPairedCoincidentProfilesForLineSegment(nLineSeg);
   // LogStream << "{ ";
   // for (int nn = 0; nn < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetNumCoincidentProfilesInLineSegment(nLineSeg); nn++)
   // LogStream << prVTmp[nn].first << "[" << prVTmp[nn].second << "] ";
   // LogStream << "} ";
   // }
   // LogStream << endl;
   //
   // for (int nPoint = 0; nPoint < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetProfileSize()-1; nPoint++)
   // {
   // CGeom2DPoint
   // Pt1 = *m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nPoint),
   // Pt2 = *m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nPoint+1);
   //
   // if (Pt1 == Pt2)
   // LogStream << m_ulIter << ": IDENTICAL POINTS before changes, in profile {" << nThisProfile << "} points = " << nPoint << " and " << nPoint+1 << endl;
   // }
   // }
   // // DEBUG CODE ******************************************************************

   return RTN_OK;
}

//===============================================================================================================================
//! Truncate a profile at the point of intersection, and do the same for all its co-incident profiles
//===============================================================================================================================
void CSimulation::TruncateProfileAndAppendNew(int const nCoast, CGeomProfile* pProfileToRetain, int const nMainProfileIntersectLineSeg, vector<CGeom2DPoint> const* pPtVProfileLastPart, vector<vector<pair<int, int>>> const* pprVLineSegLastPart)
{
   // // DEBUG CODE ****************************************************************
   // Get the index numbers of all coincident profiles for the 'main' profile for the line segment in which intersection occurred
   // vector<pair<int, int> > prVCoincidentProfilesCHECK1 = *m_VCoast[nCoast].pGetProfile(nMainProfile)->pprVGetPairedCoincidentProfilesForLineSegment(nMainProfileIntersectLineSeg);
   // int nNumCoincidentCHECK1 = prVCoincidentProfilesCHECK1.size();
   //
   // LogStream << "\tTruncating profile {" << nMainProfile << "}, intersection is at [" << dIntersectX << ", " << dIntersectY << "] in line segment " << nMainProfileIntersectLineSeg << endl;
   // for (int nn = 0; nn < nNumCoincidentCHECK1; nn++)
   // {
   // int nThisProfile = prVCoincidentProfilesCHECK1[nn].first;
   // LogStream << "\tBEFORE TruncateProfileAndAppendNew(): " << (nThisProfile == nMainProfile ? "MAIN" : "COINCIDENT") << " to-truncate profile {" << nThisProfile << "} has " << m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetProfileSize() << " points (";
   // for (int nn = 0; nn < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetProfileSize(); nn++)
   // LogStream << "[" << m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nn)->dGetX() << ", " << m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nn)->dGetY() << "] ";
   // LogStream << "), and " << m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetNumLineSegments() << " line segments, co-incident profiles are ";
   // for (int nLineSeg = 0; nLineSeg < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetNumLineSegments(); nLineSeg++)
   // {
   // vector<pair<int, int> > prVTmp = *m_VCoast[nCoast].pGetProfile(nThisProfile)->pprVGetPairedCoincidentProfilesForLineSegment(nLineSeg);
   // LogStream << "{ ";
   // for (int nn = 0; nn < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetNumCoincidentProfilesInLineSegment(nLineSeg); nn++)
   // LogStream << prVTmp[nn].first << "[" << prVTmp[nn].second << "] ";
   // LogStream << "} ";
   // }
   // LogStream << endl;
   //
   // for (int nPoint = 0; nPoint < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetProfileSize()-1; nPoint++)
   // {
   // CGeom2DPoint
   // Pt1 = *m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nPoint),
   // Pt2 = *m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nPoint+1);
   //
   // if (Pt1 == Pt2)
   // LogStream << m_ulIter << ": IDENTICAL POINTS before changes, in profile {" << nThisProfile << "} points = " << nPoint << " and " << nPoint+1 << endl;
   // }
   // }
   // LogStream << "\tPart-profile to append is ";
   // for (int mm = 0; mm < pPtVProfileLastPart->size(); mm++)
   // LogStream << "[" << pPtVProfileLastPart->at(mm).dGetX() << ", " << pPtVProfileLastPart->at(mm).dGetY() << "] ";
   // LogStream << endl;
   // LogStream << "\tPart line-segment to append is ";
   // for (int mm = 0; mm < pprVLineSegLastPart->size(); mm++)
   // {
   // vector<pair<int, int> > prVTmp = pprVLineSegLastPart->at(mm);
   // LogStream << "{ ";
   // for (int nn = 0; nn < prVTmp.size(); nn++)
   // LogStream << prVTmp[nn].first << "[" << prVTmp[nn].second << "] ";
   // LogStream << "} ";
   // }
   // LogStream << endl;
   // // DEBUG CODE ******************************************************************

   // Get the index numbers of all coincident profiles for the 'main' profile for the line segment in which intersection occurs
   vector<pair<int, int>> prVCoincidentProfiles = *pProfileToRetain->pprVGetPairedCoincidentProfilesForLineSegment(nMainProfileIntersectLineSeg);
   int const nNumCoincident = static_cast<int>(prVCoincidentProfiles.size());

   for (int nn = 0; nn < nNumCoincident; nn++)
   {
      // Do this for the main to-truncate profile, and do the same for all its co-incident profiles
      int const nThisProfile = prVCoincidentProfiles[nn].first;
      int const nThisProfileLineSeg = prVCoincidentProfiles[nn].second;
      CGeomProfile* pThisProfile = m_VCoast[nCoast].pGetProfile(nThisProfile);

      // if (nThisProfile == nMainProfile)
      // assert(nThisProfileLineSeg == nMainProfileIntersectLineSeg);

      // Truncate the profile
      // LogStream << "\tTruncating " << (nThisProfile == nMainProfile ? "MAIN" : "COINCIDENT") << " to-truncate profile {" << nThisProfile << "} at line segment " << nThisProfileLineSeg+1 << endl;
      pThisProfile->TruncateProfile(nThisProfileLineSeg + 1);

      // Reduce the number of line segments for this profile
      pThisProfile->TruncateLineSegments(nThisProfileLineSeg + 1);

      // Append the profile points from the last part of the retain-profile
      for (unsigned int mm = 0; mm < pPtVProfileLastPart->size(); mm++)
      {
         CGeom2DPoint const Pt = pPtVProfileLastPart->at(mm);
         pThisProfile->AppendPointInProfile(&Pt);
      }

      // Append the line segments, and their co-incident profile numbers, from the last part of the retain-profile
      for (unsigned int mm = 0; mm < pprVLineSegLastPart->size(); mm++)
      {
         vector<pair<int, int>> prVTmp = pprVLineSegLastPart->at(mm);

         pThisProfile->AppendLineSegment(&prVTmp);
      }

      // Fix the line seg numbers for this profile
      vector<int> nVProf;
      vector<int> nVProfsLineSeg;

      for (int nSeg = 0; nSeg < pThisProfile->nGetNumLineSegments(); nSeg++)
      {
         for (int nCoinc = 0; nCoinc < pThisProfile->nGetNumCoincidentProfilesInLineSegment(nSeg); nCoinc++)
         {
            int const nProf = pThisProfile->nGetProf(nSeg, nCoinc);
            int const nProfsLineSeg = pThisProfile->nGetProfsLineSeg(nSeg, nCoinc);

            auto it = find(nVProf.begin(), nVProf.end(), nProf);

            if (it == nVProf.end())
            {
               // Not found
               nVProf.push_back(nProf);
               nVProfsLineSeg.push_back(nProfsLineSeg);
            }

            else
            {
               // Found
               int const nPos = static_cast<int>(it - nVProf.begin());
               int nNewProfsLineSeg = nVProfsLineSeg[nPos];
               nNewProfsLineSeg++;

               nVProfsLineSeg[nPos] = nNewProfsLineSeg;
               pThisProfile->SetProfsLineSeg(nSeg, nCoinc, nNewProfsLineSeg);
            }
         }
      }

      // assert(pThisProfile->nGetProfileSize() > 1);
   }

   // // DEBUG CODE ****************************************************************
   // Get the index numbers of all coincident profiles for the 'main' to-truncate profile for the line segment in which intersection occurred
   // vector<pair<int, int> > prVToTruncateCoincidentProfilesCHECK2 = *m_VCoast[nCoast].pGetProfile(nMainProfile)->pprVGetPairedCoincidentProfilesForLineSegment(nMainProfileIntersectLineSeg);
   // int nNumToTruncateCoincidentCHECK2 = prVToTruncateCoincidentProfilesCHECK2.size();
   // for (int nn = 0; nn < nNumToTruncateCoincidentCHECK2; nn++)
   // {
   // int nThisProfile = prVToTruncateCoincidentProfilesCHECK2[nn].first;
   // LogStream << "\tAFTER TruncateProfileAndAppendNew(): " << (nThisProfile == nMainProfile ? "MAIN" : "COINCIDENT") << " to-truncate profile {" << nThisProfile << "} has " << m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetProfileSize() << " points (";
   // for (int nn = 0; nn < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetProfileSize(); nn++)
   // LogStream << "[" << m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nn)->dGetX() << ", " << m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nn)->dGetY() << "] ";
   // LogStream << "), and " << m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetNumLineSegments() << " line segments, co-incident profiles are ";
   // for (int mm = 0; mm < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetNumLineSegments(); mm++)
   // {
   // vector<pair<int, int> > prVTmp = *m_VCoast[nCoast].pGetProfile(nThisProfile)->pprVGetPairedCoincidentProfilesForLineSegment(mm);
   // LogStream << "{ ";
   // for (int nn = 0; nn < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetNumCoincidentProfilesInLineSegment(mm); nn++)
   // LogStream << prVTmp[nn].first << "[" << prVTmp[nn].second << "] ";
   // LogStream << "} ";
   // }
   // LogStream << endl;
   //
   // for (int nPoint = 0; nPoint < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetProfileSize()-1; nPoint++)
   // {
   // CGeom2DPoint
   // Pt1 = *m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nPoint),
   // Pt2 = *m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nPoint+1);
   //
   // if (Pt1 == Pt2)
   // LogStream << m_ulIter << ": IDENTICAL POINTS before changes, in profile {" << nThisProfile << "} points = " << nPoint << " and " << nPoint+1 << endl;
   // }
   // }
   // // DEBUG CODE ******************************************************************
}

