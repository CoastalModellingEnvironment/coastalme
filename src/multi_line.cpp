/*!
   \file multi_line.cpp
   \brief CGeomMultiLine routines
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

#include "multi_line.h"
#include "2d_point.h"
#include "line.h"

//! Constructor, no parameters
CGeomMultiLine::CGeomMultiLine(void)
{
}

//! Destructor
CGeomMultiLine::~CGeomMultiLine(void)
{
}

//! Returns a pointer to the points of the CGeomLine
vector<CGeom2DPoint>& CGeomMultiLine::pGetPoints(void)
{
   return CGeomLine::m_VPtPoints;
}

// //! Replaces the points of the CGeomLine
// void CGeomMultiLine::SetPoints(vector<CGeom2DPoint> const& pVPts)
// {
//    CGeomLine::m_VPtPoints = pVPts;
// }

//! Appends a new empty line segment
void CGeomMultiLine::AppendLineSegment(void)
{
   m_prVVLineSegment.push_back(vector<pair<int, int>>());
}

//! Appends a new line segment which is supplied as a parameter
void CGeomMultiLine::AppendLineSegment(vector<pair<int, int>>* pprVIn)
{
   m_prVVLineSegment.push_back(*pprVIn);
}

//! Appends a line segment which then inherits from the preceding line segments
// void CGeomMultiLine::AppendLineSegmentAndInherit(void)
// {
// vector<pair<int, int> > prVNewLineSeg;
// m_prVVLineSegment.push_back(prVNewLineSeg);
//
//    // Must inherit any profile numbers stored in earlier (i.e. coastward) line segments
// int nSize = m_prVVLineSegment.size();
// if (nSize > 1)
// {
// for (int n = 0; n < m_prVVLineSegment[nSize-2].size(); n++)
// {
// int
// nPrevProfile = m_prVVLineSegment[nSize-2][n].first,
// nPrevLineSeg = m_prVVLineSegment[nSize-2][n].second + 1;
//
// m_prVVLineSegment[nSize-1].push_back(make_pair(nPrevProfile, nPrevLineSeg));
// }
// }
// }

//! Returns the CGeomMultiLine object's number of line segments
int CGeomMultiLine::nGetNumLineSegments(void) const
{
   return static_cast<int>(m_prVVLineSegment.size());
}

//! Cuts short the number of line segments
void CGeomMultiLine::TruncateLineSegments(int const nSize)
{
   m_prVVLineSegment.resize(nSize);
}

//! Inserts a line segment, inheriting coincident pairs from preceding line segments
void CGeomMultiLine::InsertLineSegmentWithInheritance(int const nSegment)
{
   assert(nSegment < static_cast<int>(m_prVVLineSegment.size()));

   // The new vector of pairs is identical to the existing vector of pairs i.e. we inherit profile/line seg details from the previous line seg
   vector<pair<int, int>> prVPrev = m_prVVLineSegment[nSegment];

   // Store the profile numbers that are in this existing vector of pairs, these are the profiles that will be affected by this insertion
   vector<int> nVProfsAffected;
   nVProfsAffected.reserve(prVPrev.size());

   for (unsigned int i = 0; i < prVPrev.size(); i++)
      nVProfsAffected.push_back(prVPrev[i].first);

   vector<vector<pair<int, int>>>::iterator it;
   it = m_prVVLineSegment.begin();

   m_prVVLineSegment.insert(it + nSegment + 1, prVPrev);

   // Must now increment the profile's own line seg numbers, but only for those profile numbers which were affected by the insertion. Do this for the new line seg and every line seg after that
   for (unsigned int m = nSegment + 1; m < m_prVVLineSegment.size(); m++)
   {
      for (unsigned int n = 0; n < m_prVVLineSegment[m].size(); n++)
      {
         for (unsigned int i = 0; i < nVProfsAffected.size(); i++)
         {
            if (m_prVVLineSegment[m][n].first == nVProfsAffected[i])
               m_prVVLineSegment[m][n].second++;
         }
      }
   }
}

//! Returns a vector of the line segments which follow the specified line segment
vector<vector<pair<int, int>>> CGeomMultiLine::prVVGetAllLineSegmentsAfter(int const nSegment)
{
   vector<vector<pair<int, int>>> prVTmp;

   for (unsigned int n = nSegment; n < m_prVVLineSegment.size(); n++)
      prVTmp.push_back(m_prVVLineSegment[n]);

   return prVTmp;
}

// //! Removes a line segment
// void CGeomMultiLine::RemoveLineSegment(int const nSegment)
// {
//    m_prVVLineSegment.erase(m_prVVLineSegment.begin() + nSegment);
// }

//! Appends a line segment pair (profile and line segment) to the CGeomMultiLine object's final line segment, if the profiule isn't already present in the line segment
void CGeomMultiLine::AppendPairToFinalLineSegment(pair<int, int> const prIn)
{
   m_prVVLineSegment.back().push_back(prIn);
}

//! Adds a coincident pair (profile and line segment) to a pre-existing line segment of the CGeomMultiLine object, if the profile isn't already present in the line segment
void CGeomMultiLine::AddCoincidentPairToExistingLineSegmentIfNotAlready(int const nSegment, int const nProfile, int const nLineSeg)
{
   for (size_t n = 0; n < m_prVVLineSegment[nSegment].size(); n++)
   {
      if (m_prVVLineSegment[nSegment][n].first == nProfile)
         return;
   }

   m_prVVLineSegment[nSegment].push_back(make_pair(nProfile, nLineSeg));
}

//! Returns a vector of pairs (a line segment)
vector<pair<int, int>>* CGeomMultiLine::pprVGetCoincidentPairsForLineSegment(int const nSegment)
{
   // TODO 055 No check to see if nSegment < size()
   return &m_prVVLineSegment[nSegment];
}

//! Returns the numbers of coincident pairs for a given line segment, or -1 if this line segment does not exist, or -2 if there are no coincident pairs for this line segment
int CGeomMultiLine::pprVGetFirstFromCoincidentPairForLineSegment(int const nSegment, int const nCoinc) const
{
   // Safety check
   if ((nSegment < 0) || (nSegment >= static_cast<int>(m_prVVLineSegment.size())))
      return -1;

   // Safety check
   if ((nCoinc < 0) || (nCoinc >= static_cast<int>(m_prVVLineSegment[nSegment].size())))
      return -2;

   return m_prVVLineSegment[nSegment][nCoinc].first;
}

//! Returns the count of coincident profiles in a specified line segment, or -1 if the line segment does not exist
int CGeomMultiLine::nGetNumCoincidentPairsForLineSegment(int const nSegment)
{
   // Safety check
   if (nSegment > static_cast<int>(m_prVVLineSegment.size()) - 1)
      return -1;

   return static_cast<int>(m_prVVLineSegment[nSegment].size());
}

//! Returns true if the given profile number is amongst the coincident profiles of the CGeomMultiLine object's final line segment
bool CGeomMultiLine::bFindFirstFromCoincidentPairsInLastLineSegment(int const nProfile)
{
   long unsigned int const nLineSegSize = m_prVVLineSegment.size();

   // Note no check to ensure that nLineSegSize < 0
   long unsigned int const nCoincidentSize = m_prVVLineSegment[nLineSegSize - 1].size();

   for (unsigned int i = 0; i < nCoincidentSize; i++)
   {
      if (m_prVVLineSegment[nLineSegSize - 1][i].first == nProfile)
         return true;
   }

   return false;
}

//! Returns true if the given pair-first is one of the coincident pair-firsts of the specified line segment
// bool CGeomMultiLine::bFindFirstInCoincidentPairsOfLineSegment(int const nPairFirst, int const nSegment)
// {
//    // Note no check to see if nSegment < m_prVVLineSegment.size()
// int nCoincidentSize = m_prVVLineSegment[nSegment].size();
//
// for (int i = 0; i < nCoincidentSize; i++)
// if (m_prVVLineSegment[nSegment][i].first == nPairFirst)
// return true;
//
// return false;
// }

//! Returns true if the given pair-first is a coincident pair of any line segment of the CGeomMultiLine object
bool CGeomMultiLine::bFindFirstInCoincidentPairs(int const nPairFirst)
{
   int const nSegSize = static_cast<int>(m_prVVLineSegment.size());

   if (nSegSize == 0)
      return false;

   for (int i = nSegSize - 1; i >= 0; i--)
   {
      for (unsigned int j = 0; j < m_prVVLineSegment[i].size(); j++)
      {
         if (m_prVVLineSegment[i][j].first == nPairFirst)
            return true;
      }
   }

   return false;
}

//! Searches for the lowest-numbered line segment for which the given pair-first is coincident, or -1 if not coincident. If coincident, also finds the line segment of the other profile
void CGeomMultiLine::SearchForLowestNumberedCoincidentLineSegments(int const nPairFirst, int& nThisLineSegment, int& nOtherLineSegment)
{
   nThisLineSegment = -1;
   nOtherLineSegment = -1;

   long unsigned int const nSegSize = m_prVVLineSegment.size();

   if (nSegSize == 0)
      return;

   for (unsigned int i = 0; i < nSegSize; i++)
   {
      for (unsigned int j = 0; j < m_prVVLineSegment[i].size(); j++)
      {
         if (m_prVVLineSegment[i][j].first == nPairFirst)
         {
            nThisLineSegment = i;
            nOtherLineSegment = m_prVVLineSegment[i][j].second;

            return;
         }
      }
   }
}

//! Returns the pair-first, given a line segment number and the index of a co-incident pair-first
int CGeomMultiLine::nGetPairFirst(int const nSegment, int const nCoinc) const
{
   return m_prVVLineSegment[nSegment][nCoinc].first;
}

//! Returns the pair-second, given a line segment number and the index of a co-incident pair-first
int CGeomMultiLine::nGetPairSecond(int const nSegment, int const nCoinc) const
{
   return m_prVVLineSegment[nSegment][nCoinc].second;
}

//! Sets the pair-second, given a line segment number and the index of a co-incident pair-first
void CGeomMultiLine::SetPairSecond(int const nSegment, int const nCoinc, int const nLineSeg)
{
   // Note no check to see if nSegment < m_prVVLineSegment.size() or to see if nCoinc < m_prVVLineSegment[nSegment].size()
   m_prVVLineSegment[nSegment][nCoinc].second = nLineSeg;
}

// //! Returns the number of the last line segment which includes the given pair-first as a co-incident
// int CGeomMultiLine::nFindLastSegForCoincPairFirst(int const nPairFirst) const
// {
// int nSeg = -1;
// for (int i = static_cast<int>(m_prVVLineSegment.size()-1); i >= 0; i--)
// {
// for (unsigned int j = 0; j < m_prVVLineSegment[i].size(); j++)
// {
// if (m_prVVLineSegment[i][j].first == nPairFirst)
// nSeg = i;
// }
// }
//
// return nSeg;
// }
