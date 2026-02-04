/*!
   \file yaml_parser.cpp
   \brief Simple YAML parser implementation for CoastalME
   \details A lightweight YAML parser using only standard C++ library
   \author Wilf Chun
   \author David Favis-Mortlock
   \author Andres Payo
   \date 2026
   \copyright GNU General Public License
*/

/* ==============================================================================================================================
   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public  License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
==============================================================================================================================*/
#include <iostream>
using std::cerr;
using std::endl;

#include <string>
using std::getline;
using std::stod;
using std::stoi;
using std::to_string;

#include <sstream>
using std::stringstream;

#include <algorithm>

#include <cctype>
using std::isspace;

#include <exception>
using std::exception;

#include "yaml_parser.h"
#include "simulation.h"

//===============================================================================================================================
// CYamlNode implementation
//===============================================================================================================================
CYamlNode::CYamlNode() : m_bIsSequence(false)
{
}

CYamlNode::~CYamlNode()
{
}

void CYamlNode::SetValue(string const* pStrValue)
{
   m_strValue = *pStrValue;
}

void CYamlNode::AddChild(string const* pStrKey, CYamlNode const* pNode)
{
   m_mapChildren[*pStrKey] = *pNode;
}

void CYamlNode::AddSequenceItem(CYamlNode const* pNode)
{
   m_VYamlChildren.push_back(*pNode);
   m_bIsSequence = true;
}

string* CYamlNode::pstrGetValue() const
{
   return const_cast<string*>(&m_strValue);
}

bool CYamlNode::bHasChild(const char chKey[]) const
{
   return m_mapChildren.find(chKey) != m_mapChildren.end();
}

CYamlNode CYamlNode::GetChild(string const& strKey) const
{
   auto it = m_mapChildren.find(strKey);
   if (it != m_mapChildren.end())
      return it->second;
   return CYamlNode();      // Return empty node if not found
}

vector<CYamlNode> CYamlNode::GetSequence() const
{
   return m_VYamlChildren;
}

bool CYamlNode::IsSequence() const
{
   return m_bIsSequence;
}

int CYamlNode::nGetSequenceSize() const
{
   return static_cast<int>(m_VYamlChildren.size());
}

int CYamlNode::nGetIntValue(int nDefault) const
{
   try
   {
      if (! m_strValue.empty())
         return stoi(m_strValue);
   }
   catch (...)
   {
      // Return default on conversion error
   }
   return nDefault;
}

unsigned long CYamlNode::GetULongValue(unsigned long nDefault) const
{
   try
   {
      if (! m_strValue.empty())
         return stoul(m_strValue);
   }
   catch (...)
   {
      // Return default on conversion error
   }
   return nDefault;
}

double CYamlNode::dGetDoubleValue(double dDefault) const
{
   try
   {
      if (! m_strValue.empty())
         return stod(m_strValue);
   }
   catch (...)
   {
      // Return default on conversion error
   }
   return dDefault;
}

bool CYamlNode::bGetBoolValue(bool bDefault) const
{
   if (m_strValue.empty())
      return bDefault;

   string strLower = CSimulation::strToLower(&m_strValue);

   if (strLower == "true" || strLower == "yes" || strLower == "y" || strLower == "1")
      return true;
   else if (strLower == "false" || strLower == "no" || strLower == "n" || strLower == "0")
      return false;

   return bDefault;
}

vector<string> const CYamlNode::VstrGetStringSequence() const
{
   vector<string> VstrResult;
   for (auto const& node : m_VYamlChildren)
   {
      VstrResult.push_back(*node.pstrGetValue());
   }
   return VstrResult;
}

//===============================================================================================================================
// CYamlParser implementation
//===============================================================================================================================
CYamlParser::CYamlParser() : m_nCurrentLine(0)
{
}

CYamlParser::~CYamlParser()
{
}

bool CYamlParser::bParseFile(string const& strFileName)
{
   m_strFileName = strFileName;
   m_nCurrentLine = 0;
   m_strError.clear();

   ifstream fileStream(strFileName);
   if (! fileStream.is_open())
   {
      m_strError = "Cannot open file: " + strFileName;
      return false;
   }

   try
   {
      m_RootNode = ParseSection(fileStream, -1);
   }
   catch (exception const& e)
   {
      m_strError = "Parse error at line " + to_string(m_nCurrentLine) + ": " + e.what();
      return false;
   }

   fileStream.close();
   return true;
}

CYamlNode CYamlParser::GetRoot() const
{
   return m_RootNode;
}

string const* CYamlParser::pstrGetError() const
{
   return &m_strError;
}

bool CYamlParser::bHasError() const
{
   return ! m_strError.empty();
}

int CYamlParser::nGetIndentLevel(string const& strLine) const
{
   int nIndent = 0;
   for (char c : strLine)
   {
      if (c == ' ')
         nIndent++;
      else if (c == '\t')
         nIndent += 4;      // Treat tab as 4 spaces
      else
         break;
   }
   return nIndent;
}

string CYamlParser::strTrimLeft(string const& strLine) const
{
   auto it = strLine.begin();
   while (it != strLine.end() && isspace(*it))
      it++;
   return string(it, strLine.end());
}

string CYamlParser::strTrimRight(string const& strLine) const
{
   auto it = strLine.rbegin();
   while (it != strLine.rend() && isspace(*it))
      it++;
   return string(strLine.begin(), it.base());
}

string CYamlParser::strTrim(string const& strLine) const
{
   return strTrimLeft(strTrimRight(strLine));
}

bool CYamlParser::bIsComment(string const& strLine) const
{
   string strTrimmed = strTrimLeft(strLine);
   return strTrimmed.empty() || strTrimmed[0] == '#';
}

bool CYamlParser::bIsEmpty(string const& strLine) const
{
   return strTrim(strLine).empty();
}

bool CYamlParser::bParseLine(string const& strLine, string& strKey, string& strValue, bool& bIsSequence) const
{
   string strTrimmed = strTrimLeft(strLine);
   bIsSequence = false;

   // Check for sequence item
   if (strTrimmed.length() > 0 && strTrimmed[0] == '-')
   {
      bIsSequence = true;
      strKey.clear();
      strValue = strTrim(strTrimmed.substr(1));
      strValue = strRemoveQuotes(strValue);
      return true;
   }

   // Look for key-value pair
   size_t nColonPos = strTrimmed.find(':');
   if (nColonPos == string::npos)
      return false;

   strKey = strTrim(strTrimmed.substr(0, nColonPos));
   if (nColonPos + 1 < strTrimmed.length())
      strValue = strTrim(strTrimmed.substr(nColonPos + 1));
   else
      strValue.clear();

   // Remove quotes from string values
   strValue = strRemoveQuotes(strValue);

   return true;
}

string CYamlParser::strRemoveQuotes(string const& strValue) const
{
   string result = strValue;

   // Remove surrounding double quotes
   if (result.length() >= 2 && result.front() == '"' && result.back() == '"')
   {
      result = result.substr(1, result.length() - 2);
   }
   // Remove surrounding single quotes
   else if (result.length() >= 2 && result.front() == '\'' && result.back() == '\'')
   {
      result = result.substr(1, result.length() - 2);
   }

   return result;
}

CYamlNode CYamlParser::ParseSection(ifstream& fileStream, int nBaseIndent)
{
   CYamlNode currentNode;
   string strLine;

   while (getline(fileStream, strLine))
   {
      m_nCurrentLine++;

      // Skip comments and empty lines
      if (bIsComment(strLine) || bIsEmpty(strLine))
         continue;

      int nIndent = nGetIndentLevel(strLine);

      // If we've gone back to a lower indentation level, we're done with this section
      if (nBaseIndent >= 0 && nIndent <= nBaseIndent)
      {
         // Put the line back for the parent to process
         fileStream.seekg(static_cast<int>(fileStream.tellg()) - static_cast<int>(strLine.length()) - 1);
         m_nCurrentLine--;
         break;
      }

      string strKey;
      string strValue;
      bool bIsSequence;

      if (bParseLine(strLine, strKey, strValue, bIsSequence))
      {
         if (bIsSequence)
         {
            // Handle sequence item
            CYamlNode itemNode;
            if (! strValue.empty())
            {
               itemNode.SetValue(&strValue);
            }
            else
            {
               // Multi-line sequence item
               itemNode = ParseSection(fileStream, nIndent);
            }
            currentNode.AddSequenceItem(&itemNode);
         }
         else
         {
            // Handle key-value pair
            CYamlNode childNode;
            if (! strValue.empty())
            {
               childNode.SetValue(&strValue);
            }
            else
            {
               // Multi-line value
               childNode = ParseSection(fileStream, nIndent);
            }
            currentNode.AddChild(&strKey, &childNode);
         }
      }
   }

   return currentNode;
}
