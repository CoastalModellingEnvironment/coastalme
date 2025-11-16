/*!

   \file yaml_parser.h
   \brief Simple YAML parser for CoastalME configuration files
   \details A lightweight YAML parser using only standard C++ library, designed specifically for CoastalME configuration needs
   \author Wilf Chun
   \author David Favis-Mortlock
   \author Andres Payo
   \date 2025
   \copyright GNU General Public License

*/

/* ==============================================================================================================================

   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public  License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

==============================================================================================================================*/
#ifndef YAML_PARSER_H
#define YAML_PARSER_H

#include <string>
#include <map>
#include <vector>
#include <fstream>

using std::ifstream;
using std::map;
using std::string;
using std::vector;

//! Simple YAML node class to represent parsed values
class CYamlNode
{
   private:
   string m_strValue;
   map<string, CYamlNode> m_mapChildren;
   vector<CYamlNode> m_vecChildren;
   bool m_bIsSequence;

   public:
   CYamlNode();
   ~CYamlNode();

   void SetValue(string const& strValue);
   void AddChild(string const& strKey, CYamlNode const& node);
   void AddSequenceItem(CYamlNode const& node);

   string const* pstrGetValue() const;
   bool HasChild(string const& strKey) const;
   CYamlNode GetChild(string const& strKey) const;
   vector<CYamlNode> GetSequence() const;
   bool IsSequence() const;
   int GetSequenceSize() const;

   // Convenience methods for common types
   int GetIntValue(int nDefault = 0) const;
   unsigned long GetULongValue(unsigned long nDefault = 0) const;
   double dGetDoubleValue(double dDefault = 0.0) const;
   bool bGetBoolValue(bool bDefault = false) const;
   vector<string> VstrGetStringSequence() const;
};

//! Simple YAML parser class
class CYamlParser
{
   private:
   CYamlNode m_RootNode;
   string m_strFileName;
   int m_nCurrentLine;
   string m_strError;

   // Helper methods
   int nGetIndentLevel(string const& strLine) const;
   string strTrimLeft(string const& strLine) const;
   string strTrimRight(string const& strLine) const;
   string strTrim(string const& strLine) const;
   bool bIsComment(string const& strLine) const;
   bool bIsEmpty(string const& strLine) const;
   bool bParseLine(string const& strLine, string& strKey, string& strValue, bool& bIsSequence) const;
   string strRemoveQuotes(string const& strValue) const;
   CYamlNode ParseSection(ifstream& fileStream, int nBaseIndent);

   public:
   CYamlParser();
   ~CYamlParser();

   bool bParseFile(string const& strFileName);
   CYamlNode GetRoot() const;
   string const* pstrGetError() const;
   bool bHasError() const;
};

#endif      // YAML_PARSER_H
