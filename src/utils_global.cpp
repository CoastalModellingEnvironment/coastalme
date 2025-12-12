/*!
   \file utils_global.cpp
   \brief Globally-available utility routines
   \details TODO 001 A more detailed description of these routines.
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
#include <cmath>

#include <cstdio>
using std::size_t;

#include <sstream>
using std::stringstream;
using std::istringstream;

#include <ios>
using std::fixed;
using std::right;

#include "cme.h"

//===============================================================================================================================
//! Correctly rounds doubles
//===============================================================================================================================
double dRound(double const d)
{
   // Rounds positive or negative doubles correctly
   return ((d < 0.0) ? ceil(d - 0.5) : floor(d + 0.5));
}

//===============================================================================================================================
//! Correctly rounds doubles, returns an int
//===============================================================================================================================
int nRound(double const d)
{
   // Rounds positive or negative doubles correctly
   return static_cast<int>((d < 0.0) ? ceil(d - 0.5) : floor(d + 0.5));
}

// bool bIsWhole(double d)
// {
//    // From http://answers.yahoo.com/question/index?qid=20110320132617AAMdb7u
// return (static_cast<int>(d) == d);
// }

//===============================================================================================================================
//! Checks to see if a string can be read as a valid double number. Does not find trailing (i.e.post-number) rubbish, but then neither does strtod(). From https://stackoverflow.com/questions/392981/how-can-i-convert-string-to-double-in-c
//===============================================================================================================================
bool bIsStringValidDouble(string& str)
{
   istringstream iStr(str);
   double dDummy;

   if (! (iStr >> dDummy))
      return false;

   return true;
}

//===============================================================================================================================
//! Checks to see if a string can be read as a valid integer, from https://stackoverflow.com/questions/2844817/how-do-i-check-if-a-c-string-is-an-int
//===============================================================================================================================
bool bIsStringValidInt(string& str)
{
   // Trim leading whitespace
   size_t const nPos = str.find_first_not_of(" \t");

   if (nPos != string::npos)
      str = str.substr(nPos);

   // If the first character is the sign, remove it
   if ((str[0] == '-') || (str[0] == '+'))
      str.erase(0, 1);

   // Now check that the string contains only numbers
   return (str.find_first_not_of("0123456789") == string::npos);
}

//===============================================================================================================================
//! Operator that inserts a given fill character, to a given width, into an output stream. From http://stackoverflow.com/questions/2839592/equivalent-of-02d-with-stdstringstream
//===============================================================================================================================
ostream& operator<<(ostream& ostr, const FillToWidth& args)
{
   ostr.fill(args.chFill);
   ostr.width(args.nWidth);

   return ostr;
}

// //===============================================================================================================================
// //! Converts double to string with specified number of places after the decimal. From https://stackoverflow.com/questions/14765155/how-can-i-easily-format-my-data-table-in-c
// //===============================================================================================================================
// string strDbl(double const dX, int const nDigits)
// {
//    stringstream ss;
//    ss << fixed;
//    ss.precision(nDigits);     // Set the number of places after decimal
//    ss << dX;
//    return ss.str();
// }

//===============================================================================================================================
//! Converts double to string with specified number of decimal places, within a field of given width, pads with blank spaces to enforce right alignment. Modified from https://stackoverflow.com/questions/14765155/how-can-i-easily-format-my-data-table-in-c
//===============================================================================================================================
string strDblRight(double const dX, int const nDigits, int const nWidth, bool const bShowDash)
{
   stringstream ss;
   ss << fixed << right;
   ss.fill(' ');
   ss.width(nWidth - 1);

   if (bFPIsEqual(dX, 0.0, TOLERANCE))
   {
      if (bShowDash)
         ss << "-";

      else
         ss << SPACE;
   }

   else
   {
      ss.precision(nDigits);  // Set number of places after decimal
      ss << dX;
   }

   ss << " ";                 // Add a final space
   return ss.str();
}

//===============================================================================================================================
//! Converts int to string within a field of given width, pads with blank spaces to enforce alignment.. From https://stackoverflow.com/questions/14765155/how-can-i-easily-format-my-data-table-in-c
//===============================================================================================================================
string strIntRight(int const nX, int const nWidth)
{
   stringstream ss;
   ss << fixed << right;
   ss.fill(' ');              // Fill space around displayed number
   ss.width(nWidth - 1);      // Set width around displayed number
   ss << nX;
   ss << " "; // Add a final space
   return ss.str();
}

//===============================================================================================================================
//! Centre-aligns char array within a field of given width, pads with blank spaces to enforce alignment. From https://stackoverflow.com/questions/14765155/how-can-i-easily-format-my-data-table-in-c
//===============================================================================================================================
string strCentre(const char* pchIn, int const nWidth)
{
   string const strIn(pchIn);
   stringstream ss, spaces;
   int const nPadding = nWidth - static_cast<int>(strIn.size());

   for (int i = 0; i < nPadding / 2; ++i)
      spaces << " ";

   ss << spaces.str() << strIn << spaces.str();

   if (nPadding > 0 && nPadding % 2 != 0) // If odd number, add one space
      ss << " ";

   return ss.str();
}

//===============================================================================================================================
//! Centre-aligns string within a field of given width, pads with blank spaces to enforce alignment. From https://stackoverflow.com/questions/14765155/how-can-i-easily-format-my-data-table-in-c
//===============================================================================================================================
string strCentre(const string& strIn, int const nWidth)
{
   stringstream ss, spaces;
   int const nPadding = nWidth - static_cast<int>(strIn.size());

   for (int i = 0; i < nPadding / 2; ++i)
      spaces << " ";

   ss << spaces.str() << strIn << spaces.str();

   if (nPadding > 0 && nPadding % 2 != 0) // If odd number, add one space
      ss << " ";

   return ss.str();
}

//===============================================================================================================================
//! Right-aligns string within a field of given width, pads with blank spaces to enforce alignment. From https://stackoverflow.com/questions/14765155/how-can-i-easily-format-my-data-table-in-c
//===============================================================================================================================
string strRight(const string& strIn, int const nWidth)
{
   stringstream ss, spaces;
   int const nPadding = nWidth - static_cast<int>(strIn.size()) - 1;

   for (int i = 0; i < nPadding; ++i)
      spaces << " ";

   ss << spaces.str() << strIn;
   ss << " ";
   return ss.str();
}

//===============================================================================================================================
//! Right-aligns char array within a field of given width, pads with blank spaces to enforce alignment. From https://stackoverflow.com/questions/14765155/how-can-i-easily-format-my-data-table-in-c
//===============================================================================================================================
string strRight(const char* pchIn, int const nWidth)
{
   string const strIn(pchIn);
   stringstream ss, spaces;
   int const nPadding = nWidth - static_cast<int>(strIn.size()) - 1;

   for (int i = 0; i < nPadding; ++i)
      spaces << " ";

   ss << spaces.str() << strIn;
   ss << " ";
   return ss.str();
}

//===============================================================================================================================
//! Left-aligns string within a field of given width, pads with blank spaces to enforce alignment. From https://stackoverflow.com/questions/14765155/how-can-i-easily-format-my-data-table-in-c
//===============================================================================================================================
string strLeft(const string& strIn, int const nWidth)
{
   stringstream ss, spaces;
   int const nPadding = nWidth - static_cast<int>(strIn.size());

   for (int i = 0; i < nPadding; ++i)
      spaces << " ";

   ss << strIn << spaces.str();
   return ss.str();
}

//===============================================================================================================================
//! Left-aligns char array within a field of given width, pads with blank spaces to enforce alignment. From https://stackoverflow.com/questions/14765155/how-can-i-easily-format-my-data-table-in-c
//===============================================================================================================================
string strLeft(const char* pchIn, int const nWidth)
{
   string const strIn(pchIn);
   stringstream ss, spaces;
   int const nPadding = nWidth - static_cast<int>(strIn.size());

   for (int i = 0; i < nPadding; ++i)
      spaces << " ";

   ss << strIn << spaces.str();
   return ss.str();
}

//===============================================================================================================================
//! Calculates a percentage from two numbers then, if the result is non-zero, right-aligns the result as a string within a field of given width, pads with blank spaces to enforce alignment. Modified from https://stackoverflow.com/questions/14765155/how-can-i-easily-format-my-data-table-in-c
//===============================================================================================================================
string strRightPerCent(double const d1, double const d2, int const nWidth, int const nDigits, bool const bShowDash)
{
   stringstream ss;
   ss << fixed << right;

   // Are either of the inputs zero?
   if ((bFPIsEqual(d1, 0.0, TOLERANCE)) || (bFPIsEqual(d2, 0.0, TOLERANCE)))
   {
      ss.fill(' ');
      ss.width(nWidth - 1);

      if (bShowDash)
         ss << "-";

      else
         ss << SPACE;
   }
   else
   {
      // Non-zero, so calculate the percentage
      double const dResult = 100 * d1 / d2;

      stringstream ssResult;
      ssResult << fixed << right;
      ssResult.precision(nDigits);
      ssResult << "(" << dResult << "%)";

      long int const nResultWidth = ssResult.str().size();

      for (int i = 0; i < (nWidth - nResultWidth - 1); i++)
         ss << SPACE;

      ss << ssResult.str();
   }

   // Append a final space
   ss << " ";

   return ss.str();
}
