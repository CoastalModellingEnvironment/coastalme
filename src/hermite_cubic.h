/*!
   \brief Definitions of some routines from the hermite_cubic library
   \details TODO 001 This is a more detailed description of the hermite_cubic routines.
   \author John Burkardt
   \author Modified by David Favis-Mortlock and Andres Payo
   \date 2026
   \copyright GNU Lesser General Public License
   \file hermite_cubic.h
   \brief Contains definitions of hermite-cubic routines
*/

#ifndef HERMITE_H
#define HERMITE_H
/* ===============================================================================================================================
   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
===============================================================================================================================*/
void r8vec_bracket3(int const, double const*, double const, int*);
void hermite_cubic_value(double const, double const, double const, double const, double const, double const, int const, double const*, double* const, double* const, double* const, double* const);
void hermite_cubic_spline_value(int const, double const*, double const*, double const*, int const, double const*, double*, double*, double*, double*);
#endif // HERMITE_H
