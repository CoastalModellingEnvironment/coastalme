/*!
   \file cme.cpp
   \brief The start-up routine for CoastalME
   \details TODO 001 A more detailed description of this routine
   \author David Favis-Mortlock
   \author Andres Payo
   \author Wilf Chun
   \date 2025
   \copyright GNU General Public License
*/

/* ===============================================================================================================================
   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
===============================================================================================================================*/
#include <clocale>

#include "simulation.h"

// #include <fenv.h>    // Include this to check for first appearance in NaN when debugging (comment out, otherwise)

//===============================================================================================================================
//! CoastalME's main function
//===============================================================================================================================
int main(int argc, char const* argv[])
{
   // This is to check for first appearance of NaN when debugging (comment out, otherwise)
   // #ifdef __APPLE__
   // #else
   //    feenableexcept(FE_INVALID | FE_OVERFLOW);
   // #endif

   // Enable the use of UTF-8 symbols in CoastalME output
   setlocale(LC_ALL, "en_GB.UTF-8");

   // Create a CSimulation object
   CSimulation* pSimulation = new CSimulation;

   // Run the simulation and then check how it ends
   int const nRtn = pSimulation->nDoSimulation(argc, argv);
   pSimulation->DoSimulationEnd(nRtn);

   // Get rid of the CSimulation object and close files
   delete pSimulation;

   // Go back to the OS
   return nRtn;
}
