/*!
   \class CRWIntervention
   \brief Real-world class used to represent the 'intervention' category of coastal landform objects
   \details TODO 001 This is a more detailed description of the CRWIntervention class.
   \author David Favis-Mortlock
   \author Andres Payo
   \author Wilf Chun
   \date 2025
   \copyright GNU General Public License
   \file intervention.h
   \brief Contains CRWIntervention definitions
*/

#ifndef INTERVENTION_H
#define INTERVENTION_H
/* ===============================================================================================================================
   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
===============================================================================================================================*/
#include "coast.h"
#include "coast_landform.h"

class CRWIntervention : public CACoastLandform
{
 private:
 public:
   CRWIntervention(CRWCoast*, int const, int const, int const);
   ~CRWIntervention(void) override;

   void Display(void) override;
};
#endif // INTERVENTION_H
