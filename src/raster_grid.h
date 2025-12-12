/*!
   \class CGeomRasterGrid
   \brief Geometry cass used to represent the raster grid of cell objects
   \details TODO 001 This is a more detailed description of the CGeomRasterGrid class.
   \author David Favis-Mortlock
   \author Andres Payo
   \author Wilf Chun
   \date 2025
   \copyright GNU General Public License
   \file raster_grid.h
   \brief Contains CGeomRasterGrid definitions
*/

#ifndef RASTERGRID_H
#define RASTERGRID_H
/* ===============================================================================================================================
   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
===============================================================================================================================*/

class CGeomCell;   // Forward declaration
class CSimulation; // Ditto

class CGeomRasterGrid
{
   //! The CSimulation class is a friend of the CGeomRasterGrid class
   friend class CSimulation;

   //! The CGeomProfile class is a friend of the CGeomRasterGrid class
   friend class CGeomProfile;

 private:
   //! The d50 of fine-sized  sediment
   double m_dD50Fine;

   //! The d50 of sand-sized sediment
   double m_dD50Sand;

   //! The d50 of coarse-sized sediment
   double m_dD50Coarse;

   //! A pointer to the CSimulation object
   CSimulation* m_pSim;

   //! The 2D array of m_Cell objects. A c-style 2D array seems to be faster than using 2D STL vectors
   CGeomCell** m_Cell;

 protected:
 public:
   explicit CGeomRasterGrid(CSimulation*);
   ~CGeomRasterGrid(void);

   CSimulation* pGetSim(void);
   // CGeomCell* pGetCell(int const, int const);
   int nCreateGrid(void);
};
#endif // RASTERGRID_H
