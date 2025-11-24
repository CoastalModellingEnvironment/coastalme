/*!
   \file raster_grid_inline.h
   \brief Inline implementations for CGeomRasterGrid
   \details This file contains inline method implementations that require the full
            CGeomCell definition. It must be included after both raster_grid.h and cell.h
            to avoid circular dependency issues.
   \author David Favis-Mortlock
   \author Andres Payo
   \date 2025
   \copyright GNU General Public License
*/

#ifndef RASTERGRID_INLINE_H
#define RASTERGRID_INLINE_H

// Ensure required headers are included before this file
#ifndef RASTERGRID_H
#error "raster_grid.h must be included before raster_grid_inline.h"
#endif

#ifndef CELL_H
#error "cell.h must be included before raster_grid_inline.h"
#endif

// Inline accessor implementations - eliminates function call overhead in hot loops
inline CGeomCell& CGeomRasterGrid::Cell(int const nX, int const nY)
{
   return m_CellData[nX * m_nYSize + nY];
}

inline CGeomCell const& CGeomRasterGrid::Cell(int const nX, int const nY) const
{
   return m_CellData[nX * m_nYSize + nY];
}

#endif // RASTERGRID_INLINE_H
