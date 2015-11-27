/***************************************************************************
 *   Copyright (C) 2013 by Daniel Iglesias                                 *
 *   http://code.google.com/p/mknix                                        *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "boundarygroup.h"
#include "node.h"
#include "cellboundary.h"

namespace mknix {

BoundaryGroup::BoundaryGroup()
{
}

// /**
//  * @brief Constructor with 1 parameter
//  *
//  * @param title_in Name of boundary in the body.
//  **/
// BoundaryGroup::BoundaryGroup( std::string title_in )
//     : title( title_in )
// {
// }


BoundaryGroup::~BoundaryGroup()
{
  std::map<int,CellBoundary*>::iterator it_cells;
  for(it_cells=cells.begin();
      it_cells!=cells.end();
      ++it_cells){
    delete(it_cells->second);
  }
}

/**
 * @brief Cascade initialization funtion. Calls the initialize methods for each of the Cells
 *        and tells them to compute their shapefunction values. Both loops are parallelized.
 *
 * @return void
 **/
void BoundaryGroup::initialize()
{
    int end_int = this->cells.size();
    
//     #pragma omp parallel for
    for (int i=0;
            i < end_int;
            ++i)
    {
        this->cells[i]->initialize( this->nodes );
    }

//     #pragma omp parallel for
    for (int i=0;
            i < end_int;
            ++i)
    {
        this->cells[i]->computeShapeFunctions( );
    }

//  // Checking the output of a shapefunction:
//   int mid_int = this->cells.size()/2;
//   // Initialize individual output files
//   std::ofstream cell_data(std::string("cell_data_"+title+".dat").c_str());
//   std::ofstream gpoint_data(std::string("cell_gpoint_data_"+title+".dat").c_str());
//   this->cells[mid_int]->gnuplotOut(cell_data, gpoint_data); // Bus error
}

/**
 * @brief Computes the local volumetric heat vector of the material body by calling each cell's cascade function.
 *
 * @return void
 **/
void BoundaryGroup::calcExternalHeat( )
{
    int end_int = this->cells.size();
//     #pragma omp parallel for
    for (int i=0;
            i < end_int;
            ++i)
    {
        this->cells[i]->computeQextGaussPoints( this->loadThermalBoundaryGroup );
    }

}

/**
 * @brief Assembles the local volumetric heat into the global heat load vector by calling each cell's cascade function.
 *
 * @return void
 **/
void BoundaryGroup::assembleExternalHeat( lmx::Vector<data_type> & globalExternalHeat )
{
    int end_int = this->cells.size();
//     #pragma omp parallel for
    for (int i=0;
            i < end_int;
            ++i)
    {
        this->cells[i]->assembleQextGaussPoints( globalExternalHeat );
    }
}

}
