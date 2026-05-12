/******************************************************************************
 *  Copyright (C) 2015 by Daniel Iglesias                                     *
 *                                                                            *
 *  This file is part of MkniX.                                             *
 *                                                                            *
 *  MkniX is free software: you can redistribute it and/or modify           *
 *  it under the terms of the GNU Lesser General Public License as            *
 *  published by the Free Software Foundation, either version 3 of the        *
 *  License, or (at your option) any later version.                           *
 *                                                                            *
 *  MkniX is distributed in the hope that it will be useful,                *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU Lesser General Public License for more details.                       *
 *                                                                            *
 *  You should have received a copy of the GNU Lesser General Public          *
 *  License along with MkniX.  If not, see <http://www.gnu.org/licenses/>.  *
 *****************************************************************************/

#include "analysisthermalstatic.h"
#include "simulation.h"

namespace mknix
{

AnalysisThermalStatic::AnalysisThermalStatic()
    : Analysis()
{
}


AnalysisThermalStatic::AnalysisThermalStatic( Simulation * simulation_in, double time_in )
    :   Analysis( simulation_in )
    , time( time_in )
{
    theProblem.setSystem( *theSimulation );
//   theProblem.setOutputFile("dis.dat", 0);
    theProblem.setResidue( &Simulation::staticThermalResidue );
    theProblem.setJacobian( &Simulation::staticThermalTangent );
    theProblem.setConvergence( &Simulation::staticThermalConvergence );
}


AnalysisThermalStatic::~AnalysisThermalStatic()
{
}


void AnalysisThermalStatic::solve( lmx::Vector< data_type >* q_in,
                                   lmx::Vector< data_type >* qdot_in = 0,
                                   lmx::Vector< data_type >* not_used = 0
                                 )
{
    theProblem.setInitialConfiguration( *q_in );
    theProblem.solve( 100 );
    *q_in = theProblem.getSolution();

//  std::ofstream disp("dis.dat");
//  disp << time << " ";
//  for(int i=0; i<q_in->size(); ++i){
//    disp << q_in->readElement(i) << " ";
//  }
//  disp << endl;
//  disp.close();
}


}
