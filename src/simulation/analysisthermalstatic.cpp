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
#include "analysisthermalstatic.h"
#include "simulation.h"

namespace mknix {

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
    if(OLD_CODE)theProblem.setConvergence( static_cast<void (*)(lmx::Vector<data_type>, lmx::Vector<data_type>)>(&Simulation::staticThermalConvergence) );
    else theProblem.setConvergence( static_cast<void (*)(VectorX<data_type>, VectorX<data_type>)>(&Simulation::staticThermalConvergence) );
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

void AnalysisThermalStatic::solve( VectorX< data_type >* q_in,
                                   VectorX< data_type >* qdot_in = 0,
                                   VectorX< data_type >* not_used = 0)
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
