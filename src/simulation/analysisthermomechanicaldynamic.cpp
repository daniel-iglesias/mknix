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
#include "analysisthermomechanicaldynamic.h"
#include "simulation.h"

namespace mknix {

AnalysisThermoMechanicalDynamic::AnalysisThermoMechanicalDynamic()
    : Analysis()
{
}


AnalysisThermoMechanicalDynamic::AnalysisThermoMechanicalDynamic
( Simulation* simulation_in,
  double to_in,
  double tf_in,
  double At_in,
  char * integrator_in
) : Analysis(simulation_in)
    , integratorType(integrator_in) // an additional integrator can be added
    , to(to_in)
    , tf(tf_in)
    , At(At_in)
{
    theProblem.setDiffSystem( *theSimulation );
    theProblem.setIntegrator1( integrator_in );
    theProblem.setIntegrator2( integrator_in );
    theProblem.setTimeParameters( to_in, tf_in, At_in );

    // Setting up the THERMAL problem
//  theProblem.setOutputFile("dis.dat", 0);
    if (theProblem.isIntegratorExplicit() ) {
        theProblem.setEvaluation1( &Simulation::explicitThermalEvaluation );
    }
    else {
        theProblem.setOutputFile1("flux.dat", 1);
        theProblem.setEvaluation1( &Simulation::dynamicThermalEvaluation );
        theProblem.setResidue1( &Simulation::dynamicThermalResidue );
        theProblem.setJacobian1( &Simulation::dynamicThermalTangent );
        if (epsilon == 0.0)
            theProblem.setConvergence1( &Simulation::dynamicThermalConvergenceInThermomechanical );
        else
            theProblem.setConvergence( epsilon );
    }

    // Setting up the MECHANICAL problem
    if (theProblem.isIntegratorExplicit() ) {
        theProblem.setEvaluation2( &Simulation::explicitAcceleration );
    }
    else {
        theProblem.setOutputFile2("vel.dat", 1);
        theProblem.setOutputFile2("acc.dat", 2);
        theProblem.setEvaluation2( &Simulation::dynamicAcceleration );
        theProblem.setResidue2( &Simulation::dynamicResidue );
        theProblem.setJacobian2( &Simulation::dynamicTangent );
        if (epsilon == 0.0)
            theProblem.setConvergence2( &Simulation::dynamicConvergence );
        else
            theProblem.setConvergence( epsilon );
    }

    theProblem.setStepTriggered( &Simulation::stepTriggered );
}


AnalysisThermoMechanicalDynamic::~AnalysisThermoMechanicalDynamic()
{
}


void AnalysisThermoMechanicalDynamic::solve
( lmx::Vector< data_type > * qt_in,
  lmx::Vector< data_type > * q_in,
  lmx::Vector< data_type >* qdot_in = 0
)
{
    if( lmx::getMatrixType() == 1 ) {
//    theProblem.setSparsePatternJacobian( theSimulation->getSparsePattern() ); // TBD for 1-DOF
        theProblem.setSparsePatternJacobian2( theSimulation->getSparsePattern() );
    }
    theProblem.setInitialConfiguration1( *qt_in );
    theProblem.setInitialConfiguration2( *q_in, *qdot_in );
    theProblem.solve();
}

void AnalysisThermoMechanicalDynamic::solve
( VectorX< data_type > * qt_in,
  VectorX< data_type > * q_in,
  VectorX< data_type >* qdot_in = 0)
{
    if( lmx::getMatrixType() == 1 ) {
//    theProblem.setSparsePatternJacobian( theSimulation->getSparsePattern() ); // TBD for 1-DOF
        theProblem.setSparsePatternJacobian2( theSimulation->getSparsePattern() );
    }
    theProblem.setInitialConfiguration1( *qt_in );
    theProblem.setInitialConfiguration2( *q_in, *qdot_in );
    theProblem.solve();
}


}
