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
#include "analysisdynamic.h"
#include "simulation.h"

namespace mknix {

AnalysisDynamic::AnalysisDynamic()
        : Analysis()
{
}


AnalysisDynamic::AnalysisDynamic(Simulation * simulation_in,
                                 double to_in,
                                 double tf_in,
                                 double At_in,
                                 const char * integrator_in,
                                 double par1 = -1.,
                                 double par2 = -1.,
                                 double par3 = -1.)
        : Analysis(simulation_in)
        , integratorType(integrator_in)
        , to(to_in)
        , tf(tf_in)
        , At(At_in)
{
    theProblem.setDiffSystem(*theSimulation);
    if (par1 == -1.) {
        theProblem.setIntegrator(integrator_in);
    }
    else if (par2 == -1.) {
        theProblem.setIntegrator(integrator_in, par1);
    }
    else if (par3 == -1.) {
        theProblem.setIntegrator(integrator_in, par1, par2);
    }
    else {
        theProblem.setIntegrator(integrator_in, par1, par2, par3);
    }
    theProblem.setTimeParameters(to_in, tf_in, At_in);

//  theProblem.setOutputFile("dis.dat", 0);
    theProblem.setOutputFile("vel.dat", 1);
    theProblem.setOutputFile("acc.dat", 2);
    theProblem.setEvaluation(&Simulation::dynamicAcceleration);
    theProblem.setResidue(&Simulation::dynamicResidue);
    theProblem.setJacobian(&Simulation::dynamicTangent);
    if (epsilon == 0.0) {
        theProblem.setConvergence(&Simulation::dynamicConvergence);
    } else {
        theProblem.setConvergence(epsilon);
    }
    theProblem.setStepTriggered(&Simulation::stepTriggered);
}


AnalysisDynamic::~AnalysisDynamic()
{
}


void AnalysisDynamic::solve(lmx::Vector<data_type> * q_in,
                            lmx::Vector<data_type> * qdot_in,
                            lmx::Vector<data_type> *
)
{
    if (lmx::getMatrixType() == 1) {
        theProblem.setSparsePatternJacobian(theSimulation->getSparsePattern());
    }
    theProblem.setInitialConfiguration(*q_in, *qdot_in);
    theProblem.solve();
}


void AnalysisDynamic::nextStep()
{
    theProblem.stepSolve();
}

void AnalysisDynamic::init(lmx::Vector<data_type> * qt_in, lmx::Vector<data_type> * qdot_in)
{
    theProblem.setInitialConfiguration(*qt_in, *qdot_in);
    theProblem.initialize();
}

}
