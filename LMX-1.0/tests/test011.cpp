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


// #define HAVE_GMM
// #define HAVE_SUPERLU

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

// #include <iostream>
// #include <cstdlib>

#include"LMX/lmx.h"
#include"LMX/lmx_diff_problem_second.h"

using namespace std;

class MyDiffSystem
{
public:
    MyDiffSystem() {}

    ~MyDiffSystem() {}

    void myEvaluation( const lmx::Vector<double>& q,
                       const lmx::Vector<double>& qdot,
                       lmx::Vector<double>& qddot,
                       double time
                     )
    {
//      qddot.fillIdentity( time );
        qddot(0) = time;
        qddot(1) = time;
    }
};

int main(int argc, char* argv[])
{
    lmx::setMatrixType( 0 );
    lmx::setVectorType( 0 );
    lmx::setLinSolverType( 0 );

    lmx::DiffProblemSecond< MyDiffSystem > theProblem;
    MyDiffSystem theSystem;
    lmx::Vector<double> q0(2);
    lmx::Vector<double> qdot0(2);

    theProblem.setDiffSystem( theSystem );
    theProblem.setIntegrator( "AB-1" );
    theProblem.setInitialConfiguration( q0, qdot0 );
    theProblem.setTimeParameters( 0, 1, 0.1 );
    theProblem.setEvaluation( &MyDiffSystem::myEvaluation );
    theProblem.solve();

    cout << "End configuration: " << theProblem.getConfiguration( 0, 0);

    return EXIT_SUCCESS;
}
