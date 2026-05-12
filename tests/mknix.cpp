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



#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>

#include <simulation/simulation.h>

using namespace std;

int main(int argc, char * argv[])
{
//   try{
    if (argc > 2)
    {
        lmx::setMatrixType(atoi(argv[2]));
        if (argc == 4)
        {
            lmx::setLinSolverType(atoi(argv[3]));
        }
    }
    if (argc >= 2 && argc < 5)
    {
//     cout << "EXECUTING: $ mknix " << argv[1] << endl;
//     system("echo -n '1. Current Directory is '; pwd");
        mknix::Simulation mySimulation;
        cout << argv[1];
        mySimulation.inputFromFile(argv[1]);
        mySimulation.run();
//       mySimulation.outputFile("mknix.msg");
//     mySimulation.geometryFile("geometry.mbr");
//     mySimulation.resultsFile("res.mbr");
//       mySimulation.run();

    }
    else if (argc != 2)
    {
        cout << "Need at least one parameter with file to read from." << endl
             << "Usage: $./mknix \"file_name\"" << endl;
    }
    return EXIT_SUCCESS;
}
