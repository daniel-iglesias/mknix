/******************************************************************************
 *  Copyright (C) 2015 by Daniel Iglesias                                     *
 *                                                                            *
 *  This file is part of Nemesis.                                             *
 *                                                                            *
 *  Nemesis is free software: you can redistribute it and/or modify           *
 *  it under the terms of the GNU Lesser General Public License as            *
 *  published by the Free Software Foundation, either version 3 of the        *
 *  License, or (at your option) any later version.                           *
 *                                                                            *
 *  Nemesis is distributed in the hope that it will be useful,                *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU Lesser General Public License for more details.                       *
 *                                                                            *
 *  You should have received a copy of the GNU Lesser General Public          *
 *  License along with Nemesis.  If not, see <http://www.gnu.org/licenses/>.  *
 *****************************************************************************/

#ifndef MKNIXRIGIDBODY_H
#define MKNIXRIGIDBODY_H

#include "body.h"
// #include "node.h"

namespace mknix {

class Node;

/**
	@author AUTHORS <MAILS>
*/
class RigidBody : public Body
{
public:
    RigidBody();

    RigidBody(std::string);

    virtual ~RigidBody();

    void assembleMassMatrix(lmx::Matrix<data_type>&) override;

    void assembleExternalForces(lmx::Vector<data_type>&) override;

//     virtual Node* getDomainNode( std::string ); // Redefined by RigidBar

    void setMass(double mass_in)
    {
        mass = mass_in;
    }

    void setDensityFactor(double density_in)
    {
        densityFactor = density_in;
    }

    virtual void setInertia(double, int) = 0;

    virtual void setPosition(std::vector<double>&) = 0;

    void setOutput(std::string) override;

    Node * getNode(int) override;

    void outputStep(const lmx::Vector<data_type>&, const lmx::Vector<data_type>&) override;

    void outputStep(const lmx::Vector<data_type>&) override;

    void outputToFile(std::ofstream *) override;

    void writeBodyInfo(std::ofstream *) override;

    virtual void writeBoundaryNodes(std::vector<Point *>&) override;

    virtual void writeBoundaryConnectivity(std::vector<std::vector<Point *> >&) override;


protected:
    bool computeEnergy;
    int dim;
    double mass, densityFactor;
    std::vector<Node *> frameNodes;
    std::vector<lmx::Vector<data_type> *> domainConf;
    lmx::Vector<data_type> externalForces;
    lmx::DenseMatrix<data_type> localMassMatrix;
    std::vector<lmx::Vector<data_type> *> energy;
};

}

#endif
