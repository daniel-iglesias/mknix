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
class RigidBody : public Body {
public:
    RigidBody();

    RigidBody( std::string );

    virtual ~RigidBody();

    void assembleMassMatrix( SparseMatrix<data_type> & );

    void assembleExternalForces( VectorX<data_type> & );

//     virtual Node* getDomainNode( std::string ); // Redefined by RigidBar

    void setMass( double mass_in ) {
        mass = mass_in;
    }
    void setDensityFactor( double density_in ) {
        densityFactor = density_in;
    }
    virtual void setInertia( double, int ) = 0;
    virtual void setPosition( std::vector<double>& ) = 0;

    void setOutput( std::string );

    Node* getNode( int );

    void outputStep( const VectorX<data_type>&, const VectorX<data_type>& );

    void outputStep( const VectorX<data_type>& );

    void outputToFile( std::ofstream* );

    void writeBodyInfo( std::ofstream* );

    virtual void writeBoundaryNodes( std::vector<Point*>& );

    virtual void writeBoundaryConnectivity(std::vector< std::vector<Point*> >&);


protected:
    bool computeEnergy;
    int dim;
    double mass, densityFactor;
    std::vector<Node*> frameNodes;
    std::vector< VectorX<data_type>* > domainConf;
    VectorX<data_type> externalForces;
    lmx::DenseMatrix< data_type > localMassMatrix;
    std::vector< VectorX<data_type>* > energy;
};

}

#endif
