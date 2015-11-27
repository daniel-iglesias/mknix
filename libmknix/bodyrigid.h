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

    void assembleMassMatrix( lmx::Matrix<data_type> & );

    void assembleExternalForces( lmx::Vector<data_type> & );

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

    void outputStep( const lmx::Vector<data_type>&, const lmx::Vector<data_type>& );

    void outputStep( const lmx::Vector<data_type>& );

    void outputToFile( std::ofstream* );

    void writeBodyInfo( std::ofstream* );

    virtual void writeBoundaryNodes( std::vector<Point*>& );

    virtual void writeBoundaryConnectivity(std::vector< std::vector<Point*> >&);


protected:
    int dim;
    double mass, densityFactor;
    std::vector<Node*> frameNodes;
    std::vector< lmx::Vector<data_type>* > domainConf;
    lmx::Vector<data_type> externalForces;
    lmx::DenseMatrix< data_type > localMassMatrix;
    std::vector< lmx::Vector<data_type>* > energy;
    bool computeEnergy;
};

}

#endif
