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

#ifndef MKNIXCONSTRAINT_H
#define MKNIXCONSTRAINT_H

#include "common.h"
#include "LMX/lmx.h"


namespace mknix {
class Node;
/**
	@author AUTHORS <MAILS>
*/
class Constraint {
public:
    Constraint();

    Constraint(double&, std::string&);

    Constraint(double&, std::string&, int);

    virtual ~Constraint();

    void setTitle( std::string& title_in )
    { title = title_in; }

    VectorX<data_type>& getInternalForces( ) {
        return this->internalForces;
    }

    lmx::DenseMatrix<data_type>& getStiffnessMatrix( ) {
        return this->stiffnessMatrix;
    }

    virtual void writeJointInfo( std::ofstream* );

    virtual void calcPhi( ) = 0 ;

    virtual void calcPhiq( ) = 0 ;

    virtual void calcPhiqq( ) = 0 ;

    virtual void calcInternalForces( );

    virtual void calcTangentMatrix( );

    virtual void assembleInternalForces( VectorX<data_type> & );

    virtual void assembleTangentMatrix( SparseMatrix<data_type> & );

    virtual bool checkAugmented();

    virtual void clearAugmented();

    virtual Node* getNode( int nodeNumber )
    {
        return nodes[nodeNumber];
    }

    void outputStep( const VectorX<data_type>&, const VectorX<data_type>& );

    void outputStep( const VectorX<data_type>& );

    void outputToFile( std::ofstream* );

protected:
    int dim, iter_augmented;
    double alpha;
    std::string method;
    std::string title;
    std::vector<Node*> nodes;
    VectorX<data_type> internalForces;
    std::vector< VectorX<data_type> > internalForcesOutput;
    lmx::DenseMatrix<data_type> stiffnessMatrix;
    std::vector< double > lambda;
    std::vector< double > phi;
    std::vector< VectorX<data_type> > phi_q;
    std::vector< lmx::DenseMatrix<data_type> > phi_qq;

};

}

#endif
