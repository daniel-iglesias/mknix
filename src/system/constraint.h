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
  
    virtual void writeJointInfo( std::ofstream* );

    virtual void calcPhi( ) = 0 ;

    virtual void calcPhiq( ) = 0 ;

    virtual void calcPhiqq( ) = 0 ;

    virtual void calcInternalForces( );

    virtual void calcTangentMatrix( );

    virtual void assembleInternalForces( lmx::Vector<data_type> & );

    virtual void assembleTangentMatrix( lmx::Matrix<data_type> & );

    virtual bool checkAugmented();

    virtual void clearAugmented();

    virtual Node* getNode( int nodeNumber )
    {
        return nodes[nodeNumber];
    }

    void outputStep( const lmx::Vector<data_type>&, const lmx::Vector<data_type>& );

    void outputStep( const lmx::Vector<data_type>& );

    void outputToFile( std::ofstream* );

protected:
    int dim;
    double alpha;
    std::string method;
    std::string title;
    std::vector<Node*> nodes;
    lmx::Vector<data_type> internalForces;
    std::vector< lmx::Vector<data_type> > internalForcesOutput;
    lmx::DenseMatrix<data_type> stiffnessMatrix;
    std::vector< double > lambda;
    std::vector< double > phi;
    std::vector< lmx::Vector<data_type> > phi_q;
    std::vector< lmx::DenseMatrix<data_type> > phi_qq;

};

}

#endif
