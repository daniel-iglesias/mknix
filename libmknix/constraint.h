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
    std::string title;
    std::string method;
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
