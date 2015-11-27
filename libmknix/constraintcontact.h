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
#ifndef MKNIXCONSTRAINTCONTACT_H
#define MKNIXCONSTRAINTCONTACT_H

#include <constraint.h>

namespace mknix {

/**
	@author AUTHORS <MAILS>
*/
class ConstraintContact : public Constraint
{
public:
    ConstraintContact();

    ~ConstraintContact();

    ConstraintContact( Node* , Node*, Node* , double&, std::string& );

    void calcPhi( ) ;

    void calcPhiq( ) ;

    void calcPhiqq( ) ;

    lmx::Vector<data_type>& getInternalForces( ) {
        return this->internalForces;
    }

    lmx::DenseMatrix<data_type>& getStiffnessMatrix( ) {
        return this->stiffnessMatrix;
    }

    double getGap() {
        return rt;
    }

protected:
    std::vector<double> normal;
    double rh;
    double rt;
};

}

#endif
