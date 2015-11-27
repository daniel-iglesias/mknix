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
#ifndef MKNIXCONSTRAINTTHERMALFIXED_H
#define MKNIXCONSTRAINTTHERMALFIXED_H

#include "constraintthermal.h"

namespace mknix {

/**
	@author AUTHORS <MAILS>
*/
class ConstraintThermalFixed : public ConstraintThermal
{
public:
    ConstraintThermalFixed();

    ConstraintThermalFixed( Node* , Node*, double&, std::string& );

    ~ConstraintThermalFixed();

    void calcPhi( ) ;

    void calcPhiq( ) ;

    void calcPhiqq( ) ;

    lmx::Vector<data_type>& getInternalForces( ) {
        return this->internalForces;
    }

    lmx::DenseMatrix<data_type>& getStiffnessMatrix( ) {
        return this->stiffnessMatrix;
    }

protected:
    double To;
    double Tt;

};

}

#endif
