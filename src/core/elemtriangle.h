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

#ifndef MKNIXELEMTRIANGLE_H
#define MKNIXELEMTRIANGLE_H

#include "celltriang.h"

namespace mknix {

/**
	@author AUTHORS <MAILS>
*/
class ElemTriangle : public CellTriang
{
public:
    ElemTriangle();

//    ElemTriangle( Material&, double, int,
//                  double, double,
//                  double, double,
//                  double, double );

    ElemTriangle( Material&, double, int,
                  Node*,
                  Node*,
                  Node*
                );

    ~ElemTriangle();

    void initialize( std::vector<Node*> & );

    void computeShapeFunctions(  );

protected:
    void createGaussPoints_MC( );

};

}

#endif
