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
#ifndef MKNIXELEMTETRAHEDRON_H
#define MKNIXELEMTETRAHEDRON_H

#include <celltetrahedron.h>

namespace mknix {

/**
	@author AUTHORS <MAILS>
*/
class ElemTetrahedron : public CellTetrahedron
{
public:
    ElemTetrahedron();

//    ElemTetrahedron( Material&, double, int,
//                  double, double,
//                  double, double,
//                  double, double );

    ElemTetrahedron( Material&, double, int,
                     Node*,
                     Node*,
                     Node*,
                     Node*
                   );

    ~ElemTetrahedron();

    void initialize( std::vector<Node*> & );

    void computeShapeFunctions(  );

    void createGaussPoints_MC();
};

}

#endif
