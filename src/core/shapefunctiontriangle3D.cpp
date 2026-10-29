/******************************************************************************
 *  Copyright (C) 2015 by Daniel Iglesias                                     *
 *                                                                            *
 *  This file is part of MkniX.                                             *
 *                                                                            *
 *  MkniX is free software: you can redistribute it and/or modify           *
 *  it under the terms of the GNU Lesser General Public License as            *
 *  published by the Free Software Foundation, either version 3 of the        *
 *  License, or (at your option) any later version.                           *
 *                                                                            *
 *  MkniX is distributed in the hope that it will be useful,                *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU Lesser General Public License for more details.                       *
 *                                                                            *
 *  You should have received a copy of the GNU Lesser General Public          *
 *  License along with MkniX.  If not, see <http://www.gnu.org/licenses/>.  *
 *****************************************************************************/

#include "shapefunctiontriangle3D.h"
#include "point.h"
#include "node.h"

namespace mknix
{

/**
 * @brief Default constructor for ShapeFunctionTriangleSigned.
 */
ShapeFunctionTriangleSigned::ShapeFunctionTriangleSigned()
    : ShapeFunction()
{
//   this->phi.resize(6, 3);
}


/**
 * @brief Constructs a signed triangular shape function for a triangle embedded in 3D space.
 * @param gp_in Pointer to the evaluation Point (expects exactly 3 support nodes).
 */
ShapeFunctionTriangleSigned::ShapeFunctionTriangleSigned( Point* gp_in )
    : ShapeFunction(gp_in)
{
    this->phi.resize(3, 3);
}


/**
 * @brief Destructor for ShapeFunctionTriangleSigned.
 */
ShapeFunctionTriangleSigned::~ShapeFunctionTriangleSigned()
{
}


/**
 * @brief Computes barycentric (signed-area) shape functions for a triangle embedded in 3D space,
 *        using cross products to obtain area ratios relative to the total triangle area.
 */
void ShapeFunctionTriangleSigned::calc()
{
    const Point* first = gp->getSupportNodes()[0];
    const Point* second = gp->getSupportNodes()[1];
    const Point* third = gp->getSupportNodes()[2];

    const double v0x = second->getX() - first->getX();
    const double v0y = second->getY() - first->getY();
    const double v0z = second->getZ() - first->getZ();
    const double v1x = third->getX() - first->getX();
    const double v1y = third->getY() - first->getY();
    const double v1z = third->getZ() - first->getZ();
    const double v2x = gp->getX() - first->getX();
    const double v2y = gp->getY() - first->getY();
    const double v2z = gp->getZ() - first->getZ();

    const double d00 = v0x * v0x + v0y * v0y + v0z * v0z;
    const double d01 = v0x * v1x + v0y * v1y + v0z * v1z;
    const double d11 = v1x * v1x + v1y * v1y + v1z * v1z;
    const double d20 = v2x * v0x + v2y * v0y + v2z * v0z;
    const double d21 = v2x * v1x + v2y * v1y + v2z * v1z;
    const double denominator = d00 * d11 - d01 * d01;
    const double secondValue = (d11 * d20 - d01 * d21) / denominator;
    const double thirdValue = (d00 * d21 - d01 * d20) / denominator;

    phi.writeElement(1.0 - secondValue - thirdValue, 0, 0);
    phi.writeElement(secondValue, 0, 1);
    phi.writeElement(thirdValue, 0, 2);
}


}
