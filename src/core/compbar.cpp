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


#ifdef HAVE_VTK

#include "node.h"
#include "compbar.h"

#include <vtkLineSource.h>
#include <vtkTubeFilter.h>
#include <vtkActor.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>

#include <cmath>

namespace mknix
{

/**
 * @brief Default constructor for CompBar.
 */
CompBar::CompBar()
{
}


/**
 * @brief Constructs a CompBar between two nodes and creates the corresponding VTK visualization objects.
 * @param mat_in Material index for this bar.
 * @param nodeA_in Pointer to the first end node.
 * @param nodeB_in Pointer to the second end node.
 */
CompBar::CompBar(int mat_in, Node * nodeA_in, Node *nodeB_in)
    : mat(mat_in)
    , nodeA(nodeA_in)
    , nodeB(nodeB_in)
{
    double length;
    length = std::sqrt( std::pow(nodeA->getqx(0)-nodeB->getqx(0), 2) +
                        std::pow(nodeA->getqx(1)-nodeB->getqx(1), 2) +
                        std::pow(nodeA->getqx(2)-nodeB->getqx(2), 2) );
    line = vtkLineSource::New();
    line->SetResolution(10);

    lineTubes = vtkTubeFilter::New();
    lineTubes->SetInputConnection(line->GetOutputPort());
    lineTubes->SetRadius(length/20.);
    lineTubes->SetNumberOfSides(8);

    lineMapper = vtkPolyDataMapper::New();
    lineMapper->SetInputConnection( lineTubes->GetOutputPort() );

    lineActor = vtkActor::New();
    lineActor->SetMapper( lineMapper );
    lineActor->GetProperty()->SetColor(0.4235,0.6667,0.000);

}


/**
 * @brief Destructor for CompBar.
 */
CompBar::~CompBar()
{
}


/**
 * @brief Updates the VTK line source endpoints to reflect the current node positions.
 */
void CompBar::updatePoints()
{
    line->SetPoint1(nodeA->getqx(0),
                    nodeA->getqx(1),
                    nodeA->getqx(2)
                   );
    line->SetPoint2(nodeB->getqx(0),
                    nodeB->getqx(1),
                    nodeB->getqx(2)
                   );
}

/**
 * @brief Adds this bar's VTK actor to the given renderer.
 * @param renderer_in VTK renderer to which the actor is added.
 */
void CompBar::addToRender(vtkRenderer * renderer_in)
{
    renderer_in->AddActor( lineActor );

}

/**
 * @brief Removes this bar's VTK actor from the given renderer.
 * @param renderer_in VTK renderer from which the actor is removed.
 */
void CompBar::removeFromRender(vtkRenderer * renderer_in)
{
    renderer_in->RemoveActor( lineActor );
}

}

#endif // HAVE_VTK