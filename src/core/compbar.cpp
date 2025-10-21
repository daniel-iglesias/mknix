/***************************************************************************
 *   Copyright (C) 2007 by Daniel Iglesias   *
 *   daniel@extremo   *
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

CompBar::CompBar()
{
}


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


CompBar::~CompBar()
{
}


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

void CompBar::addToRender(vtkRenderer * renderer_in)
{
    renderer_in->AddActor( lineActor );

}

void CompBar::removeFromRender(vtkRenderer * renderer_in)
{
    renderer_in->RemoveActor( lineActor );
}

}

#endif // HAVE_VTK