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

#ifndef COMPBAR_H
#define COMPBAR_H


class vtkRenderer;
class vtkLineSource;
class vtkTubeFilter;
class vtkPolyDataMapper;
class vtkActor;

/**
	@author Daniel Iglesias <daniel@extremo>
*/
namespace mknix {

class Node;

class CompBar {
public:
    CompBar();

    CompBar(int, Node*, Node*);

    ~CompBar();

    void addToRender( vtkRenderer* );

    void removeFromRender( vtkRenderer* );

    void updatePoints();

private:
    int mat;
    Node* nodeA;
    Node* nodeB;
    vtkLineSource *line;
    vtkTubeFilter *lineTubes;
    vtkPolyDataMapper *lineMapper;
    vtkActor *lineActor;
};

}

#endif
