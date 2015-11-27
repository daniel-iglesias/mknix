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
#ifndef MKNIXGENERALCONTACT_H
#define MKNIXGENERALCONTACT_H

#include <vector>
#include <map>

class vtkPoints;
class vtkCellArray;
class vtkPolyData;
class vtkDelaunay2D;
class vtkLineSource;
class vtkPolyDataMapper;
class vtkActor;
class vtkExtractEdges;
class vtkTubeFilter;
class vtkPolyDataMapper;
class vtkActor;
class vtkRenderer;
class vtkRenderWindow;
class vtkRenderWindowInteractor;
class vtkWindowToImageFilter;
class vtkJPEGWriter;

namespace mknix {

class Simulation;
class Node;
class Constraint;
class CompBar;

class Contact {
public:
    Contact( );
    Contact( Simulation*, double );
    ~Contact() {}

    void createPoints();

    void updatePoints();

    void createPolys();

    void updateLines();

    void createDelaunay();

    void updateDelaunay();

    void createDrawingObjects();

    void createDrawingContactObjects();

    void drawObjects();

private:
    void readDelaunay();

    bool orderTriangleNodes();

    void updateContactElements();

private:
    Simulation* theSimulation;
    double delta;
    double alpha;
    int image_number;
    vtkPoints* points;
    vtkCellArray* polys;
    vtkPolyData* polyData;
    vtkDelaunay2D* delny;
    vtkPolyData* delnyData;
    vtkPolyDataMapper* mapMesh;
    vtkActor* meshActor;
    vtkExtractEdges* extract;
    vtkTubeFilter* tubes;
    vtkPolyDataMapper* mapEdges;
    vtkActor* edgeActor;
    vtkRenderer* ren;
    vtkRenderWindow* renWin;
    vtkRenderWindowInteractor *iren;
    vtkWindowToImageFilter* w2if;
    vtkJPEGWriter* psw;

    std::vector< Node* > nodes;
    std::map< int,int > mapNodes;
    std::vector< std::vector< Node* > > boundaries;
    int boundary_first_number;
    std::vector<int> vTriangles;
    std::vector< Constraint* > constraints;
    std::vector<CompBar*> bars;
};


}

#endif
