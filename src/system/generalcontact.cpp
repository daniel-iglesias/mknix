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

#ifdef USE_VTK

#include "generalcontact.h"
#include "constraintcontact.h"
#include "system.h"

#include <core/compbar.h>
#include <simulation/simulation.h>

#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkPolyData.h"
#include "vtkDelaunay2D.h"
#include <vtkLineSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCamera.h>
#include <vtkExtractEdges.h>
#include <vtkTubeFilter.h>
#include <vtkProperty.h>
#include <vtkExtractPolyDataGeometry.h>
#include <vtkWindowToImageFilter.h>
#include <vtkJPEGWriter.h>

namespace mknix
{

Contact::Contact( )
{
}

Contact::Contact( Simulation* theSimulation_in, double alpha_in )
    : theSimulation(theSimulation_in)
    , delta(alpha_in*2.5)
    , alpha(alpha_in)
    , image_number(0)
{
}

void Contact::createPoints()
{
    points = vtkPoints::New();
    std::ofstream outfile;
//    cout<<" outfile.is_open() = "<<outfile.is_open()<<endl;
//    cout<<" outfile.fail() = "<<outfile.fail()<<endl;
//    cout<<" outfile.bad() = "<<outfile.bad()<<endl;
//    cout<<" outfile.good() = "<<outfile.good()<<endl;
    outfile.open("boundary_nodes");
    theSimulation->baseSystem->writeBoundaryNodes( nodes );
//     cout << "Contact Points: " << nodes.size() << endl;

    for(size_type i=0; i<nodes.size(); ++i)
    {
        points->InsertPoint(i,
                            nodes[i]->getx(),
                            nodes[i]->gety(),
                            nodes[i]->getz()
                           );
        mapNodes[nodes[i]->getNumber()] = i;
//      cout<<" outfile.is_open() = "<<outfile.is_open()<<endl;
//      cout<<" outfile.fail() = "<<outfile.fail()<<endl;
//      cout<<" outfile.bad() = "<<outfile.bad()<<endl;
//      cout<<" outfile.good() = "<<outfile.good()<<endl;
        outfile << nodes[i]->getNumber() << "(";
//      cout<<"nodes[i]->getNumber() = " <<nodes[i]->getNumber()<< endl;
//      cout<<"nodes[i]->getx() = " <<nodes[i]->getx()<< endl;
//      cout<<"nodes[i]->gety() = " <<nodes[i]->gety()<< endl;
        outfile << nodes[i]->getx() <<","
                << nodes[i]->gety() <<","
                //<< nodes[i]->getz() << ")"
                << endl;
//      cout << "OK! " << i << endl;
//      points->GetData()->Print( cout);
//      boundary_first_number = points->GetData()->GetMaxId() + 1;
//      cout << "BOUNDARY FIRST NUMBER = " << boundary_first_number << endl;
    }
    boundary_first_number = points->GetData()->GetMaxId() + 1;
    double bounds[6];
    points->GetBounds(bounds);
    points->InsertPoint(boundary_first_number,
                        bounds[0]-delta,
                        bounds[2]-delta,
                        bounds[4]
                       );
    points->InsertPoint(boundary_first_number+1,
                        bounds[1]+delta,
                        bounds[2]-delta,
                        bounds[4]
                       );
    points->InsertPoint(boundary_first_number+2,
                        bounds[1]+delta,
                        bounds[3]+delta,
                        bounds[4]
                       );
    points->InsertPoint(boundary_first_number+3,
                        bounds[0]-delta,
                        bounds[3]+delta,
                        bounds[4]
                       );

    cout << "BOUNDS = (" << bounds[0]
         << "," << bounds[1]
         << "," << bounds[2]
         << "," << bounds[3]
         << "," << bounds[4]
         << "," << bounds[5]
         << ")" << endl;
}

void Contact::updatePoints()
{
//    nodes.clear();
//    theSimulation->baseSystem->writeBoundaryNodes( nodes );
//    cout << "Contact Points: " << nodes.size() << endl;

    for(size_type i=0; i<nodes.size(); ++i)
    {
        points->SetPoint(   i,
                            nodes[i]->getx(),
                            nodes[i]->gety(),
                            nodes[i]->getz()
                        );
    }
//  boundary_first_number = points->GetData()->GetMaxId() + 1;
    double bounds[6];
    points->GetBounds(bounds);
    points->SetPoint(   boundary_first_number,
                        bounds[0]-delta,
                        bounds[2]-delta,
                        bounds[4]
                    );
    points->SetPoint(   boundary_first_number+1,
                        bounds[1]+delta,
                        bounds[2]-delta,
                        bounds[4]
                    );
    points->SetPoint(boundary_first_number+2,
                     bounds[1]+delta,
                     bounds[3]+delta,
                     bounds[4]
                    );
    points->SetPoint(  boundary_first_number+3,
                       bounds[0]-delta,
                       bounds[3]+delta,
                       bounds[4]
                    );

    cout << "BOUNDS = (" << bounds[0]
         << "," << bounds[1]
         << "," << bounds[2]
         << "," << bounds[3]
         << "," << bounds[4]
         << "," << bounds[5]
         << ")" << endl;
}

void Contact::createPolys()
{
    polys = vtkCellArray::New();
    theSimulation->baseSystem->writeBoundaryConnectivity( boundaries );

    //boundary:
    polys->InsertNextCell( 4 );
    polys->InsertCellPoint( boundary_first_number );
    polys->InsertCellPoint( boundary_first_number+1 );
    polys->InsertCellPoint( boundary_first_number+2 );
    polys->InsertCellPoint( boundary_first_number+3 );

    std::vector< std::vector< Node* > >::iterator it_boundaries;
    std::vector< Node* >::reverse_iterator it_nodes;
    std::vector< std::vector< Node* > >::iterator
    it_boundaries_end = boundaries.end();
    std::vector< Node* >::reverse_iterator it_nodes_end;
    std::ofstream outfile("boundary");
    for( it_boundaries = boundaries.begin();
            it_boundaries!= it_boundaries_end;
            ++it_boundaries )
    {
        polys->InsertNextCell( it_boundaries->size() );
        it_nodes_end = it_boundaries->rend();
        for( it_nodes = it_boundaries->rbegin();
                it_nodes!= it_nodes_end;
                ++it_nodes
           )
        {
//         cout << "Node A:" << (*it_nodes)->getNumber() << cout.flush();
            polys->InsertCellPoint( mapNodes[(*it_nodes)->getNumber()] );
            outfile << (*it_nodes)->getNumber() << " "
                    << (*it_nodes)->getx() << " "
                    << (*it_nodes)->gety() << " "
                    << (*it_nodes)->getz() << endl;
        }
        outfile << endl;
    }

    polyData = vtkPolyData::New();
    polyData->SetPoints(points);
    polyData->SetPolys(polys);

    int mat=0;
    //Create lines in vector:
    for( it_boundaries = boundaries.begin();
            it_boundaries!= it_boundaries_end;
            ++it_boundaries )
    {
        polys->InsertNextCell( it_boundaries->size() );
        it_nodes_end = --it_boundaries->rend();
        for( it_nodes = it_boundaries->rbegin();
                it_nodes!= it_nodes_end;
                ++it_nodes
           )
        {
//         cout << "Node 1:" << (*it_nodes)->getNumber() << cout.flush();
//         cout << ", Node 2: " << (*it_nodes+1)->getNumber() << endl;
            bars.push_back( new CompBar(mat, *it_nodes, *(it_nodes+1) ) );
            bars.back()->updatePoints( );
        }
        ++mat;
    }
    this->updateLines();
}

void Contact::updateLines()
{
    std::vector<CompBar*>::iterator it_bars;
    for( it_bars = bars.begin();
            it_bars!= bars.end();
            ++it_bars )
    {
        (*it_bars)->updatePoints( );
    }
}

void Contact::readDelaunay()
{
    int j=0;
    int h;
    vtkIdType npts=3;
    vtkIdType *pts;
    vtkCellArray* oCellArray = vtkCellArray::New();
    cout << "BOUNDARY FIRST NODE: " << boundary_first_number << endl;

    vTriangles.clear();
    oCellArray = delnyData->GetPolys();
    for(int i=0; i<delnyData->GetNumberOfPolys(); ++i)
    {
        h=oCellArray->GetNextCell(npts, pts);
        if(h==0)
            break;
        if(npts==3)
        {
            if( pts[0] < nodes.size() &&
                    pts[1] < nodes.size() &&
                    pts[2] < nodes.size() )
            {
                vTriangles.push_back(pts[0]);
                vTriangles.push_back(pts[1]);
                vTriangles.push_back(pts[2]);
                cout << "Triangle " << h << " = ("
                     << vTriangles[j] << ", "
                     << vTriangles[j+1] << ", "
                     << vTriangles[j+2] << ")" << endl;
                j+=3;
            }
        }
    }
    this->orderTriangleNodes();
    this->updateContactElements();
}

bool Contact::orderTriangleNodes()
{
    size_type i;
    int j, k; //maybe j isn't necessary
    Node* first_node, * second_node, * third_node;
    std::vector< std::vector< Node* > >::iterator it_boundaries;
    std::vector< Node* >::iterator it_nodes;
    size_type number_of_triangles = vTriangles.size()/3;

    for(i=0; i<number_of_triangles; ++i)
    {
        cout << vTriangles[3*i] << "," << std::flush;
        cout << nodes[vTriangles[3*i]]->getNumber() << "; " << std::flush;
        cout << vTriangles[3*i+1] << "," << std::flush;
        cout << nodes[vTriangles[3*i+1]]->getNumber() << "; " << std::flush;
        cout << vTriangles[3*i+2] << "," << std::flush;
        cout << nodes[vTriangles[3*i+2]]->getNumber() << endl;
        first_node  = theSimulation->nodes[ nodes[vTriangles[3*i]]->getNumber() ];
        second_node = theSimulation->nodes[ nodes[vTriangles[3*i+1]]->getNumber() ];
        third_node  = theSimulation->nodes[ nodes[vTriangles[3*i+2]]->getNumber() ];
        for( it_boundaries = boundaries.begin();
                it_boundaries!= boundaries.end();
                ++it_boundaries )
        {
            it_nodes = find( it_boundaries->begin(),
                             it_boundaries->end(),
                             first_node
                           );
            if( it_nodes != it_boundaries->end() )
            {
                j = it_boundaries - boundaries.begin(); // num. of body where the node is
                k = it_nodes - it_boundaries->begin(); // position in boundary
                break; //node found, stop searching
            }
        }
        it_nodes = find( it_boundaries->begin(),
                         it_boundaries->end(),
                         second_node
                       );
        if( it_nodes != it_boundaries->end() )
        {
            if( it_nodes - it_boundaries->begin() < k )
            {
                vTriangles[3*i] = mapNodes[second_node->getNumber()];
                vTriangles[3*i+1] = mapNodes[first_node->getNumber()];
            }
            // else : triangle is ordered, nothing to do.
//        // Now we look for three adjacent nodes in the boundary... TBD
            it_nodes = find( it_boundaries->begin(),
                             it_boundaries->end(),
                             third_node
                           );
            // this can be done more efficiently, sure! posible bug:
//        if( it_nodes != it_boundaries->end() ){ //if its found in the same boundary, closed to the other points...
//          vTriangles[3*i] = -1;
//          vTriangles[3*i+1] = -1;
//          vTriangles[3*i+2] = -1;
//        }
        }
        else   //node2 not found in body of node1... searching for node3:
        {
            it_nodes = find( it_boundaries->begin(),
                             it_boundaries->end(),
                             third_node
                           );
            if( it_nodes != it_boundaries->end() )
            {
                if( it_nodes - it_boundaries->begin() < k )
                {
                    vTriangles[3*i] = mapNodes[third_node->getNumber()];
                    vTriangles[3*i+1] = mapNodes[first_node->getNumber()];
                }
                else
                {
                    vTriangles[3*i] = mapNodes[first_node->getNumber()];
                    vTriangles[3*i+1] = mapNodes[third_node->getNumber()];
                }
                vTriangles[3*i+2] = mapNodes[second_node->getNumber()];
            }
            else   //node 2 & 3 are in the same body... searching:
            {
                for( it_boundaries = boundaries.begin();
                        it_boundaries!= boundaries.end();
                        ++it_boundaries )
                {
                    it_nodes = find( it_boundaries->begin(),
                                     it_boundaries->end(),
                                     second_node
                                   );
                    if( it_nodes != it_boundaries->end() )
                    {
                        j = it_boundaries - boundaries.begin(); // num. of body where the node is
                        k = it_nodes - it_boundaries->begin(); // position in boundary
                        break; //node found, stop searching
                    }
                }
                it_nodes = find( it_boundaries->begin(),
                                 it_boundaries->end(),
                                 third_node
                               );
                if( it_nodes != it_boundaries->end() )
                {
                    if( it_nodes - it_boundaries->begin() < k )
                    {
                        vTriangles[3*i] = mapNodes[third_node->getNumber()];
                        vTriangles[3*i+1] = mapNodes[second_node->getNumber()];
                    }
                    else
                    {
                        vTriangles[3*i] = mapNodes[second_node->getNumber()];
                        vTriangles[3*i+1] = mapNodes[third_node->getNumber()];
                    }
                    vTriangles[3*i+2] = mapNodes[first_node->getNumber()];
                }
            }
        }
        cout << "Triangle " << i << " = ("
             << vTriangles[3*i] << ", "
             << vTriangles[3*i+1] << ", "
             << vTriangles[3*i+2] << ")" << endl;
    }
    return true;
}

void Contact::updateContactElements()
{
    std::vector< Constraint* >::iterator it_constraints;
    std::map< std::string, Constraint* >::iterator it_to_constraint;
    std::vector< std::string > indexes;
    for(it_to_constraint = theSimulation->baseSystem->constraints.begin();
            it_to_constraint!= theSimulation->baseSystem->constraints.end();
            ++it_to_constraint
       )
    {
        // deleting the pointers from basesystem:
        it_constraints = find( constraints.begin(),
                               constraints.end(),
                               it_to_constraint->second );
        if( it_constraints != constraints.end() )
        {
            if( (static_cast<ConstraintContact*>(*it_constraints))->getGap() > 0 )
            {
                delete it_to_constraint->second;
                theSimulation->baseSystem->constraints.erase( it_to_constraint );
                constraints.erase( it_constraints );
            }
            else
            {
                indexes.push_back( it_to_constraint->first );
                theSimulation->baseSystem->constraints.erase( it_to_constraint );
            }
        }
        else
            cerr << "ERROR IN CONTACT CONTAINERS!!!" << endl;
    }
    //renaming keys of contacts with negative gaps:
    int i; //index number
    std::vector< std::string >::iterator it_indexes;
    for(i=0; i<constraints.size(); ++i)
    {
        theSimulation->baseSystem->constraints
        ["CONTACT.N"+i] = constraints[i];

    }
    //creating new contact elements
    size_type old_size = constraints.size();
    size_type number_of_triangles = vTriangles.size()/3;
    for(size_type j=0; j<number_of_triangles; ++j)
    {
//      if(vTriangles[3*j] >= 0){
        constraints.push_back(
            new ConstraintContact( theSimulation->nodes[ nodes[vTriangles[3*j]]->getNumber() ]
                                   , theSimulation->nodes[ nodes[vTriangles[3*j+1]]->getNumber() ]
                                   , theSimulation->nodes[ nodes[vTriangles[3*j+2]]->getNumber() ]
                                   , Simulation::alpha
                                   , Simulation::constraintMethod
                                 )
        );
        theSimulation->baseSystem->constraints
        ["CONTACT.N"+j+old_size] = constraints.back();
//      }
    }

}

void Contact::createDelaunay()
{
    delny = vtkDelaunay2D::New();
    delny->SetAlpha(alpha);
    delny->SetInput(polyData);
    delny->SetSource(polyData);
    delny->Update();
    delnyData = vtkPolyData::New();
    delnyData = delny->GetOutput();
//    delnyData->Print(cout);
//    vtkCellArray* cell_array = delnyData->GetPolys();
//    cout << *cell_array << endl;
//    cout << cell_array->GetData()->GetValue(1) << endl;
//    vtkExtractEdges* extract = vtkExtractEdges::New();
//    extract->SetInputConnection(delny->GetOutputPort());
//    extract->Print(cout);
}

void Contact::updateDelaunay()
{
    delny->Delete();
    delny = vtkDelaunay2D::New();
    delny->SetAlpha(alpha);
    delny->SetInput(polyData);
    delny->SetSource(polyData);
    delny->Update();
    delnyData = delny->GetOutput();

    this->readDelaunay();
    if(mapMesh) mapMesh->SetInputConnection(delny->GetOutputPort());
    if(extract) extract->SetInputConnection(delny->GetOutputPort());

}


void Contact::createDrawingObjects()
{
    ren = vtkRenderer::New();
    renWin = vtkRenderWindow::New();
//    iren = vtkRenderWindowInteractor::New();
//    iren->SetRenderWindow(renWin);
    w2if = vtkWindowToImageFilter::New();
    psw = vtkJPEGWriter::New();
    renWin->AddRenderer(ren);
    w2if->SetInput(renWin);
    psw->SetInput(w2if->GetOutput());
//    vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::New();
//    iren->SetRenderWindow(renWin);
    std::vector<CompBar*>::iterator it_bars;
    for( it_bars = bars.begin();
            it_bars!= bars.end();
            ++it_bars )
    {
        (*it_bars)->addToRender( ren );
    }
    ren->SetBackground(0, 0, 0);
    renWin->SetSize(450, 300);
//    renWin->MakeRenderWindowInteractor();

//    iren->Initialize();
//    iren->Start();
}

void Contact::createDrawingContactObjects()
{
    this->readDelaunay();

    mapMesh = vtkPolyDataMapper::New();
    mapMesh->SetInputConnection(delny->GetOutputPort());
    meshActor = vtkActor::New();
    meshActor->SetMapper(mapMesh);

    extract = vtkExtractEdges::New();
    extract->SetInputConnection(delny->GetOutputPort());
    tubes = vtkTubeFilter::New();
    tubes->SetInputConnection(extract->GetOutputPort());
    tubes->SetRadius(0.5);
    tubes->SetNumberOfSides(6);
    mapEdges = vtkPolyDataMapper::New();
    mapEdges->SetInputConnection(tubes->GetOutputPort());
    edgeActor = vtkActor::New();
    edgeActor->SetMapper(mapEdges);
    edgeActor->GetProperty()->SetColor(1,0,0);
    edgeActor->GetProperty()->SetSpecularColor(1, 1, 1);
    edgeActor->GetProperty()->SetSpecular(0.3);
    edgeActor->GetProperty()->SetSpecularPower(20);
    edgeActor->GetProperty()->SetAmbient(0.2);
    edgeActor->GetProperty()->SetDiffuse(0.8);

    ren->AddActor(meshActor);
    ren->AddActor(edgeActor);
}


void Contact::drawObjects()
{
    ren->ResetCamera();
//    ren->GetActiveCamera()->Zoom(2.5);

//    iren->Initialize();
    renWin->Render();
//    iren->Start();

// Commented because it was annoying having lots of images saved in the simulation folder. Waiting for ideas...
    // output image:
//     w2if->Modified();
//     std::ostringstream title;
//     if (image_number <10 )
//       title << "image00" << image_number << ".jpg" ;
//     else if (image_number <100 )
//       title << "image0" << image_number << ".jpg" ;
//     else
//       title << "image" << image_number << ".jpg" ;
//     psw->SetFileName(title.str().c_str());
//     psw->Write();
//
// //    psw->Delete();
//     ++image_number;


}

}

#endif // USE_VTK