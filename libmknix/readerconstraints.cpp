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
#include "readerconstraints.h"
#include "bodyflex.h"
#include "bodyrigid.h"
#include "constraint.h"
#include "constraintclearance.h"
#include "constraintcontact.h"
#include "constraintdistance.h"
#include "constraintfixedaxis.h"
#include "constraintfixedcoordinates.h"
#include "constraintthermalfixed.h"
#include "node.h"
#include "simulation.h"
#include "system.h"

namespace mknix {

ReaderConstraints::ReaderConstraints()
    : theSimulation(0)
    , output(0)
    , input(0)
    , p_nodeA(0)
    , p_nodeB(0)
{
}

ReaderConstraints::ReaderConstraints( Simulation* simulation_in,
                                      std::ofstream & output_in,
                                      std::ifstream & input_in
                                    )
    : theSimulation( simulation_in )
    , output(& output_in)
    , input(& input_in)
    , p_nodeA(0)
    , p_nodeB(0)
{
}


ReaderConstraints::~ReaderConstraints()
{
}


} // namespace mknix

void mknix::ReaderConstraints::readConstraints(System * system_in)
{
    char keyword[20];
    /*consTitle[30],*/
    std::string consTitle;

    while(*input >> keyword)
    {
        if(!strcmp( keyword,"ENDJOINTS") )
        {
            return;
        }
        else if( !strcmp(keyword,"PENALTY") )
        {
            Simulation::constraintMethod = "PENALTY";
            *output << "PENALTY set" << endl;
        }
        else if( !strcmp(keyword,"AUGMENTED") )
        {
            Simulation::constraintMethod = "AUGMENTED";
            *output << "AUGMENTED set" << endl;
        }
        else if( !strcmp(keyword,"ALPHA") )
        {
            *input  >> Simulation::alpha;
            *output << "ALPHA: "
                    << Simulation::getAlpha()
                    << endl;
        }
        else if(!strcmp(keyword,"SPHERICAL"))
        {
            /* Igual a una restriccion de distancia constante */
            std::string bodyTitleA, bodyTitleB;
            std::string nodeA, nodeB;

            p_nodeA=0;
            p_nodeB=0;

            *input >> consTitle;
            *output << "SPHERICAL: "
                    << system_in->getTitle()
                    << "."
                    << consTitle << std::endl;

            while( *input >> keyword )
            {
                if( !strcmp( keyword,"ENDSPHERICAL") )
                {
                    break;
                }
                else if( !strcmp(keyword,"NODEA") )
                    this->readNodeName( bodyTitleA, nodeA );
                else if( !strcmp(keyword,"NODEB") )
                    this->readNodeName( bodyTitleB, nodeB );
            }
            this->assignConstraintNodes( system_in,
                                         bodyTitleA,
                                         nodeA,
                                         bodyTitleB,
                                         nodeB
                                       );
            system_in->constraints[consTitle]
                = new ConstraintFixedCoordinates( p_nodeA,
                                                  p_nodeB,
                                                  Simulation::alpha,
                                                  Simulation::constraintMethod
                                                );
	    system_in->constraints[consTitle]->setTitle(consTitle);

            this->outputConstraintNode( system_in, consTitle, "NODEA",
                                        bodyTitleA, nodeA, 0 );
            this->outputConstraintNode( system_in, consTitle, "NODEB",
                                        bodyTitleB, nodeB, 1 );
        }

        else if(!strcmp(keyword,"DISTANCE"))
        {
            /* Igual a una restriccion de distancia constante */
            std::string bodyTitleA, bodyTitleB;
            std::string nodeA, nodeB;

            p_nodeA=0;
            p_nodeB=0;

            *input >> consTitle;
            *output << "DISTANCE: "
                    << system_in->getTitle()
                    << "."
                    << consTitle << std::endl;

            while( *input >> keyword )
            {
                if( !strcmp( keyword,"ENDDISTANCE") )
                {
                    break;
                }
                else if( !strcmp(keyword,"NODEA") )
                    this->readNodeName( bodyTitleA, nodeA );
                else if( !strcmp(keyword,"NODEB") )
                    this->readNodeName( bodyTitleB, nodeB );
            }
            this->assignConstraintNodes( system_in,
                                         bodyTitleA,
                                         nodeA,
                                         bodyTitleB,
                                         nodeB
                                       );
            system_in->constraints[consTitle]
                = new ConstraintDistance( p_nodeA,
					  p_nodeB,
					  Simulation::alpha,
					  Simulation::constraintMethod
					 );
	    system_in->constraints[consTitle]->setTitle(consTitle);

            this->outputConstraintNode( system_in, consTitle, "NODEA",
                                        bodyTitleA, nodeA, 0 );
            this->outputConstraintNode( system_in, consTitle, "NODEB",
                                        bodyTitleB, nodeB, 1 );
        }

        else if(!strcmp(keyword,"AXIS")) {
            // Igual a una restricción de distancia constante
            std::string axisName;
            std::string bodyTitleA, bodyTitleB;
            std::string nodeA, nodeB;

            p_nodeA=0;
            p_nodeB=0;

            *input >> consTitle;
            *output << "AXIS: "
                    << system_in->getTitle()
                    << "."
                    << consTitle << std::endl;
            while(*input >> keyword) {
                if(!strcmp(keyword,"ENDAXIS")) break;
                else if(!strcmp(keyword,"DIRECTION")) {
                    *input >> axisName;
                    *output << "DIRECTION: " << axisName << "-axis" << endl;
                }
                else if( !strcmp(keyword,"NODEA") )
                    this->readNodeName( bodyTitleA, nodeA );
                else if( !strcmp(keyword,"NODEB") )
                    this->readNodeName( bodyTitleB, nodeB );
            }
            this->assignConstraintNodes( system_in,
                                         bodyTitleA,
                                         nodeA,
                                         bodyTitleB,
                                         nodeB
                                       );
            system_in->constraints[consTitle]
                = new ConstraintFixedAxis( p_nodeA,
                                           p_nodeB,
                                           axisName,
                                           Simulation::alpha,
                                           Simulation::constraintMethod
                                         );
	    system_in->constraints[consTitle]->setTitle(consTitle);

            this->outputConstraintNode( system_in, consTitle, "NODEA",
                                        bodyTitleA, nodeA, 0 );
            this->outputConstraintNode( system_in, consTitle, "NODEB",
                                        bodyTitleB, nodeB, 1 );
        }
        else if(!strcmp(keyword,"CLEARANCE")) {
            // Igual a una restricción de distancia constante
            std::string bodyTitleA, bodyTitleB;
            std::string nodeA, nodeB;
            double tolerance=0.;

            p_nodeA=0;
            p_nodeB=0;

            *input >> consTitle;
            *output << "CLEARANCE: "
                    << system_in->getTitle()
                    << "."
                    << consTitle << std::endl;
            while(*input >> keyword) {
                if(!strcmp(keyword,"ENDCLEARANCE")) break;
                else if(!strcmp(keyword,"TOLERANCE")) {
                    *input >> tolerance;
                    *output << "TOLERANCE: " << tolerance << endl;
                }
                else if( !strcmp(keyword,"NODEA") )
                    this->readNodeName( bodyTitleA, nodeA );
                else if( !strcmp(keyword,"NODEB") )
                    this->readNodeName( bodyTitleB, nodeB );
            }
            this->assignConstraintNodes( system_in,
                                         bodyTitleA,
                                         nodeA,
                                         bodyTitleB,
                                         nodeB
                                       );
            system_in->constraints[consTitle]
                = new ConstraintClearance( p_nodeA,
                                           p_nodeB,
                                           tolerance,
                                           Simulation::alpha,
                                           Simulation::constraintMethod
                                         );
	    system_in->constraints[consTitle]->setTitle(consTitle);

            this->outputConstraintNode( system_in, consTitle, "NODEA",
                                        bodyTitleA, nodeA, 0 );
            this->outputConstraintNode( system_in, consTitle, "NODEB",
                                        bodyTitleB, nodeB, 1 );
        }
        else if(!strcmp(keyword,"THERMALSPHERICAL"))
        {
            /* Igual a una restriccion de distancia constante */
            std::string bodyTitleA, bodyTitleB;
            std::string nodeA, nodeB;

            p_nodeA=0;
            p_nodeB=0;

            *input >> consTitle;
            *output << "THERMALSPHERICAL: "
                    << system_in->getTitle()
                    << "."
                    << consTitle << std::endl;

            while( *input >> keyword )
            {
                if( !strcmp( keyword,"ENDTHERMALSPHERICAL") )
                {
                    break;
                }
                else if( !strcmp(keyword,"NODEA") )
                    this->readNodeName( bodyTitleA, nodeA );
                else if( !strcmp(keyword,"NODEB") )
                    this->readNodeName( bodyTitleB, nodeB );
            }
            this->assignConstraintNodes( system_in,
                                         bodyTitleA,
                                         nodeA,
                                         bodyTitleB,
                                         nodeB
                                       );
            system_in->constraintsThermal[consTitle]
                = new ConstraintThermalFixed( p_nodeA,
                                                  p_nodeB,
                                                  Simulation::alpha,
                                                  Simulation::constraintMethod
                                                );
	    system_in->constraintsThermal[consTitle]->setTitle(consTitle);

            this->outputConstraintThermalNode( system_in, consTitle, "NODEA",
                                        bodyTitleA, nodeA, 0 );
            this->outputConstraintThermalNode( system_in, consTitle, "NODEB",
                                        bodyTitleB, nodeB, 1 );
        }

        else if(!strcmp(keyword,"OTRO")) {
        }
    }
}


void mknix::ReaderConstraints::readNodeName( std::string & bodyTitle, std::string & node)
{
    char a;

    input->get(a); // get blank space
    while(input->get(a))
    {
        if (a == '.') break;
        else if (a=='\n') break;
        else
        {
            bodyTitle.push_back( a );
        }
    }
    if ( a!='\n' )
    {
        while(input->get(a))
        {
            if (a == '.') break;
            else if (a=='\n') break;
            else
            {
                node.push_back( a );;
            }
        }
    }
    cout << "NODE read: " << bodyTitle << "." << node << endl;
}


void mknix::ReaderConstraints::assignConstraintNodes( System* system_in,
        std::string & bodyTitleA,
        std::string & nodeA,
        std::string & bodyTitleB,
        std::string & nodeB
                                                    )
{
    if ( bodyTitleA=="GROUND" )
    {
        /* If the body is a system */
        if( system_in->subSystems.find(bodyTitleB) !=
                system_in->subSystems.end() )
            system_in->groundNodes.push_back
            (  new Node( *system_in->subSystems[bodyTitleB]
                         ->getNode( atoi(nodeB.c_str()) )
                       )
            );
        else {
            if( system_in->rigidBodies.find(bodyTitleB) !=
                    system_in->rigidBodies.end() )
                system_in->groundNodes.push_back
                (  new Node( *system_in->rigidBodies[bodyTitleB]
                             //                         ->getDomainNode( nodeB )
                             ->getNode( atoi(nodeB.c_str()) )
                           )
                );
            /* The body is a flexbody*/
            else
                system_in->groundNodes.push_back
                ( new Node( *system_in->flexBodies[bodyTitleB]
                            ->getNode( atoi( nodeB.c_str() ) )
                          )
                );
        }
        p_nodeA = system_in->groundNodes.back();
        p_nodeA->setNumber(-1);
        p_nodeA->setThermalNumber(-1);
    }
    if ( bodyTitleB=="GROUND" )
    {
        if( system_in->subSystems.find(bodyTitleA) !=
                system_in->subSystems.end() )
            system_in->groundNodes.push_back
            (  new Node( *system_in->subSystems[bodyTitleA]
                         ->getNode( atoi(nodeA.c_str()) )
                       )
            );
        else {
            /* If the body is a rigidbody */
            if( system_in->rigidBodies.find(bodyTitleA) !=
                    system_in->rigidBodies.end() )
            {
                system_in->groundNodes.push_back
                ( new Node( system_in->rigidBodies[bodyTitleA]
                            //                         ->getDomainNode(  nodeA )
                            ->getNode( atoi(nodeA.c_str()) )
                          )
                );
            }
            /* the body is a flexbody */
            else
                system_in->groundNodes.push_back
                ( new Node( *system_in->flexBodies[bodyTitleA]
                            ->getNode( atoi( nodeA.c_str() ) )
                          )
                );
        }
        p_nodeB = system_in->groundNodes.back();
        p_nodeB->setNumber(-1);
        p_nodeB->setThermalNumber(-1);
//         cout << system_in->groundNodes.size()
//              << "(i)x="
//              << *system_in
//                  ->groundNodes[system_in->groundNodes.size()-1]->getx()
//              << endl;
    }
    /* if it's not grounded */
    if( p_nodeA == 0 )
    {
        if( system_in->subSystems.find(bodyTitleA) !=
                system_in->subSystems.end() )
            p_nodeA = system_in->subSystems[bodyTitleA]
                      ->getNode( atoi(nodeA.c_str()) );
        else {
            if( system_in->rigidBodies.find(bodyTitleA) !=
                    system_in->rigidBodies.end() ) //if the body is a rigidbody
                p_nodeA = system_in->rigidBodies[bodyTitleA]
                          // 		  ->getDomainNode( nodeA );
                          ->getNode( atoi(nodeA.c_str()) ) ;
            else { //the body is a flexbody
                p_nodeA = system_in->flexBodies[bodyTitleA]
                          ->getNode( atoi( nodeA.c_str() ) );
            }
        }
    }
    /* if it's not grounded */
    if( p_nodeB == 0)
    {
        if( system_in->subSystems.find(bodyTitleB) !=
                system_in->subSystems.end() )
            p_nodeB = system_in->subSystems[bodyTitleB]
                      ->getNode( atoi(nodeB.c_str()) );
        else {
            if( system_in->rigidBodies.find(bodyTitleB) !=
                    system_in->rigidBodies.end() ) //if the body is a rigidbody
                p_nodeB = system_in->rigidBodies[bodyTitleB]
                          // 		  ->getDomainNode( nodeB );
                          ->getNode( atoi(nodeB.c_str()) );
            else //the body is a flexbody
                p_nodeB = system_in->flexBodies[bodyTitleB]
                          ->getNode( atoi( nodeB.c_str() ) );
        }
    }
}


void mknix::ReaderConstraints::outputConstraintNode( System* system_in,
        std::string & consTitle,
        const char * nodeName,
        std::string & bodyTitle,
        std::string & node,
        int i
                                                   )
{
    *output << "\t"
            << consTitle << "." << nodeName << ": "
            << bodyTitle << "." << node
            << ": " << system_in->constraints[consTitle]->getNode(i)->getNumber()
            << ": " << system_in->constraints[consTitle]->getNode(i)->getX()
            << ", " << system_in->constraints[consTitle]->getNode(i)->getY()
            << ", " << system_in->constraints[consTitle]->getNode(i)->getZ()
            << std::endl;

}


void mknix::ReaderConstraints::outputConstraintThermalNode( System* system_in,
        std::string & consTitle,
        const char * nodeName,
        std::string & bodyTitle,
        std::string & node,
        int i
                                                   )
{
    *output << "\t"
            << consTitle << "." << nodeName << ": "
            << bodyTitle << "." << node
            << ": " << system_in->constraintsThermal[consTitle]->getNode(i)->getNumber()
            << ": " << system_in->constraintsThermal[consTitle]->getNode(i)->getX()
            << ", " << system_in->constraintsThermal[consTitle]->getNode(i)->getY()
            << ", " << system_in->constraintsThermal[consTitle]->getNode(i)->getZ()
            << std::endl;

}

