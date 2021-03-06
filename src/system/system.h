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

#ifndef MKNIXSYSTEM_H
#define MKNIXSYSTEM_H

#include <map>
#include <string>
#include "LMX/lmx.h"
#include "common.h"

namespace mknix {

class Body;
class RigidBody;
class FlexBody;
class Constraint;
class ConstraintThermal;
class Load;
class LoadThermal;
class Node;
class Point;
class Motion;

/**
  @author AUTHORS <MAILS>
*/
class System {

    friend class Reader;
    friend class ReaderFlex;
    friend class ReaderRigid;
    friend class ReaderConstraints;
    friend class Contact;

public:
    System();

    System(const std::string& title);

    virtual ~System();

    bool outputMaxInterfaceTemp;

    std::string getTitle( )
    {
        return this->title;
    }

    int getNumberOfNodes()
    { return groundNodes.size(); }
    
    Node* getNode( int index)
    { return groundNodes[index]; }
    
    System* getSystem( std::string sysName )
    { return subSystems.at(sysName); }
    
    ConstraintThermal* getConstraintThermal( std::string constraintName )
    { return constraintsThermal.at(constraintName); }
    
    void getThermalNodes( std::vector<double>& );

    void getOutputSignalThermal(double* );

    void updateThermalLoads( double* );

    virtual void update( double );
    
    void initFlexBodies();

    void writeRigidBodies( std::ofstream* );

    void writeFlexBodies( std::ofstream* );

    void writeJoints( std::ofstream* );
  
    void calcCapacityMatrix( );

    void calcConductivityMatrix( );

    void calcInternalHeat( );

    void calcExternalHeat( );

    void calcThermalTangentMatrix( );

    void assembleCapacityMatrix( lmx::Matrix<data_type>& );

    void assembleConductivityMatrix( lmx::Matrix<data_type>& );

    void assembleExternalHeat( lmx::Vector<data_type>& );

    void assembleInternalHeat( lmx::Vector<data_type>& );

    void assembleThermalTangentMatrix( lmx::Matrix<data_type>& );

    void calcMassMatrix( );

    void calcInternalForces( );

    void calcExternalForces( );

    void calcTangentMatrix( );

    void assembleMassMatrix( lmx::Matrix<data_type>& );

    void assembleInternalForces( lmx::Vector<data_type>& );

    void assembleExternalForces( lmx::Vector<data_type>& );

    void assembleTangentMatrix( lmx::Matrix<data_type>& );

    void assembleConstraintForces( lmx::Vector<data_type>& );

    void setMechanical( );

    void outputStep( const lmx::Vector<data_type>&, const lmx::Vector<data_type>& );

    void outputStep( const lmx::Vector<data_type>& );

    void outputToFile( std::ofstream* );

    bool checkAugmented();

    void clearAugmented();

    void writeBoundaryNodes( std::vector<Point*>& );

    void writeBoundaryConnectivity( std::vector< std::vector<Point*> >& );

protected:
    std::string title;
    std::vector<Node*> groundNodes;

    std::map< std::string, System* > subSystems;
    std::map< std::string, RigidBody* > rigidBodies;
    std::map< std::string, FlexBody* > flexBodies;
    std::map< std::string, Body* > thermalBodies;
    std::map< std::string, Constraint* > constraints;
    std::map< std::string, ConstraintThermal* > constraintsThermal;
    std::vector< Load* > loads;
    std::vector< LoadThermal* > loadsThermal;
    std::vector< Node* > outputSignalThermal;
    std::vector< Motion* > motions;

};

}

#endif
