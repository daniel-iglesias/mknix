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
#include "system.h"

#include "body.h"
#include "bodyflex.h"
#include "bodyrigid.h"
#include "constraint.h"
#include "constraintthermal.h"
#include "load.h"
#include "loadthermal.h"
#include "motion.h"
#include "node.h"


namespace mknix {

System::System()
  : outputMaxInterfaceTemp(false)
{
}


System::System(const char * title_in)
  : title(title_in)
  , outputMaxInterfaceTemp(false)
{
}


System::~System()
{
    for ( itSubSystems = subSystems.begin();
            itSubSystems!= subSystems.end();
            ++itSubSystems
        )
    {
        delete(itSubSystems->second);
    }
    for ( itRigidBodies = rigidBodies.begin();
            itRigidBodies!= rigidBodies.end();
            ++itRigidBodies
        )
    {
        delete(itRigidBodies->second);
    }
    for ( itFlexBodies = flexBodies.begin();
            itFlexBodies!= flexBodies.end();
            ++itFlexBodies
        )
    {
        delete(itFlexBodies->second);
    }
//     for ( itConstraints = constraints.begin();
//             itConstraints!= constraints.end();
//             ++itConstraints
//         )
//     {
//         delete(itConstraints->second);
//     }
//     for ( itConstraintsThermal = constraintsThermal.begin();
//             itConstraintsThermal!= constraintsThermal.end();
//             ++itConstraintsThermal
//         )
//     {
//         delete(itConstraintsThermal->second);
//     }
    for ( itLoads = loads.begin();
            itLoads!= loads.end();
            ++itLoads
        )
    {
        delete(*itLoads);
    }
    for ( itLoadsThermal = loadsThermal.begin();
            itLoadsThermal!= loadsThermal.end();
            ++itLoadsThermal
        )
    {
        delete(*itLoadsThermal);
    }
}

// BUG: Specific for tiles susbsystem. That can change in input
void System::getThermalNodes(std::vector< double >& x_coordinates)
{
  for ( itLoadsThermal = subSystems["tiles"]->loadsThermal.begin();
        itLoadsThermal!= subSystems["tiles"]->loadsThermal.end();
        ++itLoadsThermal
      )
  {
    (*itLoadsThermal)->insertNodesXCoordinates( x_coordinates );
  }
}


void System::getOutputSignalThermal(double* vector_in)
{
  int counter=0;
  for ( itOutputSignalThermal = subSystems["tiles"]->outputSignalThermal.begin();
        itOutputSignalThermal!= subSystems["tiles"]->outputSignalThermal.end();
        ++itOutputSignalThermal
      )
  {
    vector_in[counter] = (*itOutputSignalThermal)->getTemp();
    ++counter;
  }
  
  if(subSystems["tiles"]->outputMaxInterfaceTemp){
    vector_in[counter] = 0;
    for ( itLoadsThermal = subSystems["tiles"]->loadsThermal.begin();
        itLoadsThermal!= subSystems["tiles"]->loadsThermal.end();
        ++itLoadsThermal
      )
    {
      (*itLoadsThermal)->getMaxTemp( vector_in[counter] );
    }
  }
}


void System::updateThermalLoads(double* vector_in)
{
  int counter=0;
  for ( itLoadsThermal = subSystems["tiles"]->loadsThermal.begin();
        itLoadsThermal!= subSystems["tiles"]->loadsThermal.end();
        ++itLoadsThermal
      )
  {
    (*itLoadsThermal)->updateLoad(vector_in[counter]);
    ++counter;
  }
}

void System::update(double time)
{
    for ( itMotions = motions.begin();
            itMotions!= motions.end();
            ++itMotions
        )
    {
        (*itMotions)->update( time );
    }
    for ( itSubSystems = subSystems.begin();
            itSubSystems!= subSystems.end();
            ++itSubSystems
        )
    {
        itSubSystems->second->update( time );
    }
}


void System::initFlexBodies()
{
    for ( itFlexBodies = flexBodies.begin();
            itFlexBodies!= flexBodies.end();
            ++itFlexBodies
        )
    {
        itFlexBodies->second->initialize( );
    }
}


void System::writeRigidBodies( std::ofstream* outFile )
{
    for ( itRigidBodies = rigidBodies.begin();
            itRigidBodies!= rigidBodies.end();
            ++itRigidBodies
        )
    {
        itRigidBodies->second->writeBodyInfo( outFile );
    }

    for ( itSubSystems = subSystems.begin();
            itSubSystems!= subSystems.end();
            ++itSubSystems
        )
    {
        itSubSystems->second->writeRigidBodies( outFile );
    }


}

void System::writeFlexBodies( std::ofstream* outFile )
{
    for ( itFlexBodies = flexBodies.begin();
            itFlexBodies!= flexBodies.end();
            ++itFlexBodies
        )
    {
        itFlexBodies->second->writeBodyInfo( outFile );
    }

    for ( itSubSystems = subSystems.begin();
            itSubSystems!= subSystems.end();
            ++itSubSystems
        )
    {
        itSubSystems->second->writeFlexBodies( outFile );
    }
}

  void System::writeJoints( std::ofstream* outFile )
  {
    for ( itConstraints = constraints.begin();
         itConstraints!= constraints.end();
         ++itConstraints
         )
    {
      itConstraints->second->writeJointInfo( outFile );
    }
    
    for ( itSubSystems = subSystems.begin();
         itSubSystems!= subSystems.end();
         ++itSubSystems
         )
    {
      itSubSystems->second->writeJoints( outFile );
    }
  }

} // Namespace mknix


void mknix::System::calcCapacityMatrix( )
{
    for ( itThermalBodies = thermalBodies.begin();
            itThermalBodies!= thermalBodies.end();
            ++itThermalBodies
        )
    {
        itThermalBodies->second->calcCapacityMatrix();
    }

    for ( itSubSystems = subSystems.begin();
            itSubSystems!= subSystems.end();
            ++itSubSystems
        )
    {
        itSubSystems->second->calcCapacityMatrix();
    }
}

void mknix::System::calcConductivityMatrix( )
{
    for ( itThermalBodies = thermalBodies.begin();
            itThermalBodies!= thermalBodies.end();
            ++itThermalBodies
        )
    {
        itThermalBodies->second->calcConductivityMatrix();
    }

//   for ( itConstraints = constraints.begin();
//         itConstraints!= constraints.end();
//         ++itConstraints
//       )
//   {
//     itConstraints->second->calcConductivityMatrix();
//   }

    for ( itSubSystems = subSystems.begin();
            itSubSystems!= subSystems.end();
            ++itSubSystems
        )
    {
        itSubSystems->second->calcConductivityMatrix();
    }
}

void mknix::System::calcExternalHeat( )
{
    for ( itThermalBodies = thermalBodies.begin();
            itThermalBodies!= thermalBodies.end();
            ++itThermalBodies
        )
    {
        itThermalBodies->second->calcExternalHeat();
    }

    for ( itSubSystems = subSystems.begin();
            itSubSystems!= subSystems.end();
            ++itSubSystems
        )
    {
        itSubSystems->second->calcExternalHeat();
    }
}

void mknix::System::calcInternalHeat( )
{
    for ( itConstraintsThermal = constraintsThermal.begin();
            itConstraintsThermal!= constraintsThermal.end();
            ++itConstraintsThermal
        )
    {
        itConstraintsThermal->second->calcInternalForces();
    }


    for ( itSubSystems = subSystems.begin();
            itSubSystems!= subSystems.end();
            ++itSubSystems
        )
    {
        itSubSystems->second->calcInternalHeat();
    }
}

void mknix::System::calcThermalTangentMatrix()
{
    for ( itConstraintsThermal = constraintsThermal.begin();
            itConstraintsThermal!= constraintsThermal.end();
            ++itConstraintsThermal
        )
    {
        itConstraintsThermal->second->calcTangentMatrix();
    }

    for ( itSubSystems = subSystems.begin();
            itSubSystems!= subSystems.end();
            ++itSubSystems
        )
    {
        itSubSystems->second->calcThermalTangentMatrix();
    }

}

void mknix::System::assembleCapacityMatrix( lmx::Matrix<data_type>& globalCapacity_in)
{
    for ( itThermalBodies = thermalBodies.begin();
            itThermalBodies!= thermalBodies.end();
            ++itThermalBodies
        )
    {
        itThermalBodies->second->assembleCapacityMatrix(globalCapacity_in);
    }

    for ( itSubSystems = subSystems.begin();
            itSubSystems!= subSystems.end();
            ++itSubSystems
        )
    {
        itSubSystems->second->assembleCapacityMatrix(globalCapacity_in);
    }

}

void mknix::System::assembleConductivityMatrix( lmx::Matrix<data_type>& globalConductivity_in)
{
    for ( itThermalBodies = thermalBodies.begin();
            itThermalBodies!= thermalBodies.end();
            ++itThermalBodies
        )
    {
        itThermalBodies->second->assembleConductivityMatrix(globalConductivity_in);
    }

//   for ( itConstraints = constraints.begin();
//         itConstraints!= constraints.end();
//         ++itConstraints
//       )
//   {
//     itConstraints->second->assembleConductivityMatrix(globalConductivity_in);
//   }

    for ( itSubSystems = subSystems.begin();
            itSubSystems!= subSystems.end();
            ++itSubSystems
        )
    {
        itSubSystems->second->assembleConductivityMatrix(globalConductivity_in);
    }

}

void mknix::System::assembleExternalHeat( lmx::Vector<data_type>& externalHeat_in)
{
    for ( itThermalBodies = thermalBodies.begin();
            itThermalBodies!= thermalBodies.end();
            ++itThermalBodies
        )
    {
        itThermalBodies->second->assembleExternalHeat(externalHeat_in);
    }

  for ( itLoadsThermal = loadsThermal.begin();
        itLoadsThermal!= loadsThermal.end();
        ++itLoadsThermal
      )
  {
    (*itLoadsThermal)->assembleExternalHeat(externalHeat_in);
  }

//  cout << "External heat in System (1) = " << externalHeat_in;
    for ( itSubSystems = subSystems.begin();
            itSubSystems!= subSystems.end();
            ++itSubSystems
        )
    {
        itSubSystems->second->assembleExternalHeat(externalHeat_in);
    }
}


void mknix::System::assembleInternalHeat( lmx::Vector<data_type>& internalHeat_in)
{
    for ( itSubSystems = subSystems.begin();
            itSubSystems!= subSystems.end();
            ++itSubSystems
        )
    {
        itSubSystems->second->assembleInternalHeat(internalHeat_in);
    }
    for ( itConstraintsThermal = constraintsThermal.begin();
            itConstraintsThermal!= constraintsThermal.end();
            ++itConstraintsThermal
        )
    {
        itConstraintsThermal->second->assembleInternalForces(internalHeat_in);
    }

}


void mknix::System::assembleThermalTangentMatrix(lmx::Matrix< data_type > & globalTangent_in)
{
    for ( itConstraintsThermal = constraintsThermal.begin();
            itConstraintsThermal!= constraintsThermal.end();
            ++itConstraintsThermal
        )
    {
        itConstraintsThermal->second->assembleTangentMatrix(globalTangent_in);
    }

    for ( itSubSystems = subSystems.begin();
            itSubSystems!= subSystems.end();
            ++itSubSystems
        )
    {
        itSubSystems->second->assembleThermalTangentMatrix(globalTangent_in);
    }
}


void mknix::System::calcMassMatrix()
{
    for ( itFlexBodies = flexBodies.begin();
            itFlexBodies!= flexBodies.end();
            ++itFlexBodies
        )
    {
        itFlexBodies->second->calcMassMatrix();
    }

    for ( itRigidBodies = rigidBodies.begin();
            itRigidBodies!= rigidBodies.end();
            ++itRigidBodies
        )
    {
        itRigidBodies->second->calcMassMatrix();
    }

    for ( itSubSystems = subSystems.begin();
            itSubSystems!= subSystems.end();
            ++itSubSystems
        )
    {
        itSubSystems->second->calcMassMatrix();
    }
}

void mknix::System::calcInternalForces()
{
    for ( itFlexBodies = flexBodies.begin();
            itFlexBodies!= flexBodies.end();
            ++itFlexBodies
        )
    {
        itFlexBodies->second->calcInternalForces();
    }

    for ( itConstraints = constraints.begin();
            itConstraints!= constraints.end();
            ++itConstraints
        )
    {
        itConstraints->second->calcInternalForces();
    }

    for ( itSubSystems = subSystems.begin();
            itSubSystems!= subSystems.end();
            ++itSubSystems
        )
    {
        itSubSystems->second->calcInternalForces();
    }

}

void mknix::System::calcExternalForces()
{
    for ( itFlexBodies = flexBodies.begin();
            itFlexBodies!= flexBodies.end();
            ++itFlexBodies
        )
    {
        itFlexBodies->second->calcExternalForces();
    }

    for ( itRigidBodies = rigidBodies.begin();
            itRigidBodies!= rigidBodies.end();
            ++itRigidBodies
        )
    {
        itRigidBodies->second->calcExternalForces();
    }

    for ( itSubSystems = subSystems.begin();
            itSubSystems!= subSystems.end();
            ++itSubSystems
        )
    {
        itSubSystems->second->calcExternalForces();
    }

}

void mknix::System::calcTangentMatrix()
{
    for ( itFlexBodies = flexBodies.begin();
            itFlexBodies!= flexBodies.end();
            ++itFlexBodies
        )
    {
        itFlexBodies->second->calcTangentMatrix();
    }


    for ( itConstraints = constraints.begin();
            itConstraints!= constraints.end();
            ++itConstraints
        )
    {
        itConstraints->second->calcTangentMatrix();
    }

    for ( itSubSystems = subSystems.begin();
            itSubSystems!= subSystems.end();
            ++itSubSystems
        )
    {
        itSubSystems->second->calcTangentMatrix();
    }

}

void mknix::System::assembleMassMatrix(lmx::Matrix< data_type > & globalMass_in)
{
    for ( itFlexBodies = flexBodies.begin();
            itFlexBodies!= flexBodies.end();
            ++itFlexBodies
        )
    {
        itFlexBodies->second->assembleMassMatrix(globalMass_in);
    }

    for ( itRigidBodies = rigidBodies.begin();
            itRigidBodies!= rigidBodies.end();
            ++itRigidBodies
        )
    {
        itRigidBodies->second->assembleMassMatrix(globalMass_in);
    }

    for ( itSubSystems = subSystems.begin();
            itSubSystems!= subSystems.end();
            ++itSubSystems
        )
    {
        itSubSystems->second->assembleMassMatrix(globalMass_in);
    }

}

void mknix::System::assembleInternalForces(lmx::Vector< data_type > & internalForces_in)
{
    for ( itFlexBodies = flexBodies.begin();
            itFlexBodies!= flexBodies.end();
            ++itFlexBodies
        )
    {
        itFlexBodies->second->assembleInternalForces(internalForces_in);
    }

    for ( itConstraints = constraints.begin();
            itConstraints!= constraints.end();
            ++itConstraints
        )
    {
        itConstraints->second->assembleInternalForces(internalForces_in);
    }

    for ( itSubSystems = subSystems.begin();
            itSubSystems!= subSystems.end();
            ++itSubSystems
        )
    {
        itSubSystems->second->assembleInternalForces(internalForces_in);
    }

}

void mknix::System::assembleExternalForces(lmx::Vector< data_type > & externalForces_in)
{

    for ( itFlexBodies = flexBodies.begin();
            itFlexBodies!= flexBodies.end();
            ++itFlexBodies
        )
    {
        itFlexBodies->second->assembleExternalForces(externalForces_in);
    }

    for ( itRigidBodies = rigidBodies.begin();
            itRigidBodies!= rigidBodies.end();
            ++itRigidBodies
        )
    {
        itRigidBodies->second->assembleExternalForces(externalForces_in);
    }

    for ( itLoads = loads.begin();
            itLoads!= loads.end();
            ++itLoads
        )
    {
        (*itLoads)->assembleExternalForces(externalForces_in);
    }

//  cout << "External in System (1) = " << externalForces_in;
    for ( itSubSystems = subSystems.begin();
            itSubSystems!= subSystems.end();
            ++itSubSystems
        )
    {
        itSubSystems->second->assembleExternalForces(externalForces_in);
    }

//  cout << "External in System (2) = " << externalForces_in;
}

void mknix::System::assembleTangentMatrix(lmx::Matrix< data_type > & globalTangent_in)
{
    for ( itFlexBodies = flexBodies.begin();
            itFlexBodies!= flexBodies.end();
            ++itFlexBodies
        )
    {
        itFlexBodies->second->assembleTangentMatrix(globalTangent_in);
    }

    for ( itConstraints = constraints.begin();
            itConstraints!= constraints.end();
            ++itConstraints
        )
    {
        itConstraints->second->assembleTangentMatrix(globalTangent_in);
    }

    for ( itSubSystems = subSystems.begin();
            itSubSystems!= subSystems.end();
            ++itSubSystems
        )
    {
        itSubSystems->second->assembleTangentMatrix(globalTangent_in);
    }
}


void mknix::System::assembleConstraintForces(lmx::Vector< data_type > & internalForces_in)
{
    for ( itConstraints = constraints.begin();
            itConstraints!= constraints.end();
            ++itConstraints
        )
    {
        itConstraints->second->calcInternalForces();
        itConstraints->second->assembleInternalForces(internalForces_in);
    }

//   for ( itFlexBodies = flexBodies.begin();
//         itFlexBodies!= flexBodies.end();
//         ++itFlexBodies
//       )
//   {
//     itFlexBodies->second->calcInternalForces();
//     itFlexBodies->second->assembleInternalForces(internalForces_in);
//   }

    for ( itSubSystems = subSystems.begin();
            itSubSystems!= subSystems.end();
            ++itSubSystems
        )
    {
        itSubSystems->second->assembleConstraintForces(internalForces_in);
    }
}


void mknix::System::setMechanical(  )
{
    for ( itRigidBodies = rigidBodies.begin();
            itRigidBodies!= rigidBodies.end();
            ++itRigidBodies
        )
    {
        itRigidBodies->second->setMechanical( );
    }

    for ( itFlexBodies = flexBodies.begin();
            itFlexBodies!= flexBodies.end();
            ++itFlexBodies
        )
    {
        itFlexBodies->second->setMechanical( );
    }

    for ( itSubSystems = subSystems.begin();
            itSubSystems!= subSystems.end();
            ++itSubSystems
        )
    {
        itSubSystems->second->setMechanical( );
    }
}

void mknix::System::outputStep( const lmx::Vector<data_type>& q, const lmx::Vector<data_type>& qdot )
{
    for ( itRigidBodies = rigidBodies.begin();
            itRigidBodies!= rigidBodies.end();
            ++itRigidBodies
        )
    {
        itRigidBodies->second->outputStep( q, qdot );
    }

    for ( itFlexBodies = flexBodies.begin();
            itFlexBodies!= flexBodies.end();
            ++itFlexBodies
        )
    {
        itFlexBodies->second->outputStep( q, qdot );
    }

    for ( itConstraints = constraints.begin();
            itConstraints!= constraints.end();
            ++itConstraints
        )
    {
        itConstraints->second->outputStep(q, qdot);
    }

    for ( itSubSystems = subSystems.begin();
            itSubSystems!= subSystems.end();
            ++itSubSystems
        )
    {
        itSubSystems->second->outputStep( q, qdot );
    }
}


void mknix::System::outputStep( const lmx::Vector<data_type>& q )
{
    for ( itRigidBodies = rigidBodies.begin();
            itRigidBodies!= rigidBodies.end();
            ++itRigidBodies
        )
    {
        itRigidBodies->second->outputStep( q );
    }

    for ( itFlexBodies = flexBodies.begin();
            itFlexBodies!= flexBodies.end();
            ++itFlexBodies
        )
    {
        itFlexBodies->second->outputStep( q );
    }

    for ( itConstraints = constraints.begin();
            itConstraints!= constraints.end();
            ++itConstraints
        )
    {
        itConstraints->second->outputStep( q );
    }

    for ( itSubSystems = subSystems.begin();
            itSubSystems!= subSystems.end();
            ++itSubSystems
        )
    {
        itSubSystems->second->outputStep( q );
    }
}


void mknix::System::outputToFile(std::ofstream * outFile)
{
    for ( itRigidBodies = rigidBodies.begin();
            itRigidBodies!= rigidBodies.end();
            ++itRigidBodies
        )
    {
        itRigidBodies->second->outputToFile( outFile );
    }

    for ( itFlexBodies = flexBodies.begin();
            itFlexBodies!= flexBodies.end();
            ++itFlexBodies
        )
    {
        itFlexBodies->second->outputToFile( outFile );
    }

    for ( itLoads = loads.begin();
            itLoads!= loads.end();
            ++itLoads
        )
    {
        (*itLoads)->outputToFile( outFile );
    }

    for ( itConstraints = constraints.begin();
            itConstraints!= constraints.end();
            ++itConstraints
        )
    {
        itConstraints->second->outputToFile( outFile );
    }

    for ( itSubSystems = subSystems.begin();
            itSubSystems!= subSystems.end();
            ++itSubSystems
        )
    {
        itSubSystems->second->outputToFile( outFile );
    }
}


bool mknix::System::checkAugmented()
{
    bool convergence = 1;

    for ( itConstraints = constraints.begin();
            itConstraints!= constraints.end();
            ++itConstraints
        )
    {
        if( !itConstraints->second->checkAugmented() )
            convergence = 0;
    }

    for ( itSubSystems = subSystems.begin();
            itSubSystems!= subSystems.end();
            ++itSubSystems
        )
    {
        if( !itSubSystems->second->checkAugmented() )
            convergence = 0;
    }

    return convergence;
}

void mknix::System::clearAugmented()
{
    for ( itConstraints = constraints.begin();
            itConstraints!= constraints.end();
            ++itConstraints
        )
    {
        itConstraints->second->clearAugmented();
    }
    for ( itSubSystems = subSystems.begin();
            itSubSystems!= subSystems.end();
            ++itSubSystems
        )
    {
        itSubSystems->second->clearAugmented();
    }

}


void mknix::System::writeBoundaryNodes( std::vector<Point*>& boundary_nodes )
{
    for ( itRigidBodies = rigidBodies.begin();
            itRigidBodies!= rigidBodies.end();
            ++itRigidBodies
        )
    {
        itRigidBodies->second->writeBoundaryNodes( boundary_nodes );
    }

    for ( itFlexBodies = flexBodies.begin();
            itFlexBodies!= flexBodies.end();
            ++itFlexBodies
        )
    {
        itFlexBodies->second->writeBoundaryNodes( boundary_nodes );
    }

    for ( itSubSystems = subSystems.begin();
            itSubSystems!= subSystems.end();
            ++itSubSystems
        )
    {
        itSubSystems->second->writeBoundaryNodes( boundary_nodes );
    }
}


void mknix::System::writeBoundaryConnectivity( std::vector< std::vector<Point*> >& connectivity_nodes )
{
    for ( itRigidBodies = rigidBodies.begin();
            itRigidBodies!= rigidBodies.end();
            ++itRigidBodies
        )
    {
        itRigidBodies->second->writeBoundaryConnectivity( connectivity_nodes );
    }

    for ( itFlexBodies = flexBodies.begin();
            itFlexBodies!= flexBodies.end();
            ++itFlexBodies
        )
    {
        itFlexBodies->second->writeBoundaryConnectivity( connectivity_nodes );
    }

    for ( itSubSystems = subSystems.begin();
            itSubSystems!= subSystems.end();
            ++itSubSystems
        )
    {
        itSubSystems->second->writeBoundaryConnectivity( connectivity_nodes );
    }
}
