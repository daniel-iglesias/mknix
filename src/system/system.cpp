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
//#include "gpu/functions_cpu.h"

namespace mknix {

System::System()
        : outputMaxInterfaceTemp(false)
{
}


System::System(const std::string& title_in)
        : outputMaxInterfaceTemp(false)
        , title(title_in)
{
}


System::~System()
{
    for (auto& system : subSystems) {
        delete system.second;
    }
    for (auto& body : rigidBodies) {
        delete body.second;
    }
    for (auto& body : flexBodies) {
        delete body.second;
    }
//     for ( itConstraints = constraints.begin();
//             itConstraints!= constraints.end();
//             ++itConstraints
//         )
//     {
//         delete(constraint.second);
//     }
//     for ( itConstraintsThermal = constraintsThermal.begin();
//             itConstraintsThermal!= constraintsThermal.end();
//             ++itConstraintsThermal
//         )
//     {
//         delete(itConstraintsThermal->second);
//     }
    for (auto& load : loads) {
        delete load;
    }
    for (auto& load : loadsThermal) {
        delete load;
    }
}

// FIXME: Specific for tiles subsystem. That can change in input
void System::getThermalNodes(std::vector<double>& x_coordinates)
{
    for (auto& load : subSystems["tiles"]->loadsThermal) {
        load->insertNodesXCoordinates(x_coordinates);
    }
}

void System::getOutputSignalThermal(double * vector_in)
{
    int counter = 0;
    for (auto& signal : subSystems["tiles"]->outputSignalThermal) {
        vector_in[counter] = signal->getTemp();
        ++counter;
    }

    if (subSystems["tiles"]->outputMaxInterfaceTemp) {
        vector_in[counter] = 0;
        for (auto& load : subSystems["tiles"]->loadsThermal) {
            load->getMaxTemp(vector_in[counter]);
        }
    }
}

void System::updateThermalLoads(double * vector_in)
{
    int counter = 0;
    for (auto& load : subSystems["tiles"]->loadsThermal) {
        load->updateLoad(vector_in[counter]);
        ++counter;
    }
}

void System::update(double time)
{
    for (auto& motion : motions) {
        motion->update(time);
    }
    for (auto& subSystem : subSystems) {
        subSystem.second->update(time);
    }
}


void System::initFlexBodies()
{
    for (auto& flexBody : flexBodies) {
        flexBody.second->initialize();
    }
}


void System::writeRigidBodies(std::ofstream * outFile)
{
    for (auto& rigidBody : rigidBodies) {
        rigidBody.second->writeBodyInfo(outFile);
    }

    for (auto& subSystem : subSystems) {
        subSystem.second->writeRigidBodies(outFile);
    }
}

void System::writeFlexBodies(std::ofstream * outFile)
{
    for (auto& body : flexBodies) {
        body.second->writeBodyInfo(outFile);
    }

    for (auto& system : subSystems) {
        system.second->writeFlexBodies(outFile);
    }
}

void System::writeJoints(std::ofstream * outFile)
{
    for (auto& constraint : constraints) {
        constraint.second->writeJointInfo(outFile);
    }

    for (auto& system : subSystems) {
        system.second->writeJoints(outFile);
    }
}

} // Namespace mknix

void mknix::System::setMaterialTable(MaterialTable* mt_ptr)
{
  for (auto& body : thermalBodies) {
      body.second->setMaterialTable(mt_ptr);
  }

  for (auto& system : subSystems) {
      system.second->setMaterialTable(mt_ptr);
  }
}

void mknix::System::setTemperatureVector(lmx::Vector<data_type>& q)
{
  for (auto& body : thermalBodies) {
      body.second->setTemperatureVector(q);
  }

  for (auto& system : subSystems) {
      system.second->setTemperatureVector(q);
  }
}

void mknix::System::calcFactors()
{
  for (auto& body : thermalBodies) {
      body.second->calcFactors();
  }

  for (auto& system : subSystems) {
      system.second->calcFactors();
  }
}

void mknix::System::calcCapacityMatrix()
{
    for (auto& body : thermalBodies) {
        body.second->calcCapacityMatrix();
    }

    for (auto& system : subSystems) {
        system.second->calcCapacityMatrix();
    }
}

void mknix::System::calcConductivityMatrix()
{
    for (auto& body : thermalBodies) {
        body.second->calcConductivityMatrix();
    }

//   for ( itConstraints = constraints.begin();
//         itConstraints!= constraints.end();
//         ++itConstraints
//       )
//   {
//     constraint.second->calcConductivityMatrix();
//   }

    for (auto& system : subSystems) {
        system.second->calcConductivityMatrix();
    }
}

void mknix::System::calcExternalHeat()
{
    for (auto& body : thermalBodies) {
        body.second->calcExternalHeat();
    }

    for (auto& system : subSystems) {
        system.second->calcExternalHeat();
    }
}

void mknix::System::calcInternalHeat()
{
    for (auto& constraint : constraintsThermal) {
        constraint.second->calcInternalForces();
    }


    for (auto& system : subSystems) {
        system.second->calcInternalHeat();
    }
}

void mknix::System::calcThermalTangentMatrix()
{
    for (auto& constraint : constraintsThermal) {
        constraint.second->calcTangentMatrix();
    }

    for (auto& system : subSystems) {
        system.second->calcThermalTangentMatrix();
    }

}

void mknix::System::assembleCapacityMatrix(lmx::Matrix<data_type>& globalCapacity_in)
{
    for (auto& body : thermalBodies) {
        body.second->assembleCapacityMatrix(globalCapacity_in);
    }

    for (auto& system : subSystems) {
        system.second->assembleCapacityMatrix(globalCapacity_in);
    }

}

void mknix::System::assembleConductivityMatrix(lmx::Matrix<data_type>& globalConductivity_in)
{
    for (auto& body : thermalBodies) {
        body.second->assembleConductivityMatrix(globalConductivity_in);
    }

//   for ( itConstraints = constraints.begin();
//         itConstraints!= constraints.end();
//         ++itConstraints
//       )
//   {
//     constraint.second->assembleConductivityMatrix(globalConductivity_in);
//   }

    for (auto& system : subSystems) {
        system.second->assembleConductivityMatrix(globalConductivity_in);
    }

}

void mknix::System::assembleExternalHeat(lmx::Vector<data_type>& externalHeat_in)
{
    for (auto& body : thermalBodies) {
        body.second->assembleExternalHeat(externalHeat_in);
    }

    for (auto& load : loadsThermal) {
        load->assembleExternalHeat(externalHeat_in);
    }

//  cout << "External heat in System (1) = " << externalHeat_in;
    for (auto& system : subSystems) {
        system.second->assembleExternalHeat(externalHeat_in);
    }
}


void mknix::System::assembleInternalHeat(lmx::Vector<data_type>& internalHeat_in)
{
    for (auto& system : subSystems) {
        system.second->assembleInternalHeat(internalHeat_in);
    }
    for (auto& constraint : constraintsThermal) {
        constraint.second->assembleInternalForces(internalHeat_in);
    }

}


void mknix::System::assembleThermalTangentMatrix(lmx::Matrix<data_type>& globalTangent_in)
{
    for (auto& constraint : constraintsThermal) {
        constraint.second->assembleTangentMatrix(globalTangent_in);
    }

    for (auto& system : subSystems) {
        system.second->assembleThermalTangentMatrix(globalTangent_in);
    }
}


void mknix::System::calcMassMatrix()
{
    for (auto& body : flexBodies) {
        body.second->calcMassMatrix();
    }

    for (auto& body : rigidBodies) {
        body.second->calcMassMatrix();
    }

    for (auto& system : subSystems) {
        system.second->calcMassMatrix();
    }
}

void mknix::System::calcInternalForces()
{
    for (auto& body : flexBodies) {
        body.second->calcInternalForces();
    }

    for (auto& constraint : constraints) {
        constraint.second->calcInternalForces();
    }

    for (auto& system : subSystems) {
        system.second->calcInternalForces();
    }

}

void mknix::System::calcExternalForces()
{
    for (auto& body : flexBodies) {
        body.second->calcExternalForces();
    }

    for (auto& body : rigidBodies) {
        body.second->calcExternalForces();
    }

    for (auto& system : subSystems) {
        system.second->calcExternalForces();
    }

}

void mknix::System::calcTangentMatrix()
{
    for (auto& body : flexBodies) {
        body.second->calcTangentMatrix();
    }


    for (auto& constraint : constraints) {
        constraint.second->calcTangentMatrix();
    }

    for (auto& system : subSystems) {
        system.second->calcTangentMatrix();
    }

}

void mknix::System::assembleMassMatrix(lmx::Matrix<data_type>& globalMass_in)
{
    for (auto& body : flexBodies) {
        body.second->assembleMassMatrix(globalMass_in);
    }

    for (auto& body : rigidBodies) {
        body.second->assembleMassMatrix(globalMass_in);
    }

    for (auto& system : subSystems) {
        system.second->assembleMassMatrix(globalMass_in);
    }

}

void mknix::System::assembleInternalForces(lmx::Vector<data_type>& internalForces_in)
{
    for (auto& body : flexBodies) {
        body.second->assembleInternalForces(internalForces_in);
    }

    for (auto& constraint : constraints) {
        constraint.second->assembleInternalForces(internalForces_in);
    }

    for (auto& system : subSystems) {
        system.second->assembleInternalForces(internalForces_in);
    }

}

void mknix::System::assembleExternalForces(lmx::Vector<data_type>& externalForces_in)
{

    for (auto& body : flexBodies) {
        body.second->assembleExternalForces(externalForces_in);
    }

    for (auto& body : rigidBodies) {
        body.second->assembleExternalForces(externalForces_in);
    }

    for (auto& load : loads) {
        load->assembleExternalForces(externalForces_in);
    }

//  cout << "External in System (1) = " << externalForces_in;
    for (auto& system : subSystems) {
        system.second->assembleExternalForces(externalForces_in);
    }

//  cout << "External in System (2) = " << externalForces_in;
}

void mknix::System::assembleTangentMatrix(lmx::Matrix<data_type>& globalTangent_in)
{
    for (auto& body : flexBodies) {
        body.second->assembleTangentMatrix(globalTangent_in);
    }

    for (auto& constraint : constraints) {
        constraint.second->assembleTangentMatrix(globalTangent_in);
    }

    for (auto& system : subSystems) {
        system.second->assembleTangentMatrix(globalTangent_in);
    }
}


void mknix::System::assembleConstraintForces(lmx::Vector<data_type>& internalForces_in)
{
    for (auto& constraint : constraints) {
        constraint.second->calcInternalForces();
        constraint.second->assembleInternalForces(internalForces_in);
    }

//   for ( itFlexBodies = flexBodies.begin();
//         itFlexBodies!= flexBodies.end();
//         ++itFlexBodies
//       )
//   {
//     body.second->calcInternalForces();
//     body.second->assembleInternalForces(internalForces_in);
//   }

    for (auto& system : subSystems) {
        system.second->assembleConstraintForces(internalForces_in);
    }
}


void mknix::System::setMechanical()
{
    for (auto& body : rigidBodies) {
        body.second->setMechanical();
    }

    for (auto& body : flexBodies) {
        body.second->setMechanical();
    }

    for (auto& system : subSystems) {
        system.second->setMechanical();
    }
}

void mknix::System::outputStep(const lmx::Vector<data_type>& q, const lmx::Vector<data_type>& qdot)
{
    for (auto& body : rigidBodies) {
        body.second->outputStep(q, qdot);
    }

    for (auto& body : flexBodies) {
        body.second->outputStep(q, qdot);
    }

    for (auto& constraint : constraints) {
        constraint.second->outputStep(q, qdot);
    }

    for (auto& system : subSystems) {
        system.second->outputStep(q, qdot);
    }
}


void mknix::System::outputStep(const lmx::Vector<data_type>& q)
{
    for (auto& body : rigidBodies) {
        body.second->outputStep(q);
    }

    for (auto& body : flexBodies) {
        body.second->outputStep(q);
    }

    for (auto& constraint : constraints) {
        constraint.second->outputStep(q);
    }

    for (auto& system : subSystems) {
        system.second->outputStep(q);
    }
}


void mknix::System::outputToFile(std::ofstream * outFile)
{
    for (auto& body : rigidBodies) {
        body.second->outputToFile(outFile);
    }

    for (auto& body : flexBodies) {
        body.second->outputToFile(outFile);
    }

    for (auto& load : loads) {
        load->outputToFile(outFile);
    }

    for (auto& constraint : constraints) {
        constraint.second->outputToFile(outFile);
    }

    for (auto& system : subSystems) {
        system.second->outputToFile(outFile);
    }
}


bool mknix::System::checkAugmented()
{
    bool convergence = 1;

    for (auto& constraint : constraints) {
//         bool convergence = 1;
        if (!constraint.second->checkAugmented()) {
            convergence = 0;
        }
    }

    for (auto& constraintThermal : constraintsThermal) {
//         bool convergence = 1;
        if (!constraintThermal.second->checkAugmented()) {
            convergence = 0;
        }
    }

    for (auto& system : subSystems) {
        if (!system.second->checkAugmented()) {
            convergence = 0;
        }
    }

    return convergence;
}

void mknix::System::clearAugmented()
{
    for (auto& constraint : constraints) {
        constraint.second->clearAugmented();
    }
    for (auto& constraint : constraintsThermal) {
        constraint.second->clearAugmented();
    }
    for (auto& system : subSystems) {
        system.second->clearAugmented();
    }

}


void mknix::System::writeBoundaryNodes(std::vector<Point *>& boundary_nodes)
{
    for (auto& body : rigidBodies) {
        body.second->writeBoundaryNodes(boundary_nodes);
    }

    for (auto& body : flexBodies) {
        body.second->writeBoundaryNodes(boundary_nodes);
    }

    for (auto& system : subSystems) {
        system.second->writeBoundaryNodes(boundary_nodes);
    }
}


void mknix::System::writeBoundaryConnectivity(std::vector<std::vector<Point *> >& connectivity_nodes)
{
    for (auto& body : rigidBodies) {
        body.second->writeBoundaryConnectivity(connectivity_nodes);
    }

    for (auto& body : flexBodies) {
        body.second->writeBoundaryConnectivity(connectivity_nodes);
    }

    for (auto& system : subSystems) {
        system.second->writeBoundaryConnectivity(connectivity_nodes);
    }
}
