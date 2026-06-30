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

#include "system.h"

#include "body.h"
#include "bodyflex.h"
#include "bodyrigid.h"
#include "constraint.h"
#include "constraintthermal.h"
#include "load.h"
#include "loadthermal.h"
#include "motion.h"

namespace mknix
{

/**
 * @brief Default constructor.
 */
System::System()
    : outputMaxInterfaceTemp(false)
{
}


/**
 * @brief Constructor with title.
 * @param title_in Name identifier for this system.
 */
System::System(const std::string& title_in)
    : outputMaxInterfaceTemp(false)
    , title(title_in)
{
}


/**
 * @brief Destructor. Releases memory for sub-systems, bodies, loads and thermal loads.
 */
System::~System()
{
    for (auto& system : subSystems)
    {
        delete system.second;
    }
    for (auto& body : rigidBodies)
    {
        delete body.second;
    }
    for (auto& body : flexBodies)
    {
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
    for (auto& load : loads)
    {
        delete load;
    }
    for (auto& load : loadsThermal)
    {
        delete load;
    }
}

// FIXME: Specific for tiles subsystem. That can change in input
/**
 * @brief Collects the X coordinates of all thermal-load nodes in the "tiles" sub-system.
 * @param x_coordinates Vector to which the node X coordinates are appended.
 */
void System::getThermalNodes(std::vector<double>& x_coordinates)
{
    for (auto& load : subSystems["tiles"]->loadsThermal)
    {
        load->insertNodesXCoordinates(x_coordinates);
    }
}

/**
 * @brief Writes the current temperature of each output-signal node in the "tiles" sub-system
 *        (and optionally the maximum interface temperature) into the provided array.
 * @param vector_in Output array filled with temperature values.
 */
void System::getOutputSignalThermal(double * vector_in)
{
    int counter = 0;
    for (auto& signal : subSystems["tiles"]->outputSignalThermal)
    {
        vector_in[counter] = signal->getTemp();
        ++counter;
    }

    if (subSystems["tiles"]->outputMaxInterfaceTemp)
    {
        vector_in[counter] = 0;
        for (auto& load : subSystems["tiles"]->loadsThermal)
        {
            load->getMaxTemp(vector_in[counter]);
        }
    }
}

/**
 * @brief Updates the thermal loads in the "tiles" sub-system from an external array.
 * @param vector_in Array of new load values, one per thermal load.
 */
void System::updateThermalLoads(double * vector_in)
{
    int counter = 0;
    for (auto& load : subSystems["tiles"]->loadsThermal)
    {
        load->updateLoad(vector_in[counter]);
        ++counter;
    }
}

/**
 * @brief Propagates the current time to all motions and recursively to all sub-systems.
 * @param time Current simulation time.
 */
void System::update(double time)
{
    for (auto& motion : motions)
    {
        motion->update(time);
    }
    for (auto& subSystem : subSystems)
    {
        subSystem.second->update(time);
    }
}


/**
 * @brief Initializes all flexible bodies in this system.
 */
void System::initFlexBodies()
{
    for (auto& flexBody : flexBodies)
    {
        flexBody.second->initialize();
    }
}


/**
 * @brief Writes rigid-body descriptions to the output file for this system and all sub-systems.
 * @param outFile Pointer to the output file stream.
 */
void System::writeRigidBodies(std::ofstream * outFile)
{
    for (auto& rigidBody : rigidBodies)
    {
        rigidBody.second->writeBodyInfo(outFile);
    }

    for (auto& subSystem : subSystems)
    {
        subSystem.second->writeRigidBodies(outFile);
    }
}

/**
 * @brief Writes flexible-body descriptions to the output file for this system and all sub-systems.
 * @param outFile Pointer to the output file stream.
 */
void System::writeFlexBodies(std::ofstream * outFile)
{
    for (auto& body : flexBodies)
    {
        body.second->writeBodyInfo(outFile);
    }

    for (auto& system : subSystems)
    {
        system.second->writeFlexBodies(outFile);
    }
}

/**
 * @brief Writes joint/constraint descriptions to the output file for this system and all sub-systems.
 * @param outFile Pointer to the output file stream.
 */
void System::writeJoints(std::ofstream * outFile)
{
    for (auto& constraint : constraints)
    {
        constraint.second->writeJointInfo(outFile);
    }

    for (auto& system : subSystems)
    {
        system.second->writeJoints(outFile);
    }
}

} // Namespace mknix


/**
 * @brief Computes the thermal capacity matrices for all thermal bodies and sub-systems.
 */
void mknix::System::calcCapacityMatrix()
{
    for (auto& body : thermalBodies)
    {
        body.second->calcCapacityMatrix();
    }

    for (auto& system : subSystems)
    {
        system.second->calcCapacityMatrix();
    }
}

/**
 * @brief Computes the thermal conductivity matrices for all thermal bodies and sub-systems.
 */
void mknix::System::calcConductivityMatrix()
{
    for (auto& body : thermalBodies)
    {
        body.second->calcConductivityMatrix();
    }

//   for ( itConstraints = constraints.begin();
//         itConstraints!= constraints.end();
//         ++itConstraints
//       )
//   {
//     constraint.second->calcConductivityMatrix();
//   }

    for (auto& system : subSystems)
    {
        system.second->calcConductivityMatrix();
    }
}

/**
 * @brief Computes the external heat vectors for all thermal bodies and sub-systems.
 */
void mknix::System::calcExternalHeat()
{
    for (auto& body : thermalBodies)
    {
        body.second->calcExternalHeat();
    }

    for (auto& system : subSystems)
    {
        system.second->calcExternalHeat();
    }
}

/**
 * @brief Computes the internal heat vectors for all thermal constraints and sub-systems.
 */
void mknix::System::calcInternalHeat()
{
    for (auto& constraint : constraintsThermal)
    {
        constraint.second->calcInternalForces();
    }


    for (auto& system : subSystems)
    {
        system.second->calcInternalHeat();
    }
}

/**
 * @brief Computes the thermal tangent matrices for all thermal constraints and sub-systems.
 */
void mknix::System::calcThermalTangentMatrix()
{
    for (auto& constraint : constraintsThermal)
    {
        constraint.second->calcTangentMatrix();
    }

    for (auto& system : subSystems)
    {
        system.second->calcThermalTangentMatrix();
    }

}

/**
 * @brief Assembles the capacity matrices of all thermal bodies and sub-systems into the global matrix.
 * @param globalCapacity_in Reference to the global thermal capacity matrix.
 */
void mknix::System::assembleCapacityMatrix(lmx::Matrix<data_type>& globalCapacity_in)
{
    for (auto& body : thermalBodies)
    {
        body.second->assembleCapacityMatrix(globalCapacity_in);
    }

    for (auto& system : subSystems)
    {
        system.second->assembleCapacityMatrix(globalCapacity_in);
    }

}

/**
 * @brief Assembles the conductivity matrices of all thermal bodies and sub-systems into the global matrix.
 * @param globalConductivity_in Reference to the global thermal conductivity matrix.
 */
void mknix::System::assembleConductivityMatrix(lmx::Matrix<data_type>& globalConductivity_in)
{
    for (auto& body : thermalBodies)
    {
        body.second->assembleConductivityMatrix(globalConductivity_in);
    }

//   for ( itConstraints = constraints.begin();
//         itConstraints!= constraints.end();
//         ++itConstraints
//       )
//   {
//     constraint.second->assembleConductivityMatrix(globalConductivity_in);
//   }

    for (auto& system : subSystems)
    {
        system.second->assembleConductivityMatrix(globalConductivity_in);
    }

}

/**
 * @brief Assembles external heat contributions from all thermal bodies, thermal loads, and sub-systems.
 * @param externalHeat_in Reference to the global external heat vector.
 */
void mknix::System::assembleExternalHeat(lmx::Vector<data_type>& externalHeat_in)
{
    for (auto& body : thermalBodies)
    {
        body.second->assembleExternalHeat(externalHeat_in);
    }

    for (auto& load : loadsThermal)
    {
        load->assembleExternalHeat(externalHeat_in);
    }

//  cout << "External heat in System (1) = " << externalHeat_in;
    for (auto& system : subSystems)
    {
        system.second->assembleExternalHeat(externalHeat_in);
    }
}


/**
 * @brief Assembles internal heat contributions from all thermal constraints and sub-systems.
 * @param internalHeat_in Reference to the global internal heat vector.
 */
void mknix::System::assembleInternalHeat(lmx::Vector<data_type>& internalHeat_in)
{
    for (auto& system : subSystems)
    {
        system.second->assembleInternalHeat(internalHeat_in);
    }
    for (auto& constraint : constraintsThermal)
    {
        constraint.second->assembleInternalForces(internalHeat_in);
    }

}


/**
 * @brief Assembles the thermal tangent matrix from all thermal constraints and sub-systems.
 * @param globalTangent_in Reference to the global thermal tangent matrix.
 */
void mknix::System::assembleThermalTangentMatrix(lmx::Matrix<data_type>& globalTangent_in)
{
    for (auto& constraint : constraintsThermal)
    {
        constraint.second->assembleTangentMatrix(globalTangent_in);
    }

    for (auto& system : subSystems)
    {
        system.second->assembleThermalTangentMatrix(globalTangent_in);
    }
}


/**
 * @brief Computes the mass matrices for all flexible and rigid bodies and sub-systems.
 */
void mknix::System::calcMassMatrix()
{
    for (auto& body : flexBodies)
    {
        body.second->calcMassMatrix();
    }

    for (auto& body : rigidBodies)
    {
        body.second->calcMassMatrix();
    }

    for (auto& system : subSystems)
    {
        system.second->calcMassMatrix();
    }
}

/**
 * @brief Computes internal forces for all flexible bodies, constraints, and sub-systems.
 */
void mknix::System::calcInternalForces()
{
    for (auto& body : flexBodies)
    {
        body.second->calcInternalForces();
    }

    for (auto& constraint : constraints)
    {
        constraint.second->calcInternalForces();
    }

    for (auto& system : subSystems)
    {
        system.second->calcInternalForces();
    }

}

/**
 * @brief Computes external forces for all flexible bodies, rigid bodies, and sub-systems.
 */
void mknix::System::calcExternalForces()
{
    for (auto& body : flexBodies)
    {
        body.second->calcExternalForces();
    }

    for (auto& body : rigidBodies)
    {
        body.second->calcExternalForces();
    }

    for (auto& system : subSystems)
    {
        system.second->calcExternalForces();
    }

}

/**
 * @brief Computes tangent stiffness matrices for all flexible bodies, constraints, and sub-systems.
 */
void mknix::System::calcTangentMatrix()
{
    for (auto& body : flexBodies)
    {
        body.second->calcTangentMatrix();
    }


    for (auto& constraint : constraints)
    {
        constraint.second->calcTangentMatrix();
    }

    for (auto& system : subSystems)
    {
        system.second->calcTangentMatrix();
    }

}

/**
 * @brief Assembles mass matrices from all bodies and sub-systems into the global matrix.
 * @param globalMass_in Reference to the global mass matrix.
 */
void mknix::System::assembleMassMatrix(lmx::Matrix<data_type>& globalMass_in)
{
    for (auto& body : flexBodies)
    {
        body.second->assembleMassMatrix(globalMass_in);
    }

    for (auto& body : rigidBodies)
    {
        body.second->assembleMassMatrix(globalMass_in);
    }

    for (auto& system : subSystems)
    {
        system.second->assembleMassMatrix(globalMass_in);
    }

}

/**
 * @brief Assembles internal forces from all bodies, constraints, and sub-systems into the global vector.
 * @param internalForces_in Reference to the global internal force vector.
 */
void mknix::System::assembleInternalForces(lmx::Vector<data_type>& internalForces_in)
{
    for (auto& body : flexBodies)
    {
        body.second->assembleInternalForces(internalForces_in);
    }

    for (auto& constraint : constraints)
    {
        constraint.second->assembleInternalForces(internalForces_in);
    }

    for (auto& system : subSystems)
    {
        system.second->assembleInternalForces(internalForces_in);
    }

}

/**
 * @brief Assembles external forces from all bodies, loads, and sub-systems into the global vector.
 * @param externalForces_in Reference to the global external force vector.
 */
void mknix::System::assembleExternalForces(lmx::Vector<data_type>& externalForces_in)
{

    for (auto& body : flexBodies)
    {
        body.second->assembleExternalForces(externalForces_in);
    }

    for (auto& body : rigidBodies)
    {
        body.second->assembleExternalForces(externalForces_in);
    }

    for (auto& load : loads)
    {
        load->assembleExternalForces(externalForces_in);
    }

//  cout << "External in System (1) = " << externalForces_in;
    for (auto& system : subSystems)
    {
        system.second->assembleExternalForces(externalForces_in);
    }

//  cout << "External in System (2) = " << externalForces_in;
}

/**
 * @brief Assembles the tangent stiffness matrix from all bodies, constraints, and sub-systems.
 * @param globalTangent_in Reference to the global tangent matrix.
 */
void mknix::System::assembleTangentMatrix(lmx::Matrix<data_type>& globalTangent_in)
{
    for (auto& body : flexBodies)
    {
        body.second->assembleTangentMatrix(globalTangent_in);
    }

    for (auto& constraint : constraints)
    {
        constraint.second->assembleTangentMatrix(globalTangent_in);
    }

    for (auto& system : subSystems)
    {
        system.second->assembleTangentMatrix(globalTangent_in);
    }
}


/**
 * @brief Recomputes and assembles constraint internal forces into the global force vector.
 * @param internalForces_in Reference to the global internal force vector.
 */
void mknix::System::assembleConstraintForces(lmx::Vector<data_type>& internalForces_in)
{
    for (auto& constraint : constraints)
    {
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

    for (auto& system : subSystems)
    {
        system.second->assembleConstraintForces(internalForces_in);
    }
}


/**
 * @brief Marks all rigid and flexible bodies as mechanical (non-thermal), propagating to sub-systems.
 */
void mknix::System::setMechanical()
{
    for (auto& body : rigidBodies)
    {
        body.second->setMechanical();
    }

    for (auto& body : flexBodies)
    {
        body.second->setMechanical();
    }

    for (auto& system : subSystems)
    {
        system.second->setMechanical();
    }
}

/**
 * @brief Collects step output data for a dynamic step from all bodies, constraints, and sub-systems.
 * @param q    Global configuration vector.
 * @param qdot Global velocity vector.
 */
void mknix::System::outputStep(const lmx::Vector<data_type>& q, const lmx::Vector<data_type>& qdot)
{
    for (auto& body : rigidBodies)
    {
        body.second->outputStep(q, qdot);
    }

    for (auto& body : flexBodies)
    {
        body.second->outputStep(q, qdot);
    }

    for (auto& constraint : constraints)
    {
        constraint.second->outputStep(q, qdot);
    }

    for (auto& system : subSystems)
    {
        system.second->outputStep(q, qdot);
    }
}


/**
 * @brief Collects step output data for a static step from all bodies, constraints, and sub-systems.
 * @param q Global configuration vector.
 */
void mknix::System::outputStep(const lmx::Vector<data_type>& q)
{
    for (auto& body : rigidBodies)
    {
        body.second->outputStep(q);
    }

    for (auto& body : flexBodies)
    {
        body.second->outputStep(q);
    }

    for (auto& constraint : constraints)
    {
        constraint.second->outputStep(q);
    }

    for (auto& system : subSystems)
    {
        system.second->outputStep(q);
    }
}


/**
 * @brief Streams all stored results for bodies, loads, and constraints to the output file.
 * @param outFile Pointer to the output file stream.
 */
void mknix::System::outputToFile(std::ofstream * outFile)
{
    for (auto& body : rigidBodies)
    {
        body.second->outputToFile(outFile);
    }

    for (auto& body : flexBodies)
    {
        body.second->outputToFile(outFile);
    }

    for (auto& load : loads)
    {
        load->outputToFile(outFile);
    }

    for (auto& constraint : constraints)
    {
        constraint.second->outputToFile(outFile);
    }

    for (auto& system : subSystems)
    {
        system.second->outputToFile(outFile);
    }
}


/**
 * @brief Checks whether all constraints (mechanical and thermal) in this system and sub-systems
 *        have satisfied the augmented-Lagrangian convergence criterion.
 * @return True if all constraints converged; false otherwise.
 */
bool mknix::System::checkAugmented()
{
    bool convergence = 1;

    for (auto& constraint : constraints)
    {
//         bool convergence = 1;
        if (!constraint.second->checkAugmented())
        {
            convergence = 0;
        }
    }

    for (auto& constraintThermal : constraintsThermal)
    {
//         bool convergence = 1;
        if (!constraintThermal.second->checkAugmented())
        {
            convergence = 0;
        }
    }

    for (auto& system : subSystems)
    {
        if (!system.second->checkAugmented())
        {
            convergence = 0;
        }
    }

    return convergence;
}

/**
 * @brief Resets the Lagrange multipliers of all constraints and sub-systems to zero.
 */
void mknix::System::clearAugmented()
{
    for (auto& constraint : constraints)
    {
        constraint.second->clearAugmented();
    }
    for (auto& constraint : constraintsThermal)
    {
        constraint.second->clearAugmented();
    }
    for (auto& system : subSystems)
    {
        system.second->clearAugmented();
    }

}


/**
 * @brief Collects boundary-node pointers from all bodies and sub-systems.
 * @param boundary_nodes Vector to which the boundary node pointers are appended.
 */
void mknix::System::writeBoundaryNodes(std::vector<Point *>& boundary_nodes)
{
    for (auto& body : rigidBodies)
    {
        body.second->writeBoundaryNodes(boundary_nodes);
    }

    for (auto& body : flexBodies)
    {
        body.second->writeBoundaryNodes(boundary_nodes);
    }

    for (auto& system : subSystems)
    {
        system.second->writeBoundaryNodes(boundary_nodes);
    }
}


/**
 * @brief Builds the ordered boundary connectivity from all bodies and sub-systems.
 * @param connectivity_nodes Vector of node chains to which each body's boundary is appended.
 */
void mknix::System::writeBoundaryConnectivity(std::vector<std::vector<Point *> >& connectivity_nodes)
{
    for (auto& body : rigidBodies)
    {
        body.second->writeBoundaryConnectivity(connectivity_nodes);
    }

    for (auto& body : flexBodies)
    {
        body.second->writeBoundaryConnectivity(connectivity_nodes);
    }

    for (auto& system : subSystems)
    {
        system.second->writeBoundaryConnectivity(connectivity_nodes);
    }
}

/**
 * @brief Finds and returns a body by its system and body names.
 * @param system_name Name of the sub-system containing the body.
 * @param body_name   Name of the body to retrieve.
 * @return Pointer to the found Body.
 * @throws std::out_of_range if the system or body is not found.
 */
mknix::Body * mknix::System::getBody(const std::string& system_name, const std::string& body_name)
{
    auto it = std::find_if(subSystems.begin(), subSystems.end(),
                           [&system_name](std::pair<std::string, System *> el)
    {
        return el.first == system_name;
    });

    if (it == subSystems.end())
    {
        throw std::out_of_range("system " + system_name + " not found");
    }

    auto system = it->second;

    if (system->rigidBodies.count(body_name))
    {
        return system->rigidBodies[body_name];
    }
    else if (system->flexBodies.count(body_name))
    {
        return system->flexBodies[body_name];
    }

    throw std::out_of_range("body " + body_name + " not found");
}
