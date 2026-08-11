/******************************************************************************
 *  Copyright (C) 2015 by Daniel Iglesias                                     *
 *                                                                            *
 *  This file is part of MkniX.                                               *
 *                                                                            *
 *  MkniX is free software: you can redistribute it and/or modify             *
 *  it under the terms of the GNU Lesser General Public License as            *
 *  published by the Free Software Foundation, either version 3 of the        *
 *  License, or (at your option) any later version.                           *
 *                                                                            *
 *  MkniX is distributed in the hope that it will be useful,                  *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU Lesser General Public License for more details.                       *
 *                                                                            *
 *  You should have received a copy of the GNU Lesser General Public          *
 *  License along with MkniX.  If not, see <http://www.gnu.org/licenses/>.    *
 *****************************************************************************/

#include "version.h" 
#include "simulation.h"
#include "analysisdynamic.h"

#include <system/system.h>
#include <system/constraintthermal.h>
#include <reader/reader.h>
#include <system/generalcontact.h>
#include <system/body.h>

namespace mknix
{


// Static variables:
double Simulation::stepTime = 0;
double Simulation::oldClockTime = 0;
lmx::Vector<double> Simulation::gravity = lmx::Vector<double>(3);
double Simulation::alpha = 1E4;
double Simulation::augmentedTolerance = 5.0;
int Simulation::dimension = 2;
std::string Simulation::contact = "NONE";
bool Simulation::visualization = 0;
bool Simulation::outputMatrices = 0;
std::string Simulation::constraintMethod = "PENALTY";
double Simulation::epsilon = 1E-5;
std::string Simulation::smoothingType = "GLOBAL";

/**
 * @brief Returns a component of the global gravity vector.
 * @param component Index of the gravity vector component (0, 1, or 2).
 * @return The gravity value for the specified component.
 */
double Simulation::getGravity(int component)
{
    return gravity.readElement(component);
}

/**
 * @brief Returns the penalty parameter alpha used for constraint enforcement.
 * @return The alpha penalty parameter.
 */
double Simulation::getAlpha()
{
    return alpha;
}

/**
 * @brief Returns the convergence tolerance for the augmented Lagrangian method.
 * @return The augmented Lagrangian tolerance value.
 */
double Simulation::getAugmentedTolerance()
{
    return augmentedTolerance;
}

/**
 * @brief Sets the convergence tolerance for the augmented Lagrangian method.
 * @param tolerance The new tolerance value.
 */
void Simulation::setAugmentedTolerance(double tolerance)
{
    augmentedTolerance = tolerance;
}

/**
 * @brief Returns the current simulation step time.
 * @return The current value of stepTime.
 */
double Simulation::getTime()
{
    return Simulation::stepTime;
}

/**
 * @brief Returns the spatial dimension of the simulation (2 or 3).
 * @return The number of spatial dimensions.
 */
int Simulation::getDim()
{
    return dimension;
}

/**
 * @brief Returns the name of the constraint enforcement method in use.
 * @return A string identifying the constraint method (e.g. "PENALTY").
 */
std::string Simulation::getConstraintMethod()
{
    return constraintMethod;
}

/**
 * @brief Returns the smoothing type used in the simulation.
 * @return A string identifying the smoothing type (e.g. "GLOBAL").
 */
std::string Simulation::getSmoothingType()
{
    return smoothingType;
}

/**
 * @brief Constructs a Simulation object and initializes output files and the global timer.
 */
Simulation::Simulation()
    : baseSystem(nullptr)
//  , stepTime(0.)
    , timerFile(nullptr)
    , configurationFile(nullptr)
    , iterationsNLSolver(0)
    , outputFilesDetail(2)
//     , outFile(0)
    , initialTemperature(0)
{
    globalTimer = new lmx::ExactStopwatch;
    globalTimer->setQuiet();
    if (outputFilesDetail > 0)
    {
        timerFile = new std::ofstream("simulation_times.dat");
        if (outputFilesDetail > 1)
        {
            configurationFile = new std::ofstream("dis.dat");
        }
    }
}

/**
 * @brief Destroys the Simulation object and frees allocated timers and output file streams.
 */
Simulation::~Simulation()
{
    if (globalTimer) delete globalTimer;
    if (timerFile) delete timerFile;
    if (configurationFile) delete configurationFile;
}


/**
 * @brief Reads and parses the simulation configuration from an input file.
 * @param FileIn Path to the input file to read.
 */
void Simulation::inputFromFile(const std::string& FileIn)
{
    auto reader = make_unique<Reader>(this);
    if (!baseSystem)
    {
        baseSystem = make_unique<System>("baseSystem");
    }
    reader->inputFromFile(FileIn);
}

/**
 * @brief Returns the number of nodes in a named subsystem interface.
 * @param name Name of the subsystem.
 * @return Number of nodes in the specified subsystem.
 */
size_t Simulation::getInterfaceNumberOfNodes(const std::string& name) const
{
    return baseSystem->getSystem(name)->getNumberOfNodes();
}

/**
 * @brief Returns a node from a named subsystem by index.
 * @param system_name Name of the subsystem.
 * @param num Zero-based index of the node.
 * @return Pointer to the requested Node.
 */
Node* Simulation::getInterfaceNode(const std::string& system_name, size_t num) const
{
    return baseSystem->getSystem(system_name)->getNode(num);
}

/**
 * @brief Returns the signal nodes associated with a named signal in a subsystem.
 * @param system_name Name of the subsystem.
 * @param name Name of the signal.
 * @return Vector of pointers to the signal nodes.
 */
std::vector<Node*> Simulation::getSignalNodes(const std::string& system_name, const std::string& name) const
{
    return baseSystem->getSystem(system_name)->getSignalNodes(name);
}

/**
 * @brief Returns an output node identified by name from a subsystem.
 * @param system_name Name of the subsystem.
 * @param name Name of the output node.
 * @return Pointer to the requested output Node.
 */
Node* Simulation::getOuputNode(const std::string& system_name, const std::string& name) const
{
    return baseSystem->getSystem(system_name)->getOutputNode(name);
}

/**
 * @brief Returns the coordinates of all thermal interface nodes.
 * @return A vector of coordinate values for all thermal nodes.
 */
std::vector<double> Simulation::getInterfaceNodesCoords()
{
    // Loads have access to the nodes, and are part of the system.
    std::vector<double> temp_x_coordinates;
    baseSystem->getThermalNodes(temp_x_coordinates);

    return temp_x_coordinates;
}


/**
 * @brief Returns the names of all constraints defined in a subsystem.
 * @param systemName Name of the subsystem.
 * @return Vector of constraint name strings.
 */
std::vector<std::string> Simulation::getConstraintNames(const std::string& systemName) const
{
    auto system = baseSystem->getSystem(systemName);
    return system->getConstraintNames();
}

/**
 * @brief Returns a constraint by name from a subsystem, searching both mechanical and thermal constraints.
 * @param constraintName Name of the constraint.
 * @param systemName Name of the subsystem.
 * @return Pointer to the matching Constraint, or nullptr if not found.
 */
Constraint* Simulation::getConstraint(const std::string& constraintName, const std::string& systemName) const
{
    auto system = baseSystem->getSystem(systemName);
    auto constraint = system->getConstraint(constraintName);
    if (constraint == nullptr)
    {
        constraint = system->getConstraintThermal(constraintName);
    }
    return constraint;
}

/**
 * @brief Returns the internal reaction force of a constraint for a given component.
 * @param constraintName Name of the constraint.
 * @param systemName Name of the subsystem containing the constraint.
 * @param component Index of the force component to retrieve.
 * @return The negated internal force value for the specified component.
 */
double Simulation::getConstraintOutput(const std::string& constraintName, const std::string& systemName,
                                       size_t component)
{
    auto constraint = getConstraint(constraintName, systemName);
    return -(constraint->getInternalForces().readElement(component));
}


/**
 * @brief Sets the uniform initial temperature applied to all thermal nodes.
 * @param temp_in Initial temperature value.
 */
void Simulation::setInitialTemperatures(double temp_in)
{
    initialTemperature = temp_in;
}

/**
 * @brief Sets the initial temperature for a specific thermal node, overriding the uniform value.
 * @param node Pointer to the node whose temperature is to be set.
 * @param temperature Initial temperature value for this node.
 */
void Simulation::setThermalNodeInitialTemperature(Node* node, double temperature)
{
    if (node == nullptr)
    {
        return;
    }

    if (node->getThermalNumber() >= 0)
    {
        thermalNodeInitialOverrides[node->getThermalNumber()] = temperature;
    }
    else
    {
        groundedNodeInitialOverrides[node] = temperature;
    }

    node->setqt(temperature);
}

// Part copy of run(), limited to preparation and thermal dynamic analysis
/**
 * @brief Initializes all analyses without executing the full simulation loop.
 * @param verbosity Verbosity level controlling diagnostic output.
 */
void Simulation::init(int verbosity)
{
    if (outputFilesDetail > 1)
    {
        writeSystem();
        outFile << "ANALYSIS " << endl;
        outFile << "CONFIGURATION" << endl;
    }
    if (outputFilesDetail > 0)
    {
        *timerFile << stepTime << "\t" << globalTimer->getTime() << std::endl;
    }

    for (auto& analysis : analyses)
    {
        if (analysis->type() == "THERMAL" || analysis->type() == "THERMALSTATIC")
        {
            initThermalSimulation(analysis.get(), verbosity);
        }
        else
        {
            initMechanicalSimulation(analysis.get(), verbosity);
        }
    }
}

/**
 * @brief Initializes the global matrices and state vector for a thermal analysis.
 * @param theAnalysis_in Pointer to the Analysis object to initialize.
 * @param verbosity Verbosity level controlling diagnostic output.
 * @param init If true, calls the analysis init routine after setup.
 * @return The initial temperature state vector.
 */
lmx::Vector<data_type> Simulation::initThermalSimulation(Analysis* theAnalysis_in, int verbosity, bool init)
{
    theAnalysis = theAnalysis_in;
    auto gdlSize = nodes.size();
    lmx::Vector<data_type> q(gdlSize);
    q.fillIdentity(initialTemperature);

    for (const auto& thermalOverride : thermalNodeInitialOverrides)
    {
        if (thermalOverride.first >= 0
            && static_cast<std::size_t>(thermalOverride.first) < static_cast<std::size_t>(q.size()))
        {
            q(thermalOverride.first) = thermalOverride.second;
        }
    }

    globalConductivity.resize(gdlSize, gdlSize);
    globalCapacity.resize(gdlSize, gdlSize);
    globalRHSHeat.resize(gdlSize);
    globalExternalHeat.resize(gdlSize);
    globalInternalHeat.resize(gdlSize);

    baseSystem->calcConductivityMatrix();
    baseSystem->assembleConductivityMatrix(globalConductivity);
    baseSystem->calcCapacityMatrix();
    baseSystem->assembleCapacityMatrix(globalCapacity);

    // Set initial temperatures in nodes:
    for (auto& node : thermalNodes)
    {
        node.second->setqt(q);
    }

    for (const auto& groundedOverride : groundedNodeInitialOverrides)
    {
        if (groundedOverride.first != nullptr)
        {
            groundedOverride.first->setqt(groundedOverride.second);
        }
    }

    writeConfStep();

    if (outputFilesDetail > 1 && theAnalysis->type() == "THERMAL")
    {
        systemOuputStep(q);
    }

    if (init)
    {
        theAnalysis->init(&q, nullptr, verbosity);
    }

    return q;
}

/**
 * @brief Initializes the global matrices and state vector for a mechanical analysis.
 * @param analysis Pointer to the Analysis object to initialize.
 * @param verbosity Verbosity level controlling diagnostic output.
 * @param init If true, calls the analysis init routine after setup.
 * @return The initial displacement state vector.
 */
lmx::Vector<data_type> Simulation::initMechanicalSimulation(Analysis* analysis, int verbosity, bool init)
{
    theAnalysis = analysis;
    auto gdlSize = nodes.size() * Simulation::dimension;
    lmx::Vector<data_type> q(gdlSize);

    globalMass.resize(gdlSize, gdlSize);
    globalRHSForces.resize(gdlSize);
    globalInternalForces.resize(gdlSize);
    globalExternalForces.resize(gdlSize);

    auto i = 0u;
    for (auto& node : nodes)
    {
        q(Simulation::dimension * i) = node.second->getqx(0);
        q(Simulation::dimension * i + 1) = node.second->getqx(1);
        if (Simulation::getDim() == 3)
        {
            q(Simulation::dimension * i + 2) = node.second->getqx(2);
        }
        ++i;
    }

    baseSystem->calcMassMatrix();
    baseSystem->assembleMassMatrix(globalMass);
    // maybe is better to make an specific function call for the sparse
    // pattern, but this should work...
    if (lmx::getMatrixType() == 1)
    {
        globalSparsePattern.resize(gdlSize, gdlSize);
        baseSystem->calcTangentMatrix();
        baseSystem->assembleTangentMatrix(globalSparsePattern);
    }

    // Output matrices in initial configuration:
    if (outputMatrices)
    {
        lmx::Matrix<data_type> K_temp(gdlSize, gdlSize);
        baseSystem->calcTangentMatrix();
        baseSystem->assembleTangentMatrix(K_temp);
        K_temp.harwellBoeingSave((char*)"K.mat");
        globalMass.harwellBoeingSave((char*)"M.mat");
        // save raw matrices:
        std::ofstream Kfile("K");
        Kfile << K_temp;
        Kfile.close();
    }
    // Output to file the initial configuration:
    writeConfStep();

    if (outputFilesDetail > 1 && theAnalysis->type() == "DYNAMIC")
    {
        systemOuputStep(q);
    }

    if (init)
    {
        lmx::Vector<data_type> qdot(gdlSize);
        theAnalysis->init(&q, &qdot, verbosity);
    }

    return q;
}

/**
 * @brief Placeholder function for setting a signal on a node (currently a no-op).
 * @param node Name of the node.
 */
void setSignal(std::string node, std::vector<double>)
{
    return;
}

/**
 * @brief Placeholder function for retrieving a signal from a node (currently returns an empty vector).
 * @param node Name of the node.
 * @return An empty vector.
 */
std::vector<double> getSignal(const std::string& node)
{
    return { };
}

/**
 * @brief Advances the current analysis by one time step.
 */
void Simulation::solveStep()
{
    theAnalysis->nextStep();
}

/**
 * @brief Advances the current analysis by one time step, applying an input signal and capturing an output signal.
 * @param signal Pointer to the input signal array used to update thermal loads.
 * @param outputSignal Pointer to the output signal array; filled with thermal output if non-null.
 */
void Simulation::solveStep(double* signal, double* outputSignal)
{
    baseSystem->updateThermalLoads(signal);
    theAnalysis->nextStep();
    if (outputSignal)
    {
        baseSystem->getOutputSignalThermal(outputSignal);
    }
}

/**
 * @brief Finalizes the simulation, closes output files, and writes remaining results.
 */
void Simulation::endSimulation()
{
    configurationFile->close();

    analyses.clear();

    if (outputFilesDetail > 1)
    {
        std::ifstream disp("dis.dat");

        char a;
//       std::string aa;

        if (disp.is_open())
        {
            if (outFile.is_open())
            {
// 	  while (disp >> aa) {
// 	    outFile << aa;
// 	  }
                while (disp.get(a))
                {
                    outFile.put(a);
                    outFile.flush();
                }
                outFile << "ENDCONFIGURATION" << endl;

                // output extra flexible bodies data...
                baseSystem->outputToFile(&outFile);

            }
        }
    }

//       cout << "q(" << q.size()/2 << ") = " << q(q.size()/2) << endl;
//       cout << "q(" << q.size()-1 << ") = " << q(q.size()-1) << endl;

    // output f_int of constraints...
//     lmx::Vector<double> constr_forces(nodes.size()*dimension);
//     baseSystem->assembleConstraintForces( constr_forces );
//     for(size_type i=0; i< constr_forces.size(); ++i ) {
//         if(constr_forces(i) != 0. )
//             cout << "R(" << i << ") = " << constr_forces(i) << endl;
//     }
}


/**
 * @brief Executes the complete simulation by running all defined analyses in sequence.
 */
void Simulation::run()
{
#ifdef HAVE_VTK
    if(Simulation::contact == "GLOBAL" || Simulation::visualization == 1)
    {
        this->theContact = new Contact(this, 10.);
        this->theContact->createPoints();
        this->theContact->createPolys();
        this->theContact->createDrawingObjects();
        if(Simulation::contact == "GLOBAL")
        {
            this->theContact->createDelaunay();
            this->theContact->createDrawingContactObjects();
        }
        this->theContact->drawObjects();
    }
#endif
    writeSystem();
    outFile << "ANALYSIS " << endl;
    outFile << "CONFIGURATION" << endl;
//   << "FILE " << "dis.dat" << endl;
    *timerFile << stepTime << "\t" << globalTimer->getTime() << std::endl;

    for (auto& analysis : analyses)
    {
        if (analysis->type() == "THERMAL" || analysis->type() == "THERMALSTATIC")
        {
            this->runThermalAnalysis(analysis.get());
        }
        else
        {
            if (analysis->type() == "STATIC" || analysis->type() == "DYNAMIC")
            {
                this->baseSystem->setMechanical();
            }
            this->runMechanicalAnalysis(analysis.get());
        }
    }
}

/**
 * @brief Initializes and solves a thermal or thermalstatic analysis.
 * @param theAnalysis_in Pointer to the thermal Analysis object to run.
 */
void Simulation::runThermalAnalysis(Analysis* theAnalysis_in)
{
    auto q = initThermalSimulation(theAnalysis_in, 1, false);

    if (theAnalysis->type() == "THERMAL")
    {
        theAnalysis->solve(&q);
    }

    else if (theAnalysis->type() == "THERMALSTATIC")
    {
        // write initial configuration...
        outFile << "0 "; //time=0
//           systemOuputStep( q ); // produces output of temperatures, incompatible with mknixpost-static
        for (auto& point : outputPoints)
        {
            outFile << point.second->getConf(0) << " ";
            outFile << point.second->getConf(1) << " ";
            if (Simulation::getDim() == 3)
            {
                outFile << point.second->getConf(2) << " ";
            }
        }
        outFile << endl;
        theAnalysis->solve(&q);
    }

    if (outputFilesDetail > 1)
    {
        std::ifstream disp("dis.dat");
        char a;

        if (outFile.is_open())
        {
            while (disp.get(a))
            {
                outFile.put(a);
                outFile.flush();
            }
            outFile << "ENDCONFIGURATION" << endl;

            // output extra flexible bodies data...
            baseSystem->outputToFile(&outFile);

            // output material data...
            outFile << "MATERIALS data:" << endl;
            for (auto& mat : materials)
            {
                outFile << "Material " << mat.first << ": " ;
                mat.second.outputToFile(&outFile);
            }
        }
    }

    // output f_int of constraints...
    lmx::Vector<double> constr_forces(nodes.size() * dimension);
    baseSystem->assembleConstraintForces(constr_forces);
    for (size_type i = 0; i < constr_forces.size(); ++i)
    {
//         if(constr_forces(i) != 0. )
//             cout << "R(" << i << ") = " << constr_forces(i) << endl;
    }
}

/**
 * @brief Initializes and solves a mechanical (static, dynamic, or thermo-mechanical dynamic) analysis.
 * @param theAnalysis_in Pointer to the mechanical Analysis object to run.
 */
void Simulation::runMechanicalAnalysis(Analysis* theAnalysis_in)
{
    auto q = initMechanicalSimulation(theAnalysis_in, 1, false);

    auto gdlSize = nodes.size() * Simulation::dimension;

    if (theAnalysis->type() == "DYNAMIC")
    {
        lmx::Vector<data_type> qdot(gdlSize);
        // output first step data
        theAnalysis->setEpsilon(epsilon);
        theAnalysis->solve(&q, &qdot);
    }
    else if (theAnalysis->type() == "STATIC")
    {
        theAnalysis->solve(&q);
    }
    else if (theAnalysis->type() == "THERMOMECHANICALDYNAMIC")
    {
        // Init thermal conf vector and a zero velocity vector
        auto thermalSize = thermalNodes.size();
        lmx::Vector<data_type> qdot(gdlSize);
        lmx::Vector<data_type> qt(thermalSize);
        auto i = 0u;
        for (auto& node : thermalNodes)
        {
            qt(i) = node.second->getqt();
            ++i;
        }
        globalCapacity.resize(thermalSize, thermalSize);
        globalConductivity.resize(thermalSize, thermalSize);
        globalRHSHeat.resize(thermalSize);
        globalExternalHeat.resize(thermalSize);
        globalInternalHeat.resize(thermalSize);

        baseSystem->calcConductivityMatrix();
        baseSystem->assembleConductivityMatrix(globalConductivity);
        baseSystem->calcCapacityMatrix();
        baseSystem->assembleCapacityMatrix(globalCapacity);

        // output first step data
        systemOuputStep(q);
        theAnalysis->setEpsilon(epsilon);
        theAnalysis->solve(&qt, &q, &qdot);
    }

    if (outputFilesDetail > 1)
    {
        std::ifstream disp("dis.dat");
        char a;

        if (outFile.is_open())
        {
            while (disp.get(a))
            {
                outFile.put(a);
                outFile.flush();
            }
            outFile << "ENDCONFIGURATION" << endl;

            // output extra flexible bodies data...
            baseSystem->outputToFile(&outFile);

            // output material data...
            outFile << "MATERIALS data:" << endl;
            for (auto& mat : materials)
            {
                mat.second.outputToFile(&outFile);
            }
        }
    }

//       cout << "q(" << q.size()/2 << ") = " << q(q.size()/2) << endl;
//       cout << "q(" << q.size()-1 << ") = " << q(q.size()-1) << endl;

    // output f_int of constraints...
//       lmx::Vector<double> constr_forces(nodes.size()*dimension);
//       baseSystem->assembleConstraintForces( constr_forces );
//       for(size_type i=0; i< constr_forces.size(); ++i ){
//         if(constr_forces(i) != 0. )
//           cout << "R(" << i << ") = " << constr_forces(i) << endl;
//       }
}

/**
 * @brief Writes the initial system topology (nodes, rigid/flex bodies, joints) to the main output file
 *        and a nodes.dat file.
 */
void Simulation::writeSystem()
{
    std::stringstream ss;
    ss << title << ".mec";
    auto outFileName = ss.str();
    outFile.open(outFileName.c_str(), std::ofstream::out);
    if (outFile.fail())
    {
        cout << "\n\nCOULD NOT SAVE OUTPUT. Error opening file!\n";
        system("pause");
        return;
    }

    outFile << "GIT_REV is " << GIT_REV << endl;
    outFile << "GIT_REV_VAL is " << GIT_REV_VAL << endl;
    outFile << "GIT_REV_FULL is " << GIT_REV_FULL << endl;
    outFile << "GIT_TAG is " << GIT_TAG << endl;
    outFile << "GIT_BRANCH is " << GIT_BRANCH << endl << endl;


    outFile << "DIMENSION " << Simulation::dimension << endl;

    outFile << "SYSTEM" << endl;

    outFile << "NODES" << endl;
    for (auto& point : outputPoints)
    {
        outFile
                << "\t" << point.first
                << "\t" << point.second->getConf(0)
                << "\t" << point.second->getConf(1)
                << "\t" << point.second->getConf(2)
                << endl;
    }
    outFile << "ENDNODES" << endl;

    outFile << "RIGIDBODIES" << endl;
    baseSystem->writeRigidBodies(&outFile);
    outFile << "ENDRIGIDBODIES" << endl;

    outFile << "FLEXBODIES" << endl;
    baseSystem->writeFlexBodies(&outFile);
    outFile << "ENDFLEXBODIES" << endl;

    outFile << "JOINTS" << endl;
    baseSystem->writeJoints(&outFile);
    outFile << "ENDJOINTS" << endl;

    outFile << "ENDSYSTEM" << endl;

    // write a standard file for nodal info:
    std::ofstream nodeFile("nodes.dat");
    for (auto& node : nodes)
    {
        nodeFile
                << "\t" << node.first
                << "\t" << node.second->getConf(0)
                << "\t" << node.second->getConf(1)
                << "\t" << node.second->getConf(2)
                << endl;
    }

}


/**
 * @brief Computes the residue vector for the static thermal equilibrium problem.
 * @param residue Output residue vector (K*q + f_int - f_ext).
 * @param q Current temperature state vector.
 */
void Simulation::staticThermalResidue(lmx::Vector<data_type>& residue,
                                      lmx::Vector<data_type>& q
                                     )
{
    for (auto& node : thermalNodes)
    {
        node.second->setqt(q);
    }

    globalConductivity.reset();
    globalExternalHeat.reset();
    globalInternalHeat.reset();

    baseSystem->calcConductivityMatrix();
    baseSystem->calcExternalHeat();
    baseSystem->calcInternalHeat();
    baseSystem->assembleConductivityMatrix(globalConductivity);
    baseSystem->assembleExternalHeat(globalExternalHeat);
    baseSystem->assembleInternalHeat(globalInternalHeat);

    residue = globalConductivity * q;
    residue += globalInternalHeat;
    residue -= globalExternalHeat;

//   cout << endl << "RESIDUE PARTS: " << (globalConductivity*q).norm2()
//     << " " << globalInternalHeat.norm2() << " " << globalExternalHeat.norm2() << endl;

}

/**
 * @brief Computes the tangent (Jacobian) matrix for the static thermal problem.
 * @param tangent_in Output tangent matrix.
 * @param q Current temperature state vector.
 */
void Simulation::staticThermalTangent(lmx::Matrix<data_type>& tangent_in,
                                      lmx::Vector<data_type>& q
                                     )
{
    tangent_in.reset();
    baseSystem->calcThermalTangentMatrix();
    baseSystem->assembleThermalTangentMatrix(tangent_in);
//     cout << tangent_in << endl;
    tangent_in += globalConductivity;
//     cout << tangent_in << endl;
}

/**
 * @brief Checks whether the static thermal solver has converged and handles post-convergence output.
 * @param res Current residue vector.
 * @param q Current temperature state vector.
 * @return True if the step has converged, false otherwise.
 */
bool Simulation::staticThermalConvergence(lmx::Vector<data_type>& res,
        lmx::Vector<data_type>& q
                                         )
{
//   lmx::Vector<data_type> res( qddot.size() );
//   res =  globalInternalForces - globalExternalForces;
    stepConverged = (res.norm2() <= epsilon); 
    if (stepConverged) // if convergence...
    {
        stepConverged = baseSystem->checkAugmented();
        if (stepConverged)  //... and augmented system is also converged
        {
            stepTime = 1.;
            systemOuputStep(q);
            baseSystem->clearAugmented();
            stepTriggered();
        }
    }
    updateMaterials(stepConverged);
    return stepConverged;
}


/**
 * @brief Evaluates the thermal time derivative for explicit time integration.
 * @param qt Current temperature state vector.
 * @param qtdot Output temperature rate vector (solved from C*qtdot = -(K*qt + f_int - f_ext)).
 * @param time Current simulation time.
 */
void Simulation::explicitThermalEvaluation
(const lmx::Vector<data_type>& qt, lmx::Vector<data_type>& qtdot, double time
)
{
    for (auto& node : thermalNodes)
    {
        node.second->setqt(qt);
    }

//     globalConductivity.reset();
//     globalCapacity.reset();
    globalExternalHeat.reset();
    globalInternalHeat.reset();

//     baseSystem->calcConductivityMatrix();
//     baseSystem->calcCapacityMatrix();
    baseSystem->calcExternalHeat();
    baseSystem->calcInternalHeat();
//     baseSystem->assembleConductivityMatrix(globalConductivity);
//     baseSystem->assembleCapacityMatrix(globalCapacity);
    baseSystem->assembleExternalHeat(globalExternalHeat);
    baseSystem->assembleInternalHeat(globalInternalHeat);
    globalRHSHeat = globalConductivity * qt;
    globalRHSHeat += globalInternalHeat;
    globalRHSHeat -= globalExternalHeat;

//   cout << globalMass << endl;
    lmx::LinearSystem<data_type> theLSolver(globalCapacity, qtdot, globalRHSHeat);
    theLSolver.solveYourself();
//  cout << "initial_flux :" << qtdot << endl;
    if (theAnalysis->type() != "THERMOMECHANICALDYNAMIC")   // for regular THERMAL dynamic problems
    {
        stepTime = time;
        systemOuputStep(qt);
    }


}

/**
 * @brief Evaluates the thermal time derivative for implicit dynamic time integration,
 *        recomputing the conductivity and capacity matrices at each call.
 * @param qt Current temperature state vector.
 * @param qtdot Output temperature rate vector.
 * @param time Current simulation time.
 */
void Simulation::dynamicThermalEvaluation(const lmx::Vector<data_type>& qt,
        lmx::Vector<data_type>& qtdot,
        double time
                                         )
{
    globalCapacity.reset();
    globalConductivity.reset();
    globalExternalHeat.reset();
    globalInternalHeat.reset();

    baseSystem->calcConductivityMatrix();
    baseSystem->calcCapacityMatrix();
    baseSystem->calcExternalHeat();
    baseSystem->calcInternalHeat();
    baseSystem->assembleCapacityMatrix(globalCapacity);
    baseSystem->assembleConductivityMatrix(globalConductivity);
    baseSystem->assembleExternalHeat(globalExternalHeat);
    baseSystem->assembleInternalHeat(globalInternalHeat);
    globalRHSHeat = globalConductivity * qt;
    globalRHSHeat += globalInternalHeat;
    globalRHSHeat -= globalExternalHeat;

//     cout << "H = " << globalConductivity << endl;
//     cout << "C = " << globalCapacity << endl;
//     cout << globalRHSHeat << endl;
    lmx::LinearSystem<data_type> theLSolver(globalCapacity, qtdot, globalRHSHeat);
    theLSolver.solveYourself();
//    cout << "initial_flux :" << qtdot << endl;

    stepTime = time;
}

/**
 * @brief Computes the residue vector for the dynamic thermal problem (C*qdot + K*q + f_int - f_ext).
 * @param residue Output residue vector.
 * @param q Current temperature state vector.
 * @param qdot Current temperature rate vector.
 */
void Simulation::dynamicThermalResidue(lmx::Vector<data_type>& residue,
                                       const lmx::Vector<data_type>& q,
                                       const lmx::Vector<data_type>& qdot,
                                       double /*time*/
                                      )
{
    for (auto& node : thermalNodes)
    {
        node.second->setqt(q);
    }

    globalCapacity.reset();
    globalConductivity.reset();
    globalExternalHeat.reset();
    globalInternalHeat.reset();

    baseSystem->calcConductivityMatrix();
    baseSystem->calcCapacityMatrix();
    baseSystem->calcExternalHeat();
    baseSystem->calcInternalHeat();

    baseSystem->assembleCapacityMatrix(globalCapacity);
    baseSystem->assembleConductivityMatrix(globalConductivity);
    baseSystem->assembleExternalHeat(globalExternalHeat);
    baseSystem->assembleInternalHeat(globalInternalHeat);

    residue = globalCapacity * qdot;
    residue += globalConductivity * q;
    residue += globalInternalHeat;
    residue -= globalExternalHeat;

//   cout << endl << "RESIDUE PARTS: " << (globalCapacity*qdot).norm2() << " "
//   << (globalConductivity*q).norm2() << " " << globalExternalHeat.norm2() << endl;
//     cout << "H = " << globalConductivity << endl;
//     cout << "C = " << globalCapacity << endl;
//     cout << globalExternalHeat << endl;
//     cout << "q = " << q << endl;
//     cout << "qdot = " << qdot << endl;
//     cout << "residue = " << residue << endl;
//     cout << "globalRHSHeat.norm1 = " << globalRHSHeat.norm1() << endl;
//     cout << "q = " << q << endl;

}

/**
 * @brief Computes the tangent matrix for the dynamic thermal problem
 *        as partial_qdot * C + K.
 * @param tangent_in Output tangent matrix.
 * @param q Current temperature state vector.
 * @param partial_qdot Partial derivative of qdot with respect to the unknown (from the time integrator).
 */
void Simulation::dynamicThermalTangent(lmx::Matrix<data_type>& tangent_in,
                                       const lmx::Vector<data_type>& q,
                                       double partial_qdot,
                                       double /*time*/
                                      )
{
    tangent_in.reset();
//   baseSystem->calcTangentMatrix(  );
//   baseSystem->assembleTangentMatrix( tangent_in );
    tangent_in += (data_type)partial_qdot * globalCapacity;
    tangent_in += globalConductivity;
}

/**
 * @brief Checks convergence for the dynamic thermal solver and handles post-convergence output.
 * @param q Current temperature state vector.
 * @param qdot Current temperature rate vector.
 * @param time Current simulation time.
 * @return True if the step has converged, false otherwise.
 */
bool Simulation::dynamicThermalConvergence(const lmx::Vector<data_type>& q,
        const lmx::Vector<data_type>& qdot,
        double time
                                          )
{
    ++iterationsNLSolver;
    lmx::Vector<data_type> res(qdot.size());
//  double energy_max, energy_sum;
//  cout << "\n"
//       << (globalMass*qddot).norm2() << "\t"
//       << globalInternalForces.norm2() << "\t"
//       << globalExternalForces.norm2() << "\n";
    res = globalCapacity * qdot + globalConductivity * q + globalInternalHeat - globalExternalHeat;
//  energy_max = std::max( std::fabs(globalMass*qddot*q)
//                       , std::fabs(globalInternalForces*q) );
//  energy_max = std::max( energy_max, std::fabs(globalExternalForces*q) );
//  energy_sum = std::fabs(globalMass*qddot*q)
//             + std::fabs(globalInternalForces*q)
//             + std::fabs(globalExternalForces*q);
//      cout << "            : MAX_ENERGY = " << energy_max << endl
//           << "              SUM_ENERGY = " << energy_sum << endl;
    stepConverged = (res.norm2() <= epsilon); 
    if (stepConverged) // if convergence...
    {
//  if( (energy_max / energy_sum) <= epsilon ){
        stepConverged = baseSystem->checkAugmented();
        if (stepConverged)  //... and augmented system is also converged
        {
//      cout << " CONVERGENCE: MAX_ENERGY = " << energy_max << endl
//           << "              SUM_ENERGY = " << energy_sum << endl;
            stepTime = time;
            systemOuputStep(q);
//             baseSystem->clearAugmented();
        }
    }
    updateMaterials(stepConverged);
    return stepConverged;
}

/**
 * @brief Checks thermal convergence within a coupled thermo-mechanical dynamic analysis.
 * @param q Current temperature state vector.
 * @param qdot Current temperature rate vector.
 * @param time Current simulation time.
 * @return True if the thermal step has converged, false otherwise.
 */
bool Simulation::dynamicThermalConvergenceInThermomechanical(const lmx::Vector<data_type>& q,
        const lmx::Vector<data_type>& qdot,
        double time
                                                            )
{
    lmx::Vector<data_type> res(qdot.size());
    res = globalCapacity * qdot + globalConductivity * q + globalInternalHeat - globalExternalHeat;
    stepConverged = (res.norm2() <= epsilon); 
    if (stepConverged) // if convergence...
    {
        stepConverged = baseSystem->checkAugmented();
        if (stepConverged)  //... and augmented system is also converged
        {
            stepTime = time;
            baseSystem->clearAugmented();
        }
    }
    updateMaterials(stepConverged);
    return stepConverged;
}


/**
 * @brief Computes the nodal acceleration for explicit mechanical time integration
 *        by solving M*qddot = f_ext - f_int.
 * @param q Current displacement state vector.
 * @param qdot Current velocity state vector.
 * @param qddot Output acceleration vector.
 * @param time Current simulation time.
 */
void Simulation::explicitAcceleration(const lmx::Vector<data_type>& q,
                                      const lmx::Vector<data_type>& qdot,
                                      lmx::Vector<data_type>& qddot,
                                      double time
                                     )
{
    for (auto& node : nodes)
    {
        node.second->setqx(q, getDim());
    }

//   globalMass.reset();
    globalInternalForces.reset();
    globalExternalForces.reset();

//   baseSystem->calcMassMatrix();
    baseSystem->calcInternalForces();
    baseSystem->calcExternalForces();
//   baseSystem->assembleMassMatrix( globalMass );
    baseSystem->assembleInternalForces(globalInternalForces);
    baseSystem->assembleExternalForces(globalExternalForces);
    globalRHSForces = globalExternalForces;
    globalRHSForces -= globalInternalForces;

//   cout << globalMass << endl;
    lmx::LinearSystem<data_type> theLSolver(globalMass, qddot, globalRHSForces);
    theLSolver.solveYourself();

    stepTime = time;
    systemOuputStep(q, qdot);

}

/**
 * @brief Computes the nodal acceleration for implicit dynamic mechanical analysis,
 *        recomputing the mass matrix at each call.
 * @param q Current displacement state vector.
 * @param qdot Current velocity state vector.
 * @param qddot Output acceleration vector.
 */
void Simulation::dynamicAcceleration(const lmx::Vector<data_type>& q,
                                     const lmx::Vector<data_type>& qdot,
                                     lmx::Vector<data_type>& qddot,
                                     double /*time*/
                                    )
{
    globalMass.reset();
    globalInternalForces.reset();
    globalExternalForces.reset();

    for (auto& node : nodes)
    {
        node.second->setqx(q, getDim());
    }
    baseSystem->calcMassMatrix();
    baseSystem->calcInternalForces();
    baseSystem->calcExternalForces();
    baseSystem->assembleMassMatrix(globalMass);
    baseSystem->assembleInternalForces(globalInternalForces);
    baseSystem->assembleExternalForces(globalExternalForces);
    globalRHSForces = globalExternalForces;
    globalRHSForces -= globalInternalForces;

//    cout << globalMass << endl;
    lmx::LinearSystem<data_type> theLSolver(globalMass, qddot, globalRHSForces);
    theLSolver.solveYourself();
//    cout << "initial_acceleration :" << qddot << endl;
}

/**
 * @brief Computes the residue vector for the dynamic mechanical problem
 *        (M*qddot + f_int - f_ext).
 * @param residue Output residue vector.
 * @param q Current displacement state vector.
 * @param qdot Current velocity state vector.
 * @param qddot Current acceleration state vector.
 * @param time Current simulation time.
 */
void Simulation::dynamicResidue(lmx::Vector<data_type>& residue,
                                const lmx::Vector<data_type>& q,
                                const lmx::Vector<data_type>& qdot,
                                const lmx::Vector<data_type>& qddot,
                                double time
                               )
{
    for (auto& node : nodes)
    {
        node.second->setqx(q, getDim());
    }
// At this time globalMass is always the same...
// so the commented lines increment efficiency by 12% aprox.

//   globalMass.reset();
    globalInternalForces.reset();
    globalExternalForces.reset();

    baseSystem->update(time);
//   baseSystem->calcMassMatrix();
    baseSystem->calcInternalForces();
    baseSystem->calcExternalForces();
//   baseSystem->assembleMassMatrix( globalMass );
    baseSystem->assembleInternalForces(globalInternalForces);
    baseSystem->assembleExternalForces(globalExternalForces);

    residue = globalMass * qddot;
    residue += globalInternalForces;
    residue -= globalExternalForces;

//     cout << "qddot : " << qddot;
//     cout << "q : " << q;
//     cout << "globalMass*qddot : " << globalMass*qddot;
//     cout << "globalInternalForces : " << globalInternalForces;
//     cout << "globalExternalForces : " << globalExternalForces;
}

/**
 * @brief Computes the tangent matrix for the dynamic mechanical problem
 *        as K + partial_qddot * M.
 * @param tangent_in Output tangent matrix.
 * @param q Current displacement state vector.
 * @param qdot Current velocity state vector.
 * @param partial_qddot Partial derivative of qddot with respect to the unknown (from the time integrator).
 */
void Simulation::dynamicTangent(lmx::Matrix<data_type>& tangent_in,
                                const lmx::Vector<data_type>& q,
                                const lmx::Vector<data_type>& qdot,
                                double /*partial_qdot*/,
                                double partial_qddot,
                                double /*time*/
                               )
{
    tangent_in.reset();
    baseSystem->calcTangentMatrix();
    baseSystem->assembleTangentMatrix(tangent_in);
    tangent_in += (data_type)partial_qddot * globalMass;
}

/**
 * @brief Checks convergence for the dynamic mechanical solver and handles post-convergence output.
 * @param q Current displacement state vector.
 * @param qdot Current velocity state vector.
 * @param qddot Current acceleration state vector.
 * @param time Current simulation time.
 * @return True if the step has converged, false otherwise.
 */
bool Simulation::dynamicConvergence(const lmx::Vector<data_type>& q,
                                    const lmx::Vector<data_type>& qdot,
                                    const lmx::Vector<data_type>& qddot,
                                    double time
                                   )
{
    ++iterationsNLSolver;
    lmx::Vector<data_type> res(qddot.size());
//  double energy_max, energy_sum;
//  cout << "\n"
//       << "(globalMass*qddot).norm2() = " << "\t"
//       << (globalMass*qddot).norm2() << "\n"
//       << "globalInternalForces.norm2() = " << "\t"
//       << globalInternalForces.norm2() << "\n"
//       << "globalExternalForces.norm2() = " << "\t"
//       << globalExternalForces.norm2() << "\n";
    res = globalMass * qddot + globalInternalForces - globalExternalForces;
//  energy_max = std::max( std::fabs(globalMass*qddot*q)
//                       , std::fabs(globalInternalForces*q) );
//  energy_max = std::max( energy_max, std::fabs(globalExternalForces*q) );
//  energy_sum = std::fabs(globalMass*qddot*q)
//             + std::fabs(globalInternalForces*q)
//             + std::fabs(globalExternalForces*q);
//      cout << "            : MAX_ENERGY = " << energy_max << endl
//           << "              SUM_ENERGY = " << energy_sum << endl;
    stepConverged = (res.norm2() <= epsilon); 
    if (stepConverged) // if convergence...
    {
//  if( (energy_max / energy_sum) <= epsilon ){
        stepConverged = baseSystem->checkAugmented();
        if (stepConverged)  //... and augmented system is also converged
        {
//      cout << " CONVERGENCE: MAX_ENERGY = " << energy_max << endl
//           << "              SUM_ENERGY = " << energy_sum << endl;
            stepTime = time;
            systemOuputStep(q, qdot);
            baseSystem->clearAugmented();
        }
    }
    updateMaterials(stepConverged);
    return stepConverged;

}


/**
 * @brief Computes the residue vector for the static mechanical equilibrium problem
 *        (f_int - f_ext).
 * @param residue Output residue vector.
 * @param q Current displacement state vector.
 */
void Simulation::staticResidue(lmx::Vector<data_type>& residue,
                               lmx::Vector<data_type>& q
                              )
{
    for (auto& node : nodes)
    {
        node.second->setqx(q, getDim());
    }

    globalInternalForces.reset();
    globalExternalForces.reset();

    baseSystem->calcInternalForces();
    baseSystem->calcExternalForces();
    baseSystem->assembleInternalForces(globalInternalForces);
    baseSystem->assembleExternalForces(globalExternalForces);

    residue = globalInternalForces;
    residue -= globalExternalForces;

//   cout << "residue:" << residue;
//
//   cout << "q : " << q;
//   cout << "globalInternalForces : " << globalInternalForces;
//   cout << "globalExternalForces : " << globalExternalForces;
}

/**
 * @brief Computes the tangent (stiffness) matrix for the static mechanical problem.
 * @param tangent_in Output tangent matrix.
 * @param q Current displacement state vector.
 */
void Simulation::staticTangent(lmx::Matrix<data_type>& tangent_in,
                               lmx::Vector<data_type>& q
                              )
{
    tangent_in.reset();
    baseSystem->calcTangentMatrix();
    baseSystem->assembleTangentMatrix(tangent_in);
//  cout << "TANGENT:\n" << tangent_in;
}

/**
 * @brief Checks convergence for the static mechanical solver and handles post-convergence output.
 * @param res Current residue vector.
 * @param q Current displacement state vector.
 * @return True if the step has converged, false otherwise.
 */
bool Simulation::staticConvergence(lmx::Vector<data_type>& res,
                                   lmx::Vector<data_type>& q
                                  )
{
//   lmx::Vector<data_type> res( qddot.size() );
//   res =  globalInternalForces - globalExternalForces;
    stepConverged = (res.norm2() <= epsilon); 
    if (stepConverged) // if convergence...
    {
        stepConverged = baseSystem->checkAugmented();
        if (stepConverged)  //... and augmented system is also converged
        {
            stepTime = 1.;
            systemOuputStep(q);
// 	    this->storeTimeConfiguration(q);
            baseSystem->clearAugmented();
            stepTriggered();
            return 1;
        }
    }
    updateMaterials(stepConverged);
    return stepConverged;
}


/**
 * @brief Callback invoked after a successful step convergence; writes configuration and updates timing output.
 */
void Simulation::stepTriggered()
{
#ifdef HAVE_VTK
    if(contact=="GLOBAL" || visualization==1)
    {
        this->theContact->updatePoints();
        this->theContact->updateLines();
        if(contact=="GLOBAL")
            this->theContact->updateDelaunay();
        this->theContact->drawObjects();
    }
#endif
// Output configuration of time step:
    writeConfStep();

    // Output timer info:
    if (outputFilesDetail > 0)
    {
        double theTime = globalTimer->getTime();
        *timerFile << stepTime << "\t"
                   << theTime - oldClockTime << "\t"
                   << theTime << "\t"
                   << iterationsNLSolver << std::endl;
        oldClockTime = theTime;
        iterationsNLSolver = 0;
    }
}

/**
 * @brief Writes the current nodal configuration at the current step time to the displacement output file.
 */
void Simulation::writeConfStep()
{
    if (outputFilesDetail > 1)
    {
        configurationFile->setf(std::ios::scientific, std::ios::floatfield);
        configurationFile->precision(6);
        *configurationFile << stepTime << "\t";
        int i;
        for (auto& point : outputPoints)
        {
            for (i = 0; i < dimension; ++i)
            {
                *configurationFile << point.second->getConf(i) << "\t";
            }
        }
        *configurationFile << endl;
    }
//    for (auto& node : nodes) {
// 	    cout << endl;
//             cout << it_nodes->second->getqx(0) << " "
// 		 << it_nodes->second->getqx(1) << " ";
//             if (Simulation::getDim() == 3)
//               cout << it_nodes->second->getqx(2) << " ";
// 	    cout << endl;
//    }

}

/**
 * @brief Writes the system state output for the current time step using the displacement vector.
 * @param q Current state vector (displacements or temperatures).
 */
void Simulation::systemOuputStep(const lmx::Vector<data_type>& q)
{
    if (outputFilesDetail > 1)
    {
        baseSystem->outputStep(q);
    }

}

/**
 * @brief Writes the system state output for the current time step using both displacement and velocity vectors.
 * @param q Current displacement state vector.
 * @param qdot Current velocity state vector.
 */
void Simulation::systemOuputStep(const lmx::Vector<data_type>& q, const lmx::Vector<data_type>& qdot)
{
    if (outputFilesDetail > 1)
    {
        baseSystem->outputStep(q, qdot);
    }
}

/**
 * @brief Returns the names of all flexible bodies in the base system.
 * @return Vector of flexible body name strings.
 */
std::vector<std::string> Simulation::bodyNames()
{
    return baseSystem->flexBodyNames();
}

/**
 * @brief Returns the current nodal coordinates of all nodes belonging to a named body in a subsystem.
 * @param system_name Name of the subsystem containing the body.
 * @param name Name of the body.
 * @return Flat vector of (x, y, z) coordinate triples for each node.
 */
std::vector<double> Simulation::bodyPoints(const std::string& system_name, const std::string& name) const
{
    std::vector<double> points;

    auto body = baseSystem->getBody(system_name, name);

    for (const auto& node : body->getNodes())
    {
        points.push_back(node->getConf(0));
        points.push_back(node->getConf(1));
        points.push_back(node->getConf(2));
    }

    return points;
}

}

