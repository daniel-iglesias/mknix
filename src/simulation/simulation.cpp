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

#include "simulation.h"
#include "analysisdynamic.h"

#include <system/system.h>
#include <system/constraintthermal.h>
#include <core/node.h>
#include <reader/reader.h>
#include <system/generalcontact.h>

#include "gpu/cpu_run_type.h"
#include "gpu/cpu_solvers.h"

namespace mknix {


// Static variables:
double Simulation::stepTime = 0;
double Simulation::oldClockTime = 0;
lmx::Vector<double> Simulation::gravity = lmx::Vector<double>(3);
double Simulation::alpha = 1E4;
int Simulation::dimension = 2;
std::string Simulation::contact = "NONE";
bool Simulation::visualization = 0;
bool Simulation::outputMatrices = 0;
std::string Simulation::constraintMethod = "PENALTY";
double Simulation::epsilon = 1E-5;
std::string Simulation::smoothingType = "GLOBAL";

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
    if (outputFilesDetail > 0) {
        timerFile = new std::ofstream("simulation_times.dat");
        if (outputFilesDetail > 1) {
            configurationFile = new std::ofstream("dis.dat");
        }
    }
}

Simulation::~Simulation()
{
    if (globalTimer) delete globalTimer;
    if (timerFile) delete timerFile;
    if (configurationFile) delete configurationFile;
}


void Simulation::inputFromFile(const std::string& FileIn)
{
    auto reader = make_unique<Reader>(this);
    if (!baseSystem) {
        baseSystem = make_unique<System>("baseSystem");
    }
    reader->inputFromFile(FileIn);
}

int Simulation::getInterfaceNumberOfNodes(string name)
{
    return baseSystem->getSystem(name)->getNumberOfNodes();
}


Node* Simulation::getInterfaceNode(std::string name, int num ){
    return baseSystem->getSystem(name)->getNode(num);
}

std::vector<double> Simulation::getInterfaceNodesCoords()
{
    // Loads have access to the nodes, and are part of the system.
    std::vector<double> temp_x_coordinates;
    baseSystem->getThermalNodes(temp_x_coordinates);

    return temp_x_coordinates;
}

double Simulation::getConstraintOutput( std::string constraintName, std::string systemName, int component)
{
    return -( baseSystem->getSystem(systemName)->getConstraintThermal(constraintName)->getInternalForces().readElement(component) );
}


void Simulation::setInitialTemperatures(double temp_in)
{
    initialTemperature = temp_in;
}

// Part copy of run(), limited to preparation and thermal dynamic analysis
void Simulation::init(int vervosity)
{
    if (outputFilesDetail > 1) {
        writeSystem();
        outFile << "ANALYSIS " << endl;
        outFile << "CONFIGURATION" << endl;
    }
    if (outputFilesDetail > 0) {
        *timerFile << stepTime << "\t" << globalTimer->getTime() << std::endl;
    }

    for (auto& analysis : analyses) {
        if (analysis->type() == "THERMAL" || analysis->type() == "THERMALSTATIC") {
            initThermalSimulation(analysis.get(), vervosity);
        }
        else {
            initMechanicalSimulation(analysis.get());
        }
    }
}

lmx::Vector<data_type> Simulation::initThermalSimulation(Analysis * theAnalysis_in, int vervosity, bool init)
{

std::cout << "\n\n  --- void Simulation::initThermalSimulation --- \n\n" << std::endl;

    theAnalysis = theAnalysis_in;
    auto gdlSize = nodes.size();
    lmx::Vector<data_type> q(gdlSize);
    q.fillIdentity(initialTemperature);

if(OLD_CODE){
    globalConductivity.resize(gdlSize, gdlSize);
    globalCapacity.resize(gdlSize, gdlSize);
    globalRHSHeat.resize(gdlSize);
    globalExternalHeat.resize(gdlSize);
    globalInternalHeat.resize(gdlSize);

    baseSystem->calcConductivityMatrix();
    baseSystem->assembleConductivityMatrix(globalConductivity);
    baseSystem->calcCapacityMatrix();
    baseSystem->assembleCapacityMatrix(globalCapacity);
} else {

  globalConductivity.resize(gdlSize, gdlSize);
  globalCapacity.resize(gdlSize, gdlSize);
  globalRHSHeat.resize(gdlSize);
  globalExternalHeat.resize(gdlSize);
  globalInternalHeat.resize(gdlSize);

  baseSystem->setTemperatureVector(q);
  baseSystem->setMaterialTable(myMaterialTable);//setting soa structture
  baseSystem->calcConductivityMatrix();
  baseSystem->assembleConductivityMatrix(globalConductivity);
  baseSystem->calcCapacityMatrix();
  baseSystem->assembleCapacityMatrix(globalCapacity);
}
    writeConfStep();

    if (outputFilesDetail > 1 && theAnalysis->type() == "THERMAL") {
        systemOuputStep(q);
    }

    if (init) {
        theAnalysis->init(&q, vervosity);
    }

    return q;
}

lmx::Vector<data_type> Simulation::initMechanicalSimulation(Analysis * analysis, bool init)
{
    theAnalysis = analysis;
    auto gdlSize = nodes.size() * Simulation::dimension;
    lmx::Vector<data_type> q(gdlSize);

    globalMass.resize(gdlSize, gdlSize);
    globalRHSForces.resize(gdlSize);
    globalInternalForces.resize(gdlSize);
    globalExternalForces.resize(gdlSize);

    auto i = 0u;
    for (auto& node : nodes) {
        q(Simulation::dimension * i) = node.second->getqx(0);
        q(Simulation::dimension * i + 1) = node.second->getqx(1);
        if (Simulation::getDim() == 3) {
            q(Simulation::dimension * i + 2) = node.second->getqx(2);
        }
        ++i;
    }

    baseSystem->calcMassMatrix();
    baseSystem->assembleMassMatrix(globalMass);
    // maybe is better to make an specific function call for the sparse
    // pattern, but this should work...
    if (lmx::getMatrixType() == 1) {
        globalSparsePattern.resize(gdlSize, gdlSize);
        baseSystem->calcTangentMatrix();
        baseSystem->assembleTangentMatrix(globalSparsePattern);
    }

    // Output matrices in initial configuration:
    if (outputMatrices) {
        lmx::Matrix<data_type> K_temp(gdlSize, gdlSize);
        baseSystem->calcTangentMatrix();
        baseSystem->assembleTangentMatrix(K_temp);
        K_temp.harwellBoeingSave((char *) "K.mat");
        globalMass.harwellBoeingSave((char *) "M.mat");
        // save raw matrices:
        std::ofstream Kfile("K");
        Kfile << K_temp;
        Kfile.close();
    }
    // Output to file the initial configuration:
    writeConfStep();

    if (outputFilesDetail > 1 && theAnalysis->type() == "DYNAMIC") {
        systemOuputStep(q);
    }

    return q;
}

void setSignal(std::string node, std::vector<double>)
{
    return;
}


std::vector<double> getSignal(std::string node)
{
    return { };
}

void Simulation::solveStep()
{
    theAnalysis->nextStep();
}

void Simulation::solveStep(double * signal, double * outputSignal)
{
    baseSystem->updateThermalLoads(signal);
    theAnalysis->nextStep();
    if (outputSignal) {
        baseSystem->getOutputSignalThermal(outputSignal);
    }
}

void Simulation::endSimulation()
{
    configurationFile->close();

    analyses.clear();

    if (outputFilesDetail > 1) {
        std::ifstream disp("dis.dat");

        char a;
//       std::string aa;

        if (disp.is_open()) {
            if (outFile.is_open()) {
// 	  while (disp >> aa) {
// 	    outFile << aa;
// 	  }
                while (disp.get(a)) {
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


void Simulation::run()
{
#ifdef HAVE_VTK
    if(Simulation::contact == "GLOBAL" || Simulation::visualization == 1) {
        this->theContact = new Contact(this, 10.);
        this->theContact->createPoints();
        this->theContact->createPolys();
        this->theContact->createDrawingObjects();
        if(Simulation::contact == "GLOBAL") {
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

    for (auto& analysis : analyses) {
        if (analysis->type() == "THERMAL" || analysis->type() == "THERMALSTATIC") {
          std::cout << "About to run this->runThermalAnalysis(analysis.get());"<< std::endl;
            this->runThermalAnalysis(analysis.get());
        }
        else {
            if (analysis->type() == "STATIC" || analysis->type() == "DYNAMIC") {
                this->baseSystem->setMechanical();
            }
            this->runMechanicalAnalysis(analysis.get());
        }
    }
}

void Simulation::runThermalAnalysis(Analysis * theAnalysis_in)
{
    auto q = initThermalSimulation(theAnalysis_in, 2, false);

    if (theAnalysis->type() == "THERMAL") {
        theAnalysis->solve(&q);
    }

    else if (theAnalysis->type() == "THERMALSTATIC") {
        // write initial configuration...
        outFile << "0 "; //time=0
//           systemOuputStep( q ); // produces output of temperatures, incompatible with mknixpost-static
        for (auto& point : outputPoints) {
            outFile << point.second->getConf(0) << " ";
            outFile << point.second->getConf(1) << " ";
            if (Simulation::getDim() == 3) {
                outFile << point.second->getConf(2) << " ";
            }
        }
        outFile << endl;
        theAnalysis->solve(&q);
    }

    if (outputFilesDetail > 1) {
        std::ifstream disp("dis.dat");
        char a;

        if (outFile.is_open()) {
            while (disp.get(a)) {
                outFile.put(a);
                outFile.flush();
            }
            outFile << "ENDCONFIGURATION" << endl;

            // output extra flexible bodies data...
            baseSystem->outputToFile(&outFile);
        }
    }

    // output f_int of constraints...
    lmx::Vector<double> constr_forces(nodes.size() * dimension);
    baseSystem->assembleConstraintForces(constr_forces);
    for (size_type i = 0; i < constr_forces.size(); ++i) {
//         if(constr_forces(i) != 0. )
//             cout << "R(" << i << ") = " << constr_forces(i) << endl;
    }
}

void Simulation::runMechanicalAnalysis(Analysis * theAnalysis_in)
{
    auto q = initMechanicalSimulation(theAnalysis_in, false);

    auto gdlSize = nodes.size() * Simulation::dimension;

    if (theAnalysis->type() == "DYNAMIC") {
        lmx::Vector<data_type> qdot(gdlSize);
        // output first step data
        theAnalysis->setEpsilon(epsilon);
        theAnalysis->solve(&q, &qdot);
    }
    else if (theAnalysis->type() == "STATIC") {
        theAnalysis->solve(&q);
    }
    else if (theAnalysis->type() == "THERMOMECHANICALDYNAMIC") {
        // Init thermal conf vector and a zero velocity vector
        auto thermalSize = thermalNodes.size();
        lmx::Vector<data_type> qdot(gdlSize);
        lmx::Vector<data_type> qt(thermalSize);
        auto i = 0u;
        for (auto& node : thermalNodes) {
            qt(i) = node.second->getqt();
            ++i;
        }
        globalCapacity.resize(thermalSize, thermalSize);
        globalConductivity.resize(thermalSize, thermalSize);
        globalRHSHeat.resize(thermalSize);
        globalExternalHeat.resize(thermalSize);
        globalInternalHeat.resize(thermalSize);
std::cout << "\n\n  --- void Simulation::runMechanicalAnalysis with THERMOMECHANICALDYNAMIC --- \n\n" << std::endl;
        //baseSystem->calcFactors();
        baseSystem->calcConductivityMatrix();
        baseSystem->assembleConductivityMatrix(globalConductivity);
        baseSystem->calcCapacityMatrix();
        baseSystem->assembleCapacityMatrix(globalCapacity);

        // output first step data
        systemOuputStep(q);
        theAnalysis->setEpsilon(epsilon);
        theAnalysis->solve(&qt, &q, &qdot);
    }

    if (outputFilesDetail > 1) {
        std::ifstream disp("dis.dat");
        char a;

        if (outFile.is_open()) {
            while (disp.get(a)) {
                outFile.put(a);
                outFile.flush();
            }
            outFile << "ENDCONFIGURATION" << endl;

            // output extra flexible bodies data...
            baseSystem->outputToFile(&outFile);
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

void Simulation::writeSystem()

{
    std::stringstream ss;
    ss << title << ".mec";
    auto outFileName = ss.str();
    outFile.open(outFileName.c_str(), std::ofstream::out);
    outFile << "DIMENSION " << Simulation::dimension << endl;

    outFile << "SYSTEM" << endl;

    outFile << "NODES" << endl;
    for (auto& point : outputPoints) {
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
    for (auto& node : nodes) {
        nodeFile
        << "\t" << node.first
        << "\t" << node.second->getConf(0)
        << "\t" << node.second->getConf(1)
        << "\t" << node.second->getConf(2)
        << endl;
    }

}

void Simulation::staticThermalResidue(lmx::Vector<data_type>& residue,
                                      lmx::Vector<data_type>& q
)
{
    for (auto& node : thermalNodes) {
        node.second->setqt(q);
    }
std::cout << "\n\n  --- void Simulation::staticThermalResidue --- \n\n" << std::endl;
    globalConductivity.reset();
    globalExternalHeat.reset();
    globalInternalHeat.reset();

    baseSystem->setTemperatureVector(q);

    baseSystem->calcFactors();
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

void Simulation::staticThermalTangent(lmx::Matrix<data_type>& tangent_in,
                                      lmx::Vector<data_type>& q
)
{
      std::cout << "\n\n  --- void Simulation::staticThermalTangent --- \n\n" << std::endl;
    tangent_in.reset();
    baseSystem->calcThermalTangentMatrix();
    baseSystem->assembleThermalTangentMatrix(tangent_in);
//     cout << tangent_in << endl;
    tangent_in += globalConductivity;
//     cout << tangent_in << endl;
}

bool Simulation::staticThermalConvergence(lmx::Vector<data_type>& res,
                                          lmx::Vector<data_type>& q
)
{
    std::cout << "\n\n  --- void Simulation::staticThermalConvergence --- \n\n" << std::endl;
//   lmx::Vector<data_type> res( qddot.size() );
//   res =  globalInternalForces - globalExternalForces;
    if (res.norm2() <= epsilon) {
        if (baseSystem->checkAugmented()) { // if convergence...
            stepTime = 1.;
            systemOuputStep(q);
            baseSystem->clearAugmented();
            stepTriggered();
            return 1;
        }
        else { return 0; }
    }
    else { return 0; }

}


void Simulation::explicitThermalEvaluation
        (const lmx::Vector<data_type>& qt, lmx::Vector<data_type>& qtdot, double time
        )
{
    for (auto& node : thermalNodes) {
        node.second->setqt(qt);
    }
std::cout << "\n\n  --- void Simulation::explicitThermalEvaluation --- \n\n" << std::endl;
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
    if (theAnalysis->type() != "THERMOMECHANICALDYNAMIC") { // for regular THERMAL dynamic problems
        stepTime = time;
        systemOuputStep(qt);
    }


}

void Simulation::dynamicThermalEvaluation(const lmx::Vector<data_type>& qt,
                                          lmx::Vector<data_type>& qtdot,
                                          double time
)
{
 //std::cout << "\n\n  --- void Simulation::dynamicThermalEvaluation --- \n\n" << std::endl;
  if(OLD_CODE){
    globalCapacity.reset();//TODO remove in new SOA versions
    globalConductivity.reset();
    globalExternalHeat.reset();
    globalInternalHeat.reset();

    //baseSystem->calcFactors();
    baseSystem->calcConductivityMatrix();
    baseSystem->calcCapacityMatrix();
    baseSystem->calcExternalHeat();
    baseSystem->calcInternalHeat();
    baseSystem->assembleCapacityMatrix(globalCapacity);
    baseSystem->assembleConductivityMatrix(globalConductivity);
    baseSystem->assembleExternalHeat(globalExternalHeat);
    baseSystem->assembleInternalHeat(globalInternalHeat);
    globalRHSHeat = globalConductivity * qt;

  //  std::cout << std::endl << " globalRHSHeat.sumSum() = "<< globalRHSHeat.sumSum() << std::endl;
    globalRHSHeat += globalInternalHeat;
    globalRHSHeat -= globalExternalHeat;

  /*  double debugio = globalRHSHeat.sumSum();
    std::cout << std::endl << " globalRHSHeat.sumSum() = "<< debugio << std::endl;

    std::cout << std::endl << " globalCapacity.sumSum() = "<< globalCapacity.sumSum() << std::endl;
    std::cout << std::endl << " globalCapacity.trace() = "<< globalCapacity.trace() << std::endl;
    std::cout << std::endl << " globalConductivity.SumSum() = "<< globalConductivity.sumSum() << std::endl;
    std::cout << std::endl << " globalConductivity.trace() = "<< globalConductivity.trace() << std::endl;*/

    // cout << "H = " << globalConductivity << endl;
    // cout << "C = " << globalCapacity << endl;
    // cout << globalRHSHeat << endl;
    lmx::LinearSystem<data_type> theLSolver(globalCapacity, qtdot, globalRHSHeat);
    theLSolver.solveYourself();
//    cout << "initial_flux :" << qtdot << endl;
}else{
//reset not necesary now
//std::cout << "0. before resets" << std::endl;
globalCapacity.reset();
globalConductivity.reset();
globalExternalHeat.reset();
globalInternalHeat.reset();

//  std::cout << "1. before calcFactors" << std::endl;
//  baseSystem->calcFactors();
  //std::cout << "2. before calcConductivityMatrix" << std::endl;
  baseSystem->calcConductivityMatrix();
  //std::cout << "3. before calcCapacityMatrix" << std::endl;
  baseSystem->calcCapacityMatrix();
  //std::cout << "4. before calcExternalHeat" << std::endl;
  baseSystem->calcExternalHeat();
  //std::cout << "5. before calcInternalHeat" << std::endl;
  baseSystem->calcInternalHeat();
  //std::cout << "6. before assembleCapacityMatrix" << std::endl;
  baseSystem->assembleCapacityMatrix(globalCapacity);
  //std::cout << "7. before assembleConductivityMatrix" << std::endl;
  baseSystem->assembleConductivityMatrix(globalConductivity);
  //std::cout << "8. before assembleExternalHeat" << std::endl;
  baseSystem->assembleExternalHeat(globalExternalHeat);
  //std::cout << "9. before assembleInternalHeat" << std::endl;
  baseSystem->assembleInternalHeat(globalInternalHeat);
  //std::cout << "10. before globalRHSHeat = globalConductivity * qt;" << std::endl;
  globalRHSHeat = globalConductivity * qt;

  //double debugio0 = globalRHSHeat.sumSum();
  //std::cout << std::endl << " globalRHSHeat.sumSum() = "<< debugio0 << std::endl;

  //std::cout << "11. before globalRHSHeat += globalInternalHeat;" << std::endl;
  globalRHSHeat += globalInternalHeat;
  //std::cout << "12. before globalRHSHeat -= globalExternalHeat;" << std::endl;
  globalRHSHeat -= globalExternalHeat;

  //std::cout << std::endl << " globalRHSHeat.sumSum() = "<< globalRHSHeat.sumSum() << std::endl;

//     cout << "H = " << globalConductivity << endl;
//     cout << "C = " << globalCapacity << endl;
//     cout << globalRHSHeat << endl;
/*std::cout << std::endl << " globalCapacity.sumSum() = "<< globalCapacity.sumSum() << std::endl;
std::cout << std::endl << " globalCapacity.trace() = "<< globalCapacity.trace() << std::endl;
std::cout << std::endl << " globalConductivity.SumSum() = "<< globalConductivity.sumSum() << std::endl;
std::cout << std::endl << " globalConductivity.trace() = "<< globalConductivity.trace() << std::endl;*/

  //std::cout << "13. before lmx::LinearSystem<data_type> theLSolver(globalCapacity,qtdot, globalRHSHeat);" << std::endl;
  lmx::LinearSystem<data_type> theLSolver(globalCapacity,qtdot, globalRHSHeat);
  //std::cout << "14. before theLSolver.solveYourself();" << std::endl;
  theLSolver.solveYourself();
}
    stepTime = time;
}

void Simulation::dynamicThermalResidue(lmx::Vector<data_type>& residue,
                                       const lmx::Vector<data_type>& q,
                                       const lmx::Vector<data_type>& qdot,
                                       double /*time*/
)
{
    for (auto& node : thermalNodes) {//TODO:THIS ONE!
        node.second->setqt(q);
    }
std::cout << "\n\n  --- void Simulation::dynamicThermalResidue --- \n\n" << std::endl;

if(OLD_CODE){
    globalCapacity.reset();
    globalConductivity.reset();
    globalExternalHeat.reset();
    globalInternalHeat.reset();

    //baseSystem->calcFactors();
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
}else{
  //reste now is automatically done in the assembly functions

  //baseSystem->calcFactors();
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

}
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

void Simulation::dynamicThermalTangent(lmx::Matrix<data_type>& tangent_in,
                                       const lmx::Vector<data_type>& q,
                                       double partial_qdot,
                                       double /*time*/
)
{
  std::cout << "\n\n  --- void Simulation::dynamicThermalTangent --- \n\n" << std::endl;
    tangent_in.reset();
//   baseSystem->calcTangentMatrix(  );
//   baseSystem->assembleTangentMatrix( tangent_in );
if(OLD_CODE){
  cpuClock tanck;
  cpuTick(&tanck);
    tangent_in += (data_type) partial_qdot * globalCapacity;
    tangent_in += globalConductivity;
  cpuTock(&tanck, "Simulation::dynamicThermalTangent");
  }else{
    /*tangent_in += (data_type) partial_qdot * h_globalCapacity;
    tangent_in += h_globalConductivity;*/
  }
}

bool Simulation::dynamicThermalConvergence(const lmx::Vector<data_type>& q,
                                           const lmx::Vector<data_type>& qdot,
                                           double time
)
{
  std::cout << "\n\n  --- void Simulation::dynamicThermalConvergence --- \n\n" << std::endl;
    ++iterationsNLSolver;
    lmx::Vector<data_type> res(qdot.size());
//  double energy_max, energy_sum;
//  cout << "\n"
//       << (globalMass*qddot).norm2() << "\t"
//       << globalInternalForces.norm2() << "\t"
//       << globalExternalForces.norm2() << "\n";
if(OLD_CODE){
    res = globalCapacity * qdot + globalConductivity * q + globalInternalHeat - globalExternalHeat;
  } else {
  //  res = h_globalCapacity * qdot + h_globalConductivity * q + h_globalInternalHeat - h_globalExternalHeat;
  }
//  energy_max = std::max( std::fabs(globalMass*qddot*q)
//                       , std::fabs(globalInternalForces*q) );
//  energy_max = std::max( energy_max, std::fabs(globalExternalForces*q) );
//  energy_sum = std::fabs(globalMass*qddot*q)
//             + std::fabs(globalInternalForces*q)
//             + std::fabs(globalExternalForces*q);
//      cout << "            : MAX_ENERGY = " << energy_max << endl
//           << "              SUM_ENERGY = " << energy_sum << endl;
    if (res.norm2() <= epsilon) {
//  if( (energy_max / energy_sum) <= epsilon ){
        if (baseSystem->checkAugmented()) {
//      cout << " CONVERGENCE: MAX_ENERGY = " << energy_max << endl
//           << "              SUM_ENERGY = " << energy_sum << endl;
            stepTime = time;
            systemOuputStep(q);
//             baseSystem->clearAugmented();
            return 1;
        }
        else { return 0; }
    }
    else { return 0; }
}

bool Simulation::dynamicThermalConvergenceInThermomechanical(const lmx::Vector<data_type>& q,
                                                             const lmx::Vector<data_type>& qdot,
                                                             double time
)
{
  std::cout << "\n\n  --- void Simulation::dynamicThermalConvergenceInThermomechanical --- \n\n" << std::endl;
    lmx::Vector<data_type> res(qdot.size());
    res = globalCapacity * qdot + globalConductivity * q + globalInternalHeat - globalExternalHeat;
    if (res.norm2() <= epsilon) {
        if (baseSystem->checkAugmented()) {
            stepTime = time;
            baseSystem->clearAugmented();
            return 1;
        }
        else { return 0; }
    }
    else { return 0; }
}


void Simulation::explicitAcceleration(const lmx::Vector<data_type>& q,
                                      const lmx::Vector<data_type>& qdot,
                                      lmx::Vector<data_type>& qddot,
                                      double time
)
{
    std::cout << "\n\n  --- void Simulation::explicitAcceleration --- \n\n" << std::endl;
    for (auto& node : nodes) {
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

void Simulation::dynamicAcceleration(const lmx::Vector<data_type>& q,
                                     const lmx::Vector<data_type>& qdot,
                                     lmx::Vector<data_type>& qddot,
                                     double /*time*/
)
{
    globalMass.reset();
    globalInternalForces.reset();
    globalExternalForces.reset();

    for (auto& node : nodes) {
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

void Simulation::dynamicResidue(lmx::Vector<data_type>& residue,
                                const lmx::Vector<data_type>& q,
                                const lmx::Vector<data_type>& qdot,
                                const lmx::Vector<data_type>& qddot,
                                double time
)
{
    for (auto& node : nodes) {
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
    tangent_in += (data_type) partial_qddot * globalMass;
}

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
    if (res.norm2() <= epsilon) {
//  if( (energy_max / energy_sum) <= epsilon ){
        if (baseSystem->checkAugmented()) {
//      cout << " CONVERGENCE: MAX_ENERGY = " << energy_max << endl
//           << "              SUM_ENERGY = " << energy_sum << endl;
            stepTime = time;
            systemOuputStep(q, qdot);
            baseSystem->clearAugmented();
            return 1;
        }
        else { return 0; }
    }
    else { return 0; }

}


void Simulation::staticResidue(lmx::Vector<data_type>& residue,
                               lmx::Vector<data_type>& q
)
{
    for (auto& node : nodes) {
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

void Simulation::staticTangent(lmx::Matrix<data_type>& tangent_in,
                               lmx::Vector<data_type>& q
)
{
    tangent_in.reset();
    baseSystem->calcTangentMatrix();
    baseSystem->assembleTangentMatrix(tangent_in);
//  cout << "TANGENT:\n" << tangent_in;
}

bool Simulation::staticConvergence(lmx::Vector<data_type>& res,
                                   lmx::Vector<data_type>& q
)
{
//   lmx::Vector<data_type> res( qddot.size() );
//   res =  globalInternalForces - globalExternalForces;
    if (res.norm2() <= epsilon) {
        if (baseSystem->checkAugmented()) { // if convergence...
            stepTime = 1.;
            systemOuputStep(q);
// 	    this->storeTimeConfiguration(q);
            baseSystem->clearAugmented();
            stepTriggered();
            return 1;
        }
        else { return 0; }
    }
    else { return 0; }

}


void Simulation::stepTriggered()
{
#ifdef HAVE_VTK
    if(contact=="GLOBAL" || visualization==1) {
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
    if (outputFilesDetail > 0) {
        double theTime = globalTimer->getTime();
        *timerFile << stepTime << "\t"
        << theTime - oldClockTime << "\t"
        << theTime << "\t"
        << iterationsNLSolver << std::endl;
        oldClockTime = theTime;
        iterationsNLSolver = 0;
    }
}

void Simulation::writeConfStep()
{
    if (outputFilesDetail > 1) {
        configurationFile->setf(std::ios::scientific, std::ios::floatfield);
        configurationFile->precision(6);
        *configurationFile << stepTime << "\t";
        int i;
        for (auto& point : outputPoints) {
            for (i = 0; i < dimension; ++i) {
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

void Simulation::systemOuputStep(const lmx::Vector<data_type>& q)
{
    if (outputFilesDetail > 1) {
        baseSystem->outputStep(q);
    }

}

void Simulation::systemOuputStep(const lmx::Vector<data_type>& q, const lmx::Vector<data_type>& qdot)
{
    if (outputFilesDetail > 1) {
        baseSystem->outputStep(q, qdot);
    }
}

}
