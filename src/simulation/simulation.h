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

#ifndef MKNIXSIMULATION_H
#define MKNIXSIMULATION_H

#include "LMX/lmx.h"
#include "common.h"

//#include <Sparse> //eigen

#include "core/material.h"
#include "gpu/functions_cpu.h"
#include "gmm/gmm_matrix.h"
#include "gpu/cpu_run_type.h"

namespace mknix {

class Reader;

class Contact;

class System;

class Analysis;

class Node;

class Point;


/**
  @author AUTHORS <MAILS>
*/
class Simulation
{

    friend class Reader;

    friend class ReaderConstraints;

    friend class ReaderFlex;

    friend class ReaderRigid;

    friend class Contact;

    friend class SystemChain;

public:
    Simulation();

    ~Simulation();

    Simulation(const Simulation&) = delete;

    Simulation& operator=(const Simulation&) = delete;

    void inputFromFile(const std::string& fileIn);

//     void geometryFile(char*);
//     void outputFile(char*);
    int getInterfaceNumberOfNodes( std::string );
    Node* getInterfaceNode( std::string, int );
    double getConstraintOutput( std::string constraintName, std::string systemName="", int component=0 );

    std::vector<double> getInterfaceNodesCoords();

    void setOutputFilesDetail(int level_in) // 0 none, 1 only times, 2 all
    { outputFilesDetail = level_in; }

    void init(int vervosity=2);

    void setInitialTemperatures(double);

    void solveStep();

    void solveStep(double *, double * o_output = 0);

    void endSimulation();

    void run();

    void runThermalAnalysis(Analysis *);

    void runMechanicalAnalysis(Analysis *);

    void writeSystem();

    void staticThermalResidue(VectorX<data_type>& residue,
                              VectorX<data_type>& q
    );

    void staticThermalTangent(SparseMatrix<data_type>& tangent_in,
                              VectorX<data_type>& q
    );

    bool staticThermalConvergence(VectorX<data_type>& res,
                                  VectorX<data_type>& q
    );

    void explicitThermalEvaluation
            (const VectorX<data_type>& qt, VectorX<data_type>& qtdot, double time
            );

    void dynamicThermalEvaluation(const VectorX<data_type>& q,
                                  VectorX<data_type>& qdot,
                                  double time
    );

    void dynamicThermalResidue(VectorX<data_type>& residue,
                               const VectorX<data_type>& q,
                               const VectorX<data_type>& qdot,
                               double time
    );

    void dynamicThermalTangent(SparseMatrix<data_type>& tangent_in,
                               const VectorX<data_type>& q,
                               double partial_qdot,
                               double time
    );

    bool dynamicThermalConvergence(const VectorX<data_type>& q,
                                   const VectorX<data_type>& qdot,
                                   double time
    );

    bool dynamicThermalConvergenceInThermomechanical
            (const VectorX<data_type>& q,
             const VectorX<data_type>& qdot,
             double time
            );

    void explicitAcceleration(const VectorX<data_type>& q,
                              const VectorX<data_type>& qdot,
                              VectorX<data_type>& qddot,
                              double time
    );

    void dynamicAcceleration(const VectorX<data_type>& q,
                             const VectorX<data_type>& qdot,
                             VectorX<data_type>& qddot,
                             double time
    );

    void dynamicResidue(VectorX<data_type>& residue,
                        const VectorX<data_type>& q,
                        const VectorX<data_type>& qdot,
                        const VectorX<data_type>& qddot,
                        double time
    );

    void dynamicTangent(SparseMatrix<data_type>& tangent_in,
                        const VectorX<data_type>& q,
                        const VectorX<data_type>& qdot,
                        double partial_qdot,
                        double partial_qddot,
                        double time
    );

    bool dynamicConvergence(const VectorX<data_type>& q,
                            const VectorX<data_type>& qdot,
                            const VectorX<data_type>& qddot,
                            double time
    );

    void staticResidue(VectorX<data_type>& residue,
                       VectorX<data_type>& q
    );

    void staticTangent(SparseMatrix<data_type>& tangent_in,
                       VectorX<data_type>& q
    );

    bool staticConvergence(VectorX<data_type>& res,
                           VectorX<data_type>& q
    );

    void stepTriggered();

    void writeConfStep();

    lmx::DenseMatrix<data_type>& getSparsePattern()
    {
        return this->globalSparsePattern;
    }

    static double getGravity(int component)
    {
        return gravity.readElement(component);
    }

    static double getAlpha()
    {
        return alpha;
    }

    static double getTime()
    {
        return stepTime;
    }

    static int getDim()
    {
        return dimension;
    }

    static std::string getConstraintMethod()
    {
        return constraintMethod;
    }

    static std::string getSmoothingType()
    {
        return smoothingType;
    }

private:
    void storeTimeConfiguration(VectorX<data_type>& q);

private:
    void systemOuputStep(const VectorX<data_type>&, const VectorX<data_type>&);

    void systemOuputStep(const VectorX<data_type>&);

    std::string title;

    std::unique_ptr<System> baseSystem;
#ifdef HAVE_VTK
    Contact * theContact;
#endif
    Analysis * theAnalysis;

    std::vector<std::unique_ptr<Analysis>> analyses;
//     std::vector<std::vector< double > > pointsTimeConfiguration;
    std::map<int, Point *> outputPoints;
    std::map<int, Node *> nodes;
    std::map<int, Node *> thermalNodes;
    std::map<int, Material> materials;
    MaterialTable *myMaterialTable;//SOA implementation

    std::map<int, LoadThermalBoundary1D> thermalBoundaries;
    ThermalBoundaryTable *myThermalBoundary;//SOA implementation

    lmx::ExactStopwatch * globalTimer;
    std::ofstream * timerFile;
    std::ofstream * configurationFile;
    static double stepTime, oldClockTime;
    int iterationsNLSolver;
    int outputFilesDetail;

    double initialTemperature;
    SparseMatrix<data_type> globalCapacity;
    SparseMatrix<data_type> globalConductivity;

    VectorX<data_type> globalRHSHeat;
    VectorX<data_type> globalExternalHeat;
    VectorX<data_type> globalInternalHeat;

//////////////////// - - - SOA - - - ////////////////////////
    /*gmm::csc_matrix<data_type> h_globalCapacity;
    gmm::csc_matrix<data_type> h_globalConductivity;
    std::vector<data_type> h_globalRHSHeat;
    std::vector<data_type> h_globalExternalHeat;
    std::vector<data_type> h_globalInternalHeat;*/

    SparseMatrix<data_type> globalMass;
    VectorX<data_type> globalRHSForces;
    VectorX<data_type> globalInternalForces;
    VectorX<data_type> globalExternalForces;
    lmx::DenseMatrix<data_type> globalSparsePattern;

    std::ofstream outFile;

    static VectorX<double> gravity;
    static double alpha;
    static int dimension;
    static std::string contact;
    static bool visualization;
    static bool outputMatrices;
    static std::string constraintMethod;
    static double epsilon;
    static std::string smoothingType;

    VectorX<data_type> initThermalSimulation(Analysis * analysis, int, bool init = true);
    VectorX<data_type> initMechanicalSimulation(Analysis * analysis, bool init = true);
};

}

#endif
