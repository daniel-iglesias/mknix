/***************************************************************************
 *   Copyright (C) 2005 by Daniel Iglesias                                 *
 *   http://code.google.com/p/lmx                                          *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU Library General Public License as       *
 *   published by the Free Software Foundation; either version 2 of the    *
 *   License, or (at your option) any later version.                       *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU Library General Public     *
 *   License along with this program; if not, write to the                 *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef LMXDIFF_PROBLEM_H
#define LMXDIFF_PROBLEM_H


//////////////////////////////////////////// Doxygen file documentation entry:
/*!
  \file lmx_diff_problem.h

  \brief DiffProblem class implementation

  Describes an initial value for a dynamic system with an ODE or DAE description.

  This is the base file of lmx_diff systems' manipulation and solution.

  \author Daniel Iglesias

*/
//////////////////////////////////////////// Doxygen file documentation (end)
#include <map>
#include"lmx_nlsolvers.h"
#include "lmx_diff_configuration.h"
#include "lmx_diff_integrator_ab.h"
#include "lmx_diff_integrator_am.h"
#include "lmx_diff_integrator_bdf.h"
#include "lmx_diff_integrator_centraldiff.h"

namespace lmx {

/**
\class DiffProblem
\brief Template class DiffProblem.
Implementation for ODE system solvers.

This class implements methods for defining and solving initial value problems described by a TotalDiff class' derivided object, and initial conditions in the form \f$ \dot{q}(t_o) = \dot{q}_o \f$,  \f$ q(t_o) = q_o \f$.

@author Daniel Iglesias.
*/
template<typename Sys, typename T=double>
class DiffProblem
{

public:

    /** Empty constructor. */
    DiffProblem()
            : theConfiguration(0)
            , theIntegrator(0)
            , theNLSolver(0)
            , theSystem(0)
            , p_delta_q(0)
            , b_steptriggered(0)
            , vervosity(2)
            { }

    /** Destructor. */
    virtual ~DiffProblem()
    {
        if (theIntegrator) delete theIntegrator;
        if (theConfiguration) delete theConfiguration;
    }

    /**
     * @param system_in Object that defines the differential system equations.
     */
    void setDiffSystem(Sys& system_in) { theSystem = &system_in; }

    void setIntegrator(int type, int opt1 = 0, int opt2 = 0);

    void setIntegrator(const std::string& type, int opt2 = 0);

    void setInitialConfiguration(lmx::Vector<T>& q_o);

    void setInitialConfiguration(lmx::Vector<T>& q_o, lmx::Vector<T>& qdot_o);

    void setOutputFile(const char * filename, int diffOrder);

    void setTimeParameters(double to_in, double tf_in, double step_size_in);

    void iterationResidue(lmx::Vector<T>& residue, lmx::Vector<T>& q_actual);

    void setStepTriggered(void (Sys::* stepTriggered_in)());

    // needs documentation:
    void setConvergence(double eps_in) { epsilon = eps_in; }

    // needs documentation:
    const lmx::Vector<T>& getConfiguration(int order, int step = 0) { return theConfiguration->getConf(order, step); }

    bool isIntegratorExplicit()
    {
        if (theIntegrator) return this->theIntegrator->isExplicit();
        return false;
    }
    
    void setVervosity(int level){
        vervosity = level;
    }

    /**
     * Solve method to be implemented in derived classes.
     */
    virtual void initialize() = 0;

    virtual void stepSolve() = 0;

    virtual void solve() = 0;

    void advance(); ///< Advances the configuration and the integrator without solving the system.

protected:
    void writeStepFiles();

private:
    virtual void stepSolveExplicit() = 0;

    virtual void stepSolveImplicit() = 0;

protected:
    lmx::Configuration<T> * theConfiguration; ///< Pointer to the Configuration object, (auto-created).
    lmx::IntegratorBase<T> * theIntegrator; ///< Pointer to the Integrator object, (auto-created).
    lmx::NLSolver<T> * theNLSolver; ///< Pointer to the NLSolver object, (auto-created).
    Sys * theSystem; ///< Pointer to object where the differential system is defined.
    lmx::Vector<T> * p_delta_q; ///< Stores pointer to NLSolver increment.
    bool b_steptriggered; ///< 1 if stepTriggered function is set.
    double to; ///< Value of the start time stored from input.
    double tf; ///< Value of the finish time stored from input.
    double stepSize; ///< Value of the time step stored from input.
    double epsilon; ///< Value for L2 convergence.
    std::map<int, std::ofstream *> fileOutMap; ///< collection of output streams for each diff-order requested.
    void (Sys::* stepTriggered)(); ///< function called at the end of each time step
    int vervosity;
};


/////////////////////////////// Implementation of the methods defined previously

/**
 * Defines the integrator that will be used for configuration advance & actualization.
 * @param type Key of integrator family to use.
 * @param opt1 Optional value for some integrators (usually specifies the order).
 * @param opt2 Optional value for some integrators.
 */
template<typename Sys, typename T>
void DiffProblem<Sys, T>::setIntegrator(int type, int opt1, int opt2)
{
    switch (type) {
    case 0 : // integrator == 0 -> Adams-Bashford
        theIntegrator = new IntegratorAB<T>(opt1);
        break;

    case 1 : // integrator == 1 -> Adams-Moulton
        theIntegrator = new IntegratorAM<T>(opt1);
        break;

    case 2 : // integrator == 2 -> BDF
        theIntegrator = new IntegratorBDF<T>(opt1);
        break;

    case 3 : // integrator == 2 -> BDF
        theIntegrator = new IntegratorCentralDifference<T>();
        break;

    }

}

/**
 * Defines the integrator that will be used for configuration advance & actualization.
 * @param type Key of integrator to use.
 * @param opt2 Optional value for some integrators.
 */
template<typename Sys, typename T>
void DiffProblem<Sys, T>::setIntegrator(const std::string& type, int opt2)
{
    if (type == "AB-1") {
        theIntegrator = new IntegratorAB<T>(1);
    } else if (type == "AB-2") {
        theIntegrator = new IntegratorAB<T>(2);
    } else if (type == "AB-3") {
        theIntegrator = new IntegratorAB<T>(3);
    } else if (type == "AB-4") {
        theIntegrator = new IntegratorAB<T>(4);
    } else if (type == "AB-5") {
        theIntegrator = new IntegratorAB<T>(5);
    } else if (type == "AM-1") {
        theIntegrator = new IntegratorAM<T>(1);
    } else if (type == "AM-2") {
        theIntegrator = new IntegratorAM<T>(2);
    } else if (type == "AM-3") {
        theIntegrator = new IntegratorAM<T>(3);
    } else if (type == "AM-4") {
        theIntegrator = new IntegratorAM<T>(4);
    } else if (type == "AM-5") {
        theIntegrator = new IntegratorAM<T>(5);
    } else if (type == "BDF-1") {
        theIntegrator = new IntegratorBDF<T>(1);
    } else if (type == "BDF-2") {
        theIntegrator = new IntegratorBDF<T>(2);
    } else if (type == "BDF-3") {
        theIntegrator = new IntegratorBDF<T>(3);
    } else if (type == "BDF-4") {
        theIntegrator = new IntegratorBDF<T>(4);
    } else if (type == "BDF-5") {
        theIntegrator = new IntegratorBDF<T>(5);
    } else if (type == "CD") theIntegrator = new IntegratorCentralDifference<T>();
}

/**
 * Defines initial conditions for first order diff. problems.
 * @param q_o Zero-order initial configuration.
 */
template<typename Sys, typename T>
void DiffProblem<Sys, T>::setInitialConfiguration(lmx::Vector<T>& q_o)
{
    if (theConfiguration == 0) {
        theConfiguration = new Configuration<T>;
    }

    theConfiguration->setInitialCondition(0, q_o);
    if (vervosity == 0) theConfiguration->quiet();
}

/**
 * Defines initial conditions for second order diff. problems.
 * @param q_o Zero-order initial configuration.
 * @param qdot_o First-order initial configuration.
 */
template<typename Sys, typename T>
void DiffProblem<Sys, T>::setInitialConfiguration(lmx::Vector<T>& q_o, lmx::Vector<T>& qdot_o)
{
    if (theConfiguration == 0) {
        theConfiguration = new Configuration<T>;
    }

    theConfiguration->setInitialCondition(0, q_o);
    theConfiguration->setInitialCondition(1, qdot_o);
    if (vervosity == 0) theConfiguration->quiet();
}


/**
 * Defines which variables to store, specifing the file name for each diff-order.
 *
 * @param filename Name of file for storing variables.
 * @param diffOrder
 */
template<typename Sys, typename T>
void DiffProblem<Sys, T>::setOutputFile(const char * filename, int diffOrder)
{
    if (!(fileOutMap[diffOrder] == 0)) {
        cout << "WARNING: Changing opened file for diff order = " << diffOrder << endl;
        cout << "         New name: " << filename << endl;
        fileOutMap[diffOrder]->open(filename);
    }
    else {
        fileOutMap[diffOrder] = new std::ofstream(filename);
    }
}


/**
 * Defines basic time parameters.
 *
 * @param to_in Initial time.
 * @param tf_in End time
 * @param step_size_in Prefered time step.
 */
template<typename Sys, typename T>
void DiffProblem<Sys, T>::setTimeParameters(double to_in, double tf_in, double step_size_in)
{
    this->to = to_in;
    this->tf = tf_in;
    this->stepSize = step_size_in;
}

/**
 * When the configuration advances, this method is invoked for writing the requested diff-order values.
 */
template<typename Sys, typename T>
void DiffProblem<Sys, T>::writeStepFiles()
{
    std::map<int, std::ofstream *>::iterator it;
    for (it = fileOutMap.begin(); it != fileOutMap.end(); ++it) {
        it->second->setf(std::ios::scientific, std::ios::floatfield);
        it->second->precision(6);
        *(it->second) << theConfiguration->getTime() << "\t";
        for (unsigned int i = 0; i < theConfiguration->getConf(it->first, 0).size(); ++i) {
            *(it->second) << theConfiguration->getConf(it->first, 0).readElement(i) << "\t";
        }
        *(it->second) << endl;
    }
}

/**
 * Defines a function call between time steps.
 *
 */
template<typename Sys, typename T>
void DiffProblem<Sys, T>::setStepTriggered(void (Sys::* stepTriggered_in)())
{
    this->stepTriggered = stepTriggered_in;
    b_steptriggered = 1;
}

/**
 * Advances the configuration and the integrator without solving the system.
 */
template<typename Sys, typename T>
void DiffProblem<Sys, T>::advance()
{
    this->theConfiguration->nextStep(this->stepSize);
    this->theIntegrator->advance();
    if (this->b_steptriggered) (this->theSystem->*stepTriggered)();
    this->writeStepFiles();
}


}; // namespace lmx


#endif
