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

#ifndef LMXNL_SOLVERS_DOUBLE_H
#define LMXNL_SOLVERS_DOUBLE_H

#include<iostream>


//////////////////////////////////////////// Doxygen file documentation entry:
/**
 * \file lmx_nlsolvers.h
 *
 * \brief NLSolver class implementation
 *
 * Implements a typical residue "R(q) = 0" system and the methods for solving it.
 *
 * \author Daniel Iglesias
 */
//////////////////////////////////////////// Doxygen file documentation (end)



namespace lmx {

/**
 * \class NLSolver
 * \brief Template class NLSolver.
 * Non-linear systems implementation: "R(q) = 0" .
 *
 * This class permits the creation of a non-linear solver object.
 *
 * @author Daniel Iglesias.
 */

template<typename Sys, typename T=double>
class NLSolver
{

public:

    NLSolver()
            : increment(0)
            , conv1(0)
            , conv2(0)
            , conv3(0)
            , epsilon(1E-6)
            , externalConvergence1(0)
            , externalConvergence2(0)
            , externalConvergence3(0)
            , deltaInResidue(0)
            , info(1)
    /**
     * Empty constructor.
     */
    { }

    ~NLSolver()
    /**
     * Destructor.
     */
    {
        if (increment) {
            delete increment;
            increment = 0;
        }
    }

    /** Set information level.
      */
    void setInfo(int level) { info = level; }


    template<class C>
    void setInitialConfiguration(const lmx::Vector<C>& q_in)
    /**
     * Takes a lmx::Vector for initial value guess and dimensioning the problem.
     * @param q_in Initial guess values.
     */
    {
        q.resize(q_in.size());
        delta_q.resize(q_in.size());
        res_vector.resize(q_in.size());
        jac_matrix.resize(q_in.size(), q_in.size());
        q = q_in;
    }

    void setSystem(Sys& system_in)
    /**
     * Sets witch Sys object is going to be used for member function calls.
     * @param system_in The Sys object.
     */
    { theSystem = &system_in; }

    void setDeltaInResidue(bool state = 1)
    /**
     * Sets whether the parameter in Residue function corresponds to the actual variables configuration
     * or indicates the increment of those variables.
     * @param state TRUE (default) if the variable's increment is going to be passed.
     */
    { deltaInResidue = state; }


    void setResidue(void (Sys::*residue_in)(lmx::Vector<T>&, lmx::Vector<T>&))
    /**
     * Defines the member function that computes the residue.
     * @param residue_in Residue member function.
     */
    { res = residue_in; }

    void setJacobian(void (Sys::*jacobian_in)(lmx::Matrix<T>&, lmx::Vector<T>&))
    /**
     * Defines the member function that computes the tangent to the residue.
     * @param jacobian_in Jacobian member function.
     */
    { jac = jacobian_in; }

    void setConvergence(double eps_in)
    /**
    * Defines the epsilon value for the L2 norm.
    * @param eps_in Value of the maximum L2 limit.
     */
    { epsilon = eps_in; }

    void setConvergence(bool (Sys::*convergence_in)(lmx::Vector<T>&))
    /**
    * Defines the optional member function for convergence evaluation with residue parameter.
    * @param convergence_in Convergence evaluation member function.
     */
    {
        conv1 = convergence_in;
        externalConvergence1 = 1;
    }

    void setConvergence
            (bool (Sys::*convergence_in)(lmx::Vector<T>&, lmx::Vector<T>&))
    /**
      * Defines the optional member function for convergence evaluation with residue and configuration parameters.
      * @param convergence_in Convergence evaluation member function.
      */
    {
        conv2 = convergence_in;
        externalConvergence2 = 1;
    }

    void setConvergence(bool (Sys::*convergence_in)(lmx::Vector<T>&, lmx::Vector<T>&, lmx::Vector<T>&))
    /**
      * Defines the optional member function for convergence evaluation with residue, configuration and
        * increment vector parameters.
      * @param convergence_in Convergence evaluation member function.
      */
    {
        conv3 = convergence_in;
        externalConvergence3 = 1;
    }

    bool convergence();

    void solve(int max_iter = 100);

    lmx::Vector<T>& getSolution()
    /**
      * Solution vector read-write access.
      * @return Reference to the solution vector.
      */
    { return this->q; }

    void setSparse(std::vector<size_type>& rows, std::vector<size_type> columns)
    // Needs documentation
    { jac_matrix.sparsePattern(rows, columns); }

private:
    lmx::Vector<T> q;
    /**< Coordinates values for nl iterations. */
    lmx::Vector<T> delta_q;
    /**< Coordinates values for nl iterations. */
    lmx::Matrix<T> jac_matrix;
    /**< Jacobian -tangent- matrix (only used in Newton's method). */
    lmx::Vector<T> res_vector;
    /**< Residual vector. */
    lmx::LinearSystem<T> * increment; // jac_matrix*\delta q + f = 0
    // A*x = b
    Sys * theSystem;

    void (Sys::*res)(lmx::Vector<T>&, lmx::Vector<T>&);

    void (Sys::*jac)(lmx::Matrix<T>&, lmx::Vector<T>&);

    bool (Sys::*conv1)(lmx::Vector<T>&);

    bool (Sys::*conv2)(lmx::Vector<T>&, lmx::Vector<T>&);

    bool (Sys::*conv3)(lmx::Vector<T>&, lmx::Vector<T>&, lmx::Vector<T>&);

    double epsilon, energy_i, energy_0;
    bool externalConvergence1;
    bool externalConvergence2;
    bool externalConvergence3;
    bool deltaInResidue;
    int info;
    int iteration;
};

}; // namespace lmx

/////////////////////////////// Implementation of the methods defined previously

namespace lmx {


template<typename Sys, class T>
bool NLSolver<Sys, T>::convergence()
/**
 * Internal convergence criteria.
 * Is used if no external convergence function is set.
 */
{
    if (deltaInResidue == 0) { // use norm2 criteria
        if (this->res_vector.norm2() < epsilon) { return 1; }
        else { return 0; }
    }
    else { // use energetic criteria
        if (iteration == 0) {
            energy_0 = fabs(res_vector * q); // store first residual energy
//        cout << "iteration, energy = "<< iteration << " "<< energy_0 << endl;
            if (energy_0 < 1E-50) {
                return 1; // too small energy, rare
            } else { return 0; } // common return value
        }
        else { // rest of iterations
            energy_i = fabs(res_vector * q);
            energy_i /= energy_0; // dimensionless residual energy rate
            if (energy_i < epsilon) {
                return 1; // convergence!!
            } else { return 0; } // not converged
        }
    }
}

template<typename Sys, class T>
void NLSolver<Sys, T>::solve(int max_iter)
/**
 * Solve function. Initiates the nl-solver loop.
 * @param max_iter Defines the maximun number of iterations for each iteration.
 */
{
    if (res_vector.size() == 0) {
        std::stringstream message;
        message << "Error in NLSolver \"R(x) = 0\": dimension of problem not defined. \n"
        << "Use NLSystem::setInitialConfiguration( x_o ) function before the solve() function call." << endl;
        LMX_THROW(dimension_error, message.str());
    }

    switch (nl_solver_type) {

    case 0 : // select_nl_solver == 0 -> Newton's method
        if (!increment) {
            increment = new lmx::LinearSystem<T>(jac_matrix, delta_q, res_vector);
        }
        if (info > 0) {
            cout << "     iter-NL\t  | RES |\t|| RES ||\t|| Dq ||" << endl;
            cout.setf(std::ios::scientific, std::ios::floatfield);
            cout.precision(3);
        }

        for (iteration = 0; iteration < max_iter; iteration++) {
            if (info > 0) {
                cout << "\t" << iteration << "\t";
            }
            if (deltaInResidue) {
                (theSystem->*res)(res_vector, delta_q);
            } else { (theSystem->*res)(res_vector, q); }
            res_vector *= -1.;
            if (info > 0) {
                cout << res_vector.norm1() << "\t" << res_vector.norm2() << "\t";
            }

            if (externalConvergence1) {
                if ((theSystem->*conv1)(res_vector)) {
                    if (info > 0) {
                        std::cout << endl << endl;
                    }
                    break;
                }
            }
            else if (externalConvergence2) {
                if ((theSystem->*conv2)(res_vector, q)) {
                    if (info > 0) {
                        std::cout << endl << endl;
                    }
                    break;
                }
            }
            else if (externalConvergence3) {
                if ((theSystem->*conv3)(res_vector, q, increment->getSolution())) {
                    if (info > 0) {
                        std::cout << endl << endl;
                    }
                    break;
                }
            }
            else if (this->convergence()) {
                if (info > 0) {
                    std::cout << endl << endl;
                }
                break;
            }
            (theSystem->*jac)(jac_matrix, q);
            q += increment->solveYourself(); // = delta_q
            if (info > 0) {
                std::cout << increment->getSolution().norm2() << endl;
            }
        }

        break;
    }
    if (info > 0) {
        cout.unsetf(std::ios::floatfield);
    }
}

}; // namespace lmx


#endif
