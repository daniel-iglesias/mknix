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

#ifndef LMXNL_SOLVERS_H
#define LMXNL_SOLVERS_H

#include<iostream>


//////////////////////////////////////////// Doxygen file documentation entry:
/**
 * \file lmx_nlsolvers_double.h
 *
 * \brief NLSolverDoubleDouble class implementation
 *
 * Implements a partitioned residual "{R1(q), R2(q)} = 0" system and the methods for solving it.
 *
 * \author Daniel Iglesias
 */
//////////////////////////////////////////// Doxygen file documentation (end)



namespace lmx {

/**
 * \class NLSolverDouble
 * \brief Template class NLSolverDouble.
 * Non-linear systems implementation: "R(q) = 0" .
 *
 * This class permits the creation of a non-linear solver object.
 *
 * @author Daniel Iglesias.
 */

template<typename Sys, typename T=double>
class NLSolverDouble
{

public:

    NLSolverDouble()
            : increment1(0)
            , increment2(0)
            , conv1_1(0)
            , conv1_2(0)
//          , conv2(0)
//          , conv3(0)
            , epsilon(1E-6)
            , externalConvergence1_1(0)
            , externalConvergence1_2(0)
//          , externalConvergence2(0)
//          , externalConvergence3(0)
            , deltaInResidues(0)
            , info(1)
    /**
     * Empty constructor.
     */
    { }

    ~NLSolverDouble()
    /**
     * Destructor.
     */
    {
        if (increment1) {
            delete increment1;
            increment1 = 0;
        }
        if (increment2) {
            delete increment2;
            increment2 = 0;
        }
    }

    /** Set information level.
      */
    void setInfo(int level) { info = level; }


    template<class C>
    void setInitialConfiguration1(const lmx::Vector<C>& q_in)
    /**
     * Takes a lmx::Vector for initial value guess and dimensioning the problem.
     * @param q_in Initial guess values.
     */
    {
        q1.resize(q_in.size());
        delta_q1.resize(q_in.size());
        res_vector1.resize(q_in.size());
        jac_matrix1.resize(q_in.size(), q_in.size());
        q1 = q_in;
    }

    template<class C>
    void setInitialConfiguration2(const lmx::Vector<C>& q_in)
    /**
     * Takes a lmx::Vector for initial value guess and dimensioning the problem.
     * @param q_in Initial guess values.
     */
    {
        q2.resize(q_in.size());
        delta_q2.resize(q_in.size());
        res_vector2.resize(q_in.size());
        jac_matrix2.resize(q_in.size(), q_in.size());
        q2 = q_in;
    }

    void setSystem(Sys& system_in)
    /**
     * Sets witch Sys object is going to be used for member function calls.
     * @param system_in The Sys object.
     */
    { theSystem = &system_in; }

    void setDeltaInResidues(bool state = 1)
    /**
     * Sets whether the parameter in Residue function corresponds to the actual variables configuration
     * or indicates the increment of those variables.
     * @param state TRUE (default) if the variable's increment is going to be passed.
     */
    { deltaInResidues = state; }


    void setResidue1(void (Sys::*residue_in)(lmx::Vector<T>&, lmx::Vector<T>&))
    /**
     * Defines the member function that computes the residue.
     * @param residue_in Residue member function.
     */
    { res1 = residue_in; }

    void setResidue2(void (Sys::*residue_in)(lmx::Vector<T>&, lmx::Vector<T>&))
    /**
     * Defines the member function that computes the residue.
     * @param residue_in Residue member function.
     */
    { res2 = residue_in; }

    void setJacobian1(void (Sys::*jacobian_in)(lmx::Matrix<T>&, lmx::Vector<T>&))
    /**
     * Defines the member function that computes the tangent to the residue.
     * @param jacobian_in Jacobian member function.
     */
    { jac1 = jacobian_in; }

    void setJacobian2(void (Sys::*jacobian_in)(lmx::Matrix<T>&, lmx::Vector<T>&))
    /**
     * Defines the member function that computes the tangent to the residue.
     * @param jacobian_in Jacobian member function.
     */
    { jac2 = jacobian_in; }

    void setConvergence(double eps_in)
    /**
    * Defines the epsilon value for the L2 norm.
    * @param eps_in Value of the maximum L2 limit.
     */
    { epsilon = eps_in; }

    void setConvergence1(bool (Sys::*convergence_in)(lmx::Vector<T>&))
    /**
    * Defines the optional member function for convergence evaluation with residue parameter.
    * @param convergence_in Convergence evaluation member function.
     */
    {
        conv1_1 = convergence_in;
        externalConvergence1_1 = 1;
    }

    void setConvergence2(bool (Sys::*convergence_in)(lmx::Vector<T>&))
    /**
    * Defines the optional member function for convergence evaluation with residue parameter.
    * @param convergence_in Convergence evaluation member function.
     */
    {
        conv1_2 = convergence_in;
        externalConvergence1_2 = 1;
    }

//       void setConvergence
//         ( bool (Sys::*convergence_in)(lmx::Vector<T>&, lmx::Vector<T>&) )
//        /**
//          * Defines the optional member function for convergence evaluation with residue and configuration parameters.
//          * @param convergence_in Convergence evaluation member function.
//          */
//       { conv2 = convergence_in; externalConvergence2 = 1; }
// 
//       void setConvergence( bool (Sys::*convergence_in)(lmx::Vector<T>&, lmx::Vector<T>&, lmx::Vector<T>&) )
//        /**
//          * Defines the optional member function for convergence evaluation with residue, configuration and 
// 	       * increment vector parameters.
//          * @param convergence_in Convergence evaluation member function.
//          */
//       { conv3 = convergence_in; externalConvergence3 = 1; }

    bool convergence();

    void solve(int max_iter = 100);

    lmx::Vector<T>& getSolution1()
    /**
      * Solution vector read-write access.
      * @return Reference to the solution vector.
      */
    { return this->q1; }

    lmx::Vector<T>& getSolution2()
    /**
      * Solution vector read-write access.
      * @return Reference to the solution vector.
      */
    { return this->q2; }

    void setSparse1(std::vector<size_type>& rows, std::vector<size_type> columns)
    // Needs documentation
    { jac_matrix1.sparsePattern(rows, columns); }

    void setSparse2(std::vector<size_type>& rows, std::vector<size_type> columns)
    // Needs documentation
    { jac_matrix2.sparsePattern(rows, columns); }

private:
    lmx::Vector<T> q1;
    /**< Coordinates values for nl iterations. */
    lmx::Vector<T> delta_q1;
    /**< Coordinates values for nl iterations. */
    lmx::Matrix<T> jac_matrix1;
    /**< Jacobian -tangent- matrix (only used in Newton's method). */
    lmx::Vector<T> res_vector1;
    /**< Residual vector. */
    lmx::LinearSystem<T> * increment1; // jac_matrix*\delta q + f = 0
    // A*x = b
    lmx::Vector<T> q2;
    /**< Coordinates values for nl iterations. */
    lmx::Vector<T> delta_q2;
    /**< Coordinates values for nl iterations. */
    lmx::Matrix<T> jac_matrix2;
    /**< Jacobian -tangent- matrix (only used in Newton's method). */
    lmx::Vector<T> res_vector2;
    /**< Residual vector. */
    lmx::LinearSystem<T> * increment2; // jac_matrix*\delta q + f = 0
    // A*x = b
    Sys * theSystem;

    void (Sys::*res1)(lmx::Vector<T>&, lmx::Vector<T>&);

    void (Sys::*res2)(lmx::Vector<T>&, lmx::Vector<T>&);

    void (Sys::*jac1)(lmx::Matrix<T>&, lmx::Vector<T>&);

    void (Sys::*jac2)(lmx::Matrix<T>&, lmx::Vector<T>&);

    bool (Sys::*conv1_1)(lmx::Vector<T>&);

    bool (Sys::*conv1_2)(lmx::Vector<T>&);

//       bool (Sys::*conv2)(lmx::Vector<T>&, lmx::Vector<T>&);
//       bool (Sys::*conv3)(lmx::Vector<T>&, lmx::Vector<T>&, lmx::Vector<T>&);
    double epsilon, energy_i1, energy_01;
    double energy_i2, energy_02;
    bool externalConvergence1_1;
    bool externalConvergence1_2;
//       bool externalConvergence2;
//       bool externalConvergence3;
    bool deltaInResidues;
    int info;
    int iteration;
};

}; // namespace lmx

/////////////////////////////// Implementation of the methods defined previously

namespace lmx {


template<typename Sys, class T>
bool NLSolverDouble<Sys, T>::convergence()
/**
 * Internal convergence criteria.
 * Is used if no external convergence function is set.
 */
{
//     if( deltaInResidues == 0 ){ // use norm2 criteria
    if ((this->res_vector1.norm2() < epsilon) &&
        (this->res_vector2.norm2() < epsilon)) {
            return 1;
        } else { return 0; }
//     }
//     else{ // use energetic criteria
//       if( iteration == 0 ){
//         energy_01 = fabs(res_vector1 * q1); // store first residual energy
//         energy_02 = fabs(res_vector2 * q2); // store first residual energy
// //        cout << "iteration, energy = "<< iteration << " "<< energy_0 << endl;
//         if (energy_01 < 1E-50 || energy_02 < 1E-50){
// 	  cout << res_vector1 << q1 << endl;
// 	  cout << res_vector2 << q2 << endl;
// 	  return 1; // too small energy, rare
// 	}
//         else return 0; // common return value
//       }
//       else{ // rest of iterations
//         energy_i1 = fabs(res_vector1 * q1);
//         energy_i2 = fabs(res_vector2 * q2);
//         energy_i1 /= energy_01; // dimensionless residual energy rate
//         energy_i2 /= energy_02; // dimensionless residual energy rate
//         if (energy_i1 < epsilon && energy_i2 < epsilon)
// 	  return 1; // convergence!!
//         else return 0; // not converged
//       }
//     }
}

template<typename Sys, class T>
void NLSolverDouble<Sys, T>::solve(int max_iter)
/**
 * Solve function. Initiates the nl-solver loop.
 * @param max_iter Defines the maximun number of iterations for each iteration.
 */
{
    if (res_vector1.size() == 0 || res_vector2.size() == 0) {
        std::stringstream message;
        message << "Error in NLSolverDouble \"R(x) = 0\": dimension of problem not defined. \n"
        << "Use NLSystem::setInitialConfiguration( x_o ) function before the solve() function call." << endl;
        LMX_THROW(dimension_error, message.str());
    }

    switch (nl_solver_type) {

    case 0 : // select_nl_solver == 0 -> Newton's method
        if (!increment1) {
            increment1 = new lmx::LinearSystem<T>(jac_matrix1, delta_q1, res_vector1);
        }
        if (!increment2) {
            increment2 = new lmx::LinearSystem<T>(jac_matrix2, delta_q2, res_vector2);
        }
        if (info > 0) {
            cout << "     iter-NL\t|  RES1 |\t|  RES2 |\t|  Dq1  |\t|  Dq2  |" << endl;
            cout.setf(std::ios::scientific, std::ios::floatfield);
            cout.precision(3);
        }

        for (iteration = 0; iteration < max_iter; iteration++) {
            if (info > 0) {
                cout << "\t" << iteration << "\t";
            }
            if (deltaInResidues) {
                (theSystem->*res1)(res_vector1, delta_q1);
                (theSystem->*res2)(res_vector2, delta_q2);
            }
            else {
                (theSystem->*res1)(res_vector1, q1);
                (theSystem->*res2)(res_vector2, q2);
            }
            res_vector1 *= -1.;
            res_vector2 *= -1.;
            if (info > 0) {
                cout << res_vector1.norm1() << "\t"
                << res_vector2.norm1() << "\t";
            }

            if (externalConvergence1_1 &&
                externalConvergence1_2) {
                if ((theSystem->*conv1_1)(res_vector1) &&
                    (theSystem->*conv1_2)(res_vector2)
                        ) {
                    if (info > 0) {
                        std::cout << endl << endl;
                    }
                    break;
                }
            }
//           else if ( externalConvergence2 ){
//             if ( (theSystem->*conv2)(res_vector, q) ){
//               if ( info > 0 )
//                 std::cout << endl << endl;
//               break;
//             }
//           }
//           else if ( externalConvergence3 ){
//             if ( (theSystem->*conv3)(res_vector, q, increment->getSolution()) ){
//               if ( info > 0 )
//                 std::cout << endl << endl;
//               break;
//             }
//           }
            else if (this->convergence()) {
                if (info > 0) {
                    std::cout << endl << endl;
                }
                break;
            }
            (theSystem->*jac1)(jac_matrix1, q1);
            (theSystem->*jac2)(jac_matrix2, q2);
            q1 += increment1->solveYourself(); // = delta_q
            q2 += increment2->solveYourself(); // = delta_q
            if (info > 0) {
                std::cout << increment1->getSolution().norm1() << "\t";
                std::cout << increment2->getSolution().norm1() << endl;
            }
        }
        break;

    default:
        throw std::invalid_argument("Unknown solver type");
    }
    if (info > 0) {
        cout.unsetf(std::ios::floatfield);
    }
}

}; // namespace lmx


#endif
