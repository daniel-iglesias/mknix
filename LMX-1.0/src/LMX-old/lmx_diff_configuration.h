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

#ifndef LMXCONFIGURATION_H
#define LMXCONFIGURATION_H

#include"lmx_except.h"
#include <memory>

//////////////////////////////////////////// Doxygen file documentation entry:
/*!
  \file lmx_diff_configuration.h

  \brief Configuration variables and time.

  This class contains the time information at (n-m) step [m-order-integrator & n-step). Also, the variables that the integrator needs are stored in the parameters.

  \author Daniel Iglesias

*/
//////////////////////////////////////////// Doxygen file documentation (end)


namespace lmx {

/**
\class Configuration
\brief Template class Configuration.

Basic variables, derivatives and time. The data is stored in the last steps, which number depends on integrator type.

@author Daniel Iglesias .
 */
template<class T>
class Configuration
{
private:
//   public:
    int vectorSize;
    std::vector<std::vector<std::unique_ptr<lmx::Vector<T> > > > q;
//     std::vector< std::vector< lmx::Vector< T >* > > q; /**< STL vector of vectors of  coordinates with "n-steps" columns and "m-diff system order" rows. */
    std::unique_ptr<lmx::Vector<T> > temp;
//     lmx::Vector< T >* temp; /**< temporary pointer to vector for advance function */
//     std::vector< double > time; /**< Time vector. Stores all the steps... */
    double lastStepSize;
    double presentTime;
    int steps;
    bool verbose;

public:

    /** Empty constructor. */
    Configuration()
            : vectorSize(0)
            , lastStepSize(0)
            , presentTime(0)
            , steps(0)
            , verbose(true)
//     , temp(0)
    { }

    /** Standard constructor.
     * \param t_o Time at first step.
     */
    Configuration(double t_o)
            : vectorSize(0)
            , verbose(1)
            , lastStepSize(0)
            , steps(1)
            , presentTime(t_o)
//     , temp(0)
    { /*time.push_back( t_o );*/ }

    /** Destructor. */
    ~Configuration()
    {
//       if (temp) { delete temp; temp = 0;}

//       for ( int i = 0; i < q.size(); ++i ){
//         for ( int j = q[i].size(); j > 0; --j ){
//           delete q[i][j-1];
// 	  q[i][j-1] = 0;
//         }
//       }
    }

    /**
     * @brief Needed to improve the output of DiffProblemFirstSecond runs
     */
    void quiet() { verbose = 0; }

    void nextStep(double& stepSize);

    /**
     * @param step Indicates the (actual - step) time step.
     * @return The time value of the step.
     */
    double getTime(int step = 0)
//       { return this->time[time.size() - 1 - step]; }
    { return (presentTime - step * lastStepSize); }

    /**
     * Access to size of time line.
     */
    int getTimeSize()
//     { return this->time.size(); }
    { return steps; }

    /**
     * @return Value of last time increment.
     */
    double getLastStepSize() { return this->lastStepSize; }

    /**
     * @param order Differential order of configuration.
     * @param step Indicates the (actual - step) time step
     * @return The configuration of the diff-order and step specified.
     */
    const lmx::Vector<T>& getConf(int order, int step = 0) { return *(this->q[order][step]); }

    /**
     * @return Maximum differential order of stored configuration.
     */
    int getDiffOrder() { return this->q.size() - 1; }


    /**
     * Sets the time of actual (last) time step.
     * @param time_in Value of time to be set.
     */
    void setTime(double& time_in)
//       { time.push_back( time_in ); }
    {
        presentTime = time_in;
        ++steps;
    }

    void setInitialCondition(int diff_order, lmx::Vector<T>& q_o);

    void setStoredSteps(int steps_q_o, int steps_q_i, int steps_q_n);

    /**
     * @param diff_order Differential order of configuration.
     * @param values Values of configuration to be set.
     * @param time_step Indicates the (actual - step) time step
     */
    void setConf(int diff_order, Vector <T> values, int time_step = 0) { *q[diff_order][time_step] = values; }

    /**
     * @param diff_order Differential order of configuration.
     * @param time_step Indicates the (actual - step) time step
     * @return Values of configuration.
     */
    Vector <T>& setConf(int diff_order, int time_step = 0) { return *q[diff_order][time_step]; }

};

template<class T>
void Configuration<T>::setInitialCondition(int diff_order, lmx::Vector<T>& q_o)
/**
 * @param diff_order Order of differential system.
 * @param q_o Value of initial condition to be set.
 */
{
    if (vectorSize == 0) {
        vectorSize = q_o.size();
    }
    else if (vectorSize != q_o.size()) {
        std::stringstream message;
        message << "ERROR : trying to assing an initial condition vector of different size than the exising ones." <<
        endl;
        LMX_THROW(lmx::failure_error, message.str());
    }
    if (diff_order + 1 >= q.size()) {
        for (int i = q.size(); i <= diff_order + 1; ++i) {
            std::unique_ptr<Vector<T> > ptr(new Vector<T>(vectorSize));
            q.push_back(std::vector<std::unique_ptr<Vector<T> > >());
            q[i].push_back(std::move(ptr));
        }
    }
    *q[diff_order][0] = q_o; // copies values... perhaps should use input values instead.

    cout << "--------------------------------------------------------" << endl;
    cout << "An initial condition has been set:" << endl;
    cout << "Derivative order = " << diff_order;
    cout << ", size of vector = " << q_o.size() << endl;
//       << ", q_0 = " << *q[diff_order][0] << endl;
    cout << "--------------------------------------------------------" << endl;
}


/**
 * Reshapes storing vectors to the required dimension for last steps storage.
 * @param steps_q_o Number of steps to store the basic variable.
 * @param steps_q_i Number of steps to store from first to (n-1) derivatives.
 * @param steps_q_n Number of steps to store the n derivative.
 */
template<class T>
void Configuration<T>::setStoredSteps(int steps_q_o, int steps_q_i, int steps_q_n)
{
    if (vectorSize == 0) {
        std::stringstream message;
        message << "ERROR : Initial conditions must be assigned before defining step storing." << endl;
        LMX_THROW(lmx::failure_error, message.str());
    }
    // Add new columns to the vectors, starting in column 1:
    unsigned int i, j;
    for (j = 1; j < steps_q_o; ++j) {
        std::unique_ptr<Vector<T> > ptr(new Vector<T>(vectorSize));
        q[0].push_back(std::move(ptr));
    }
    for (i = 1; i < (q.size() - 1); ++i) {
        for (j = 1; j < steps_q_i; ++j) {
            std::unique_ptr<Vector<T> > ptr(new Vector<T>(vectorSize));
            q[i].push_back(std::move(ptr));
        }
    }
    for (j = 1; j < steps_q_n; ++j) {
        std::unique_ptr<Vector<T> > ptr(new Vector<T>(vectorSize));
        q[q.size() - 1].push_back(std::move(ptr));
    }

    cout << "--------------------------------------------------------" << endl;
    cout << "Configuration has been resized to the following vectors:" << endl;
    for (i = 0; i < q.size(); ++i) {
        cout << "Derivative order = " << i << ", time line size = " << q[i].size() << endl;
    }
    cout << "--------------------------------------------------------" << endl;
}

/**
 * Reshapes storing vectors to the required dimension for last steps storage.
 * @param step_size time increment between last and next step.
 */
template<class T>
void Configuration<T>::nextStep(double& step_size)
{
    unsigned int i, j;

//   cout << "             Step number " << time.size() << " solved, time = " << time.back() << endl;
//   cout << "--------------------------------------------------------" << endl;
//   for ( i=0; i<q.size(); ++i)
//     cout << "Derivative order = " << i << ", q_t = " << *q[i][0] << endl;
//   cout << "--------------------------------------------------------" << endl;

    for (i = 0; i < q.size(); ++i) {
        if (q[i].size() > 1)
        {
            temp = std::move(q[i].back()); // Save direction to last element.
            for (j = q[i].size() - 1; j > 0; --j) {
                q[i][j] = std::move(q[i][j - 1]); // Move back elements
            }
            q[i][0] = std::move(temp);
            // Optional, may improve convergence but increases step-time
            *q[i][0] = *q[i][1];
        }
    }

    lastStepSize = step_size;
    presentTime += step_size;
    ++steps;
//   time.push_back( time.back() + lastStepSize );

    if (verbose) {
        cout << "--------------------------------------------------------" << endl;
//     cout << "             Solving step number " << time.size()-1 << " time = " << time.back() << endl;
        cout << "             Solving step number " << steps << " time = " << presentTime << endl;
        cout << "--------------------------------------------------------" << endl;
    }
}

}; // namespace lmx

#endif
