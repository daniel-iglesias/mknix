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

#ifndef MKNIXANALYSISTHERMOMECHANICALDYNAMIC_H
#define MKNIXANALYSISTHERMOMECHANICALDYNAMIC_H

#include "analysis.h"

#include <LMX/lmx_diff_problem_first_second.h>

namespace mknix
{

/**
	@author AUTHORS <MAILS>
*/
class AnalysisThermoMechanicalDynamic : public Analysis
{
public:
    AnalysisThermoMechanicalDynamic();

    AnalysisThermoMechanicalDynamic(Simulation *, double, double, double, char *);

    ~AnalysisThermoMechanicalDynamic();

    std::string type()
    {
        return std::string("THERMOMECHANICALDYNAMIC");
    }

    void solve(lmx::Vector<data_type> *, lmx::Vector<data_type> *, lmx::Vector<data_type> *) override;

    void nextStep() override { };

    void init(lmx::Vector<data_type> * q_in, lmx::Vector<data_type> * qdot_in, int verbosity) override { };

private:
    char * integratorType;
    double to, tf, At;
    lmx::DiffProblemFirstSecond<Simulation, data_type> theProblem;

};

}

#endif
