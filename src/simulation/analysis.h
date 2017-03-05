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

#ifndef MKNIXANALYSIS_H
#define MKNIXANALYSIS_H

#include "common.h"
#include "LMX/lmx.h"

#include "LMX/lmx_diff_problem_first.h"
#include "LMX/lmx_diff_problem_second.h"

namespace mknix {

class Simulation;

/**
	@author AUTHORS <MAILS>
*/
class Analysis
{
public:
    Analysis();

    Analysis(Simulation *);

    virtual ~Analysis();

    virtual std::string type() = 0;

    void setEpsilon(double epsilon_in)
    {
        epsilon = epsilon_in;
    }

    virtual void solve(lmx::Vector<data_type> *,
                       lmx::Vector<data_type> * = 0,
                       lmx::Vector<data_type> * = 0
    ) = 0;

    virtual void init(lmx::Vector<data_type> * q_in, lmx::Vector<data_type> * qdot_in) = 0; // specialised in AnalysisThermalDynamic
    virtual void nextStep() = 0; // specialised in AnalysisThermalDynamic

protected:
    Simulation * theSimulation;
    double epsilon;

};

}

#endif
