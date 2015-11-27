/***************************************************************************
 *   Copyright (C) 2013 by Daniel Iglesias                                 *
 *   http://code.google.com/p/mknix                                        *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
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
class Analysis {
public:
    Analysis();

    Analysis( Simulation* );

    virtual ~Analysis();

    virtual std::string type() = 0;

    void setEpsilon( double epsilon_in) {
        epsilon = epsilon_in;
    }

    virtual void solve( lmx::Vector<data_type> *,
                        lmx::Vector<data_type> * = 0,
                        lmx::Vector<data_type> * = 0
                      ) = 0;
    
    virtual void init(lmx::Vector< data_type > *){} // speciallized in AnalysisThermalDynamic
    virtual void nextStep(){} // speciallized in AnalysisThermalDynamic

protected:
    Simulation* theSimulation;
    double epsilon;

};

}

#endif
