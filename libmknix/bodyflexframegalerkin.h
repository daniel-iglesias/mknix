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
#ifndef MKNIXFLEXFRAMEGALERKIN_H
#define MKNIXFLEXFRAMEGALERKIN_H

#include "bodyflex.h"

namespace mknix {

/**
  @author AUTHORS <MAILS>
*/
class FlexFrameGalerkin : public FlexBody
{
public:
    FlexFrameGalerkin();

    FlexFrameGalerkin( std::string );

    ~FlexFrameGalerkin();

    std::string getType() {
        return bodyType;
    }

    void setType( std::string type_in ) {
        bodyType = type_in;
    }

//     void initialize( );

    void calcMassMatrix( );

    void calcInternalForces( );

    void calcExternalForces( );

    void calcTangentMatrix( );

    void assembleMassMatrix( lmx::Matrix<data_type> & );

    void assembleInternalForces( lmx::Vector<data_type> & );

    void assembleExternalForces( lmx::Vector<data_type> & );

    void assembleTangentMatrix( lmx::Matrix<data_type> & );

    void outputStep( const lmx::Vector<data_type>&, const lmx::Vector<data_type>& );

    void outputStep( const lmx::Vector<data_type>& );

private:
    void recoverStressField( int );

private:
    std::string bodyType;

};

}

#endif
