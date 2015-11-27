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
 

#ifndef LMXDATA_VEC_H
#define LMXDATA_VEC_H

#include"lmx_mat_data.h"

//////////////////////////////////////////// Doxygen file documentation entry:
    /*!
      \file lmx_mat_data_vec.h
      
      \brief This file contains the declaration of the data_vec class' pure virtual functions.
      
      For classes derived from Data_vec pure virtual class, all these methods must be implemented. Thus, for comprobation and checking, the methods here declared are as well documented.
      
      \author Daniel Iglesias Ib�ez
      
    */
//////////////////////////////////////////// Doxygen file documentation (end)

namespace lmx {

    /**
    \class Data_vec 
    \brief Template class Data_vec.
    Container for Vector data.
    
    This class represents the skeleton for the data container used by the Vector class. No parameter nor function implementation here, just pure virtual class. See derived classes for details in implementation. Also maybe useful to see how this class is used in the Vector class.
    
    @author Daniel Iglesias Ib�ez.
    */
template <typename T> class Data_vec : public Data<T>{
public:

  /** Empty constructor. */
  Data_vec(){}

  /** Destructor. */
  virtual ~Data_vec(){}

  /** Read file */
  virtual void readDataFile(const char*) = 0;

  /** Write file */
  virtual void writeDataFile(const char*) = 0;

};

};

#endif
