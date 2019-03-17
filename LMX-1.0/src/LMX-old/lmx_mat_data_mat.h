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
 

#ifndef LMXDATA_MAT_H
#define LMXDATA_MAT_H

#include <vector>
#include"lmx_mat_data.h"

//////////////////////////////////////////// Doxygen file documentation entry:
    /*!
      \file lmx_mat_data_mat.h
      
      \brief This file contains the declaration of the data_mat class' pure virtual functions.
      
      For classes derived from Data_mat pure virtual class, all these methods must be implemented. Thus, for comprobation and checking, the methods here declared are as well documented.
      
      \author Daniel Iglesias 
      
    */
//////////////////////////////////////////// Doxygen file documentation (end)


namespace lmx {

// Forward declarations:
template <typename T> class Vector;

    /**
    \class Data_mat 
    \brief Template class Data_mat.
    Container for Matrix and Vector data.
    
    This class represents the skeleton for the data container used by Matrix and Vector classes. No parameter nor function implementation here, just pure virtual class. See derived classes for details in implementation. Also maybe useful to see how this class is used in Matrix and Vector classes.
    
    @author Daniel Iglesias .
    */
template <typename T> class Data_mat : public Data<T>{
public:

  /** Empty constructor. */
  Data_mat(){}

  /** Destructor. */
  virtual ~Data_mat(){}

  /**
   * Multiplication function.
   */
  virtual void multiply(const Data<T>* , const Data<T>* ) = 0;

  /** Read data in Matrix Market format method.
   * Opens the file specified and reads the matrix's data in it,
   * suposing it's stored in Matrix Market format. */
  virtual void read_mm_file(const char*) = 0;

  /** Read data in Harwell-Boeing format method.
   * Opens the file specified and reads the matrix's data in it, 
   * suposing it's stored in Harwell-Boeing format. */
  virtual void read_hb_file(const char*) = 0;

  /** Write data in Harwell-Boeing format method.
   * Opens the file specified and writes the matrix's data in it. */
  virtual void write_hb_file(const char*) = 0;

  /** Traspose matrix function. */
  virtual void trn() = 0;

  /**
   * Returns TRUE or FALSE depending of element existance.
   * Needed for less expensive access to CSC matrices.
   * @return TRUE if the element exists in internal storage structure.
   */
  virtual bool exists( size_type, size_type ) = 0;

  virtual void factorize() = 0;
  
  virtual void subsSolve(Vector<T>& rhs) {}
  
  /** Prepares the sparse structure of a CSC matrix. */
  virtual void setSparsePattern( Vector<size_type>&, Vector<size_type>& )
  {}

  /** Prepares the sparse structure of a CSC matrix. */
  virtual void setSparsePattern( std::vector<size_type>&,
                                 std::vector<size_type>&
                               )
  {}
  

};

};

#endif
