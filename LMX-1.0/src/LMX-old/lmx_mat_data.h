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
#ifndef LMXDATA_H
#define LMXDATA_H

#ifdef HAVE_GMM
  #include "gmm/gmm.h"
#endif

#include "lmx_mat_data_blas.h"

//////////////////////////////////////////// Doxygen file documentation entry:
    /*!
      \file lmx_mat_data.h
      
      \brief This file contains the declaration of the data class' pure virtual functions.
      
      For classes derived from Data pure virtual class, all these methods must be implemented. Thus, for comprobation and checking, the methods here declared are as well documented.
      
      \author Daniel Iglesias
      
    */
//////////////////////////////////////////// Doxygen file documentation (end)

/** Member size_type is a redefinition of the numeric type used for indexing the matrices. Actually, it's defined as size_t integer.
*   */
typedef size_t size_type;

namespace lmx {

    /**
    \class Data
    \brief Template class Data.
    Mother class for Data_mat & Data_vec.
    
    This class represents the skeleton for the data container used by the Vector and Matrix classes. No parameter nor function implementation here, just pure virtual class. See derived classes for details in implementation. Also maybe useful to see how this class is used in the Vector class.
    
    @author Daniel Iglesias.
    */
template <typename T> class Data{
public:

  /** Empty constructor.
   */
  Data(){}

  /** Destructor. 
   */
  virtual ~Data(){}

  /** Resize method for augmenting or reducing the container's dimension.
   */
  virtual void resize(size_type, size_type) = 0;

  /** Read method for accesing stored data.
   * \returns Value of the data stored in the (size_type,size_type) position.
   */
  virtual const T& readElement(const size_type&, const size_type&) const = 0;

  /** Write method for storing data.
   */
  virtual void writeElement(T, const size_type&, const size_type&) = 0;

  /** Write method for adding to data.
   */
  virtual void addElement(T, const size_type&, const size_type&) = 0;

  /** Method for knowing the number of data rows. 
   * \returns Number of rows.
   */
  virtual size_type getRows() const = 0;

  /** Method for knowing the number of data columns. 
   * \returns Number of columns.
   */
  virtual size_type getCols() const = 0;

  /** Method for equaling object's data to the parameter's data.
   * Necessary for implementing the "=" overloaded operator method in Vector and Vector class types.
   */
  virtual void equals(const Data*) = 0;

  /** Method for adding object's data to the parameter's data.
   * Necessary for implementing the "+=" overloaded operator method in Vector and Vector class types.
   */
  virtual void add(const Data*) = 0;

  /** Method for substracting object's data from the parameter's data.
   * Necessary for implementing the "-=" overloaded operator method in Vector and Vector class types.
   */
  virtual void substract(const Data*) = 0;

  /** Method for multiplying a scalar with an array and saving the result in the object's data.
   * Necessary for implementing the "*" overloaded operator method with scalars in Vector class type.
   */
  virtual void multiplyScalar(const T&) = 0;
  
  /** Method multiplying element-by-element of two arrays. One would be the object's contents and the other the parameter's contents.
   * Necessary for implementing  Vector to Vector multElements.
   */
  virtual void multiplyElements(const Data*) = 0;
  
  /**
   * Method for cleaning all elements below a given factor.
   */
  virtual void cleanBelow(const double) = 0;

  /**
   * Method for clearing all elements.
   */
  virtual void clear() = 0;

	//begin JCGO 18/03/09
  /**
   * Method for all elements to 0
   */
  virtual void reset() = 0;
  //end JCGO

};

};

#endif
