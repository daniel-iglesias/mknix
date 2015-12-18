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
#ifndef LMXELEM_REF_H
#define LMXELEM_REF_H

#include "lmx_mat_data.h"

//////////////////////////////////////////// Doxygen file documentation entry:
    /*!
      \file lmx_mat_elem_ref.h

      \brief This file contains both the declaration and implementation for Elem_ref class.

      Elem_ref class is used by Matrix and Vector objects as a parameter.

      \author Daniel Iglesias 
      
    */
//////////////////////////////////////////// Doxygen file documentation (end)

namespace lmx {

int getMatrixType();

    /**
    \class Elem_ref
    \brief Template class Elem_ref

    This class is used to gain access for read and write purposes in any data structure derived from the Data class.

    @param type Pointer to the matrix or vector structure that is a Data derived class.
    @param value Copy of the numerical value read.
    @param i Row where the original value is stored.
    @param j Column where the original value is stored.

    @author Daniel Iglesias .
    */

template <typename T> class Elem_ref{
private:
  Data<T>* type;
  T value;
  size_type i;
  size_type j;

public:
  Elem_ref();

  /**
   * Standard constructor
   * @param type_in Pointer to Matrix or Vector's data.
   */
  Elem_ref(Data<T>* type_in) :
   type(type_in)
  {
  }

  /**
   * Destructor
   */
  ~Elem_ref()
  {}

  /**
   * Cast operator. Converts the Elem_ref to a T type if its used as a right-hand-side argument in assignment.
   * @return Value pointed by type in position (i,j).
   */
  operator T() { return this->type->readElement(i,j); }

/**
 * Method for storing a given position.
 * @param m Row number.
 * @param n Column number.
 */
  void write_pos(size_type m, size_type n)
  { this->i = m;
    this->j = n;
  }

  inline void operator = (const Elem_ref& ref_in);

  inline T operator = (T data_in);

  inline T operator += (T data_in);

  inline T operator -= (T data_in);

  inline T operator *= (T data_in);

  inline T operator /= (T data_in);

  inline T operator ++ ();

  inline T operator ++ (int);

  inline T operator -- ();

  inline T operator -- (int);

};

/**
 * Empty constructor
 */
template <typename T>
    Elem_ref<T>::Elem_ref()
{}

/**
 * Overloaded operator for equaling the element in position (i,j) in
 * type structure to a given element reference.
 * @param ref_in Element reference.
 * @return Value of element.
 */
template <typename T> inline
    void Elem_ref<T>::operator = (const Elem_ref& ref_in)
{ type->writeElement( ref_in.type->readElement(ref_in.i, ref_in.j), i, j); }

/**
 * Overloaded operator for equaling the element in position (i,j) in
 * type structure to a given value.
 * @param data_in Value to equal to.
 * @return A reference to the value.
 */
template <typename T> inline
    T Elem_ref<T>::operator = (T data_in)
{
  type->writeElement( data_in, i, j );
  return data_in;
}

/** Overloaded operator for adding the element in position (i,j) in 
 * type structure the given value.
 * @param data_in Value to add.
 * @return A reference to the final value.
 */
template <typename T> inline
    T Elem_ref<T>::operator += (T data_in)
{
  type->writeElement( type->readElement(i,j) + data_in, i, j );
  return type->readElement(i,j);
}

/** Overloaded operator for substracting the element in position (i,j) in 
 * type structure the given value.
 * @param data_in Value to substract.
 * @return A reference to the final value.
 */
template <typename T> inline
    T Elem_ref<T>::operator -= (T data_in)
{
  type->writeElement( type->readElement(i,j) - data_in, i, j );
  return type->readElement(i,j);
}

/** Overloaded operator for multiplying the element in position (i,j) in 
 * type structure by the given value.
 * @param data_in Value to multiply by.
 * @return A reference to the final value.
 */
template <typename T> inline
    T Elem_ref<T>::operator *= (T data_in)
{
  type->writeElement( type->readElement(i,j) * data_in, i, j );
  return type->readElement(i,j);
}

/** Overloaded operator for dividing the element in position (i,j) in 
 * type structure by the given value.
 * @param data_in Value to divide by.
 * @return A reference to the final value.
 */
template <typename T> inline
    T Elem_ref<T>::operator /= (T data_in)
{
  type->writeElement( type->readElement(i,j) / data_in, i, j );
  return type->readElement(i,j);
}


/**
 * Left increment operator.
 * @return Incremented value.
 */
template <typename T> inline
    T Elem_ref<T>::operator ++ ()
{
  type->writeElement( type->readElement(i,j)+1, i, j);
  return type->readElement(i,j);
}

/**
 * Right increment operator.
 * @return Not incremented value.
 */
template <typename T> inline
    T Elem_ref<T>::operator ++ (int foo)
{ T temp = type->readElement(i,j);
  type->writeElement( temp+1, i, j);
  return temp;
}

/**
 * Left decrement operator.
 * @return Decremented value.
 */
template <typename T> inline
    T Elem_ref<T>::operator -- ()
{
  type->writeElement( type->readElement(i,j)-1, i, j);
  return type->readElement(i,j);
}

/**
 * Right decrement operator.
 * @return Not decremented value.
 */
template <typename T> inline
    T Elem_ref<T>::operator -- (int foo)
{ T temp = type->readElement(i,j);
  type->writeElement( temp-1, i, j);
  return temp;
}



}

#endif
