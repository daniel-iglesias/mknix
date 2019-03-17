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

#ifndef LMXDENSE_MATRIX_H
#define LMXDENSE_MATRIX_H

#include"lmx_mat_matrix.h"

//////////////////////////////////////////// Doxygen file documentation entry:
    /*!
      \file lmx_mat_dense_matrix.h
      
      \brief This file contains both the declaration and implementation for DenseMatrix class member and friend functions.
      
      \author Daniel Iglesias 

    */
//////////////////////////////////////////// Doxygen file documentation (end)
    
namespace lmx {

template <typename T> class Matrix;

    /**
    \class DenseMatrix
    \brief Template class DenseMatrix

    This class allows the creation of dense matrix objects. A DenseMatrix object owns two parameters, nrows and mcolumns, that store the dimension of the matrix container. The data is stored in an atribute (*type_matrix) that points to some class which derives from the Data_mat class.

    @param mrows The number of rows of the Data Container (Data_mat).
    @param ncolumns The number of columns of the Data Container (Data_mat).
    @param reference An Elem_ref object for r/w data access.
    @param *type_matrix The pointer to the matrix data container.

    @author Daniel Iglesias .
    */
template <typename T> class DenseMatrix : public Matrix<T>{
private:
  int type;
  Type_stdmatrix<T>* std_matrix;

public:

  DenseMatrix();

  DenseMatrix(size_type, size_type);

  DenseMatrix(const DenseMatrix&);

  ~DenseMatrix();

  inline DenseMatrix& mult ( const Matrix<T>&, const Matrix<T>&);

  inline DenseMatrix& operator = ( const Matrix<T>& );

  inline DenseMatrix& multElements ( const Matrix<T>& );

  inline DenseMatrix& multElements ( const Matrix<T>&, const Matrix<T>&);

      /** Overloaded operator for adding elements between a DenseMatrix and a Matrix object.
       *  */
  DenseMatrix operator + (const Matrix<T>& B) const
  { //Scheme: res = *this; res += B;
    DenseMatrix<T> res( *this );
    res += B;
    return res;
  }

      /** Overloaded operator for adding elements between a DenseMatrix and a Matrix object.
       *  */
  DenseMatrix operator - (const Matrix<T>& B) const
  { //Scheme: res = *this; res -= B;
    DenseMatrix<T> res( *this );
    res -= B;
    return res;
  }

  /** Overloaded operator for multiplying DenseMatrix and Matrix objects.
   *  */
  DenseMatrix<T> operator * (const Matrix<T>& B) const
  {
    DenseMatrix<T> res( this->mrows, B.cols() );
//       res.mult(*this, B);
    for (size_type i=0; i<this->mrows; ++i){
      for (size_type j=0; j<this->ncolumns; ++j){
        for (size_type k=0; k < this->ncolumns; ++k){
//             cout << "k = " << k << ", res(i,j) = " << res.readElement(i,j) << ", A(i,k) = " << this->readElement(i,k) << ", B(k,j) = " << B.readElement(k,j) << endl;
          res.writeElement(res.readElement(i,j) + this->readElement(i,k) * B.readElement(k,j), i, j);
        }
      }
    }
    return res;
  }

  /** Function that returns the value in specified position of Matrix object.
   *  */
    inline
        const T& readElement(size_type m, size_type n) const
    { return std_matrix->contents[m][n]; }

  /** Function writes the value in specified position of Matrix object.
   *  */
    inline
        void writeElement(T theValue, size_type m, size_type n) const
    { std_matrix->contents[m][n] = theValue; }

  /** Function adds the value in specified position of Matrix object.
   *  */
    inline
        void addElement(const T theValue, size_type m, size_type n) const
    { std_matrix->contents[m][n] += theValue; }

  /** Writes in two STL vectors the HB column and row indexes of the Non-zeros.
   *  */
        void writeSparsePattern( std::vector<size_type>&,
                                 std::vector<size_type>&
                               ) const;


};

}; // namespace lmx

/////////////////////////////// Implementation of the methods defined previously

namespace lmx {

  /** Empty constructor.
   *  */
template <typename T>
DenseMatrix<T>::DenseMatrix()
  : Matrix<T>(0)
{
  std_matrix = static_cast<Type_stdmatrix<T>*> (this->type_matrix);
}

  /** Standard constructor.
   *  \param rows Number of rows in DenseMatrix.
   *  \param columns Number of columns in DenseMatrix. */
template <typename T>
DenseMatrix<T>::DenseMatrix(size_type rows, size_type columns)
 : Matrix<T>(rows, columns, 0)
{
  std_matrix = static_cast<Type_stdmatrix<T>*> (this->type_matrix);
}

  /** Copy constructor.
   *  \param A DenseMatrix to copy from.
   *  */
template <typename T>
DenseMatrix<T>::DenseMatrix(const DenseMatrix& A) :
 Matrix<T>(0)
{ this->mrows = A.rows();
  this->ncolumns = A.cols();
  this->type_matrix->resize(this->mrows, this->ncolumns);
  this->type_matrix->equals(A.type_matrix);
  std_matrix = static_cast<Type_stdmatrix<T>*> (this->type_matrix);
  
}


  /** Destructor.
   *  */
template <typename T>
DenseMatrix<T>::~DenseMatrix()
{
}

/** Overloaded operator for equaling every element between a
 *  DenseMatrix and a Matrix object of the same type.
 *  \param A Matrix to be equal to.
 *  */
template <typename T>
    inline
    DenseMatrix<T>& DenseMatrix<T>::operator = (const Matrix<T>& A)
{
  this->mrows = A.rows();
  this->ncolumns = A.cols();
  for (size_type i=0; i<this->mrows; ++i){
    for (size_type j=0; j<this->ncolumns; ++j){
      this->writeElement(A.readElement(i,j), i, j);
    }
  }
  return *this;
    
}

/** Overload operator for negation.
 *  */
template <typename T>
    DenseMatrix<T> operator - (const DenseMatrix<T>& B)
{
  /// Scheme of function: res=1*-A, return res;
  DenseMatrix<T> res(B);
  res *= -1;
  return res;
}


/**
 * Matrix objects multiplication, with result storing in the DenseMatrix object.
 * @param A Left Matrix.
 * @param B Right Matrix.
 * @return Reference to multiplication result in DenseMatrix.
 */
template <typename T> inline
    DenseMatrix<T>& DenseMatrix<T>::mult(const Matrix<T>& A, const Matrix<T>& B)
{
  mat_mat_mult( A.type_matrix, B.type_matrix, this->type_matrix );
  return *this;
}

/**
 * Internal product between matrices. Multiplies each element of object to its equivalent in Matrix B.
 * @param B Matrix to multiply to.
 * @return Reference to internal product result.
 */
template <typename T> inline
    DenseMatrix<T>& DenseMatrix<T>::multElements(const Matrix<T>& B)
{
  mat_mat_multElements( B.type_matrix, this->type_matrix );
  return *(this);
}

/**
 * Internal product between two Matrix objects. Multiplies each element in Matrix A to its equivalent in Matrix B and saves the result in object from whitch the function is invoked.
 * @param A Matrix to multiply.
 * @param B Matrix to multiply.
 * @return Reference to internal product result.
 */
template <typename T> inline
    DenseMatrix<T>& DenseMatrix<T>::multElements(const Matrix<T>& A, const Matrix<T>& B)
{
  mat_mat_multElements( A.type_matrix, B.type_matrix, this->type_matrix );
  return *(this);
}

template <typename T>
        void DenseMatrix<T>::writeSparsePattern
        (std::vector<size_type>& nRows,
         std::vector<size_type>& nCols) const
{ int i,j, elem_counter(1);
    for(j=0; j < this->ncolumns; ++j){
        nCols.push_back( elem_counter );
        for(i=0; i < this->mrows; ++i){
            if (this->readElement(i,j) != T(0)){
                nRows.push_back( i+1 );
                ++elem_counter;
            }
        }
    }
    nCols.push_back( elem_counter );
}


}; // namespace lmx


#endif
