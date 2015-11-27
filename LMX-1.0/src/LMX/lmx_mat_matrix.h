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

#ifndef LMXMATRIX_H
#define LMXMATRIX_H

#include<cmath>

#include"lmx_except.h"
#include"lmx_mat_type_stdmatrix.h"
#include"lmx_mat_type_csc.h"

#ifdef HAVE_GMM
#include"lmx_mat_type_gmm_sparse1.h"
#endif

//////////////////////////////////////////// Doxygen file documentation entry:
    /*!
      \file lmx_mat_matrix.h

      \brief This file contains both the declaration and implementation
       for Matrix class member and friend functions.

      \author Daniel Iglesias

    */
//////////////////////////////////////////// Doxygen file documentation (end)

/// \cond COFE
namespace cofe{
  template <int Dim, class T> class TensorRank2;
  typedef TensorRank2<3, double> Tensor;
}

namespace lmx {
/// \endcond

template <typename T> class Vector;
template <typename T> class LinearSystem;
class LMXTester;

int setMatrixType(int);
int getMatrixType();

template <typename C>
    void latexPrint( std::ofstream& os, 
                     char* mat_name, 
                     Matrix<C>& mat, 
                     int prec
                   );

    /**
    \class Matrix 
    \brief Template class Matrix

    This class permits the creation of matrix objects. A Matrix object owns two 
    parameters, nrows and mcolumns, that store the dimension of the matrix 
    container. The data is stored in an atribute (*type_matrix) that points to 
    some class which derives from the Data_mat class.

    @param mrows The number of rows of the Data Container (Data_mat).
    @param ncolumns The number of columns of the Data Container (Data_mat).
    @param reference An Elem_ref object for r/w data access.
    @param *type_matrix The pointer to the matrix data container.

    @author Daniel Iglesias.
    */
template <typename T> class Matrix{
protected:

  size_type mrows,    /**< Number of rows in matrix object. */ 
        ncolumns; /**< Number of colums in matrix object. */
  Elem_ref<T>* reference; /**< Reference pointer to an element in type_matrix.*/
  Data_mat<T>* type_matrix; /**< Pointer to the container type. */

  Matrix(int);
  Matrix(size_type, size_type, int);


public:
  friend class Vector<T>;
  friend class DenseMatrix<T>;
  friend class LinearSystem<T>;
  friend class LMXTester;

public:

  Matrix();

  Matrix(size_type, size_type);

  Matrix(const Matrix&);

  Matrix(const DenseMatrix<T>&);

  ~Matrix();

  inline void initialize_type_matrix(int);

  /**
   * Read number of rows.
   * @return Number of rows of matrix.
   */
  inline size_type rows() const { return this->mrows; }

  /**
   * Read number of cols.
   * @return Number of cols of matrix.
   */
  inline size_type cols() const { return this->ncolumns; }

  void matrixMarketLoad(char*);

  void harwellBoeingLoad(char*);

  void harwellBoeingSave(char*);

  void fillIdentity( T factor = static_cast<T>(1) );

  void fillRandom(  T factor = static_cast<T>(1) );

  //needs documentation
  void sparsePattern( Vector<size_type>&, Vector<size_type>& );

  //needs documentation
  void sparsePattern( std::vector<size_type>& , std::vector<size_type>& );

  //needs documentation
  void sparsePattern( DenseMatrix<T>& );

  inline Elem_ref<T> operator () (size_type, size_type);

  inline Matrix& operator = (const Matrix&);

  template <typename C> inline Matrix<T>& operator = (const Matrix<C>&);

  inline Matrix& operator = (const DenseMatrix<T>&);

//   template <typename C> inline Matrix<T>& operator = (const DenseMatrix<C>&);

  inline Matrix& operator += (const Matrix&);

  inline Matrix& operator += (const DenseMatrix<T>&);

  inline Matrix& operator -= (const Matrix&);

  inline Matrix& operator -= (const DenseMatrix<T>&);

  inline Matrix& operator *= (T);

  inline Matrix& operator /= (T);

  inline Matrix& add(const Matrix& A, const Matrix& B);

  inline Matrix& add(const Matrix& A, const DenseMatrix<T>& B);

  inline Matrix& add(const DenseMatrix<T>& A, const Matrix& B);

  inline Matrix& subs(const Matrix& A, const Matrix& B);

  inline Matrix& subs(const Matrix& A, const DenseMatrix<T>& B);

  inline Matrix& subs(const DenseMatrix<T>& A, const Matrix& B);

  inline Matrix& mult(const Matrix& A, const Matrix& B);

  inline Matrix& mult(const Matrix& A, const DenseMatrix<T>& B);

  inline Matrix& mult(const DenseMatrix<T>& A, const Matrix& B);

  inline Matrix& mult(const DenseMatrix<T>& A, const DenseMatrix<T>& B);

  inline Matrix& multElements(const Matrix&);

  inline Matrix& multElements(const DenseMatrix<T>&);

  inline Matrix& multElements(const Matrix<T>&, const Matrix<T>&);

  inline Matrix& multElements(const Matrix<T>&, const DenseMatrix<T>&);

  inline Matrix& multElements(const DenseMatrix<T>&, const Matrix<T>&);

  inline Matrix& multElements(const DenseMatrix<T>&, const DenseMatrix<T>&);

  /** Function that returns the value in specified position of Matrix object.
   *  */
    inline
        const T& readElement(size_type m, size_type n) const
    { return this->type_matrix->readElement(m,n); }

  /** Function writes the value in specified position of Matrix object.
   *  */
    inline
        void writeElement(T theValue, size_type m, size_type n) const
    { this->type_matrix->writeElement(theValue, m, n); }

  /** Function adds the value in specified position of Matrix object.
   *  */
    inline
        void addElement(const T theValue, size_type m, size_type n) const
    { this->type_matrix->addElement(theValue, m, n); }

  /** Cleans all numbers below given factor.
   *  */
    inline
        void clean(T factor)
    { this->type_matrix->cleanBelow(factor); }

  //needs documentation
  /** Clears the object's contents.
   *  */
    inline
        void clear()
    { this->type_matrix->clear(); }

	//begin JCGO 18/03/09
  //needs documentation
  /** Reset contents to 0
   *  */
    inline
        void reset()
    { this->type_matrix->reset(); }
	//end JCGO

  /** Resize the Matrix with given size parameters.
   *  \param i Rows.
   *  \param j Columns.
   *  */
    inline
        void resize(size_type i, size_type j)
    { mrows = i;
      ncolumns = j;
      type_matrix->resize(mrows, ncolumns);
    }

  /** Resize the Matrix to the size of another Matrix.
   *  \param mat Constant reference to a Matrix.
   *  */
    inline
        void resize( const Matrix& mat )
    { mrows = mat.rows();
      ncolumns = mat.cols();
      type_matrix->resize(mrows, ncolumns);
    }

  /**
   * Function to transpose a matrix.
   */
    inline
        Matrix& transpose()
    {
      this->type_matrix->trn();
      mrows = this->type_matrix->getRows();
      ncolumns = this->type_matrix->getCols();
      return *this;
    }

    /** Overloaded operator for adding elements between two Matrix objects.
     *  */
    Matrix operator + (const Matrix<T>& B) const
    { //Scheme: res = *this; res += B;
      Matrix<T> res( *this );
      res += B;
      return res;
    }

    /** Overloaded operator for adding elements between a Matrix and a DenseMatrix object.
     *  */
    Matrix operator + (const DenseMatrix<T>& B) const
    { //Scheme: res = *this; res += B;
      Matrix<T> res( *this );
      res += B;
      return res;
    }

    /** Overloaded operator for substracting elements between two Matrix objects.
     *  */
    Matrix<T> operator - (const Matrix<T>& B) const
    { //Scheme: res = *this; res -= B;
      Matrix<T> res( *this );
      res -= B;
      return res;
    }

    /** Overloaded operator for substracting elements between a Matrix and a DenseMatrix object.
     *  */
    Matrix operator - (const DenseMatrix<T>& B) const
    { //Scheme: res = *this; res += B;
      Matrix<T> res( *this );
      res -= B;
      return res;
    }


  /** Overloaded operator for multiplying two Matrix objects.
    *  */
    Matrix<T> operator * (const Matrix<T>& B) const
    {
      Matrix<T> res( this->mrows, B.cols() );
      res.mult(*this, B);
      return res;
    }

  /** Overloaded operator for multiplying Matrix and DenseMatrix objects.
   *  */
    Matrix<T> operator * (const DenseMatrix<T>& B) const
    {
      Matrix<T> res( this->mrows, B.cols() );
//       res.mult(*this, B);
      for (size_type i=0; i<mrows; ++i){
        for (size_type j=0; j<ncolumns; ++j){
          for (size_type k=0; k < this->ncolumns; ++k){
            res.writeElement(res.readElement(i,j) + this->readElement(i,k) * B.readElement(k,j), i, j);
          }
        }
      }
      return res;
    }

  /** Overload operator for scaling a matrix.
   *  */
    Matrix<T> operator * (const T& scalar) const
    {
      Matrix<T> res( *this );
      res *= scalar;
      return res;
    }

  /** Retuns 1 if the element exists inside the storing structure
   * and 0 if it does not.
   *  */
    inline bool exists( size_type row, size_type col )
    { return this->type_matrix->exists(row, col); }

  /** Retuns 1 if the Matrix is symmetric
   * and 0 if it does not.
   * Needs documentation!
   *  */
    inline bool checkSymmetry( )
    {
      bool result = 1;
      for (size_type i=0; i<mrows; ++i){
        for (size_type j=0; j<i; ++j){
          if( this->readElement(i,j) != this->readElement(j,i) ){
            result = 0;
            break;
          }
        }
        if( result == 0 ) break;
      }
      return result;
    }
  
  /** DOCUMENT
   */
  void factorize(){
    this->type_matrix->factorize();
  }

  /** DOCUMENT
   */
  void subsSolve( Vector<T>& rhs ){
    return this->type_matrix->subsSolve( rhs );
  }
    
/////////////////////////////// Friend functions:

  /** Function to get a transposed copy of a matrix.
   *  \param mat_in Matrix to be transposed.
   *  \return A new Matrix object with mat_in transposed.
   *  */
  friend Matrix transposed(const Matrix& mat_in)
  { Matrix<T> transpose_out(mat_in);
    transpose_out.type_matrix->trn();
    transpose_out.mrows = transpose_out.type_matrix->getRows();
    transpose_out.ncolumns = transpose_out.type_matrix->getCols();
    return transpose_out;
  }

  /** Function to typeset a matrix to a file stream.
   *  \param os File stream to write to.
   *  \param mat_name Name that will be given to matrix in typesetting. It's not the matrix variable name.
   *  \param mat The Matrix variable to be typeset.
   *  \param prec Precision of numbers (number of digits to be used in output).
   *  */
template <typename C>
  friend void lmx::latexPrint(std::ofstream& os, char* mat_name, Matrix<C>& mat, int prec);

/// \cond COFE
  friend std::ostream& operator << (std::ostream& os, const Matrix<cofe::Tensor>& mat);
/// \endcond


};

}; // namespace lmx

/////////////////////////////// Implementation of the methods defined previously

namespace lmx {

template <typename C>
  void latexPrint(std::ofstream& os, char* mat_name, Matrix<C>& mat, int prec)
  { os.precision(prec);
    os << "\\[" << mat_name << " = " << std::endl;
    os << "\\left( \\begin{array}{";
    for(int j=0; j<mat.ncolumns; ++j) { os << "c"; }
    os << "}" << std::endl;
    for(int i=0; i<mat.mrows; ++i){
      for(int j=0; j<mat.ncolumns; ++j){
        if(j != mat.ncolumns-1) { os << mat(i,j) << " & "; }
        else if (i != mat.mrows-1) { os << mat(i,j) << " \\\\ "; }
        else { os << mat(i,j) << " \\end{array} \\right) \\]" << endl; }
      }
    }
  }


  /**
   * Protected constructor for special type matrices, used in inherited Class DenseMatrix.
   */
  template <typename T>
      Matrix<T>::Matrix(int type_in) : mrows(0), ncolumns(0)
  {
    initialize_type_matrix(type_in);
  }

  /** Protected Standard constructor for DenseMatrix inherited class.
   *  \param rows Number of rows in Matrix.
   *  \param columns Number of columns in Matrix.
   *  \param type_in Experimental, don't use, please.
   */
  template <typename T>
      Matrix<T>::Matrix(size_type rows, size_type columns, int type_in) : mrows(rows), ncolumns(columns)
  {
    initialize_type_matrix(type_in);
    type_matrix->resize(mrows, ncolumns);
  }


  /** Empty constructor.
   *  */
template <typename T>
    Matrix<T>::Matrix() : mrows(0), ncolumns(0)
{
  initialize_type_matrix(getMatrixType());
}

  /** Standard constructor.
   *  \param rows Number of rows in Matrix.
   *  \param columns Number of columns in Matrix. */
template <typename T>
    Matrix<T>::Matrix(size_type rows, size_type columns) : mrows(rows), ncolumns(columns)
{ 
  initialize_type_matrix(getMatrixType());
  type_matrix->resize(mrows, ncolumns);
}


  /** Copy constructor.
   *  \param A Matrix to copy from.
   *  */
template <typename T>
    Matrix<T>::Matrix(const Matrix& A) :
 mrows(A.mrows), ncolumns(A.ncolumns)
{
  initialize_type_matrix( getMatrixType() );
  type_matrix->resize(mrows, ncolumns);
  type_matrix->equals(A.type_matrix);
}

/** Copy constructor.
 *  \param A DenseMatrix to copy from.
 *  */
template <typename T>
    Matrix<T>::Matrix(const DenseMatrix<T>& A) :
    mrows(A.mrows), ncolumns(A.ncolumns)
{
  initialize_type_matrix( getMatrixType() );
  type_matrix->resize(mrows, ncolumns); //can be time consuming for CSC (check)
  if( getMatrixType() == 1 ){ // Type is CSC:
      std::vector<size_type> ia,ja;
      A.writeSparsePattern( ia, ja );
      type_matrix->setSparsePattern( ia, ja );
  }
  for(int i=0; i<mrows; ++i){
    for(int j=0; j<ncolumns; ++j){
      this->writeElement( A.readElement(i,j), i, j);
    }
  }
}


  /** Destructor.
   *  */
template <typename T>
    Matrix<T>::~Matrix()
{ // implementar los delete
      delete this->reference;
      this->reference = 0;

      delete this->type_matrix;
      this->type_matrix = 0;
}

/**
 * Function that creates specified matrix and a new Elem_ref with it.
 * @param type Integer with matrix type identifier.
 */
template <typename T> inline
    void Matrix<T>::initialize_type_matrix(int type)
{
  switch (type) {
    case 0 :
      type_matrix = new Type_stdmatrix< T >;
    break;

    case 1 :
      type_matrix = new Type_csc< T >;
      break;

    case 2 :
#ifdef HAVE_GMM
      type_matrix = new Type_gmm< T >;
#else
      {
          std::stringstream message;
          message << "gmm++ not defined.\nYou must set \"#define HAVE_GMM\" in your file in order to use this library." << endl;
          LMX_THROW(failure_error, message.str() );
      }
#endif
    break;

    case 3 :
#ifdef HAVE_GMM
      type_matrix = new Type_gmm_sparse< T >;
#else
      {
          std::stringstream message;
          message << "gmm++ not defined.\nYou must set \"#define HAVE_GMM\" in your file in order to use this library." << endl;
          LMX_THROW(failure_error, message.str() );
      }
#endif
    break;

  }

  reference = new Elem_ref<T>(type_matrix);


}


  /** Method for reading a matrix in Matrix-Market format from a file.
   *  \param input_file Name of the file to read.
   *  */
template <typename T>
    void Matrix<T>::matrixMarketLoad(char* input_file)
{
  this->type_matrix->read_mm_file(input_file);
  this->mrows = this->type_matrix->getRows();
  this->ncolumns = this->type_matrix->getCols();
}


  /** Method for reading a matrix in Harwell-Boeing format from a file.
   *  \param input_file Name of the file to read.
   *  */
template <typename T>
    void Matrix<T>::harwellBoeingLoad(char* input_file)
{
  this->type_matrix->read_hb_file(input_file);
  this->mrows = this->type_matrix->getRows();
  this->ncolumns = this->type_matrix->getCols();
}

  /** Method for reading a matrix in Harwell-Boeing format from a file.
 *  \param input_file Name of the file to read.
   *  */
template <typename T>
    void Matrix<T>::harwellBoeingSave(char* input_file)
{
  this->type_matrix->write_hb_file(input_file);
}


/**
 * Cleans matrix and sets diagonal terms to specified value. Matrix must be square.
 * @param factor Value of diagonal terms. Default is one.
 */
template <typename T>
    void Matrix<T>::fillIdentity( T factor )
{
  if (this->ncolumns != this->mrows){
    std::stringstream message;
    message << "Trying to make identity a non-squared matrix.\nSize of matrix(" << this->mrows << ", " << this->ncolumns << ")." << endl;
    LMX_THROW(dimension_error, message.str() );
  }

  for (size_type i=0; i<this->mrows; ++i){
    for (size_type j=0; j<this->ncolumns; ++j){
      this->type_matrix->writeElement(static_cast<T>(0),i,j);
    }
  }

  for (size_type i=0; i<this->mrows; ++i){
    this->type_matrix->writeElement(static_cast<T>( factor ),i,i);
  }

}

/**
 * Function for filling a matrix with random numbers
 * @param factor Scales the random numbers (default value is unity).
 */
template <typename T>
    void Matrix<T>::fillRandom( T factor )
{
  for (size_type i=0; i<this->mrows; ++i){
    for (size_type j=0; j<this->ncolumns; ++j){
      this->type_matrix->writeElement( factor * static_cast<T>( std::rand() ) / static_cast<T>(RAND_MAX),i,j);
    }
  }
}

/**
 * \brief Function for preparing a sparse non-zero pattern in matrix.
 * Uses a Harwell-Boeing (CSC) like vectors for describing the pattern.
 * Specially indicated for preparing the CSC matrix for efficient writing. Does nothing in the rest of matrix types.
 * @param row_index Position of elements in rows.
 * @param col_index Position of the first element in each column.
 */
template <typename T>
    void Matrix<T>::sparsePattern( Vector<size_type>& row_index,
                                   Vector<size_type>& col_index 
                                 )
{
  this->type_matrix->setSparsePattern( row_index, col_index );
}

/**
 * \brief Function for preparing a sparse non-zero pattern in matrix.
 * Uses a Harwell-Boeing (CSC) like vectors for describing the pattern.
 * Specially indicated for preparing the CSC matrix for efficient writing. Does nothing in the rest of matrix types.
 * @param row_index Position of elements in rows.
 * @param col_index Position of the first element in each column.
 */
template <typename T>
    void Matrix<T>::sparsePattern( std::vector<size_type>& row_index,
                                   std::vector<size_type>& col_index
                                 )
{
  this->type_matrix->setSparsePattern( row_index, col_index );
}

/**
 * \brief Function for preparing a sparse non-zero pattern in matrix.
 * Uses a Harwell-Boeing (CSC) like vectors for describing the pattern.
 * Specially indicated for preparing the CSC matrix for efficient writing. Does nothing in the rest of matrix types.
 * @param dense_matrix is the matrix with nonzero values.
 */
template <typename T>
    void Matrix<T>::sparsePattern( lmx::DenseMatrix<T>& aDenseMatrix )
{
  std::vector<size_type> row_index, col_index;
  aDenseMatrix.writeSparsePattern( row_index, col_index );
  this->type_matrix->setSparsePattern( row_index, col_index );
}

/** Overloaded operator for extracting elements from the Matrix object.
 *  \param m Row position of element.
 *  \param n Column position of element.
 */
template <typename T>
    inline
    Elem_ref<T> Matrix<T>::operator () (size_type m, size_type n)
{
  if (m >= this->mrows){
    std::stringstream message;
    message << "Row index exceeds matrix dimensions. Trying to access element (" << m << ", " << n << "of matrix with dimension " << this->mrows << "x" << this->ncolumns << "." << endl;
    LMX_THROW(dimension_error, message.str() );
  }
  else if (n >= this->ncolumns){
    std::stringstream message;
    message << "Column index exceeds matrix dimensions. Trying to access element (" << m << ", " << n << "of matrix with dimension " << this->mrows << "x" << this->ncolumns << "." << endl;
    LMX_THROW(dimension_error, message.str() );
  }

  reference->write_pos(m,n);

  return *reference;
}

/** Overloaded operator for equaling every element between two
 *  Matrix objects of the same type.
 *  \param A Matrix to be equal to.
 */
template <typename T>
    inline
    Matrix<T>& Matrix<T>::operator = (const Matrix<T>& A)
{ 
  mrows = A.mrows;
  ncolumns = A.ncolumns;
  type_matrix->equals(A.type_matrix);
  return *this;
}

/** Overloaded operator for equaling every element between two
 *  Matrix objects of different types.
 *  \param A Matrix to be equal to.
 */
template <typename T>
    template <typename C>
    inline
    Matrix<T>& Matrix<T>::operator = (const Matrix<C>& A)
{
  mrows = A.rows();
  ncolumns = A.cols();
  this->type_matrix->resize(mrows,ncolumns);
  for (int i=0; i<mrows ; ++i){
    for (int j=0; j<ncolumns ; ++j){
      this->operator ( )(i,j) = A.readElement(i,j);
    }
  }

  return *this;
}


/** Overloaded operator for equaling every element between a
 *  Matrix and a DenseMatrix object of the same type.
 *  \param A DenseMatrix to be equal to.
 */
template <typename T>
    inline
    Matrix<T>& Matrix<T>::operator = (const DenseMatrix<T>& A)
{
  mrows = A.mrows;
  ncolumns = A.ncolumns;

  switch ( getMatrixType() ) {
//     case 0 :
//       ;
//       break;
    case 1 :
      copy<T>( static_cast<const Type_stdmatrix<T>*>(A.type_matrix),
               static_cast<Type_csc<T>*>(this->type_matrix) );
      break;

    default :
      for (size_type i=0; i<mrows; ++i){
        for (size_type j=0; j<mrows; ++j){
          this->writeElement(A.readElement(i,j), i, j);
        }
      }
      break;
  }

  return *this;
}

/** Overloaded operator for assigning the addition of elements
 *  between two Matrix objects.
 *  \param A Matrix to be added.
 */
template <typename T>
    inline
    Matrix<T>& Matrix<T>::operator += (const Matrix& A)
{ 
// Scheme of function: (*this.type_matrix) + A.type_matrix, return *this;
  type_matrix->add(A.type_matrix);
  return *this;
}

/** Overloaded operator for assigning the addition of elements
 *  between a Matrix and a DenseMatrix object.
 *  \param A DenseMatrix to be added.
 */
template <typename T>
    inline
    Matrix<T>& Matrix<T>::operator += (const DenseMatrix<T>& A)
{
// Scheme of function: (*this.type_matrix) + A.type_matrix, return *this;
  for (size_type i=0; i<mrows; ++i){
    for (size_type j=0; j<mrows; ++j){
      this->writeElement(this->readElement(i,j) + A.readElement(i,j), i, j);
    }
  }
  return *this;
}

/** Overloaded operator for assigning the substraction of elements
 *  between two Matrix objects.
 *  \param A Matrix to be substracted.
 */
template <typename T>
    inline
    Matrix<T>& Matrix<T>::operator -= (const Matrix& A)
{ 
// Scheme of function: (*this.type_matrix) - A.type_matrix, return *this;
  type_matrix->substract(A.type_matrix);
  return *this;
}

/** Overloaded operator for assigning the substraction of elements
 *  between a Matrix and a DenseMatrix object.
 *  \param A DenseMatrix to be substracted.
 */
template <typename T>
    inline
    Matrix<T>& Matrix<T>::operator -= (const DenseMatrix<T>& A)
{
// Scheme of function: (*this.type_matrix) + A.type_matrix, return *this;
  for (size_type i=0; i<mrows; ++i){
    for (size_type j=0; j<mrows; ++j){
      this->writeElement(this->readElement(i,j) - A.readElement(i,j), i, j);
    }
  }
  return *this;
}

/**
 * Scalar multiplication operator.
 * @param scalar Factor of scaling.
 * @return Reference to scaled matrix.
 */
template <typename T> inline
    Matrix<T>& Matrix<T>::operator *= (T scalar)
{
  type_matrix->multiplyScalar(scalar);
  return *this;
}

/**
 * Scalar term-by-term division operator.
 * @param scalar Factor of division.
 * @return Reference to scaled matrix.
 */
template <typename T> inline
  Matrix<T>& Matrix<T>::operator /= (T scalar)
{
  T div = (T)(1./div);
  type_matrix->multiplyScalar(div);
  return *this;
}

/**
 * Matrix addition.
 * @param A Left Matrix.
 * @param B Right Matrix.
 * @return Reference to addition result.
 */
template <typename T> inline
    Matrix<T>& Matrix<T>::add(const Matrix<T>& A, const Matrix<T>& B)
{
  this->type_matrix->equals(A.type_matrix);
  this->type_matrix->add(B.type_matrix);
  return *this;
}

/**
 * Matrix and DenseMatrix addition.
 * @param A Left Matrix.
 * @param B Right DenseMatrix.
 * @return Reference to addition result.
 */
template <typename T> inline
    Matrix<T>& Matrix<T>::add(const Matrix<T>& A, const DenseMatrix<T>& B)
{
  mat_mat_add( A.type_matrix, B.type_matrix, this->type_matrix );
  return *this;
}

/**
 * DenseMatrix and Matrix addition.
 * @param A Left DenseMatrix.
 * @param B Right Matrix.
 * @return Reference to addition result.
 */
template <typename T> inline
    Matrix<T>& Matrix<T>::add(const DenseMatrix<T>& A, const Matrix<T>& B)
{
  mat_mat_add( A.type_matrix, B.type_matrix, this->type_matrix );
  return *this;
}


/**
 * Matrix substraction.
 * @param A Left Matrix.
 * @param B Right Matrix.
 * @return Reference to substraction result.
 */
template <typename T> inline
    Matrix<T>& Matrix<T>::subs(const Matrix<T>& A, const Matrix<T>& B)
{
  this->type_matrix->equals(A.type_matrix);
  this->type_matrix->substract(B.type_matrix);
  return *this;
}

/**
 * Matrix and DenseMatrix substraction.
 * @param A Left Matrix.
 * @param B Right DenseMatrix.
 * @return Reference to substraction result.
 */
template <typename T> inline
    Matrix<T>& Matrix<T>::subs(const Matrix<T>& A, const DenseMatrix<T>& B)
{
  mat_mat_subs( A.type_matrix, B.type_matrix, this->type_matrix );
  return *this;
}

/**
 * DenseMatrix and Matrix substraction.
 * @param A Left DenseMatrix.
 * @param B Right Matrix.
 * @return Reference to substraction result.
 */
template <typename T> inline
    Matrix<T>& Matrix<T>::subs(const DenseMatrix<T>& A, const Matrix<T>& B)
{
  mat_mat_subs( A.type_matrix, B.type_matrix, this->type_matrix );
  return *this;
}

/**
 * Matrix multiplication.
 * @param A Left Matrix.
 * @param B Right Matrix.
 * @return Reference to multiplication result.
 */
template <typename T> inline
    Matrix<T>& Matrix<T>::mult(const Matrix<T>& A, const Matrix<T>& B)
{
  this->type_matrix->multiply(A.type_matrix, B.type_matrix);
  return *this;
}

/**
 * Matrix and DenseMatrix multiplication.
 * @param A Left Matrix.
 * @param B Right DenseMatrix.
 * @return Reference to multiplication result.
 */
template <typename T> inline
    Matrix<T>& Matrix<T>::mult(const Matrix<T>& A, const DenseMatrix<T>& B)
{
  mat_mat_mult( A.type_matrix, B.type_matrix, this->type_matrix );
  return *this;
}

/**
 * DenseMatrix and Matrix multiplication.
 * @param A Left DenseMatrix.
 * @param B Right Matrix.
 * @return Reference to multiplication result.
 */
template <typename T> inline
    Matrix<T>& Matrix<T>::mult(const DenseMatrix<T>& A, const Matrix<T>& B)
{
  mat_mat_mult( A.type_matrix, B.type_matrix, this->type_matrix );
  return *this;
}

/**
 * DenseMatrix objects multiplication, with result storing in the Matrix object.
 * @param A Left DenseMatrix.
 * @param B Right DenseMatrix.
 * @return Reference to multiplication result in Matrix.
 */
template <typename T> inline
    Matrix<T>& Matrix<T>::mult(const DenseMatrix<T>& A, const DenseMatrix<T>& B)
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
    Matrix<T>& Matrix<T>::multElements(const Matrix<T>& B)
{
  this->type_matrix->multiplyElements(B.type_matrix);
  return *(this);
}

/**
 * Internal product between matrices. Multiplies each element of object to its equivalent in Matrix B.
 * @param B Matrix to multiply to.
 * @return Reference to internal product result.
 */
template <typename T> inline
    Matrix<T>& Matrix<T>::multElements(const DenseMatrix<T>& B)
{
  mat_mat_multElements( B.type_matrix, this->type_matrix );
  return *(this);
}


/**
 * Internal product between matrices. Multiplies each element in Matrix A to its equivalent in Matrix B and saves the result in object from whitch the function is invoked.
 * @param A Matrix to multiply.
 * @param B Matrix to multiply.
 * @return Reference to internal product result.
 */
template <typename T> inline
    Matrix<T>& Matrix<T>::multElements(const Matrix<T>& A, const Matrix<T>& B)
{
  // Must select the matrix to multiply to because can be the same object! so...
  if(this->type_matrix == A.type_matrix)
    this->type_matrix->multiplyElements(B.type_matrix);
  else if(this->type_matrix == B.type_matrix)
    this->type_matrix->multiplyElements(A.type_matrix);
  else{ // no problem different matrices...
    this->type_matrix->equals(A.type_matrix);
    this->type_matrix->multiplyElements(B.type_matrix);
  }
  return *(this);
}

/**
 * Internal product between a Matrix and a DenseMatrix object.
 * Multiplies each element in Matrix A to its equivalent in DenseMatrix B and saves the result in object from whitch the function is invoked.
 * @param A Matrix to multiply.
 * @param B DenseMatrix to multiply.
 * @return Reference to internal product result.
 */
template <typename T> inline
    Matrix<T>& Matrix<T>::multElements(const Matrix<T>& A, const DenseMatrix<T>& B)
{
  mat_mat_multElements( A.type_matrix, B.type_matrix, this->type_matrix );
  return *(this);
}

/**
 * Internal product between a DenseMatrix and a Matrix object.
 * Multiplies each element in DenseMatrix A to its equivalent in Matrix B and saves the result in object from whitch the function is invoked.
 * @param A DenseMatrix to multiply.
 * @param B Matrix to multiply.
 * @return Reference to internal product result.
 */
template <typename T> inline
    Matrix<T>& Matrix<T>::multElements(const DenseMatrix<T>& A, const Matrix<T>& B)
{
  mat_mat_multElements( A.type_matrix, B.type_matrix, this->type_matrix );
  return *(this);
}

/**
 * Internal product between two DenseMatrix objects.
 * Multiplies each element in DenseMatrix A to its equivalent in DenseMatrix B and saves the result in object from whitch the function is invoked.
 * @param A DenseMatrix to multiply.
 * @param B DenseMatrix to multiply.
 * @return Reference to internal product result.
 */
template <typename T> inline
    Matrix<T>& Matrix<T>::multElements(const DenseMatrix<T>& A, const DenseMatrix<T>& B)
{
  mat_mat_multElements( A.type_matrix, B.type_matrix, this->type_matrix );
  return *(this);
}



 //////////////////////////////////////////////////////////////////
 // Overloaded operators with Matrix as rvalue:  //////////////////
 //////////////////////////////////////////////////////////////////

  /**
   * Overload operator for negation.
   */
  template <typename T>
      Matrix<T> operator - (const Matrix<T>& B)
  {
    /// Scheme of function: res=1*-A, return res;
    Matrix<T> res(B);
    res *= -1;
    return res;
  }

  /**
   * Overload operator for multiplying a scalar and a Matrix object.
   */
  template <typename T>
      Matrix<T> operator * (const T& scalar, const Matrix<T>& B)
  {
      /// Scheme of function: mult=A*B, return B;
    Matrix<T> mult(B);
    mult *= scalar;
    return mult;
  }

  /** Overloaded operator that prints a Matrix<T> object into the output stream selected.
   * \param mat Matrix<T> reference to the object that is going to be printed.
   * \param os An output std stream.
   * \return (\a os ) The stream that was passed as a parameter.
   */
  template <class T>
  std::ostream& operator << (std::ostream& os, const Matrix<T>& mat)
  {
    os << "Matrix (" << mat.rows() << "," << mat.cols() << ") = " ;
    if ( !mat.rows() && !mat.cols() ) os << "void";
    for (size_type i=0; i<mat.rows(); ++i)
    {
      os << endl;
      for (size_type j=0; j<mat.cols(); ++j)
        os << mat.readElement(i,j) << " ";
    }
    os << endl;
    return os;
  }

  /**
   * Overloaded operator that prints a Matrix<T> object into the output stream selected. Specialized for cofe::TensorRank2<3> order 2 Tensors.
   * \param mat Matrix<T> reference to the object that is going to be printed.
   * \param os An output std stream.
   * \return (\a os ) The stream that was passed as a parameter.
   */
  template <class T>
  std::ostream& operator << (std::ostream& os, const Matrix< cofe::TensorRank2<3, T> >& mat)
  {
    os << "Matriz (" << mat.rows() << "," << mat.cols() << ") = " ;
    if ( !mat.rows() && !mat.cols() ) os << "void";
      os << endl;
    for (size_type i=0; i<mat.rows(); ++i)
    {
      os << endl;
      for (size_type j=0; j<mat.cols(); ++j)
        os << mat.readElement(i,j) << " ";
    }
    os << endl;
    return os;
  }

}; // namespace lmx


#endif
