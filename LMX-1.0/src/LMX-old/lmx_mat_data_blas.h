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
#ifndef LMXDATA_BLAS_H
#define LMXDATA_BLAS_H

#include<algorithm>

//////////////////////////////////////////// Doxygen file documentation entry:
    /*!
      \file lmx_mat_data_blas.h
      
      \brief Basic Linear Algebra Methods for Data_mat and Data_vec derived classes.
      
      Basic linear algebra numerical methodso including specialization for some Matrix and Vector data types.
      
      \author Daniel Iglesias 
      
    */
//////////////////////////////////////////// Doxygen file documentation (end)

/**
 * Integer for indexing matrices and vectors.
 */
typedef size_t size_type;


namespace lmx{

template <typename T> class Data_mat;
template <typename T> class Data_vec;

template <typename T> class Type_stdmatrix;
template <typename T> class Type_csc;
template <typename T> class Type_stdVector;
#ifdef HAVE_GMM
template <typename T> class Type_gmm;
template <typename T> class Type_gmm_sparse;
template <typename T> class Type_gmmVector_sparse;
#endif


// Matrix generalized methods (not very efficient but robust):

/**
 * Matrix and DenseMatrix compatible addition.
 * Calculates the operation A+B = C.
 * @param A LHS Data_mat<T>.
 * @param B RHS Data_mat<T>.
 * @param C Result Data_mat<T>, C=A+B.
 */
template <typename T>
    void mat_mat_add( const Data_mat<T>* A,
                       const Data_mat<T>* B,
                       Data_mat<T>* C)
{
  for (size_type i=0; i < C->getRows(); ++i){
    for (size_type j=0; j < C->getCols(); ++j){
      C->writeElement( A->readElement(i,j) + B->readElement(i,j) , i, j );
    }
  }
}

/**
 * Matrix and DenseMatrix compatible substraction.
 * Calculates the operation A-B = C.
 * @param A LHS Data_mat<T>.
 * @param B RHS Data_mat<T>.
 * @param C Result Data_mat<T>, C=A-B.
 */
template <typename T>
    void mat_mat_subs( const Data_mat<T>* A,
                       const Data_mat<T>* B,
                       Data_mat<T>* C)
{
  for (size_type i=0; i < C->getRows(); ++i){
    for (size_type j=0; j < C->getCols(); ++j){
      C->writeElement( A->readElement(i,j) - B->readElement(i,j) , i, j );
    }
  }
}

/**
 * Matrix and DenseMatrix compatible multiplication.
 * Calculates the product A*B = C.
 * @param A LHS Data_mat<T>.
 * @param B RHS Data_mat<T>.
 * @param C Result Data_mat<T>, C=A*B.
 */
template <typename T>
    void mat_mat_mult( const Data_mat<T>* A,
                       const Data_mat<T>* B,
                       Data_mat<T>* C)
{
  // Emmit an error if C is A or B...
    if(C == A || C == B){
      std::stringstream message;
      message << "Trying to multiply and save results on same data at the same time."
	    << endl << "  This cannot be done." << endl;
      LMX_THROW(failure_error, message.str() );
    }
  size_type i, j, k;
  C->resize( A->getRows(), B->getCols() );
  for (i=0; i< C->getRows(); ++i){
    for (j=0; j< C->getCols(); ++j){
      C->writeElement( T(0), i, j );
	}
  }
  // This can have a good optimization...
  for (i=0; i < C->getRows(); ++i){
    for (j=0; j < C->getCols(); ++j){
      C->writeElement( 0, i, j );
      for (k=0; k < A->getCols(); ++k){
        C->writeElement( C->readElement(i,j) + A->readElement(i,k) * B->readElement(k,j) , i, j );
      }
    }
  }
}

/**
 * Matrix and DenseMatrix compatible element-by-element multiplication with one operand.
 * Calculates the operation C(i,j) *= A(i,j).
 * @param A Data_mat<T>.
 * @param C Result Data_mat<T>.
 */
template <typename T>
    void mat_mat_multElements( const Data_mat<T>* A,
                           Data_mat<T>* C)
{
  for (size_type i=0; i < C->getRows(); ++i){
    for (size_type j=0; j < C->getCols(); ++j){
      C->writeElement( C->readElement(i,j) * A->readElement(i,j) , i, j );
    }
  }
}

/**
 * Matrix and DenseMatrix compatible element-by-element multiplication.
 * Calculates the operation A(i,j) * B(i,j) = C(i,j).
 * @param A LHS Data_mat<T>.
 * @param B RHS Data_mat<T>.
 * @param C Result Data_mat<T>.
 */
template <typename T>
    void mat_mat_multElements( const Data_mat<T>* A,
                           const Data_mat<T>* B,
                           Data_mat<T>* C)
{
  for (size_type i=0; i < C->getRows(); ++i){
    for (size_type j=0; j < C->getCols(); ++j){
      C->writeElement( A->readElement(i,j) * B->readElement(i,j) , i, j );
    }
  }
}


// Multiplication specialized methods:

/**
 * Matrix vector (pre)multiplication, specialized for Type_csc Data_mat (matrix) and Type_stdVector (STL vector) formats.
 * Calculates the product A*b = c using pointers.
 * @param matrix_in Type_csc *Matrix A.
 * @param vector_in Type_stdVector *Vector b.
 * @param vector_out Type_stdVector *Vector c = A*b.
 */
template <typename T>
    void mat_vec_mult
    ( const Type_csc<T>* matrix_in,
      const Type_stdVector<T>* vector_in,
      Type_stdVector<T>* vector_out
    )
{
  size_type a;
  size_type b;
  T c;
  size_type d;

 //rutina que multiplica la  matriz A (en formato Harwell-Boeing) con el vector X

  for(unsigned int i=0; i < (matrix_in->getRows()); ++i)
    vector_out->contents[i] = 0;

  //std::fill(vector_out->contents.begin(), vector_out->contents.end(), 0);

  for(unsigned int i=0; i < (matrix_in->getCols()); ++i) {
    a = matrix_in->ja[i];
    b = matrix_in->ja[i+1];

    for(unsigned int j=a; j<b; ++j) {
      c = matrix_in->aa[j-1];
      d = matrix_in->ia[j-1]-1;
      vector_out->contents[d] += c * vector_in->contents[i];
    }
  }
}


#ifdef HAVE_GMM
/**
 * Matrix vector (pre)multiplication, specialized for Type_gmm Data_mat (dense matrix) and Type_stdVector (STL vector) formats.
 * Calculates the product A*b = c using pointers.
 * @param matrix_in Type_gmm *Matrix A.
 * @param vector_in Type_stdVector *Vector b.
 * @param vector_out Type_stdVector *Vector c = A*b.
 */
template <typename T>
void mat_vec_mult( Type_gmm<T>* matrix_in,
                        Type_stdVector<T>* vector_in,
                        Type_stdVector<T>* vector_out) {
  gmm::mult( *(matrix_in->data_pointer() ), *(vector_in->data_pointer() ), *(vector_out->data_pointer() ) );
}

/**
 * Matrix vector (pre)multiplication, specialized for Type_gmm_sparse Data_mat (sparse matrix) and Type_stdVector (STL vector) formats.
 * Calculates the product A*b = c using pointers.
 * @param matrix_in Type_gmm_sparse *Matrix A.
 * @param vector_in Type_stdVector *Vector b.
 * @param vector_out Type_stdVector *Vector c = A*b.
 */
template <typename T>
void mat_vec_mult( Type_gmm_sparse<T>* matrix_in,
                        Type_stdVector<T>* vector_in,
                        Type_stdVector<T>* vector_out) {
  gmm::mult( *(matrix_in->data_pointer() ), *(vector_in->data_pointer() ), *(vector_out->data_pointer() )  );
}

/**
 * Matrix vector (pre)multiplication, specialized for Type_gmm Data_mat (dense matrix) and Type_gmmVector_sparse (STL vector) formats.
 * Calculates the product A*b = c using pointers.
 * @param matrix_in Type_gmm *Matrix A.
 * @param vector_in Type_gmmVector_sparse *Vector b.
 * @param vector_out Type_gmmVector_sparse *Vector c = A*b.
 */
template <typename T>
void mat_vec_mult( Type_gmm<T>* matrix_in,
                   Type_gmmVector_sparse<T>* vector_in,
                   Type_gmmVector_sparse<T>* vector_out) {
  gmm::mult( *(matrix_in->data_pointer() ), *(vector_in->data_pointer() ), *(vector_out->data_pointer() ) );
}

/**
 * Matrix vector (pre)multiplication, specialized for Type_gmm_sparse Data_mat (sparse matrix) and Type_gmmVector_sparse (STL vector) formats.
 * Calculates the product A*b = c using pointers.
 * @param matrix_in Type_gmm_sparse *Matrix A.
 * @param vector_in Type_gmmVector_sparse *Vector b.
 * @param vector_out Type_gmmVector_sparse *Vector c = A*b.
 */
template <typename T>
void mat_vec_mult( Type_gmm_sparse<T>* matrix_in,
                   Type_gmmVector_sparse<T>* vector_in,
                   Type_gmmVector_sparse<T>* vector_out) {
  gmm::mult( *(matrix_in->data_pointer() ), *(vector_in->data_pointer() ), *(vector_out->data_pointer() )  );
}
#endif

// Copy specialized:

/**
 * Matrix vector (pre)multiplication, specialized for Type_csc Data_mat (matrix) and Type_stdVector (STL vector) formats.
 * Calculates the product A*b = c using pointers.
 * @param dense_matrix_in Type_stdMatrix *Matrix B.
 * @param matrix_in Type_csc *Matrix A.
 */
template <typename T>
    void copy( const Type_stdmatrix<T>* dense_matrix_in,
               Type_csc<T>* matrix_in )
{
  size_type i,j, col_counter=1;
  size_type rows = matrix_in->getRows();
  size_type cols = matrix_in->getCols();
  T value;

  matrix_in->aa.clear();
  matrix_in->ia.clear();
  matrix_in->ja.clear();

  for (j=0; j<cols; ++j){
    matrix_in->ja.push_back(col_counter);
    for (i=0; i<rows; ++i){
      value = dense_matrix_in->contents[i][j];
      if ( value != T(0) ){
        matrix_in->ia.push_back(i+1);
        matrix_in->aa.push_back(value);
        ++col_counter;
      }
    }
  }
  matrix_in->ja.push_back(col_counter);

  matrix_in->Nrow = rows;
  matrix_in->Ncol = cols;
  matrix_in->Nnze = matrix_in->aa.size();
}                           


} //namespace lmx

#endif
