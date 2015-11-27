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
#ifndef LMX_LINSOLVERS_LAPACK_H
#define LMX_LINSOLVERS_LAPACK_H

//////////////////////////////////////////// Doxygen file documentation entry:
    /*!
      \file lmx_linsolvers_lapack.h

      \brief Lapack interface for lmx::Matrix

      Implements the fortran interface for using the linear solvers from lapack.

      \author Daniel Iglesias

     */
//////////////////////////////////////////// Doxygen file documentation (end)

/**
 * Declaration of external function. Defined in any Lapack compatible library.
 */
extern "C" void   dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, 
                         double *b, int *ldb, int *info );

namespace lmx{

/**
 *
 * \class Gesv
 * \brief Template class for lapack ?gesv routine.
 *
 * @author Daniel Iglesias .
 */
template <typename T> class Gesv{
  private:
    Vector<T>* x;
    int nrhs;
    int n;
    int lda;
    int ldb;
    T* lb;
    int* ipiv;
    int info;
    int i,j;
    T** la;

  public:
    /**
     * Empty constructor.
     */
    Gesv(){}

    Gesv( Matrix<T>*, Vector<T>*, Vector<T>* );

    ~Gesv();

    void solve();

};

template <typename T>
    /**
 * Standard constructor.
 * @param a_in Pointer to Matrix.
 * @param x_in Pointer to solution Vector.
 * @param b_in Pointer to rhs Vector.
     */
Gesv<T>::Gesv( Matrix<T>* a_in, Vector<T>* x_in, Vector<T>* b_in ) :
    x( x_in ),
    nrhs(1),
    n( a_in->rows() ),
    lda( n ),
    ldb( n ),
    lb( new T[n] ),
    ipiv( new int[n] )
{
  if( a_in->rows() != a_in->cols() ){
    std::stringstream message;
    message << "Trying to build a Gesv object with a non-squared matrix.\nSize of matrix(" << a_in->rows() << ", " << a_in->cols() << ")." << endl;
    LMX_THROW(dimension_error, message.str() );
  }
  la = new double*[n];
  double* lla = new double[n*n];
  for (i=0; i<n; ++i){
    *(la+i) = lla;
    lla += n;
  }
  for(i = 0; i < n; ++i){
    for(j = 0; j < n; ++j)
      la[j][i] = a_in->readElement(i,j);
    lb[i] = b_in->readElement(i);
  }
}

/**
 * Destructor
 */
template <typename T>
    Gesv<T>::~Gesv( )
{
  if (lb){
    delete [] lb;
    lb = 0;
  }
  if (la){
    delete [] *la;
    delete [] la;
    la = 0;
  }
} 

/**
 * Solve system
 * @return Reference to solution vector.
 */
template <typename T>
    void Gesv<T>::solve()
{
  std::stringstream message;
  message << "ERROR: Solver not implemented for this data type."
          << endl;
  LMX_THROW(to_be_done_error, message.str() );

}

/**
 * Solve system, specialized for double data type.
 * @return Reference to solution vector.
 */
template <>
    inline void Gesv<double>::solve()
{
  dgesv_(&n, &nrhs, la[0], &lda, ipiv, lb, &ldb, &info);
  for(i = 0; i < n; ++i)
    x->writeElement( lb[i], i );
}

}

#endif
