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
#ifndef LU_SOLVER_H
#define LU_SOLVER_H

//////////////////////////////////////////// Doxygen file documentation entry:
    /*!
      \file lmx_linsolvers_lu.h

      \brief LU solver class implementation

      Implements LU decomposition with partial pivoting

      \author Juan Carlos García Orden

     */
//////////////////////////////////////////// Doxygen file documentation (end)


namespace lmx{

/**
 *
 * \class LU
 * \brief Template class LU for solving linear systems.
 *
 * @author Juan Carlos García Orden.
 */
template <typename T> class LU{
  private:
    DenseMatrix<T> mat;
    Vector<T> vec;
    DenseMatrix<T> mat_rhs;
    size_type dim;
    size_type nrhs;
    Vector<T> b;	// Auxiliary vector
	Vector<int> p;	// Permutation vector
	int pivotFlag;	// 0 means no pivoting

  public:
    /**
   * Empty constructor.
     */
    LU(){}

    LU( Matrix<T>*, Vector<T>* );
    LU( Matrix<T>*, Matrix<T>* );

    /**
     * Destructor
     */
    ~LU(){}

	void lu();					// LU decomposition with pivoting
 	void forSub( Vector<T> & );	// Forward substitution
	void backSub( Vector<T>& );	// Backward substitution
    Vector<T>& solve();
    Matrix<T>& solve_nrhs();
};

template <typename T>
/**
 * Standard constructor.
 * @param mat_in Pointer to Matrix.
 * @param vec_in Pointer to Vector.
 */
LU<T>::LU( Matrix<T>* mat_in, Vector<T>* vec_in )
  : dim( mat_in->rows() )
{
  if( mat_in->rows() != mat_in->cols() ){
    std::stringstream message;
    message << "Trying to build a LU object with a non-squared matrix.\nSize of matrix(" << mat_in->rows() << ", " << mat_in->cols() << ")." << endl;
    LMX_THROW(dimension_error, message.str() );
  }

  mat.resize( dim , dim );
  for (size_type i=0; i<dim; ++i){
    for (size_type j=0; j<dim; ++j){
      mat.writeElement( mat_in->readElement(i,j), i, j );
    }
  }
  vec.resize( *vec_in );
  vec = *vec_in;
  
  pivotFlag = 0;	// ********************* 0 means no pivoting
  if(pivotFlag!=0) p.resize(dim);

}

template <typename T>
/**
 * Standard constructor for nrhs Matrix.
 * @param mat_in Pointer to Matrix.
 * @param mat_rhs_in Pointer to Matrix.
 */
LU<T>::LU( Matrix<T>* mat_in, Matrix<T>* mat_rhs_in )
  : dim( mat_in->rows() )
  , nrhs( mat_rhs_in->cols() )
{
  if( mat_in->rows() != mat_in->cols() ){
    std::stringstream message;
    message << "Trying to build a LU object with a non-squared matrix.\nSize of matrix(" << mat_in->rows() << ", " << mat_in->cols() << ")." << endl;
    LMX_THROW(dimension_error, message.str() );
  }
  if( mat_in->rows() != mat_rhs_in->rows() ){
    std::stringstream message;
    message << "Trying to build a nrhs LU object with different matrix sizes.\nSize of matrices: A(" 
            << mat_in->rows() << ", " << mat_in->cols() << "); B(" 
            << mat_rhs_in->rows() << ", " << mat_rhs_in->cols() << ")." 
            << endl;
    LMX_THROW(dimension_error, message.str() );
  }

  mat.resize( dim , dim );
  for (size_type i=0; i<dim; ++i){
    for (size_type j=0; j<dim; ++j){
      mat.writeElement( mat_in->readElement(i,j), i, j );
    }
  }
  mat_rhs.resize( dim , nrhs );
  for (size_type i=0; i<dim; ++i){
    for (size_type j=0; j<nrhs; ++j){
      mat_rhs.writeElement( mat_rhs_in->readElement(i,j), i, j );
    }
  }
  b.resize( dim );
  
  pivotFlag = 0;	// ********************* 0 means no pivoting
  if(pivotFlag!=0) p.resize(dim);

}


template <typename T>
    /**
 * Perform LU-decomposition
 * Upon return the coefficients of L and U replace those
 * of the input dim-by-dim nonsingular matrix mat.
	*/
	void LU<T>::lu()
{
	int i;
	int j;
	int k;
	T mult;
	T aux1, aux2;
	int numR = 0;
	int numRaux;
	Vector<T> vAux;
	vAux.resize(dim);
	
	if (pivotFlag!=0)
	{
		for (j = 0; j<dim; j++)	p(j) = j;	// Initialize permutation vector
	}
	
	for (k = 0; k<dim-1; k++)
	{
		if (pivotFlag!=0)	// Pivoting
		{
//      aux1 = abs(mat.readElement(k,k));
      aux1 = mat.readElement(k,k);
      if(aux1 < T(0)) aux1 = -aux1; // templatized abs();
			numR = k;
			for (i = k+1; i<dim; i++)
			{
//				aux2 = abs(mat.readElement(i,k));
        aux2 = mat.readElement(i,k);
        if(aux2 < T(0)) aux2 = -aux2; // templatized abs();
				if ( aux2 > aux1 ) 
				{
					numR = i;
					aux1 = aux2;
				}
			}
			// Interchange rows
			for (j = 0; j<dim; j++)	vAux(j) = mat.readElement(numR,j);
			for (j = 0; j<dim; j++)	mat.writeElement( mat.readElement(k,j),numR,j);
			for (j = 0; j<dim; j++)	mat.writeElement( vAux(j),k,j);
			numRaux = p.readElement(numR);
			p.writeElement(k,numR);
			p.writeElement(numRaux,k);
		}

		// This is not necessary if pivoting is employed
		if ( mat.readElement(k,k) == 0 )
		{
			std::stringstream message;
			message << "Null pivot diagonal term!!. Term position: " << k << ". Squared Matrix dimension: " << dim << "." << endl;
			LMX_THROW(internal_error, message.str() );
		}
		//
		for (i = k+1; i<dim; i++)
		{
			if( mat.readElement(i,k) != 0 )
			{
				mult = mat.readElement(i,k) / mat.readElement(k,k);
				mat.writeElement(mult, i, k);
				for (j = k+1; j<dim; j++)
				{
					mat(i,j) -= mult*mat.readElement(k,j);
				}
			}
		}
	}
}

template <typename T>
    /**
 * Forward substitution
 * Given a unit lower triangular, nonsingular dim by dim matrix mat,
 * and dim-vector c
 * obtains vector y which solves mat*y = c
 * vector y replaces c
	*/
	void LU<T>::forSub( Vector<T>& c)
{	
	int i;
	int j;

	if (pivotFlag!=0)
	{
		Vector<T> vAux = c;
		// Permute b according to p
		for (i = 0; i<dim; i++) c(i) = vAux(p(i));
	}
	for (i = 1; i<dim; i++)
		for(j = 0; j < i; j++)	c(i) -= mat.readElement(i,j) * c.readElement(j);
}

template <typename T>
    /**
 * Backward substitution
 * Given an upper triangular, nonsingular dim by dim matrix mat and
 * an dim-vector c, obtains vector y which solves mat*y=c
 * vector y replaces c
	*/
	void LU<T>::backSub( Vector<T>& c)
{

	int j;
	int k;

	c(dim-1) /= mat.readElement(dim-1,dim-1);

	for (k=dim-2; k>=0; --k)
	{
		for (j = k+1; j<dim; ++j)	c(k) -= mat.readElement(k,j)*c.readElement(j) ;
		c(k) /= mat.readElement(k,k);
	}
}


template <typename T>
    /**
 * Solve system
 * @return Reference to solution vector.
     */
    Vector<T>& LU<T>::solve()
{

//	cout << "Matriz antes de permutar " << mat << "\n";
	lu();

//	cout << "permutation vector is " << p << "\n";
	
//	cout << "Matriz antes de permutar " << mat << "\n";
	forSub (vec);
	backSub(vec);
	
	return vec;
}

template <typename T>
/**
 * Solve n-rhs system
 * @return Reference to solution vector.
 */
Matrix<T>& LU<T>::solve_nrhs()
{
	int i;
	int j;
	
	lu();
	
	// Loop on columns of rhs matrix
	for(j=0; j<nrhs; ++j)
	{
		for (i=0; i<dim; ++i)	b.writeElement( mat_rhs(i,j), i );
		forSub (b);
		backSub(b);
		for (i=0; i<dim; ++i)	mat_rhs.writeElement( b(i), i, j );
	}
	return mat_rhs;
}

}

#endif
