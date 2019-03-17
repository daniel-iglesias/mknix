/***************************************************************************
 *   Copyright (C) 2005 by Roberto Ortega                                  *
 *   rortega@mecanica.upm.es                                               *
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
#ifndef CG_SOLVER_H
#define CG_SOLVER_H


#include "lmx_mat_vector.h"


//////////////////////////////////////////// Doxygen file documentation entry:
    /*!
      \file lmx_linsolvers_cg.h

      \brief Conjugate gradient class implementation

      Implements CG linear solver method.

      \author Roberto Ortega Aguilera, Juan C. García Orden.

     */
//////////////////////////////////////////// Doxygen file documentation (end)


namespace lmx{

/**
  *
  * \class Cg
  * \brief Template class Cg.
  * This class implements the conjugate gradient iterative linear solver,
  * with and without preconditioner. It uses LMX matrix and vector facilities
  * so any kind of linked object can be used.
  *
  * @author Roberto Ortega Aguilera, Juan C. García Orden.
*/
template <typename T> class Cg{

typedef double numType;

private:
	Matrix<T> A;	//+Matrix A;
	Vector<T> b;	//+Vector b;
	Vector<T> x;	//+Unknown;
	T tol;			//+Relative residual tolerance;
	size_type nrow;	//+Matrix size;
	size_type kmax;	//+Max iterations;
	Vector<T> mp;	//+Preconditioner;

public:
	Cg(Matrix<T>*, Vector<T>*);	//+Constructor;
	~Cg();						//+Destructor;
	Vector<T> solve( int );		//+Solver;
	void precond();
};

} //namespace lmx

/////////////////////////////// Implementation of the methods defined previously


namespace lmx{

/**
 * Standard constructor.
 * @param A_in LHS Matrix
 * @param b_in RHS Vector
 */
template <typename T> Cg<T>::Cg(Matrix<T>* A_in, Vector<T>* b_in)
:	A(*A_in),
	b(*b_in),
	x(),
	tol(1.0e-6),
	nrow(0), 
	kmax(0),
	mp(*b_in)
{
	nrow = A_in->rows();
	kmax = 3*nrow;
	x.resize( nrow );
}


/**
 * Destructor
 */
template <typename T> Cg<T>::~Cg(){

}


template <typename T>
/**
 * Preconditioner function.
 */
void Cg<T>::precond(){
  
	// Diagonal preconditioner
	for (size_type i=0; i < nrow; ++i) mp(i) = 1. / A.readElement(i,i);
}


//+rutina para el metodo de los gradientes conjugados
template <typename T>
/**
 * Solver function
 * @param info output iteration information
 * @return 
 */
Vector<T> Cg<T>::solve(int info){

///////// If vector b=0 -> x=0 and end
  int i;
  for( i=0; i<b.size(); ++i ){
    if ( b.readElement(i) != T(0) ) break;
  }
  if ( i == b.size() ){
    x=b;
    cout<<":::System solved:::"<<endl;
    return x;
  }

/////////

	Vector<T> r(nrow);	//+Residual;
	Vector<T> p(nrow);	//+Search direction;
	Vector<T> h(nrow);
	Vector<T> s(nrow);
	T delta;
	T deltaPrev;
	T alpha;
	T tole;
	T bdelta;
	
	Vector<T> vAux(nrow);
	T aux;

	size_type k=1;

// -----
		
	vAux.mult(A,x);
	r.subs(b,vAux);		// r = b - A*x
	
	h.multElements(mp,r);	// h = mp*r;
	delta = r*h;			// delta = r*h;
	p = h;
	
	vAux.multElements(mp,b);
	bdelta = b*vAux;		// bdelta = b*mp*b
	
	tole = tol*tol*bdelta;
	
	while (k < kmax && delta > tole )
	{
		if (info > 0) cout<<"iteration :"<<k<<"\t";
		
		s.mult(A,p);		// s = A*p
		
		aux = p*s;
		alpha = delta/aux;	// alpha = delta/(p*s)

		vAux.mult(alpha,p);
		x += vAux;			// x = x + alpha*p
		
	    if(k%50 != 0){
			vAux.mult(alpha,s);
			r -= vAux;		// r = r - alpha*s
		}
		else
		{
			vAux.mult(A,x);
			r.subs(b,vAux);	// r = b - A*x
		}
		
		h.multElements(mp,r);	// h = mp*r

		deltaPrev = delta;
		delta =	r*h;		// delta = r*h

		aux = delta/deltaPrev;
		vAux.mult(aux,p);
		p.add(h,vAux);		// p = h + (delta/deltaPrev)*p
		
		if (info > 0) cout << "residual = " << delta << "\t";
		if (info > 0) cout << "tolerance = " << tole << endl;

		++k;
	}

	if (k==kmax)
	{
		cout<<":::WARNING:::" << endl
		<< ":::Convegence was not achieved after " << k << " iterations.:::" << endl
		<< ":::END WARNING:::" << endl ;
	}
	else
    {
		cout<<":::System solved after " << k << " iterations:::"<<endl;
	}

return x;

}

} //namespace lmx


#endif
