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

      \author Roberto Ortega Aguilera.

     */
//////////////////////////////////////////// Doxygen file documentation (end)


namespace lmx{

//Metodo de los gradientes conjugados

/**
  *
  * \class Cg
  * \brief Template class Cg.
  * This class implements the conjugate gradient iterative linear solver,
  * with and without preconditioner. It uses LMX matrix and vector facilities
  * so any kind of linked object can be used.
  *
  * @author Roberto Ortega.
*/
template <typename T> class Cg{

typedef double numType;

private:
  //+Matriz:
  Matrix<T> A;
  //+Vector b;
  Vector<T> b;
  //+Vector solucion inicial  x={0};
  Vector<T> x;
  //+Vector residuo inicial  r1={b};
  Vector<T> r;
  Vector<T> d;
  Vector<T> mp;
  Vector<T> s;
  size_type nrow;
  //+Tolerancia inicial para el residuo
  //tole=sqrt(h1/h0)
//   numType tole;
  T tole;
  T resi;
  //Definicion de error aceptable;
  T epsi;
//   numType epsi;
  //Número máximo de iteraciones:
  size_type kmax;

public:
  // Constructor:
  Cg(Matrix<T>*, Vector<T>*);

  ~Cg();

//+rutina para el precondicionador (diagonal)
  void precond();

//+rutina para el metodo de los gradientes conjugados
  Vector<T> solve( int );
};

} //namespace lmx

/////////////////////////////// Implementation of the methods defined previously


namespace lmx{

/**
 * Standard constructor.
 * @param A_in LHS Matrix
 * @param b_in RHS Vector
 */
  template <typename T> Cg<T>::Cg(Matrix<T>* A_in, Vector<T>* b_in) : A(*A_in), b(*b_in), r(*b_in), d(*b_in), tole(1), epsi(1E-6)
{
  nrow = A_in->rows();
  kmax = nrow+20;
  x.resize( nrow );
  s.resize( nrow );
  mp.resize( nrow );
  for (size_type i=0; i<nrow; ++i) mp(i) = 1;
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
  
  double temp;

  for (size_type i=0; i < nrow; ++i){
    temp = A(i,i);
    mp(i) = 1. / temp;
//     if (d(i) != 0){
      d(i) /= temp;
//     }
  }
}


//+rutina para el metodo de los gradientes conjugados
template <typename T>
/**
 * 
 * @param info 
 * @return 
 */
Vector<T> Cg<T>::solve(int info){

  //valor h y h0
  T hnew;
  hnew = r*d;

  T h0 = hnew;
  Vector<T> q(nrow);
  Vector<T> temp_vec(nrow);
  T alfa;
//   numType alfa;
  size_type k=1;

  T hold;
  T temp;
  T beta;
//   numType beta;

  while(k < kmax && tole > epsi ){
    if (info > 0) cout<<"iteracion :"<<k<<"\t";

//     q = A * d;
    q.mult( A , d );

    //for(size_type i=1;i<Nrow+1;++i)cout<<"w "<<w[i]<<endl;
    temp = d * q;
    alfa = hnew / temp;
//     x = x + alfa * d;
    temp_vec.mult(alfa, d);
    x += temp_vec;
    //for(size_type i=1;i<Nrow+1;++i)cout<<"x "<<x[i]<<endl;cout<<endl;

    if(k%50 != 0){
//       r = r - alfa * q;
      temp_vec.mult(alfa, q);
      r -= temp_vec;
    }
    else{
//       r = b - A*x;
      temp_vec.mult(A,x);
      r.subs(b, temp_vec);
    }
//       s(i) = mp(i) * r(i);
    s.multElem(mp,r);

    hold = hnew;
    hnew = r*s;
    beta = hnew / hold;
//     d = s + beta*d;
    temp_vec.mult(beta,d);
    d.add(s, temp_vec);

    resi = sqrt(hnew);
    if (info > 0) cout << "res = " << resi << "\t";
    tole = sqrt(hnew / h0);
    if (info > 0) cout << "tole = " << tole << endl;
    ++k;
  }


  if (k==kmax)
    cout<<":::WARNING:::" << endl
     << ":::Convegence was not achieved after " << k << " iterations.:::" << endl
     << ":::Check if the matrix is symmetric.:::" << endl
     << ":::END WARNING:::" << endl ;

  else
    cout<<":::System solved:::"<<endl;

  return x;

}

} //namespace lmx


#endif
