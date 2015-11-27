/***************************************************************************
 *   Copyright (C) 2005 by Daniel Iglesias                                 *
 *   diglesiasib@mecanica.upm.es                                           *
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

#ifndef LINEAR_SYSTEM_H
#define LINEAR_SYSTEM_H

#include <iostream>

#include "lmx_mat_dense_matrix.h"
#include "lmx_linsolvers_cg.h"
#include "lmx_linsolvers_gauss.h"

#ifdef HAVE_LAPACK
#include "lmx_linsolvers_lapack.h"
#endif

#ifdef HAVE_SUPERLU
#include "lmx_linsolvers_superlu_interface.h"
#endif

//////////////////////////////////////////// Doxygen file documentation entry:
    /*!
      \file lmx_linsolvers_system.h

      \brief LinearSystem class implementation

      Implements a typical "A*x = b" system

      \author Daniel Iglesias Ib��ez

    */
//////////////////////////////////////////// Doxygen file documentation (end)


// typedef size_t size_type;

namespace lmx {

int setLinSolverType(int);
int getLinSolverType();


    /**
    \class LinearSystem 
    \brief Template class LinearSystem.
    Linear systems implementation: "A*x = b" .

    This class permits the creation of a linear system object. Each object has three parameters, corresponding to each of the matrices or vectors base of the problem. The basic methods solve the problem and add functionality to control the solution procedure. Not only one solver can be used as well as the number data type (class) may be differ between the input and the one used to solve the system.

    @author Daniel Iglesias Ib��ez.
    */
template <class T> class LinearSystem{
private:
  Matrix<T>* A;
  DenseMatrix<T>* dA;
  Vector<T>* x;
  Vector<T>* b;
  bool A_new, x_new, b_new;
  int info; /**< sets level of information in std output **/
#ifdef HAVE_SUPERLU
  Superlu<T>* S;
#endif

public:

  /** Empty constructor. */
  LinearSystem() : A(0), dA(0),x(0), b(0), A_new(0), x_new(0), b_new(0)
  { 
    #ifdef HAVE_SUPERLU
        S = 0;
    #endif
  }

  /**
   * Standard constructor with two parameters.
   * @param A_in LHS Matrix.
   * @param b_in RHS Vector.
   */
  LinearSystem(Matrix<T>& A_in, Vector<T>& b_in) : A(&A_in), dA(0), b(&b_in), x(0), A_new(0), x_new(1), b_new(0)
  {
    x = new Vector<T>( b_in.size() );
//     x->resize( b_in.size() );
    *x = b_in;

#ifdef HAVE_SUPERLU
        S = 0;
#endif
  }


  /**
   * Standard constructor with two parameters (DenseMatrix).
   * @param dA_in LHS DenseMatrix.
   * @param b_in RHS Vector.
   */
  LinearSystem(DenseMatrix<T>& dA_in, Vector<T>& b_in) : A(0), dA(&dA_in), b(&b_in), x(0), A_new(0), x_new(1), b_new(0)
  {
    x = new Vector<T>;
    x->resize( b_in.size() );
    *x = b_in;

#ifdef HAVE_SUPERLU
    S = 0;
#endif
  }


  template <class C>
  /**
   * Standard constructor with two parameters with different data type as that of LinearSystem.
   * @param A_in LHS Matrix.
   * @param b_in RHS Vector.
   */
      LinearSystem(Matrix<C>& A_in, Vector<C>& b_in) : dA(0), A_new(1), x_new(1), b_new(1)
  {
    A = new Matrix<T>;
    *A = A_in;

    x = new Vector<T>;
    x->resize( rows(b_in) );
    *x = b_in;

    b = new Vector<T>;
    b->resize( rows(b_in) );
    *b = b_in;

#ifdef HAVE_SUPERLU
    S = 0;
#endif
  }


  /**
   * Standard constructor with three parameters.
   * @param A_in LHS Matrix.
   * @param x_in Solution vector.
   * @param b_in RHS Vector.
   * @return
   */
  LinearSystem(Matrix<T>& A_in, Vector<T>& x_in, Vector<T>& b_in) : A(&A_in), dA(0), x(&x_in), b(&b_in), A_new(0), x_new(0), b_new(0)
  {
#ifdef HAVE_SUPERLU
    S = 0;
#endif
  }

  /**
   * Standard constructor with three parameters (DenseMatrix).
   * @param dA_in LHS DenseMatrix.
   * @param x_in Solution vector.
   * @param b_in RHS Vector.
   */
  LinearSystem(DenseMatrix<T>& dA_in, Vector<T>& x_in, Vector<T>& b_in) : A(0), dA(&dA_in), x(&x_in), b(&b_in), A_new(0), x_new(0), b_new(0)
  {
#ifdef HAVE_SUPERLU
    S = 0;
#endif
  }

  template <class C>
  /**
   * Standard constructor with three parameters, two of them with different data types.
   * @param A_in LHS Matrix.
   * @param x_in Solution vector.
   * @param b_in RHS Vector.
   * @return
   */
  LinearSystem(Matrix<C>& A_in, Vector<T>& x_in, Vector<C>& b_in) : dA(0), x(&x_in), A_new(1), x_new(0), b_new(1)
  {
    A = new Matrix<T>;
    *A = A_in;

    b = new Vector<T>;
    b->resize( rows(b_in) );
    *b = b_in;

#ifdef HAVE_SUPERLU
    S = 0;
#endif
  }

  /**
   * Destructor.
   */
  ~LinearSystem()
  {
     if (A_new){
       delete A;
       A = 0;
     }

     if (x_new){
       delete x;
       x = 0;
     }

     if (b_new){
       delete b;
       b = 0;
     }

#ifdef HAVE_SUPERLU
     delete S;
     S = 0;
#endif
  }

  /** Set information level. 
    */
  void setInfo(int level)
  { info = level; }

  Vector<T>& solveYourself(bool);

  /**
   * Solution access.
   * @return A vector with solution values.
   */
  Vector<T>& getSolution()
  { return *(this->x); }


};

}; // namespace lmx

/////////////////////////////// Implementation of the methods defined previously

namespace lmx {


  /**
   * \brief Solve function.
   * 
   * Depending on Matrix and lin_solver types selected, the following combinations are possible:
   *
   * <table> <tr> <td>getLinSolverType()</td>    <td>getMatrixType()</td>    <td>Solver used:</td> </tr>
   *  <tr> <td> 0 </td>    <td> 0 </td>    <td> Gauss</td> </tr>
   *  <tr> <td> 0 </td>    <td> 1 </td>    <td> SuperLU </td> </tr>
   *  <tr> <td> 0 </td>    <td> 2 </td>    <td> gmm::lu_solve (SuperLU in the future)</td>    </tr>
   *  <tr> <td> 0 </td>    <td> 3 </td>    <td> gmm::lu_solve (SuperLU in the future)</td> </tr>
   *
   *  <tr> <td> 1 </td>    <td> 0 </td>    <td> Gauss</td> </tr>
   *  <tr> <td> 1 </td>    <td> 1 </td>    <td> SuperLU </td> </tr>
   *  <tr> <td> 1 </td>    <td> 2 </td>    <td> gmm::lu_solve (SuperLU in the future)</td>    </tr>
   *  <tr> <td> 1 </td>    <td> 3 </td>    <td> gmm::lu_solve (SuperLU in the future)</td> </tr>
   *
   *  <tr> <td> 2 </td>    <td> 0 </td>    <td> lmx::Cg </td> </tr>
   *  <tr> <td> 2 </td>    <td> 1 </td>    <td> lmx::Cg </td> </tr>
   *  <tr> <td> 2 </td>    <td> 2 </td>    <td> lmx::Cg (and gmm::cg possible if uncommented)</td> </tr>
   *  <tr> <td> 2 </td>    <td> 3 </td>    <td> lmx::Cg (and gmm::cg possible if uncommented)</td> </tr>
   *
   *  <tr> <td> 3 </td>    <td> 0 </td>    <td> - </td> </tr>
   *  <tr> <td> 3 </td>    <td> 1 </td>    <td> - </td> </tr>
   *  <tr> <td> 3 </td>    <td> 2 </td>    <td> - </td> </tr>
   *  <tr> <td> 3 </td>    <td> 3 </td>    <td> gmm::gmres </td> </tr>
   *  </table>
   *
   * IMPORTANT: gmm::lu_solve is not implemented when the vector_type
   * is set to code 1 (Type_gmmVector_sparse). To Be Done.
   *
   * When a solver is not available, an error will be thrown.
   *
   * @param recalc For SuperLU switches between refactoring (FALSE) or use old factoring (TRUE).
   * @return reference of solution Vector.
   */
  template <class T>
      Vector<T>& LinearSystem<T>::solveYourself(bool recalc = 0)
  {
    // Routine for DenseMatrix:
    if (A==0 && dA!=0){
      if ( ( dA->cols() != x->size() ) || (dA->rows() != b->size() ) ){
        std::stringstream message;
        message << "Error in LinearSystem \"A*x = b\": Dimensions mismatch. \n"
            << "DenseMatrix \"A\" dimension: (" << dA->rows() << "," << dA->cols() << ")" << endl
            << "LHS Vector \"x\" dimension: (" << x->size() << ")" << endl
            << "RHS Vector \"b\" dimension: (" << b->size() << ")" << endl;
        LMX_THROW(dimension_error, message.str() );
      }
      // Using built-in gauss elimination procedure:
#ifdef HAVE_LAPACK
      Gesv<T> solver( dA, x, b );
      solver.solve();
#else
      Gauss<T> solver( dA, b );
      *x = solver.solve();
#endif
      return *x;
    }
    
    else{
      if ( ( A->cols() != x->size() ) || (A->rows() != b->size() ) ){
        std::stringstream message;
        message << "Error in LinearSystem \"A*x = b\": Dimensions mismatch. \n"
            << "Matrix \"A\" dimension: (" << A->rows() << "," << A->cols() << ")" << endl
            << "LHS Vector \"x\" dimension: (" << x->size() << ")" << endl
            << "RHS Vector \"b\" dimension: (" << b->size() << ")" << endl;
        LMX_THROW(dimension_error, message.str() );
      }

      switch (getLinSolverType()) {
        case 0 : // solver_type == 0 -> directos para sistemas simetrico
          switch (getMatrixType()) {
            case 0 :
            {  // Using built-in gauss elimination procedure:
#ifdef HAVE_LAPACK
              Gesv<T> solver( A, x, b );
              solver.solve();
#else
              Gauss<T> solver( A, b );
              *x = solver.solve();
#endif
              return *x;
            }
            break;

            case 1 :
#ifdef HAVE_SUPERLU
              if (recalc == FALSE){
                cout << "S: " << S << endl;
                if (S == 0)
                  S = new Superlu<T>( static_cast<Type_csc<T>*>(A->type_matrix)-> Nrow,
                                      static_cast<Type_csc<T>*>(A->type_matrix)-> Ncol,
                                      static_cast<Type_csc<T>*>(A->type_matrix)-> Nnze,
                                      static_cast<Type_csc<T>*>(A->type_matrix)-> ia,
                                      static_cast<Type_csc<T>*>(A->type_matrix)-> ja,
                                      static_cast<Type_csc<T>*>(A->type_matrix)-> aa,
                                      *(static_cast<Type_stdVector<T>*>(x->type_vector)->data_pointer() ),
                                      *(static_cast<Type_stdVector<T>*>(b->type_vector)->data_pointer() ) );
                S->init();
                S->calc(info);
                S->get_solution( *(static_cast<Type_stdVector<T>*>(x->type_vector)->data_pointer() ) );
              }
              else{

                S->recalc(info);
                S->get_solution( *(static_cast<Type_stdVector<T>*>(x->type_vector)->data_pointer() ) );
              }
#else
              {
                  std::stringstream message;
                  message << "SuperLU not defined.\nYou must set \"#define HAVE_SUPERLU\" in your file in order to use this library." << endl;
                  LMX_THROW(failure_error, message.str() );
              }
#endif
              return *x;
            break;

            case 2 :
#ifdef HAVE_GMM
              switch (getVectorType()) {
                case 0:
                // Posible sustituci�n con SuperLU con flag de simetrico
                  gmm::lu_solve(*(static_cast<Type_gmm<T>*>(A->type_matrix)->data_pointer()), *(static_cast<Type_stdVector<T>*>(x->type_vector)->data_pointer()), *(static_cast<Type_stdVector<T>*>(b->type_vector)->data_pointer()));
    
                  break;

                case 1:
                // Posible sustituci�n con SuperLU con flag de simetrico
                  // This doesn't work, will try to fix it in the future.

  //                 gmm::lu_solve(*(static_cast<Type_gmm<T>*>(A->type_matrix)->data_pointer()), *(static_cast<Type_gmmVector_sparse<T>*>(x->type_vector)->data_pointer()), *(static_cast<Type_gmmVector_sparse<T>*>(b->type_vector)->data_pointer()));

                  std::stringstream message;
                  message << "Solver not implemented.\nLinear solver type = " << getLinSolverType() << ", Matrix type = " << getMatrixType() << ", Vector type = " << getVectorType() << "." << endl;
                  LMX_THROW(to_be_done_error, message.str() );

                  break;
              }
  #else
              {
                  std::stringstream message;
                  message << "gmm++ not defined.\nYou must set \"#define HAVE_GMM\" in your file in order to use this library." << endl;
                  LMX_THROW(failure_error, message.str() );
              }
#endif
              return *x;
              break;

            case 3 :
#ifdef HAVE_GMM
              switch (getVectorType()) {
                case 0:
              // Posible sustituci�n con SuperLU con flag de simetrico
                gmm::lu_solve(*(static_cast<Type_gmm_sparse<T>*>(A->type_matrix)->data_pointer()), *(static_cast<Type_stdVector<T>*>(x->type_vector)->data_pointer()), *(static_cast<Type_stdVector<T>*>(b->type_vector)->data_pointer()));

                case 1:
                // Posible sustituci�n con SuperLU con flag de simetrico
                  // This doesn't work, will try to fix it in the future.

                  //                                 gmm::lu_solve(*(static_cast<Type_gmm_sparse<T>*>(A->type_matrix)->data_pointer()), *(static_cast<Type_gmmVector_sparse<T>*>(x->type_vector)->data_pointer()), *(static_cast<Type_gmmVector_sparse<T>*>(b->type_vector)->data_pointer()));

                  break;
              }
#else
              {
                  std::stringstream message;
                  message << "gmm++ not defined.\nYou must set \"#define HAVE_GMM\" in your file in order to use this library." << endl;
                  LMX_THROW(failure_error, message.str() );
              }
#endif
              return *x;
              break;

          }
          break;


          case 1 : // solver_type == 1 -> directos para sistemas simetricos y no simetricos

            switch (getMatrixType()) {
            case 0 :
            {  // Using built-in gauss elimination procedure:
              Gauss<T> solver( A, b );
              *x = solver.solve();
              return *x;
            }
            break;

            case 1 :
#ifdef HAVE_SUPERLU
              if (recalc == FALSE){
                cout << "S: " << S << endl;
                if (S == 0)
                  S = new Superlu<T>( static_cast<Type_csc<T>*>(A->type_matrix)-> Nrow,
                                      static_cast<Type_csc<T>*>(A->type_matrix)-> Ncol,
                                      static_cast<Type_csc<T>*>(A->type_matrix)-> Nnze,
                                      static_cast<Type_csc<T>*>(A->type_matrix)-> ia,
                                      static_cast<Type_csc<T>*>(A->type_matrix)-> ja,
                                      static_cast<Type_csc<T>*>(A->type_matrix)-> aa,
                                      *(static_cast<Type_stdVector<T>*>(x->type_vector)->data_pointer() ),
                                      *(static_cast<Type_stdVector<T>*>(b->type_vector)->data_pointer() ) );
                S->init();
                S->calc(info);
                S->get_solution( *(static_cast<Type_stdVector<T>*>(x->type_vector)->data_pointer() ) );
              }
              else{
                S->recalc(info);
                S->get_solution( *(static_cast<Type_stdVector<T>*>(x->type_vector)->data_pointer() ) );
              }
#else
              {
                  std::stringstream message;
                  message << "SuperLU not defined.\nYou must set \"#define HAVE_SUPERLU\" in your file in order to use this library." << endl;
                  LMX_THROW(failure_error, message.str() );
              }
#endif
              return *x;
            break;

            case 2 :
#ifdef HAVE_GMM
              switch (getVectorType()) {
                case 0:
                // Posible sustituci�n con SuperLU
                  gmm::lu_solve(*(static_cast<Type_gmm<T>*>(A->type_matrix)->data_pointer()), *(static_cast<Type_stdVector<T>*>(x->type_vector)->data_pointer()), *(static_cast<Type_stdVector<T>*>(b->type_vector)->data_pointer()));
    
                  break;

                case 1:
                // Posible sustituci�n con SuperLU con flag de simetrico
                  // This doesn't work, will try to fix it in the future.

                  //                 gmm::lu_solve(*(static_cast<Type_gmm<T>*>(A->type_matrix)->data_pointer()), *(static_cast<Type_gmmVector_sparse<T>*>(x->type_vector)->data_pointer()), *(static_cast<Type_gmmVector_sparse<T>*>(b->type_vector)->data_pointer()));

                  break;
              }
#else
              {
                  std::stringstream message;
                  message << "gmm++ not defined.\nYou must set \"#define HAVE_GMM\" in your file in order to use this library." << endl;
                  LMX_THROW(failure_error, message.str() );
              }
#endif
              return *x;
            break;

            case 3 :
#ifdef HAVE_GMM
              switch (getVectorType()) {
                case 0:
              // Posible sustituci�n con SuperLU
                  gmm::lu_solve(*(static_cast<Type_gmm_sparse<T>*>(A->type_matrix)->data_pointer()), *(static_cast<Type_stdVector<T>*>(x->type_vector)->data_pointer()), *(static_cast<Type_stdVector<T>*>(b->type_vector)->data_pointer()));

                case 1:
                // Posible sustituci�n con SuperLU con flag de simetrico
                  // This doesn't work, will try to fix it in the future.

                  //                 gmm::lu_solve(*(static_cast<Type_gmm_sparse<T>*>(A->type_matrix)->data_pointer()), *(static_cast<Type_gmmVector_sparse<T>*>(x->type_vector)->data_pointer()), *(static_cast<Type_gmmVector_sparse<T>*>(b->type_vector)->data_pointer()));

                  break;
              }
#else
              {
                  std::stringstream message;
                  message << "gmm++ not defined.\nYou must set \"#define HAVE_GMM\" in your file in order to use this library." << endl;
                  LMX_THROW(failure_error, message.str() );
              }
#endif
              return *x;
              break;
            }
          break;

        case 2 : // solver_type == 2 -> iterativos para sistemas simetricos
          switch (getMatrixType()) {
            case 0 :
            {
              Cg<T> cg_solver(A, b);
              cg_solver.precond();
              *x = cg_solver.solve(info);

              return *x;
            }
            break;

            case 1 :
            {
              Cg<T> cg_solver(A, b);
              cg_solver.precond();
              *x = cg_solver.solve(info);

              return *x;
            }
            break;

            case 2 :
            { 
              // CG de lmx
              Cg<T> cg_solver(A, b);
              cg_solver.precond();
              *x = cg_solver.solve(info);

              // CG de gmm SOLO PARA VECTORES STL
  // #ifdef HAVE_GMM
  //             gmm::diagonal_precond< gmm::row_matrix< gmm::rsvector<T> > > PR(*(static_cast<Type_gmm_sparse<T>*>( A->type_matrix)->data_pointer() ) );

  //             gmm::iteration iter(1E-6);
  //             iter.set_noisy(info);

  //             gmm::identity_matrix PS;   // Optional scalar product for cg

  //             gmm::cg(*(static_cast<Type_gmm_sparse<T>*>(A->type_matrix)->data_pointer()),
  //                     *(static_cast<Type_stdVector<T>*>(x->type_vector)->data_pointer()), *(static_cast<Type_stdVector<T>*>(b->type_vector)->data_pointer()), PS, PR, iter);
  // #endif
              return *x;
            }
  //////////////////////////////////////// Just another form of calling gmm's CG, BEGIN:
  //           { gmm::csc_matrix<T> GM;
  //             gmm::copy(*(static_cast<Type_gmm_sparse<T>*>(A->type_matrix)->data_pointer() ), GM);
  // 
  //             gmm::ilut_precond< gmm::csc_matrix<T> > PR(GM, 10, 1e-4);
  // 
  //             gmm::iteration iter(1E-8);
  //             iter.set_noisy(1);
  // 
  //             gmm::identity_matrix PS;   // Optional scalar product for cg
  // 
  //             gmm::cg(GM,  *(static_cast<Type_stdVector<T>*>(x->type_vector)->data_pointer()), *(static_cast<Type_stdVector<T>*>(b->type_vector)->data_pointer()), PS, PR, iter);
  // //             gmm::cg(*(static_cast<Type_gmm_sparse<T>*>(A->type_matrix)->data_pointer()),  *(static_cast<Type_stdVector<T>*>(x->type_vector)->data_pointer()), *(static_cast<Type_stdVector<T>*>(b->type_vector)->data_pointer()), PS, P, iter);
  // 
  //             return *x;
  //            }
  //////////////////////////////////////// Just another form of calling gmm's CG, END:

            break;

            case 3 :
            {
//               Cg<T> cg_solver(A, b);
//               cg_solver.precond();
//               *x = cg_solver.solve(info);

              // CG de gmm SOLO PARA VECTORES STL
  #ifdef HAVE_GMM
              gmm::diagonal_precond< gmm::row_matrix< gmm::rsvector<T> > > PR(*(static_cast<Type_gmm_sparse<T>*>( A->type_matrix)->data_pointer() ) );

              gmm::iteration iter(1E-6);
              iter.set_noisy(info);

              gmm::identity_matrix PS;   // Optional scalar product for cg

              gmm::cg(*(static_cast<Type_gmm_sparse<T>*>(A->type_matrix)->data_pointer()),
                      *(static_cast<Type_stdVector<T>*>(x->type_vector)->data_pointer()), *(static_cast<Type_stdVector<T>*>(b->type_vector)->data_pointer()), PS, PR, iter);
  #endif
            return *x;
          }
  //////////////////////////////////////// Just another form of calling gmm's CG, BEGIN:
  //           { gmm::csc_matrix<T> GM;
  //             gmm::copy(*(static_cast<Type_gmm_sparse<T>*>(A->type_matrix)->data_pointer() ), GM);
          //
  //             gmm::ilut_precond< gmm::csc_matrix<T> > PR(GM, 10, 1e-4);
          //
  //             gmm::iteration iter(1E-8);
  //             iter.set_noisy(1);
          //
  //             gmm::identity_matrix PS;   // Optional scalar product for cg
          //
  //             gmm::cg(GM,  *(static_cast<Type_stdVector<T>*>(x->type_vector)->data_pointer()), *(static_cast<Type_stdVector<T>*>(b->type_vector)->data_pointer()), PS, PR, iter);
  // //             gmm::cg(*(static_cast<Type_gmm_sparse<T>*>(A->type_matrix)->data_pointer()),  *(static_cast<Type_stdVector<T>*>(x->type_vector)->data_pointer()), *(static_cast<Type_stdVector<T>*>(b->type_vector)->data_pointer()), PS, P, iter);
          //
  //             return *x;
  //            }
  //////////////////////////////////////// Just another form of calling gmm's CG, END:
          break;

          }
          break;

        case 3 : // solver_type == 3 -> iterativos para sistemas simetrico y no simetrico
          switch (getMatrixType()) {
            case 0 :
              // Not available (error)

            break;

            case 1 :
              // TBI

              break;

            case 2 :
              // Not available (error)

            break;

            case 3 :
              // gmres
#ifdef HAVE_GMM
              switch (getVectorType()) {
                case 0:
                {
                  gmm::ilut_precond< gmm::row_matrix< gmm::rsvector<T> > > P(*(static_cast<Type_gmm_sparse<T>*>(A->type_matrix)->data_pointer()), 10, 1e-4);

                  gmm::iteration iter(1E-6);
                  iter.set_noisy(info);

                  gmm::gmres(*(static_cast<Type_gmm_sparse<T>*>(A->type_matrix)->data_pointer()),  *(static_cast<Type_stdVector<T>*>(x->type_vector)->data_pointer()), *(static_cast<Type_stdVector<T>*>(b->type_vector)->data_pointer()), P, 50, iter);
                }
                break;

                case 1:
                {
                  gmm::ilut_precond< gmm::row_matrix< gmm::rsvector<T> > > P(*(static_cast<Type_gmm_sparse<T>*>(A->type_matrix)->data_pointer()), 10, 1e-4);

                  gmm::iteration iter(1E-6);
                  iter.set_noisy(1);

                  gmm::gmres(*(static_cast<Type_gmm_sparse<T>*>(A->type_matrix)->data_pointer()),  *(static_cast<Type_gmmVector_sparse<T>*>(x->type_vector)->data_pointer()), *(static_cast<Type_gmmVector_sparse<T>*>(b->type_vector)->data_pointer()), P, 50, iter);
                }
                break;
    
            }
#else
              {
                  std::stringstream message;
                  message << "gmm++ not defined.\nYou must set \"#define HAVE_GMM\" in your file in order to use this library." << endl;
                  LMX_THROW(failure_error, message.str() );
              }
#endif
            return *x;

            break;

          }
          break;
      }
    }
    return *x;
  }

}; // namespace lmx


#endif
