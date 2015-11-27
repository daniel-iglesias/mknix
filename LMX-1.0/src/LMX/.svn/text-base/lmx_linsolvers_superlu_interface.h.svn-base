/***************************************************************************
 *   Copyright (C) 2005 by Daniel Iglesias                                 *
 *   dani@localhost.localdomain                                            *
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

#ifndef SUPERLU_INTERFACE_H
#define SUPERLU_INTERFACE_H

#include <vector>
// #include "dsp_defs.h"
#include "slu_ddefs.h"

//////////////////////////////////////////// Doxygen file documentation entry:
    /*!
      \file lmx_linsolvers_superlu_interface.h

      \brief Implementation of class Superlu for using this library linear solver.

      \author Daniel Iglesias Ibáñez

    */
//////////////////////////////////////////// Doxygen file documentation (end)


namespace lmx {

    /**
    \class Superlu
    \brief Template class LinearSystem.
    Linear systems implementation: "A*x = b" .

    This class permits the creation of a linear system object. Each object has three parameters, corresponding to each of the matrices or vectors base of the problem. The basic methods solve the problem and add functionality to control the solution procedure. Not only one solver can be used as well as the number data type (class) may be differ between the input and the one used to solve the system.

    \author Daniel Iglesias Ibáñez.
    */
template <typename T>
class Superlu{

private:
    char           equed[1];
    yes_no_t       equil;
    trans_t        trans;
    SuperMatrix    A, A1, L, U;
    SuperMatrix    B, B1, X;
    NCformat       *Astore;
    NCformat       *Ustore;
    SCformat       *Lstore;
    double         *a, *a1;
    int            *asub, *xa, *asub1, *xa1;
    int            *perm_r; /* row permutations from partial pivoting */
    int            *perm_c; /* column permutation vector */
    int            *etree;
    void           *work;
    int            info, lwork, nrhs, ldx;
    int            i, j, m, n, nnz;
    double         *rhsb, *rhsb1, *rhsx, *xact;
    double         *R, *C;
    double         *ferr, *berr;
    double         u, rpg, rcond;
    mem_usage_t    mem_usage;
    superlu_options_t options;
    SuperLUStat_t stat;


public:
  Superlu(size_type&, size_type&, size_type&, size_type *, size_type *, T *, T *, T * );

  Superlu(size_type&, size_type&, size_type&, size_type *, size_type *, T *, std::vector<T>&, std::vector<T>& );

  Superlu(size_type&, size_type&, size_type&, std::vector<size_type>&, std::vector<size_type>&, std::vector<T>&, std::vector<T>&, std::vector<T>& );

  Superlu(size_type&, size_type&, size_type&, std::vector<size_type>&, std::vector<size_type>&, std::vector<T>&, T *, T * );

  ~Superlu();

  void init();

//   void reinit(size_type&, size_type&, size_type&, std::vector<size_type>&, std::vector<size_type>&, std::vector<T>&, std::vector<T>&);

  void calc(int);

  void recalc(int);

  void get_solution(std::vector<T>&);

}; // class superlu

template<typename T>
/**
 * Constructor from c-array style vectors.
 * @param m_in Number of rows of matrix.
 * @param n_in Number of columns of matrix.
 * @param nnz_in Number of non-zeros in sparse matrix.
 * @param asub_in Array of integers.
 * @param xa_in Array of integers.
 * @param a_in Array of matrix's values.
 * @param x_in Array of solution values.
 * @param b_in Array of RHS vector's values.
 */
    Superlu<T>::Superlu(size_type & m_in,
                        size_type & n_in,
                        size_type & nnz_in,
                        size_type * asub_in,
                        size_type * xa_in,
                        T * a_in,
                        T * x_in,
                        T * b_in)
  : m(m_in), n(n_in), nnz(nnz_in)
{
    if ( !(a = doubleMalloc(nnz)) ) ABORT("Malloc fails for a1[].");
    if ( !(asub = intMalloc(nnz)) ) ABORT("Malloc fails for asub1[].");
    if ( !(xa = intMalloc(n+1)) ) ABORT("Malloc fails for xa1[].");
    for (i = 0; i < nnz; ++i) {
      a[i] = a_in[i];
      asub[i] = asub_in[i];
    }
    for (i = 0; i < n+1; ++i)
      xa[i] = xa_in[i];

// Guardar b_in y x_in
}


template<typename T>
/**
 * Constructor from c-array style sparse matrix and STL style vectors.
 * @param m_in Number of rows of matrix.
 * @param n_in Number of columns of matrix.
 * @param nnz_in Number of non-zeros in sparse matrix.
 * @param asub_in Array of integers.
 * @param xa_in Array of integers.
 * @param a_in Array of matrix's values.
 * @param x_in Vector of solution values.
 * @param b_in Vector of RHS vector's values.
 */
    Superlu<T>::Superlu(size_type & m_in,
                        size_type & n_in,
                        size_type & nnz_in,
                        size_type * asub_in,
                        size_type * xa_in,
                        T * a_in,
                        std::vector<T>& x_in,
                        std::vector<T>& b_in)
  : m(m_in), n(n_in), nnz(nnz_in)
{
  if ( !(a = doubleMalloc(nnz)) ) ABORT("Malloc fails for a1[].");
  if ( !(asub = intMalloc(nnz)) ) ABORT("Malloc fails for asub1[].");
  if ( !(xa = intMalloc(n+1)) ) ABORT("Malloc fails for xa1[].");
  for (i = 0; i < nnz; ++i) {
    a[i] = a_in[i];
    asub[i] = asub_in[i];
  }
  for (i = 0; i < n+1; ++i)
    xa[i] = xa_in[i];

// Guardar b_in y x_in
}


template<typename T>
/**
 * Constructor for STL style sparse matrix and vectors.
 * @param m_in Number of rows of matrix.
 * @param n_in Number of columns of matrix.
 * @param nnz_in Number of non-zeros in sparse matrix.
 * @param asub_in Vector of integers.
 * @param xa_in Vector of integers.
 * @param a_in Vector of matrix's values.
 * @param x_in Vector of solution values.
 * @param b_in Vector of RHS vector's values.
 */
    Superlu<T>::Superlu(size_type& m_in,
                        size_type& n_in,
                        size_type& nnz_in,
                        std::vector<size_type>& asub_in,
                        std::vector<size_type>& xa_in,
                        std::vector<T>& a_in,
                        std::vector<T>& x_in,
                        std::vector<T>& b_in)
    : m(m_in), n(n_in), nnz(nnz_in)
{
    if ( !(a = doubleMalloc(nnz)) ) ABORT("Malloc fails for a1[].");
    if ( !(asub = intMalloc(nnz)) ) ABORT("Malloc fails for asub1[].");
    if ( !(xa = intMalloc(n+1)) ) ABORT("Malloc fails for xa1[].");
    for (i = 0; i < nnz; ++i) {
      a[i] = a_in[i];
      asub[i] = asub_in[i] - 1;
    }
    for (i = 0; i < n+1; ++i)
      xa[i] = xa_in[i] - 1;

    if ( !(rhsx = doubleMalloc(m)) ) ABORT("Malloc fails for rhsx[].");
    if ( !(rhsb = doubleMalloc(m)) ) ABORT("Malloc fails for rhsb[].");
    for (i = 0; i < m; ++i) {
      rhsx[i] = x_in[i];
      rhsb[i] = b_in[i];

    }

}


template<typename T>
/**
 * Constructor for STL style sparse matrix and STL vectors.
 * @param m_in Number of rows of matrix.
 * @param n_in Number of columns of matrix.
 * @param nnz_in Number of non-zeros in sparse matrix.
 * @param asub_in Vector of integers.
 * @param xa_in Vector of integers.
 * @param a_in Vector of matrix's values.
 * @param x_in Array of solution values.
 * @param b_in Array of RHS vector's values.
 */
    Superlu<T>::Superlu(size_type& m_in,
                        size_type& n_in,
                        size_type& nnz_in,
                        std::vector<size_type>& asub_in,
                        std::vector<size_type>& xa_in,
                        std::vector<T>& a_in,
                        T* x_in,
                        T* b_in)
  : m(m_in), n(n_in), nnz(nnz_in)
{
  if ( !(a = doubleMalloc(nnz)) ) ABORT("Malloc fails for a1[].");
  if ( !(asub = intMalloc(nnz)) ) ABORT("Malloc fails for asub1[].");
  if ( !(xa = intMalloc(n+1)) ) ABORT("Malloc fails for xa1[].");
  for (i = 0; i < nnz; ++i) {
    a[i] = a_in[i];
    asub[i] = asub_in[i] - 1;
  }
  for (i = 0; i < n+1; ++i)
    xa[i] = xa_in[i] - 1;

  if ( !(rhsx = doubleMalloc(m)) ) ABORT("Malloc fails for rhsx[].");
  if ( !(rhsb = doubleMalloc(m)) ) ABORT("Malloc fails for rhsb[].");
  for (i = 0; i < m; ++i) {
    rhsx[i] = x_in[i];
    rhsb[i] = b_in[i];

  }

}


template<typename T>
/**
 * Destructor
 */
Superlu<T>::~Superlu()
  {
    StatFree(&stat);

    SUPERLU_FREE (xact);
    SUPERLU_FREE (etree);
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    SUPERLU_FREE (R);
    SUPERLU_FREE (C);
    SUPERLU_FREE (ferr);
    SUPERLU_FREE (berr);
    Destroy_CompCol_Matrix(&A);
    Destroy_Dense_Matrix(&B);
    Destroy_Dense_Matrix(&X);
    if ( lwork >= 0 ) { /* Deallocate storage associated with L and U. */
        Destroy_SuperNode_Matrix(&L);
        Destroy_CompCol_Matrix(&U);
    }

  }

template<typename T>
/**
 * Initialize function.
 */
void Superlu<T>::init()
{
// #if ( DEBUGlevel>=1 )
//     CHECK_MALLOC("Enter main()");
// #endif


    /* Defaults */
    lwork = 0;
    nrhs  = 1;
    equil = YES;  
    u     = 1.0;
    trans = NOTRANS;

//      Set the default input options:
    options.Fact = DOFACT;
    options.Equil = YES;
    options.ColPerm = COLAMD;
    options.DiagPivotThresh = 1.0;
    options.Trans = NOTRANS;
    options.IterRefine = NOREFINE;
    options.SymmetricMode = NO;
    options.PivotGrowth = NO;
    options.ConditionNumber = NO;
    options.PrintStat = YES;

    set_default_options(&options);

//     options.ColPerm = MMD_ATA;
//     options.Equil = equil;
//     options.DiagPivotThresh = u;
//     options.Trans = trans;

    if ( lwork > 0 ) {
      work = SUPERLU_MALLOC(lwork);
      if ( !work ) {
        ABORT("DLINSOLX: cannot allocate work[]");
      }
    }

    dCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
//     Astore = A.Store;
//     printf("Dimension %dx%d; # nonzeros %d\n", A.nrow, A.ncol, Astore->nnz);

    dCreate_Dense_Matrix(&B, m, nrhs, rhsb, m, SLU_DN, SLU_D, SLU_GE);
    dCreate_Dense_Matrix(&X, m, nrhs, rhsx, m, SLU_DN, SLU_D, SLU_GE);
    xact = doubleMalloc(n * nrhs);
    ldx = n;
    dGenXtrue(n, nrhs, xact, ldx);

    if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");
    if ( !(perm_r = intMalloc(m)) ) ABORT("Malloc fails for perm_r[].");
    if ( !(etree = intMalloc(n)) ) ABORT("Malloc fails for etree[].");
    if ( !(R = (double *) SUPERLU_MALLOC(A.nrow * sizeof(double))) ) 
        ABORT("SUPERLU_MALLOC fails for R[].");
    if ( !(C = (double *) SUPERLU_MALLOC(A.ncol * sizeof(double))) )
        ABORT("SUPERLU_MALLOC fails for C[].");
    if ( !(ferr = (double *) SUPERLU_MALLOC(nrhs * sizeof(double))) )
        ABORT("SUPERLU_MALLOC fails for ferr[].");
    if ( !(berr = (double *) SUPERLU_MALLOC(nrhs * sizeof(double))) ) 
        ABORT("SUPERLU_MALLOC fails for berr[].");

    /* Initialize the statistics variables. */
    StatInit(&stat);


}


template<typename T>
/**
 * Function to compute solution vector.
 * @param noisy Sets the level of output information.
 */
void Superlu<T>::calc(int noisy) {

    /* ------------------------------------------------------------
       WE SOLVE THE LINEAR SYSTEM FOR THE FIRST TIME: AX = B
       ------------------------------------------------------------*/

//     dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);


  dgssvx(&options, &A, perm_c, perm_r, etree, equed, R, C,
         &L, &U, work, lwork, &B, &X, &rpg, &rcond, ferr, berr,
         &mem_usage, &stat, &info);

//     printf("First system: dgssvx() returns info %d\n", info);

  if (noisy > 0){

    if ( info == 0 || info == n+1 ) {
  
    if ( options.PivotGrowth ) printf("Recip. pivot growth = %e\n", rpg);
    if ( options.ConditionNumber )
      printf("Recip. condition number = %e\n", rcond);
      Lstore = (SCformat *) L.Store;
      Ustore = (NCformat *) U.Store;
      printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
      printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
      printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz - n);
      printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions %d\n",
              mem_usage.for_lu/1e6, mem_usage.total_needed/1e6,
              mem_usage.expansions);
      if ( options.IterRefine ) {
        printf("Iterative Refinement:\n");
        printf("%8s%8s%16s%16s\n", "rhsb", "Steps", "FERR", "BERR");
        for (i = 0; i < nrhs; ++i)
          printf("%8d%8d%16e%16e\n", i+1, stat.RefineSteps, ferr[i], berr[i]);
      }
      fflush(stdout);
  
    }
    else if ( info > 0 && lwork == -1 ) {
      printf("** Estimated memory: %d bytes\n", info - n);
    }

    if ( options.PrintStat ) StatPrint(&stat);

  }

}

template<typename T>
/**
 * Function to compute solution vector using previous factorization.
 * @param noisy Sets the level of output information.
 */
void Superlu<T>::recalc(int noisy) {

    /* ------------------------------------------------------------
       NOW WE SOLVE ANOTHER LINEAR SYSTEM: A1*X = B1
       ONLY THE SPARSITY PATTERN OF A1 IS THE SAME AS THAT OF A.
       ------------------------------------------------------------*/
    for (i = 0; i < nnz; ++i) {
      a[i] += 1.5E-2;
    }

//     options.Fact = SamePattern;
    options.Fact = SamePattern_SameRowPerm;
    StatInit(&stat); /* Initialize the statistics variables. */

//     dCreate_CompCol_Matrix(&A1, m, n, nnz, a1, asub1, xa1,
//                            SLU_NC, SLU_D, SLU_GE);
//     dCreate_Dense_Matrix(&B1, m, nrhs, rhsb1, m, SLU_DN, SLU_D, SLU_GE);

    dgssvx(&options, &A, perm_c, perm_r, etree, equed, R, C,
           &L, &U, work, lwork, &B, &X, &rpg, &rcond, ferr, berr,
           &mem_usage, &stat, &info);

  if (noisy > 0){

    printf("\n Resolving system: dgssvx() returns info %d\n", info);

    if ( info == 0 || info == n+1 ) {

        /* This is how you could access the solution matrix. */
//         double *sol = (double*) ((DNformat*) X.Store)->nzval; 

    if ( options.PivotGrowth ) printf("Recip. pivot growth = %e\n", rpg);
    if ( options.ConditionNumber )
        printf("Recip. condition number = %e\n", rcond);
          Lstore = (SCformat *) L.Store;
          Ustore = (NCformat *) U.Store;
    printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
        printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
        printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz - n);
    printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions %d\n",
          mem_usage.for_lu/1e6, mem_usage.total_needed/1e6,
          mem_usage.expansions);
    if ( options.IterRefine ) {
              printf("Iterative Refinement:\n");
        printf("%8s%8s%16s%16s\n", "rhs", "Steps", "FERR", "BERR");
        for (i = 0; i < nrhs; ++i)
          printf("%8d%8d%16e%16e\n", i+1, stat.RefineSteps, ferr[i], berr[i]);
    }
    fflush(stdout);
      } else if ( info > 0 && lwork == -1 ) {
          printf("** Estimated memory: %d bytes\n", info - n);
      }
    if ( options.PrintStat ) StatPrint(&stat);
  }



// #if ( DEBUGlevel>=1 )
//     CHECK_MALLOC("Exit main()");
// #endif

}

template<typename T>
/**
 * Writes solution in vector.
 * @param x_out STL-style vector for writing the computed solution.
 */
void Superlu<T>::get_solution(std::vector<T>& x_out)
{ 
    rhsx = (T*) ((DNformat*) X.Store)->nzval;

    for (i = 0; i < m; ++i) x_out[i] = rhsx[i];

}


} // namespace lmx;

#endif
