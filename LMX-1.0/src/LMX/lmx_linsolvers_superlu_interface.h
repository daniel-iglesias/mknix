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

#ifndef SUPERLU_INTERFACE_H
#define SUPERLU_INTERFACE_H

#include <vector>
// #include "dsp_defs.h"
#include "slu_ddefs.h"

//////////////////////////////////////////// Doxygen file documentation entry:
    /*!
      \file lmx_linsolvers_superlu_interface.h

      \brief Implementation of class Superlu for using this library linear solver.

      \author Daniel Iglesias

    */
//////////////////////////////////////////// Doxygen file documentation (end)


namespace lmx {

    /**
    \class Superlu
    \brief Template class LinearSystem.
    Linear systems implementation: "A*x = b" .

    This class permits the creation of a linear system object. Each object has three parameters, corresponding to each of the matrices or vectors base of the problem. The basic methods solve the problem and add functionality to control the solution procedure. Not only one solver can be used as well as the number data type (class) may be differ between the input and the one used to solve the system.

    \author Daniel Iglesias .
    */
template <typename T>
class Superlu{

private:
    char           equed[1];
    yes_no_t       equil;
    trans_t        trans;
//begin JCGO
//    SuperMatrix    A, A1, L, U;
//    SuperMatrix    B, B1, X;
    SuperMatrix    A, L, U;
    SuperMatrix    B, X;
//end JCGO
    NCformat       *Astore;
    NCformat       *Ustore;
    SCformat       *Lstore;
//begin JCGO
//    double         *a, *a1;
//    int            *asub, *xa, *asub1, *xa1;
    double         *a;
    int            *asub, *xa;
//end JCGO
    int            *perm_r; /* row permutations from partial pivoting */
    int            *perm_c; /* column permutation vector */
    int            *etree;
    void           *work;
    int            info, lwork, nrhs, ldx;
    int            i, j, m, n, nnz;
//begin JCGO
//    double         *rhsb, *rhsb1, *rhsx, *xact;
    double         *rhsb, *rhsx, *xact;
//end JCGO
    double         *R, *C;
    double         *ferr, *berr;
    double         u, rpg, rcond;
    mem_usage_t    mem_usage;
    superlu_options_t options;
    SuperLUStat_t stat;


public:
  Superlu(size_type&, size_type&, size_type&, std::vector<size_type>&, std::vector<size_type>&, std::vector<T>& );

  Superlu(size_type&, size_type&, size_type&, std::vector<size_type>&, std::vector<size_type>&, std::vector<T>&, std::vector<T>&, std::vector<T>& );
  
  Superlu(size_type&, size_type&, size_type&, std::vector<size_type>&, std::vector<size_type>&, std::vector<T>&, T *, T * );
  
  Superlu(size_type&, size_type&, size_type&, size_type *, size_type *, T *, T *, T * );

  Superlu(size_type&, size_type&, size_type&, size_type *, size_type *, T *, std::vector<T>&, std::vector<T>& );

  ~Superlu();

  void init(){ initMatrix(); initVectors(); }
  void initMatrix();
  void initVectors();

//   void reinit(size_type&, size_type&, size_type&, std::vector<size_type>&, std::vector<size_type>&, std::vector<T>&, std::vector<T>&);

  void calc(int);

  void recalc(int);

  void get_solution(std::vector<T>&);

  void get_solution(Vector<T>&);
  
  void setVectors(Vector<T>& b_in);

//begin JCGO 01/04/09
	void setb(T *);
	void setb(std::vector<T>&);
	void setA(size_type& m_in,
              size_type& n_in,
              size_type& nnz_in,
              std::vector<size_type>& asub_in,
              std::vector<size_type>& xa_in,
              std::vector<T>& a_in );
    void recalc1(int);
    void recalc2(int);
    void factorize(void);
    void subsSolve(void);
    void get_solutionB(std::vector<T>&);
//end JCGO

}; // class superlu

  template<typename T>
  /**
   * Constructor for STL style sparse matrix.
   * @param m_in Number of rows of matrix.
   * @param n_in Number of columns of matrix.
   * @param nnz_in Number of non-zeros in sparse matrix.
   * @param asub_in Vector of integers.
   * @param xa_in Vector of integers.
   * @param a_in Vector of matrix's values.
   */
  Superlu<T>::Superlu(size_type& m_in,
                      size_type& n_in,
                      size_type& nnz_in,
                      std::vector<size_type>& asub_in,
                      std::vector<size_type>& xa_in,
                      std::vector<T>& a_in)
  : m(m_in), n(n_in), nnz(nnz_in)
  {
    cout << m_in << ", " << n_in << ", " << nnz_in << endl;
    if ( !(a = doubleMalloc(nnz)) ) ABORT("Malloc fails for a[].");
    if ( !(asub = intMalloc(nnz)) ) ABORT("Malloc fails for asub[].");
    if ( !(xa = intMalloc(n+1)) ) ABORT("Malloc fails for xa[].");
    for (i = 0; i < nnz; ++i) {
      a[i] = a_in[i];
      asub[i] = asub_in[i] - 1;
    }
    for (i = 0; i < n+1; ++i)
      xa[i] = xa_in[i] - 1;
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
    cout << m_in << ", " << n_in << ", " << nnz_in << endl;
    if ( !(a = doubleMalloc(nnz)) ) ABORT("Malloc fails for a[].");
    if ( !(asub = intMalloc(nnz)) ) ABORT("Malloc fails for asub[].");
    if ( !(xa = intMalloc(n+1)) ) ABORT("Malloc fails for xa[].");
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
  
  
  // DEPRECATED: TO DELETE WITH CAUTION...
  template<typename T>
  /**
   * Constructor for STL style sparse matrix and c-array vectors.
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
    if ( !(a = doubleMalloc(nnz)) ) ABORT("Malloc fails for a[].");
    if ( !(asub = intMalloc(nnz)) ) ABORT("Malloc fails for asub[].");
    if ( !(xa = intMalloc(n+1)) ) ABORT("Malloc fails for xa[].");
    for (i = 0; i < nnz; ++i) {
      a[i] = a_in[i];
      asub[i] = asub_in[i] - 1;
    }
    for (i = 0; i < n+1; ++i)  xa[i] = xa_in[i] - 1;
    
    if ( !(rhsx = doubleMalloc(m)) ) ABORT("Malloc fails for rhsx[].");
    if ( !(rhsb = doubleMalloc(m)) ) ABORT("Malloc fails for rhsb[].");
    for (i = 0; i < m; ++i) {
      rhsx[i] = x_in[i];
      rhsb[i] = b_in[i];
      
    }
  }
    
    // DEPRECATED: TO DELETE
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
    if ( !(a = doubleMalloc(nnz)) ) ABORT("Malloc fails for a[].");
    if ( !(asub = intMalloc(nnz)) ) ABORT("Malloc fails for asub[].");
    if ( !(xa = intMalloc(n+1)) ) ABORT("Malloc fails for xa[].");
    for (i = 0; i < nnz; ++i) {
      a[i] = a_in[i];
      asub[i] = asub_in[i];
    }
    for (i = 0; i < n+1; ++i)
      xa[i] = xa_in[i];

// TODO: Guardar b_in y x_in
}


  // DEPRECATED: TO DELETE
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
  if ( !(a = doubleMalloc(nnz)) ) ABORT("Malloc fails for a[].");
  if ( !(asub = intMalloc(nnz)) ) ABORT("Malloc fails for asub[].");
  if ( !(xa = intMalloc(n+1)) ) ABORT("Malloc fails for xa[].");
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
 * Initialize A SuperMatrix.
 */
void Superlu<T>::initMatrix()
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

  /* Initialize the statistics variables. */
  StatInit(&stat);
}
  
  template<typename T>
  /**
   * Initialize B,X SuperMatrices.
   */
  void Superlu<T>::initVectors()
  {
    dCreate_Dense_Matrix(&B, m, nrhs, rhsb, m, SLU_DN, SLU_D, SLU_GE);
    dCreate_Dense_Matrix(&X, m, nrhs, rhsx, m, SLU_DN, SLU_D, SLU_GE);

    if ( !(ferr = (double *) SUPERLU_MALLOC(nrhs * sizeof(double))) )
        ABORT("SUPERLU_MALLOC fails for ferr[].");
    if ( !(berr = (double *) SUPERLU_MALLOC(nrhs * sizeof(double))) ) 
        ABORT("SUPERLU_MALLOC fails for berr[].");

    /* Initialize the statistics variables. */
//    StatInit(&stat);
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

//begin JCGO
template<typename T>
void Superlu<T>::recalc1(int noisy)
{
	options.Fact = SamePattern_SameRowPerm;
//	options.Fact = SamePattern;
	recalc(noisy);
}

template<typename T>
void Superlu<T>::recalc2(int noisy)
{
	options.Fact = FACTORED;
	recalc(noisy);
}

template<typename T>
void Superlu<T>::factorize()
{
	SuperMatrix AC;
	set_default_options(&options);
	StatInit(&stat);
	double drop_tol = 0.0;
	int relax, panel_size;
	
/*
* Get column permutation vector perm_c[], according to permc_spec: 
* permc_spec = 0: natural ordering 
* permc_spec = 1: minimum degree on structure of A’*A 
* permc_spec = 2: minimum degree on structure of A’+A 
* permc_spec = 3: approximate minimum degree for unsymmetric matrices 
*/

	int permc_spec = 3; 
	get_perm_c(permc_spec, &A, perm_c);


	sp_preorder(&options, &A, perm_c, etree, &AC); 
	panel_size = sp_ienv(1); 
	relax = sp_ienv(2); 

	dgstrf(&options, &AC, drop_tol, relax, panel_size, etree, NULL, 0, perm_c, perm_r, &L, &U, &stat,&info);

}

template<typename T>
void Superlu<T>::subsSolve()
{
	StatInit(&stat);
	/*
	typedef struct { 
		SuperMatrix *L; 
		SuperMatrix *U; 
		int *perm_c; 
		int *perm_r;
	} factors_t; 	
	
	factors_t *LUfactors;
	
	LUfactors = (factors_t*) factors[0]; 
	L = LUfactors->L; 
	U = LUfactors->U; 
	perm_c = LUfactors->perm_c; 
	perm_r = LUfactors->perm_r; 
	*/
	
	dgstrs (trans, &L, &U, perm_c, perm_r, &B, &stat, &info);
	
//	rhsx = rhsb;
//	dCreate_Dense_Matrix(&X, m, nrhs, rhsx, m, SLU_DN, SLU_D, SLU_GE);
}

//end JCGO

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

//     options.Fact = SamePattern;
//begin JCGO
//    options.Fact = SamePattern_SameRowPerm;
//end JCGO
//    B.ncol = nrhs;  /* Set the number of right-hand side */

    StatInit(&stat); /* Initialize the statistics variables. */

//     dCreate_CompCol_Matrix(&A1, m, n, nnz, a1, asub1, xa1,
//                            SLU_NC, SLU_D, SLU_GE);
//     dCreate_Dense_Matrix(&B1, m, nrhs, rhsb1, m, SLU_DN, SLU_D, SLU_GE);

    dgssvx(&options, &A, perm_c, perm_r, etree, equed, R, C,
           &L, &U, work, lwork, &B, &X, &rpg, &rcond, ferr, berr,
           &mem_usage, &stat, &info);

  if (noisy > 0){

//begin JCGO
//    printf("\n Resolving system: dgssvx() returns info %d\n", info);
      printf("\n Solving system: dgssvx() returns info %d\n", info);
//end JCGO

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


  template<typename T>
  void Superlu<T>::get_solution(Vector<T>& x_out)
  { 
    rhsx = (T*) ((DNformat*) B.Store)->nzval;
    
    for (i = 0; i < m; ++i) x_out.writeElement(rhsx[i], i);
    
  }
  
  //begin  JCGO 01/04/09

template<typename T>
void Superlu<T>::get_solutionB(std::vector<T>& x_out)
{ 
    rhsx = (T*) ((DNformat*) B.Store)->nzval;

    for (i = 0; i < m; ++i) x_out[i] = rhsx[i];

}

  
  template<typename T>
  void Superlu<T>::setVectors(Vector<T>& b_in)
  {
    if ( !(rhsx = doubleMalloc(m)) ) ABORT("Malloc fails for rhsx[].");
    if ( !(rhsb = doubleMalloc(m)) ) ABORT("Malloc fails for rhsb[].");
    for (i = 0; i < m; ++i) {
      rhsx[i] = b_in.readElement(i);
      rhsb[i] = b_in.readElement(i);
    }
    
  }
  
template<typename T>
void Superlu<T>::setb(T * b_in)
{
}

template<typename T>
void Superlu<T>::setb(std::vector<T>& b_in)
{
	/*
    SuperMatrix	Baux;
    double  *rhsbAux;
	for (i = 0; i < m; ++i)
	{
		rhsbAux[i] = b_in[i];
	}
    dCreate_Dense_Matrix(&Baux, m, nrhs, rhsbAux, m, SLU_DN, SLU_D, SLU_GE);
    B = Baux;
    */
    
    //It is assumed that m is correct
    for (i = 0; i < m; ++i)
	{
		rhsb[i] = b_in[i];
	}
	dCreate_Dense_Matrix(&B, m, nrhs, rhsb, m, SLU_DN, SLU_D, SLU_GE);
}

template<typename T>
void Superlu<T>::setA(size_type& m_in,
              	size_type& n_in,
              	size_type& nnz_in,
              	std::vector<size_type>& asub_in,
              	std::vector<size_type>& xa_in,
              	std::vector<T>& a_in )
{
/*
	SuperMatrix AAux;
	int mAux = m_in; 
	int nAux = n_in;
	int nnzAux = nnz_in;
	double *aAux;
    int  *asubAux, *xaAux;

	if ( !(aAux = doubleMalloc(nnzAux)) ) ABORT("Malloc fails for aAux[].");
  	if ( !(asubAux = intMalloc(nnzAux)) ) ABORT("Malloc fails for asubAux[].");
  	if ( !(xaAux = intMalloc(nAux+1)) ) ABORT("Malloc fails for xaAux[].");

    for (i = 0; i < nnzAux; ++i) {
      aAux[i] = a_in[i];
      asubAux[i] = asub_in[i]-1;
    }
 
    for (i = 0; i < nAux+1; ++i)	xaAux[i] = xa_in[i]-1;

    dCreate_CompCol_Matrix(&AAux, mAux, nAux, nnzAux, aAux, asubAux, xaAux,SLU_NC, SLU_D, SLU_GE); 
  
    A = AAux;
*/
	//It is assumed that m, n and nnz are correct
	for (i = 0; i < nnz; ++i) {
      a[i] = a_in[i];
      asub[i] = asub_in[i]-1;
    }
    for (i = 0; i < n+1; ++i)	xa[i] = xa_in[i]-1;

    dCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa,SLU_NC, SLU_D, SLU_GE);
}	
//end JCGO

} // namespace lmx;

#endif
