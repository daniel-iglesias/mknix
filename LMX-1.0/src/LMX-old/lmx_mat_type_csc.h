/***************************************************************************
 *   Copyright (C) 2005 by Roberto Ortega                                 *
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
#ifndef LMXTYPE_CSC_H
#define LMXTYPE_CSC_H

#include "lmx_mat_data_mat.h"
#include "lmx_base_iohb.h"

#ifdef HAVE_SUPERLU
#include "lmx_linsolvers_superlu_interface.h"
#endif

#ifdef HAVE_GMM
#include"lmx_mat_type_gmm.h"
#endif

//////////////////////////////////////////// Doxygen file documentation entry:
    /*!
      \file lmx_mat_type_csc.h

      \brief This file contains both the declaration and implementation for Type_csc class member functions.

      \author Roberto Ortega Aguilera.

     */
//////////////////////////////////////////// Doxygen file documentation (end)


namespace lmx {

  /**
    \class Type_csc

    \brief Compressed Storaged by Columns matrix.

    \author Roberto Ortega
  */
template <typename T> class Type_csc : public Data_mat<T>
{
  // Should be changed to private but needs comprobation of compiler compliance using
  // template friend funtions (needed for lmx_mat_data_blas.h functions).
public:
  std::vector<T>     aa; /**< Matrix data contents. */
  std::vector<size_type> ia; /**< Row indexer. */
  std::vector<size_type> ja; /**< Column indexer. */
  size_type Nnze; /**< Number of non-zero values. */
  size_type Nrow; /**< Number of rows. */
  size_type Ncol; /**< Number of columns. */
#ifdef HAVE_SUPERLU
  Superlu<T>* S;
#endif
  

private:
  T zero;

public:
  Type_csc();

  Type_csc(size_type rows, size_type columns);

  ~Type_csc();

  void resize(size_type, size_type);

  const T& readElement(const size_type& mrows, const size_type& ncolumns) const;

  void writeElement(T value, const size_type& mrows, const size_type& ncolumns);

  inline void addElement(T value, const size_type& mrows, const size_type& ncolumns);

  T* create_element(size_type, size_type);

  inline size_type getRows() const;

  inline size_type getCols() const;

  void equals(const Data<T>* matrix_in);

  void add(const Data<T>* matrix_in_1);

  void substract(const Data<T>* matrix_in_1);

  void multiply(const Data<T>* matrix_in_1, const Data<T>* matrix_in_2);

  void multiplyScalar(const T&);

  void multiplyElements(const Data<T>* matrix_in);

  void trn();

  void cleanBelow(const double factor);

  void clear();
  
  //begin JCGO 18/03/09
  void reset();  
  //end JCGO

  void read_mm_file(const char* input_file);

  void read_hb_file(const char* input_file);

  void write_hb_file(const char* input_file);

  /** Function for testing purposes. */
  void print()
  {
    std::cout << "aa:" << endl;
    for(int i=0; i<Nnze; ++i) std::cout << aa[i] << "; ";
    std::cout << endl;
    std::cout << "ia:" << endl;
    for(int i=0; i<Nnze; ++i) std::cout << ia[i] << "; ";
    std::cout << endl;
    std::cout << "ja:" << endl;
    for(int i=0; i<ja.size(); ++i) std::cout << ja[i] << "; ";
    std::cout << endl;
  }

  bool exists( size_type, size_type );

  void setSparsePattern( Vector<size_type>&, Vector<size_type>& );

  void setSparsePattern( std::vector<size_type>&, std::vector<size_type>& );

#ifdef HAVE_SUPERLU
  void initSLU();
  
  void factorize();
  
  void subsSolve(Vector<T>& rhs);  
#else
  void initSLU(){}
  
  void factorize(){}
  
  void subsSolve(Vector<T>& rhs){}  
  
#endif

  friend void mat_vec_mult<>( const Type_csc<T>*,
                              const Type_stdVector<T>*,
                              Type_stdVector<T>*);
};


/// Empty constructor.
template <typename T> 
    Type_csc<T>::Type_csc()
#ifdef HAVE_SUPERLU
  : S(0)
#endif
{
    Nrow = 0;
    Ncol = 0;
    Nnze = 0;
    ja.push_back( 1 );
    zero = 0;
}

/**
 * Standard constructor.
 * Creates a new object with parameter contents resized to (rows, columns) dimension.
 * \param rows Rows of H-B matrix.
 * \param columns Columns of H-B matrix.
 */
template <typename T>
    Type_csc<T>::Type_csc(size_type rows, size_type columns)
  : Data_mat<T>()
#ifdef HAVE_SUPERLU
  , S(0)
#endif
{
  Nrow = rows;
  Ncol = columns;
  Nnze = 0;
  for(size_type i=0;i<columns+1;++i){
  ja.push_back( 1 );
  }
  this->resize(rows,columns);
  zero = 0;
}

/// Destructor.
template <typename T>
    Type_csc<T>::~Type_csc()
{
}

/**
 * Resize method.
 * Changes the size of the sparse matrix.
 * \param mrows New value for rows of matrix.
 * \param ncolumns New value for columns of matrix.
 */
template <typename T>
    void Type_csc<T>::resize(size_type mrows, size_type ncolumns)
{
  size_type end_val_ja;
  size_type pos_end_col;
  size_type row_end_col;
  size_type num_col_del;
  size_type num_val_del;
  typename std::vector<T>::iterator it_aa;
  typename std::vector<size_type>::iterator it_ia;

  if( mrows < Nrow ){

    for( size_type i=0 ; i < Ncol ; ++i ){

      if( ja[i] != ja [i+1] ) {

        pos_end_col = ja[i+1]-2;
        row_end_col = ia[pos_end_col];

        while( row_end_col > mrows ){
          it_aa = aa.begin() + pos_end_col;
          aa.erase( it_aa );
          it_ia = ia.begin() + pos_end_col;
          ia.erase( it_ia );
          --Nnze;
          for ( size_type j=i+1 ; j < Ncol+1 ; ++j ) {
            --ja[j];
          }
          if( ja[i] != ja [i+1] ) {
            pos_end_col = ja[i+1]-2;
            row_end_col = ia[pos_end_col];
          }
          else{ break;}
        }
      }
    }
  }
  else if( mrows > Nrow ){

  }
  Nrow = mrows;

  if( ncolumns > Ncol ){

    end_val_ja = ja[Ncol];
    ja.resize( ncolumns + 1 );

    for( size_type i = Ncol+1 ; i < ncolumns+1 ; ++i){
      ja[i] = end_val_ja;
    }
  }
  else if( ncolumns < Ncol ){

    num_col_del = Ncol - ncolumns ;
    num_val_del = ja[Ncol] - ja[ncolumns];

    aa.resize( Nnze-num_val_del );
    ia.resize( Nnze-num_val_del );
    ja.resize( Ncol+1-num_col_del );
    Nnze = Nnze - num_val_del;
  }
  Ncol = ncolumns;
}


/**
 * Read element method.
 * Implements a method for reading data of the Harwell-Boeing matrix.
 * \param mrows Row position in Harwell-Boeing matrix.
 * \param ncolumns Column position in Harwell-Boeing matrix.
 * \return Value of the element in the position given by the parameters.
 */
template <typename T>
    const T& Type_csc<T>::readElement(const size_type& mrows, const size_type& ncolumns) const
{
  if(ja[ncolumns]==ja[ncolumns+1]){
    return zero;
  }
  else{
    for (size_type i=ja[ncolumns];i<ja[ncolumns+1];++i){
      if (mrows+1==ia[i-1]){
        return aa[i-1];
        break;
      }
    }
    return zero;
  }
}

/**
 * Write element method.
 * Implements a method for writing data on the Harwell-Boeing matrix.
 * \param mrows Row position in Harwell-Boeing matrix.
 * \param ncolumns Column position in Harwell-Boeing matrix.
 * \param value Numerical type value.
 */
template <typename T>
    void Type_csc<T>::writeElement(T value, const size_type& mrows, const size_type& ncolumns)
{
  if(ja[ncolumns]==ja[ncolumns+1]){
    *(this->create_element(mrows, ncolumns) ) = value;
  }
  else{
    if(mrows+1>ia[ja[ncolumns+1]-2]){
      *(this->create_element(mrows, ncolumns) ) = value;
    }
    else{
      for (size_type i=ja[ncolumns];i<ja[ncolumns+1];++i){
        if (mrows+1==ia[i-1]){
          aa[i-1] = value;
          return;
        }
      }
      *(this->create_element(mrows, ncolumns) ) = value;
    }
  }
}

/**
 * Add element method.
 * Implements a method for adding data on the Harwell-Boeing matrix.
 * Copy-pasted from writeElement.
 * \param mrows Row position in Harwell-Boeing matrix.
 * \param ncolumns Column position in Harwell-Boeing matrix.
 * \param value Numerical type value.
 */
template <typename T>
    void Type_csc<T>::addElement(T value, const size_type& mrows, const size_type& ncolumns)
{
  if(ja[ncolumns]==ja[ncolumns+1]){
    *(this->create_element(mrows, ncolumns) ) += value;
  }
  else{
    if(mrows+1>ia[ja[ncolumns+1]-2]){
      *(this->create_element(mrows, ncolumns) ) += value;
    }
    else{
      for (size_type i=ja[ncolumns];i<ja[ncolumns+1];++i){
        if (mrows+1==ia[i-1]){
          aa[i-1] += value;
          return;
        }
      }
      *(this->create_element(mrows, ncolumns) ) = value;
    }
  }
}

/**
 * Create element method.
 * Implements a method for creating a new data position in matrix.
 * \param mrows Row position in Harwell-Boeing matrix.
 * \param ncolumns Column position in Harwell-Boeing matrix.
 * \return A pointer to the element in the position given by the parameters.
 */
template <typename T>
    T* Type_csc<T>::create_element( size_type mrows, size_type ncolumns)
{
  if(ja[ncolumns]==ja[ncolumns+1]){
    //caso 1
    aa.insert(aa.begin()+ja[ncolumns]-1,0);
    ia.insert(ia.begin()+ja[ncolumns]-1,mrows+1);
    for(size_type i=ncolumns+1;i<ja.size();++i){
      ++ja[i];
    }
    ++Nnze;
    return &( aa[ja[ncolumns]-1] );
  }
  else{
    if(mrows+1>=ia[ja[ncolumns+1]-2]){//caso 2
      aa.insert(aa.begin()+ja[ncolumns+1]-1,0); 
      ia.insert(ia.begin()+ja[ncolumns+1]-1,mrows+1);
      for(size_type i=ncolumns+1;i<ja.size();++i){
        ++ja[i];
      }
      ++Nnze;
      return &( aa[ja[ncolumns+1]-2] );
    }
    else{
      for (size_type i=ja[ncolumns];i<ja[ncolumns+1];++i){//
        if (mrows+1<ia[i-1]){
          aa.insert(aa.begin()+i-1,0);
          ia.insert(ia.begin()+i-1,mrows+1);
          for(size_type k=ncolumns+1;k<ja.size();++k){
            ++ja[k];
          }
          ++Nnze;
          return &( aa[i-1] );
          break;
        }
      }
    }
  }
  return 0;
}

/**
 * Method for knowing the number of data rows.
 * \return Number of rows.
 */
template <typename T> inline
    size_type Type_csc<T>::getRows() const
{
  return this->Nrow;
}

/**
 * Method for knowing the number of data columns.
 * \returns Number of columns.
 */
template <typename T> inline
    size_type Type_csc<T>::getCols() const
{
  return this->Ncol;
}


/**
 * Copy method.
 * Equals the data in the object's contents to those given by the input matrix parameter.
 * \param matrix_in pointer to an object that belongs to a class derived from Data_mat.
 */
template <typename T>
    void Type_csc<T>::equals(const Data<T>* matrix_in)
{
  aa = static_cast<const Type_csc*>(matrix_in)-> aa;
  ia = static_cast<const Type_csc*>(matrix_in)-> ia;
  ja = static_cast<const Type_csc*>(matrix_in)-> ja;
  Nrow    = static_cast<const Type_csc*>(matrix_in) -> Nrow;
  Ncol    = static_cast<const Type_csc*>(matrix_in) -> Ncol;
  Nnze = static_cast<const Type_csc*>(matrix_in) -> Nnze;
}

/**
 * Add method.
 * Adds the the input matrix parameter's elements to the object's contents.
 * Necessary for overloading the "+=" operator.
 * \param matrix_in_1 pointer to an object that belongs to a class derived from Data_mat.
 */
template <typename T>
    void Type_csc<T>::add(const Data<T>* matrix_in_1)
{
  // this += matrix_in_1
  for( size_type i=0 ; i < Nrow ; ++i ){
    for( size_type j=0 ; j < Ncol ; ++j ){
      this->writeElement( this->readElement( i , j ) + matrix_in_1->readElement( i , j ) , i , j );
    }
  }
}

/**
 * Substract method.
 * Substracts the the input matrix parameter's elements to the object's contents.
 * Necessary for overloading the "-=" operator.
 * \param matrix_in_1 pointer to an object that belongs to a class derived from Data_mat.
 */
template <typename T> 
    void Type_csc<T>::substract(const Data<T>* matrix_in_1)
{
// this -= matrix_in_1
  for( size_type i=0 ; i < Nrow ; ++i ){
    for( size_type j=0 ; j < Ncol ; ++j ){
      this->writeElement( this->readElement( i , j ) - matrix_in_1->readElement( i , j ) , i , j );
    }
  }
}

/**
 * Multiply method.
 * Multiplies the input matrices and saves the result into the object's contents.
 * Necessary for overloading the "*" operator.
 * \param matrix_in_1 pointer to an object that belongs to a class derived from Data_mat.
 * \param matrix_in_2 pointer to an object that belongs to a class derived from Data_mat.
 */
template <typename T>  
    void Type_csc<T>::multiply(const Data<T>* matrix_in_1, const Data<T>* matrix_in_2)
{
// this = matrix_in_1 * matrix_in_2
  T val_mat1;
  T val_mat2;
  T val_mult;

  for( size_type i=0 ; i < matrix_in_1->getRows() ; ++i ){
    for( size_type j=0 ; j < matrix_in_2->getCols() ; ++j ){
      val_mult = 0;
      for( size_type k=0 ; k < matrix_in_1->getCols() ; ++k ){
        val_mat1 = matrix_in_1->readElement( i , k );
        val_mat2 = matrix_in_2->readElement( k , j );
        val_mult += val_mat1 * val_mat2;
      }
      this->writeElement( val_mult , i , j );
    }
  }
}

/**
 * Multiply scalar method.
 * Multiplies the object's matrix (contents) with a scalar.
 * Necessary for overloading the "*" operator.
 * \param scalar A scalar factor of template's class.
 */
template <typename T>
    void Type_csc<T>::multiplyScalar(const T& scalar)
{
  typename std::vector<T>::iterator aa_iterator;
  for(aa_iterator = aa.begin();aa_iterator != aa.end();++aa_iterator){
    *(aa_iterator) *= scalar;
  }
}

/**
 * Method multiplying element-by-element of two matrices.
 * One would be the object's contents and the other the parameter's contents.
 *
 * Necessary for implementing  Vector to Vector multElements.
 *
 * \param matrix_in pointer to an object that belongs to a class derived from Data.
 */
template <typename T>
    void Type_csc<T>::multiplyElements(const Data<T>* matrix_in)
{
  // NEEDS OPTIMIZATION FOR MULTIPLYING EXCLUSIVELY THE ELEMENTS OF aa[].
  T temp;
  for (size_type i=0; i<Nrow; ++i){
    for (size_type j=0; j<Ncol; ++j){
      temp = this->readElement(i,j);
      if( temp != static_cast<T>(0) )
        this->writeElement( temp * matrix_in->readElement(i,j), i,j);
    }
  }
}


/**
 * Traspose method.
 * Swaps elements with respect to the diagonal: A(i,j) = A(j,i)
 */
template <typename T>
    void Type_csc<T>::trn()
{

  std::vector<size_type> ja_starO(Nnze);

  size_type NrowT = Ncol;
  size_type NcolT = Nrow;

  size_type k=0;
  for( size_type i=0 ; i < Ncol ; ++i){
    for( size_type j=ja[i] ; j < ja[i+1] ; ++j){
      ja_starO[k] = i+1;
      ++k;
    }
  }

//rutina para ordenar (metodo de la burbuja)
//ordena de acuerdo a ia[] y modifica aa[] and jatmp[]
  T aux_aa;
  size_type aux_ia;
  size_type aux_ja;

  for( size_type j=Nnze-1 ; j > 0 ; j-- ){
    for( size_type i=0 ; i < j ; ++i ){
      if ( ia[i] > ia[i+1] ){
        aux_ia = ia[i];
        aux_aa = aa[i];
        aux_ja = ja_starO[i];

        ia[i] = ia[i+1];
        aa[i] = aa[i+1];
        ja_starO[i] = ja_starO[i+1];

        ia[i+1]    = aux_ia;
        aa[i+1]    = aux_aa;
        ja_starO[i+1] = aux_ja;
      }
    }
  }

  ja.resize( Nrow + 1 );
  ja[0] = 1;
  size_type  j = 0;
  size_type  l = 1;
  for( size_type i=0 ; i < Nrow+1 ; ++i){
    while( i+1 == ia[j] ){
      ++j;
    }
    ja[l] = j+1;
    ++l;
  }

  ia.resize( Nnze );
  for( size_type i=0 ; i < Nnze ; ++i ){
    ia[i] = ja_starO[i];
  }

  Nrow = NrowT;
  Ncol = NcolT;

}

/**
 * Clean below method.
 * Makes equal to zero every element below given factor.
 * \param factor Reference value for cleaning.
 */
template <typename T>
    void Type_csc<T>::cleanBelow(const double factor)
{
  for (size_type i=0; i<this->Nnze; ++i){
    if (std::abs(this->aa[i]) < factor) this->aa[i] = static_cast<T>(0);
  }
}

/**
 * Clear method.
 * Wipes all data, clearing column, row and data vectors.
 */
template <typename T>
    void Type_csc<T>::clear()
{
  aa.clear(); /**< Matrix data contents. */
  ia.clear(); /**< Row indexer. */
  ja.clear(); /**< Column indexer. */
  Nnze = 0; /**< Number of non-zero values. */
}

//begin JCGO 18/03/09
/**
 * Reset method.
 */
template <typename T>
    void Type_csc<T>::reset()
{
	for (size_type i=0; i<this->Nnze; ++i)	this->aa[i] = static_cast<T>(0); 
}
//end JCGO

/**
 * Read data in Matrix Market format method.
 * Opens the file specified and reads the matrix's data in it,
 * suposing it's stored in Matrix Market format.
 *
 * \param input_file Name of the file to be read.
 */
template <typename T> 
void Type_csc<T>::read_mm_file(const char* input_file)
{
#ifndef HAVE_GMM
  std::stringstream message;
  message
      << "---------------------------------------------\n"
      << "We are very very sorry but this feature (read_mm_file) requires at the moment the use of \nthe gmm++ library.\n"
      << "If you really feel that it is important to have the Matrix Market file reading \nas a built in feature, you can post it into our web page:\n"
      << "w3.mecanica.upm.es/lmx\n"
      << "Thank you very much for your interest and help in setting our priorities in the \ndevelopment schedule.\n"
      << "---------------------------------------------\n"
      << "gmm++ not defined.\nYou must set \"#define HAVE_GMM\" in your file in order to use this library." << endl;
  LMX_THROW(failure_error, message.str() );
#else
  gmm::csc_matrix<double> cscmat;
  gmm::col_matrix< gmm::wsvector<double> > matrix_loaded;
  gmm::MatrixMarket_load(input_file, matrix_loaded);
//   gmm::resize(contents, gmm::mat_nrows(matrix_loaded), gmm::mat_ncols(matrix_loaded));
  gmm::copy(matrix_loaded, cscmat);
  
  aa.clear();
  ja.clear();
  ia.clear();

  Ncol = cscmat.nc;
  Nrow = cscmat.nr;
  for (size_type j = 0; j<Ncol+1; ++j){
    ja.push_back( cscmat.jc[j]+1 );
  }
  cout << ja.size();
  Nnze = ja[Ncol];
  for (size_type i = 0; i<Nnze; ++i){
    ia.push_back( cscmat.ir[i]+1 );
    aa.push_back( cscmat.pr[i] );
  }  
#endif
}

/**
 * Read data in Harwell-Boeing format method.
 * Opens the file specified and reads the matrix's data in it,
 * suposing it's stored in Harwell-Boeing format.
 *
 * \param input_file Name of the file to be read.
 */
template <typename T> 
    void Type_csc<T>::read_hb_file(const char* input_file)
{
  int M;
  int N;
  int nonzeros;
  int *colptr;
  int *rowind;
  double *val;


  HarwellBoeing_IO h(input_file);
  h.read(M,N,nonzeros,colptr,rowind,val);

  Nrow = M;
  Ncol = N;
  Nnze = nonzeros;

  aa.clear();
  ja.clear();
  ia.clear();
  
  for(int i = 0 ; i < nonzeros ; ++i){
    aa.push_back(val[i]);
    ia.push_back(rowind[i]);
  }

  for(int i = 0 ; i < N+1 ; ++i){
    ja.push_back(colptr[i]);
  }
}

/**
 * Write data in Harwell-Boeing format method.
 * Opens the file specified and writes the matrix's data in it.
 *
 * \param input_file Name of the file to be read.
 */
template <typename T> 
    void Type_csc<T>::write_hb_file(const char* input_file)
{
#ifndef HAVE_GMM
  std::stringstream message;
  message
      << "---------------------------------------------\n"
      << "We are very very sorry but this feature (write_hb_file) requires at the moment the use of \nthe gmm++ library.\n"
      << "If you really feel that it is important to have the Matrix Market file reading \nas a built in feature, you can post it into our web page:\n"
      << "w3.mecanica.upm.es/lmx\n"
      << "Thank you very much for your interest and help in setting our priorities in the \ndevelopment schedule.\n"
      << "---------------------------------------------\n"
      << "gmm++ not defined.\nYou must set \"#define HAVE_GMM\" in your file in order to use this library." << endl;
  LMX_THROW(failure_error, message.str() );
#else
  typedef unsigned int IND_TYPE;
  gmm::csc_matrix<T> cscmat;
  
  if (cscmat.pr) { delete[] cscmat.pr; delete[] cscmat.ir; delete[] cscmat.jc; }

  cscmat.nr = this->Nrow;
  cscmat.nc = this->Ncol;
  
  cscmat.jc = new IND_TYPE[cscmat.nc+1];

  for (size_type j = 0; j<Ncol+1; ++j){
    cscmat.jc[j] = ja[j]-1;
  }
  
  cscmat.pr = new T[cscmat.jc[cscmat.nc]];
  cscmat.ir = new IND_TYPE[cscmat.jc[cscmat.nc]];
  for (size_type i = 0; i<Nnze; ++i){
    cscmat.ir[i] = ia[i]-1;
    cscmat.pr[i] = aa[i];
  }

  gmm::Harwell_Boeing_save(input_file, cscmat);
#endif
}


/**
 * Returns TRUE or FALSE depending of element existance.
 * 
 * \param mrows Row index.
 * \param ncolumns Column index.
 * @return TRUE if the element exists in internal storage structure.
 */
template <typename T> 
    bool Type_csc<T>::exists( size_type mrows, size_type ncolumns )
{
  if(ja[ncolumns]==ja[ncolumns+1]){
    return 0;
  }
  else{
    if(mrows+1>ia[ja[ncolumns+1]-2]){
      return 0;
    }
    else{
      for (size_type i=ja[ncolumns];i<ja[ncolumns+1];++i){
        if (mrows+1==ia[i-1]){
          return 1;
        }
      }
      return 0;
    }
  }
}

/**
 * Prepares the sparse structure of a CSC matrix.
 * @param row_index CSC row indices.
 * @param col_index CSC columns indices.
 */
template <typename T> 
    void Type_csc<T>::setSparsePattern( Vector<size_type>& row_index,
                                        Vector<size_type>& col_index 
                                      )
{
  int i;
  ia.clear();
  ja.clear();
  aa.clear();
  for (i=0; i<row_index.size(); ++i ){
    ia.push_back( row_index(i) );
    aa.push_back( T(1) );
  }
  for (i=0; i<col_index.size(); ++i ) ja.push_back( col_index(i) );
  Nnze = aa.size();
}

/**
 * Prepares the sparse structure of a CSC matrix.
 * @param row_index CSC row indices.
 * @param col_index CSC columns indices.
 */
template <typename T> 
    void Type_csc<T>::setSparsePattern( std::vector<size_type>& row_index, 
                                        std::vector<size_type>& col_index 
                                      )
{
  int i;
  ia.clear();
  ja.clear();
  aa.clear();
  ia = row_index;
  for (i=0; i<row_index.size(); ++i ){
    aa.push_back( T(0) );
  }
  ja = col_index;
  Nnze = aa.size();
}

  
#ifdef HAVE_SUPERLU
  template <class T>
  void Type_csc<T>::initSLU()
  {
    S = new Superlu<T>( Nrow,
                        Ncol,
                        Nnze,
                        ia,
                        ja,
                        aa );
    S->initMatrix();
  }
  
  template <class T>
  void Type_csc<T>::factorize()
  {
    if (S == 0) initSLU();
    S->factorize();
  }

  template <class T>
  void Type_csc<T>::subsSolve(Vector<T>& rhs_in){
    S->setVectors( rhs_in ); // hace x, b.
    S->initVectors();
    S->subsSolve();
    
    S->get_solution(rhs_in);
  }
#endif

};


#endif
