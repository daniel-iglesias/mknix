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
#ifndef LMXDATA_TYPE_GMM_CSC_H
#define LMXDATA_TYPE_GMM_CSC_H

#ifdef HAVE_GMM
  #include "gmm/gmm.h"
#endif

#include "lmx_mat_data_mat.h"


/** Member size_type is a redefinition of the numeric type used for indexing the matrices. Actually, it's defined as size_t integer.
*   */
typedef size_t size_type;

namespace lmx {

    /**
    \class Data
    \brief Template class Data.
    Mother class for Data_mat & Data_vec.

    This class represents the skeleton for the data container used by the Vector and Matrix classes. No parameter nor function implementation here, just pure virtual class. See derived classes for details in implementation. Also maybe useful to see how this class is used in the Vector class.

    @author Daniel Iglesias.
    */
template <typename T> class Type_gmm_csc : public Data_mat<T>{

private:
  /** Matrix data contents */
    gmm::csc_matrix<T> gmm_csc_matrix;
//
public:

  /** Empty constructor.
   */
  Type_gmm_csc(){}

  Type_gmm_csc(size_type, size_type);

  /** Destructor.
   */
  ~Type_gmm_csc(){}

  /** Resize method for augmenting or reducing the container's dimension.
   */
    void resize(size_type, size_type);

    void factorize(){}

  /** Read method for accesing stored data.
   * \returns Value of the data stored in the (size_type,size_type) position.
   */
  const T& readElement(const size_type& row, const size_type& column) const
  { static T reference = 0;
    return reference = gmm_csc_matrix(row, column);
  }

  /** Write method for storing data.
   */
  void writeElement(T value, const size_type& row, const size_type& column)
  {
    gmm_csc_matrix(row, column) = value;
  }

  /** Write method for adding to data.
   */
  void addElement(T value, const size_type& row, const size_type& column)
  {
    gmm_csc_matrix(row, column) += value;
  }

  /** Method for knowing the number of data rows.
   * \returns Number of rows.
   */
  size_type getRows()const { return gmm_csc_matrix.nrows(); }

  /** Method for knowing the number of data columns.
   * \returns Number of columns.
   */
  size_type getCols() const { return gmm_csc_matrix.ncols();}

  /** Method for equaling object's data to the parameter's data.
   * Necessary for implementing the "=" overloaded operator method in Vector and Vector class types.
   */
  void equals(const Data<T>* matrix_in)
{ gmm_csc_matrix = static_cast<const Type_gmm_csc*>(matrix_in)->gmm_csc_matrix; }

  /** Method for adding object's data to the parameter's data.
   * Necessary for implementing the "+=" overloaded operator method in Vector and Vector class types.
   */
  void add(const Data<T>* matrix_in_1)
  { gmm::add(static_cast<const Type_gmm_csc*>(matrix_in_1)->gmm_csc_matrix, gmm_csc_matrix); }

  /** Method for substracting object's data from the parameter's data.
   * Necessary for implementing the "-=" overloaded operator method in Vector and Vector class types.
   */
  void substract(const Data<T>* matrix_in_1)
  { gmm::scale(gmm_csc_matrix, -1.);
    gmm::add(static_cast<const Type_gmm_csc*>(matrix_in_1)->gmm_csc_matrix, gmm_csc_matrix);
    gmm::scale(gmm_csc_matrix, -1.);
  }

  /** Multiply method.
      * Multiplies the input matrices and saves the result into the object's contents.
      * Necessary for overloading the "*" operator.
      * \param matrix_in_1 pointer to an object that belongs to a class derived from Data.
      * \param matrix_in_2 pointer to an object that belongs to a class derived from Data. */
   void multiply(const Data<T>* matrix_in_1, const Data<T>* matrix_in_2)
   { gmm::mult(static_cast<const Type_gmm_csc*>(matrix_in_1)->gmm_csc_matrix,
               static_cast<const  Type_gmm_csc*>(matrix_in_2)->gmm_csc_matrix,
               gmm_csc_matrix);
   }


  /** Method for multiplying a scalar with an array and saving the result in the object's data.
   * Necessary for implementing the "*" overloaded operator method with scalars in Vector class type.
   */
  void multiplyScalar(const T& scalar)
  {
    gmm::scale(gmm_csc_matrix, scalar);
  }

  /** Method multiplying element-by-element of two arrays. One would be the object's contents and the other the parameter's contents.
   * Necessary for implementing  Vector to Vector multElements.
   */
  void multiplyElements(const Data<T>* matrix_in)
   {
     for (size_type i=0; i<gmm_csc_matrix.nrows(); ++i){
       for (size_type j=0; j<gmm_csc_matrix.ncols(); ++j)
         gmm_csc_matrix(i,j) *= matrix_in->readElement(i,j);
     }
   }

   /** Read data in Matrix Market format method.
      * Opens the file specified and reads the matrix's data in it,
      * suposing it's stored in Matrix Market format.
      * \param input_file Name of the file to be read.
      *  */
     void read_mm_file(const char* input_file)
     { gmm::col_matrix< gmm::wsvector<double> > matrix_loaded;
       gmm::MatrixMarket_load(input_file, matrix_loaded);
       gmm::resize(gmm_csc_matrix, gmm::mat_nrows(matrix_loaded), gmm::mat_ncols(matrix_loaded));
       gmm::copy(matrix_loaded, gmm_csc_matrix);
     }

     /** Read data in Harwell-Boeing format method.
      * Opens the file specified and reads the matrix's data in it,
      * suposing it's stored in Harwell-Boeing format.
      *
      * \param input_file Name of the file to be read.
      */
     void read_hb_file(const char* input_file)
     {
       gmm::csc_matrix<double> matrix_loaded;
       gmm::Harwell_Boeing_load(input_file, matrix_loaded);
       gmm::resize(this->gmm_csc_matrix, gmm::mat_nrows(matrix_loaded), gmm::mat_ncols(matrix_loaded));
       gmm::copy(matrix_loaded, gmm_csc_matrix);
     }

   /**
      * Write data in Harwell-Boeing format method.
      * Opens the file specified and writes the matrix's data in it.
      *
      * \param input_file Name of the file to be read.
    */
     void write_hb_file(const char* input_file)
     {
       gmm::csc_matrix<double> cscmat( gmm::mat_nrows(this->gmm_csc_matrix), gmm::mat_ncols(this->gmm_csc_matrix) );
       gmm::copy(this->gmm_csc_matrix, cscmat);
       gmm::Harwell_Boeing_save(input_file, cscmat);
     }

     ///added by vicen april 2017
       /** Cast a csc matrix format method.,
        * suposing it's compressed Sparse Column.
        * \param reference to the matrix to cast.
        *  */
     void cast_csc_matrix(gmm::csc_matrix<T> &matrix_to_cast)
     {
       std::cout << "GOOD CAST :-) type_gmm_csc" << std::endl;
        gmm::resize(this->gmm_csc_matrix, gmm::mat_nrows(matrix_to_cast), gmm::mat_ncols(matrix_to_cast));
        gmm::copy(matrix_to_cast, gmm_csc_matrix);
     }

       /** Traspose method.
         * Swaps elements with respect to the diagonal: A(i,j) = A(j,i) */
      void trn()
      {
        gmm::dense_matrix<T> temp( gmm::mat_ncols(gmm_csc_matrix), gmm::mat_nrows(gmm_csc_matrix) );

        copy( transposed(gmm_csc_matrix), temp );
        gmm::resize(gmm_csc_matrix, gmm::mat_nrows(temp), gmm::mat_ncols(temp));
        copy(temp, gmm_csc_matrix);
      }

  /**
   * Method for cleaning all elements below a given factor.
   */
  void cleanBelow(const double factor)
  {
    gmm::clean(gmm_csc_matrix, factor);
  }

  /**
   * Method for clearing all elements.
   */
  void clear(){
    gmm::clear(gmm_csc_matrix);
  }

	//begin JCGO 18/03/09
  /**
   * Method for all elements to 0
   */
  void reset() {
    gmm::clear(gmm_csc_matrix);
  }
  //end JCGO

  /** Data pointer method
     * Gives the direction in memory of (pointer to) the object.
     * @return A pointer to the matrix's contents (Type_gmm).
     */
    gmm::dense_matrix<T>* data_pointer()
    { return &(this->gmm_csc_matrix.pr); }

    bool exists( size_type, size_type )
    { return 1; }

  }; // class Type_gmm definitions.

    ////////////////////////////////////////////////// Implementing the methods defined previously:

      /** Standard constructor.
        * Creates a new object with parameter contents resized to (rows, columns) dimension.
        * \param rows Rows of dense matrix.
        * \param columns Columns of dense matrix. */
    template <typename T>
    Type_gmm_csc<T>::Type_gmm_csc(size_type rows, size_type columns) :
     Data_mat<T>()
      {
        gmm::resize(gmm_csc_matrix, rows, columns);
      }

      /** Resize method.
        * Changes the size of the contents parameter.
        * \param mrows New value for rows of dense matrix.
        * \param ncolumns New value for columns of dense matrix. */
    template <typename T>
    void Type_gmm_csc<T>:: resize(size_type mrows, size_type ncolumns)
      { gmm::resize(gmm_csc_matrix, mrows, ncolumns); }



};// namespace lmx

#endif//LMXDATA_TYPE_GMM_CSC_H
