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
#ifndef LMXTYPE_GMM_SPARSE1_H
#define LMXTYPE_GMM_SPARSE1_H

#include "lmx_mat_data_mat.h"

//////////////////////////////////////////// Doxygen file documentation entry:
    /*!
      \file lmx_mat_type_gmm_sparse1.h
      
      \brief This file contains both the declaration and implementation for type_gmm_sparse (sparse gmm Matrix) class member functions.
      
      \author Daniel Iglesias 
      
    */
//////////////////////////////////////////// Doxygen file documentation (end)


namespace lmx {

    /**
    \class Type_gmm_sparse 
    \brief Template class Type_gmm_sparse
    
    This class implements the methods defined in virtual class data_mat so the matrix type "gmm::row_matrix< gmm::rsvector<T> >" can be used in lmx. For details about the caracteristics of this matrix type, see the gmm library manual and code.
    
    @param contents Corresponds to a gmm::row_matrix< gmm::rsvector<T> > and it's the base of the methods implemented for this class.
    
    @author Daniel Iglesias .
    */
template <typename T> class Type_gmm_sparse : public Data_mat<T>
{
private:
  /** Matrix data contents */
    gmm::row_matrix< gmm::rsvector<T> > contents;
//     T reference;
    
public:
    
    /// Empty constructor.
    Type_gmm_sparse(){}

    Type_gmm_sparse(size_type, size_type);

    /// Destructor.
    ~Type_gmm_sparse(){}

  void resize(size_type, size_type);

  void factorize(){}
  
  /** Read element method.
    * Implements a method for reading data of the sparse matrix.
    * \param mrows Row position in sparse matrix.
    * \param ncolumns Column position in sparse matrix.
    * \return Value of the element in the position given by the parameters.
    */
    const T& readElement(const size_type& mrows, const size_type& ncolumns) const
  { static T reference = 0;
    return reference = contents(mrows, ncolumns);
  }

    /** Write element method.
      * Implements a method for writing data on the sparse matrix.
      * \param mrows Row position in sparse matrix.
      * \param ncolumns Column position in sparse matrix.
      * \param value Numerical type value. */
  void writeElement(T value, const size_type& mrows, const size_type& ncolumns)
  { contents(mrows, ncolumns) = value; }

    /** Add element method.
      * Implements a method for adding data on the Harwell-Boeing matrix.
      * Copy-pasted from writeElement.
      * \param mrows Row position in sparse matrix.
      * \param ncolumns Column position in sparse matrix.
      * \param value Numerical type value. */
  inline void addElement(T value, const size_type& mrows, const size_type& ncolumns)
  { contents(mrows, ncolumns) += value; }

   /** Method for getting the number of data rows. 
    * \returns Number of rows.
    */
  size_type getRows() const
  { return gmm::mat_nrows(contents); }

   /** Method for knowing the number of data columns. 
    * \returns Number of columns.
    */
  size_type getCols() const
  { return gmm::mat_ncols(contents); }

    /** Copy method.
      * Equals the data in the object's contents to those given by the input matrix parameter.
      * \param matrix_in pointer to an object that belongs to a class derived from Data. */
  void equals(const Data<T>* matrix_in)
  { contents = static_cast<const Type_gmm_sparse*>(matrix_in)->contents; }

    /** Add method.
      * Adds the the input matrix parameter's elements to the object's contents.
      * Necessary for overloading the "+=" operator.
      * \param matrix_in_1 pointer to an object that belongs to a class derived from Data. */
  void add(const Data<T>* matrix_in_1)
  { gmm::add(static_cast<const Type_gmm_sparse*>(matrix_in_1)->contents, contents); }

    /** Substract method.
      * Substracts the the input matrix parameter's elements to the object's contents.
      * Necessary for overloading the "-=" operator.
      * \param matrix_in_1 pointer to an object that belongs to a class derived from Data. */
  void substract(const Data<T>* matrix_in_1)
  { gmm::scale(contents, -1.);
    gmm::add(static_cast<const Type_gmm_sparse*>(matrix_in_1)->contents, contents);
    gmm::scale(contents, -1.);
  }

    /** Multiply method.
      * Multiplies the input matrices and saves the result into the object's contents.
      * Necessary for overloading the "*" operator.
      * \param matrix_in_1 pointer to an object that belongs to a class derived from Data.
      * \param matrix_in_2 pointer to an object that belongs to a class derived from Data. */
  void multiply(const Data<T>* matrix_in_1, const Data<T>* matrix_in_2)
  { gmm::mult(static_cast<const Type_gmm_sparse*>(matrix_in_1)->contents,
              static_cast<const  Type_gmm_sparse*>(matrix_in_2)->contents,
              contents);
  }
   
    /** Multiply scalar method.
      * Multiplies the object's matrix (contents) with a scalar.
      * Necessary for overloading the "*" operator.
      * \param scalar A scalar factor of template's class. */
  void multiplyScalar(const T& scalar)
   { gmm::scale(contents, scalar);
   }
   
  /** Method multiplying element-by-element of two matrices. One would be the object's contents and the other the parameter's contents.
    * Necessary for implementing  Vector to Vector multElements.
    * \param matrix_in pointer to an object that belongs to a class derived from Data.
   */
   void multiplyElements(const Data<T>* matrix_in)
   {
     for (size_type i=0; i<gmm::mat_nrows(contents); ++i){
       for (size_type j=0; j<gmm::mat_ncols(contents); ++j)
         contents(i,j) *= matrix_in->readElement(i,j);
     }
   }


    /** Read data in Matrix Market format method.
     * Opens the file specified and reads the matrix's data in it, 
     * suposing it's stored in Matrix Market format.
     * \param input_file Name of the file to be read.
     *  */
  void read_mm_file(const char* input_file)
  { gmm::col_matrix< gmm::rsvector<double> > matrix_loaded;
    gmm::MatrixMarket_load(input_file, matrix_loaded);
    gmm::resize(contents, gmm::mat_nrows(matrix_loaded), gmm::mat_ncols(matrix_loaded));
    gmm::copy(matrix_loaded, contents);
  }

    /** Read data in Harwell-Boeing format method.
     * Opens the file specified and reads the matrix's data in it, 
     * suposing it's stored in Harwell-Boeing format.
     * \param input_file Name of the file to be read.
     *  */
  void read_hb_file(const char* input_file)
  { gmm::csc_matrix<double> matrix_loaded;
    gmm::Harwell_Boeing_load(input_file, matrix_loaded);
    gmm::resize(contents, gmm::mat_nrows(matrix_loaded), gmm::mat_ncols(matrix_loaded));
    gmm::copy(matrix_loaded, contents);
  }

/**
   * Write data in Harwell-Boeing format method.
   * Opens the file specified and writes the matrix's data in it.
   *
   * \param input_file Name of the file to be read.
 */
  void write_hb_file(const char* input_file)
  {
    gmm::csc_matrix<double> cscmat( gmm::mat_nrows(this->contents), gmm::mat_ncols(this->contents) );
    gmm::copy(this->contents, cscmat);
    gmm::Harwell_Boeing_save(input_file, cscmat);
  }

    /** Traspose method.
      * Swaps elements with respect to the diagonal: A(i,j) = A(j,i) 
      */
   void trn()
   {
     gmm::row_matrix< gmm::rsvector<T> > temp( gmm::mat_ncols(contents), gmm::mat_nrows(contents) );

     copy( transposed(contents), temp );
     gmm::resize(contents, gmm::mat_nrows(temp), gmm::mat_ncols(temp));
     copy(temp, contents);
   }

    /** Clean below method.
      * Makes equal to zero every element below given factor.
      * \param factor Reference value for cleaning. */
   void cleanBelow(const double factor)
   { gmm::clean(contents, factor);
   }

   /**
    * Clear method.
    * Wipes all data. 
    */
   //TODO: Eliminate all information and size
   void clear()
   {
     gmm::clear(contents); /**< Matrix data contents. */
   }

	//begin JCGO 18/03/09
   /**
    * Reset method.
    * Resets all elements to 0.
    */
   void reset()
   {
   	gmm::clear(contents);
   }
   //end JCGO

  /** Data pointer method
   * Gives the direction in memory of (pointer to) the object.
   * @return A pointer to the matrix's contents (Type_gmm_sparse).
   */
  gmm::row_matrix< gmm::rsvector<T> >* data_pointer()
  {return &(this->contents);}

  bool exists( size_type, size_type )
  { return 1; }

}; // class Type_gmm_sparse definitions.

  ////////////////////////////////////////////////// Implementing the methods defined previously:

    /** Standard constructor.
      * Creates a new object with parameter contents resized to (rows, columns) dimension.
      * \param rows Rows of sparse matrix.
      * \param columns Columns of sparse matrix. */
  template <typename T>
  Type_gmm_sparse<T>::Type_gmm_sparse(size_type rows, size_type columns) :
   Data_mat<T>()
  {
    gmm::resize(contents, rows, columns);
  }

    /** Resize method.
      * Changes the size of the contents parameter.
      * \param mrows New value for rows of sparse matrix.
      * \param ncolumns New value for columns of sparse matrix. */
  template <typename T>
  void Type_gmm_sparse<T>::resize(size_type mrows, size_type ncolumns)
  {
    gmm::resize(contents, mrows, ncolumns);
  }

}; // namespace lmx


#endif
