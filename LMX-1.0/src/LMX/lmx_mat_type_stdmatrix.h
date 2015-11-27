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
#ifndef LMXTYPE_STDMATRIX_H
#define LMXTYPE_STDMATRIX_H

#include "lmx_mat_data_mat.h"
#include "lmx_mat_vector.h"
#include "lmx_base_iohb.h"

//////////////////////////////////////////// Doxygen file documentation entry:
/*!
  \file lmx_mat_type_stdmatrix.h

  \brief This file contains both the declaration and implementation for type_stdmatrix (std::vector type) class member functions.

  \author Daniel Iglesias

*/
//////////////////////////////////////////// Doxygen file documentation (end)


namespace lmx {

/**
\class Type_stdmatrix
\brief Template class Type_stdmatrix

This class implements the methods defined in virtual class data_mat for dense matrix built with std vectors.

@param contents Corresponds to a type std::vector< std::vector<T> >.

@author Daniel Iglesias.
*/
template<typename T>
class Type_stdmatrix : public Data_mat<T>
{
private:
    /** Matrix data contents */
    size_type rows, cols;
    std::vector<std::vector<T> > contents;

public:

    friend void copy<>(const Type_stdmatrix<T> *, Type_csc<T> *);

    friend class DenseMatrix<T>;

public:
    /// Empty constructor.
    Type_stdmatrix()
    {
        rows = 0;
        cols = 0;
    }

    Type_stdmatrix(size_type, size_type);

    /// Destructor.
    ~Type_stdmatrix()
    {
    }

    void factorize() { }

    void resize(size_type, size_type);

    /** Read element method.
      * Implements a method for reading data of the dense matrix.
      * \param mrows Row position in dense matrix.
      * \param ncolumns Column position in dense matrix.
      * \return Value of the element in the position given by the parameters. */
    const T& readElement(const size_type& mrows, const size_type& ncolumns) const { return contents[mrows][ncolumns]; }

    /** Write element method.
      * Implements a method for writing data on the dense matrix.
      * \param mrows Row position in dense matrix.
      * \param ncolumns Column position in dense matrix.
      * \param value Numerical type value. */
    void writeElement(T value, const size_type& mrows, const size_type& ncolumns) { contents[mrows][ncolumns] = value; }

    /** Add element method.
      * Implements a method for adding data on the Harwell-Boeing matrix.
      * Copy-pasted from writeElement.
      * \param mrows Row position in dense matrix.
      * \param ncolumns Column position in dense matrix.
      * \param value Numerical type value. */
    inline void addElement(T value, const size_type& mrows,
                           const size_type& ncolumns) { contents[mrows][ncolumns] += value; }

    /** Method for knowing the number of data rows.
     * \returns Number of rows.
     */
    size_type getRows() const { return rows; }

    /** Method for knowing the number of data columns.
     * \returns Number of columns.
     */
    size_type getCols() const { return cols; }

    /** Copy method.
      * Equals the data in the object's contents to those given by the input matrix parameter.
      * \param matrix_in pointer to an object that belongs to a class derived from Data. */
    void equals(const Data<T> * matrix_in) { contents = static_cast<const Type_stdmatrix *>(matrix_in)->contents; }

    /** Add method.
      * Adds the the input matrix parameter's elements to the object's contents.
      * Necessary for overloading the "+=" operator.
      * \param matrix_in_1 pointer to an object that belongs to a class derived from Data. */
    void add(const Data<T> * matrix_in_1)
    {
        for (size_type i = 0; i < rows; ++i) {
            for (size_type j = 0; j < cols; ++j) {
                contents[i][j] += matrix_in_1->readElement(i, j);
            }
        }
    }

    /** Substract method.
      * Substracts the the input matrix parameter's elements to the object's contents.
      * Necessary for overloading the "-=" operator.
      * \param matrix_in_1 pointer to an object that belongs to a class derived from Data. */
    void substract(const Data<T> * matrix_in_1)
    {
        for (size_type i = 0; i < rows; ++i) {
            for (size_type j = 0; j < cols; ++j) {
                contents[i][j] -= matrix_in_1->readElement(i, j);
            }
        }
    }


    /** Multiply method.
      * Multiplies the input matrices and saves the result into the object's contents.
      * Necessary for overloading the "*" operator.
      * \param matrix_in_1 pointer to an object that belongs to a class derived from Data.
      * \param matrix_in_2 pointer to an object that belongs to a class derived from Data. */
    void multiply(const Data<T> * matrix_in_1, const Data<T> * matrix_in_2)
    {
        // Emmit an error if called over self data...
        if (this == matrix_in_1 || this == matrix_in_2) {
            std::stringstream message;
            message << "Trying to multiply and save results on same data at the same time."
            << endl << "  This cannot be done." << endl;
            LMX_THROW(failure_error, message.str());
        }
        // This can have a good optimization...
        this->resize(matrix_in_1->getRows(), matrix_in_2->getCols());
        for (size_type i = 0; i < rows; ++i) {
            for (size_type j = 0; j < cols; ++j) {
                contents[i][j] = T(0);
            }
        }
        for (size_type i = 0; i < rows; ++i) {
            for (size_type j = 0; j < matrix_in_2->getRows(); ++j) {
                for (size_type k = 0; k < cols; ++k) {
                    contents[i][k] += matrix_in_1->readElement(i, j)
                                      * matrix_in_2->readElement(j, k);
                }
            }
        }
    }

    /** Multiply scalar method.
      * Multiplies the object's matrix (contents) with a scalar.
      * Necessary for overloading the "*" operator.
      * \param scalar A scalar factor of template's class. */
    void multiplyScalar(const T& scalar)
    {
        for (size_type i = 0; i < rows; ++i) {
            for (size_type j = 0; j < cols; ++j) {
                contents[i][j] *= scalar;
            }
        }
    }

    /** Method multiplying element-by-element of two matrices. One would be the object's contents and the other the parameter's contents.
     * Necessary for implementing  Vector to Vector multElements.
     * \param matrix_in pointer to an object that belongs to a class derived from Data.
       */
    void multiplyElements(const Data<T> * matrix_in)
    {
        for (size_type i = 0; i < rows; ++i) {
            for (size_type j = 0; j < cols; ++j) {
                contents[i][j] *= matrix_in->readElement(i, j);
            }
        }
    }


    /** Read data in Matrix Market format method.
     * Opens the file specified and reads the matrix's data in it,
     * suposing it's stored in Matrix Market format.
     * \param input_file Name of the file to be read.
     *  */
    void read_mm_file(const char * input_file)
    {
#ifndef HAVE_GMM
        std::stringstream message;
        message
        << "---------------------------------------------\n"
        <<
        "We are very very sorry but this feature (read_mm_file) requires at the moment the use of \nthe gmm++ library.\n"
        <<
        "If you really feel that it is important to have the Matrix Market file reading \nas a built in feature, you can post it into our web page:\n"
        << "w3.mecanica.upm.es/lmx\n"
        << "Thank you very much for your interest and help in setting our priorities in the \ndevelopment schedule.\n"
        << "---------------------------------------------\n"
        << "gmm++ not defined.\nYou must set \"#define HAVE_GMM\" in your file in order to use this library." << endl;
        LMX_THROW(failure_error, message.str());
#else
        gmm::csc_matrix<double> cscmat;
        gmm::col_matrix< gmm::wsvector<double> > matrix_loaded;
        gmm::MatrixMarket_load(input_file, matrix_loaded);
    //   gmm::resize(contents, gmm::mat_nrows(matrix_loaded), gmm::mat_ncols(matrix_loaded));
        gmm::copy(matrix_loaded, cscmat);

        this->resize(cscmat.nr, cscmat.nc);

        size_type col_count;
        size_type elem_count=0;

        for(int j = 0 ; j < rows+1 ; ++j){
          col_count = cscmat.jc[j];
          for(int m = cscmat.jc[j]; m < cscmat.jc[j+1]; ++m){
            this->contents [ cscmat.ir[elem_count] ] [ j ] = cscmat.pr[elem_count];
            ++elem_count;
          }
        }
        cout << col_count << endl;
#endif
    }

    /** Read data in Harwell-Boeing format method.
     * Opens the file specified and reads the matrix's data in it,
     * suposing it's stored in Harwell-Boeing format.
     *
     * \param input_file Name of the file to be read.
     */
    void read_hb_file(const char * input_file)
    {
        int M;
        int N;
        int nonzeros;
        int * colptr;
        int * rowind;
        double * val;


        HarwellBoeing_IO h(input_file);
        h.read(M, N, nonzeros, colptr, rowind, val);

        this->resize(M, N);

        size_type elem_count = 0;

        for (int j = 0; j < N + 1; ++j) {
            for (int m = colptr[j]; m < colptr[j + 1]; ++m) {
                this->contents[rowind[elem_count] - 1][j] = val[elem_count];
                ++elem_count;
            }
        }

    }

/**
   * Write data in Harwell-Boeing format method.
   * Opens the file specified and writes the matrix's data in it.
   *
   * \param input_file Name of the file to be read.
 */
    void write_hb_file(const char * input_file)
    {
#ifndef HAVE_GMM
        std::stringstream message;
        message
        << "---------------------------------------------\n"
        <<
        "We are very very sorry but this feature (write_hb_file) requires at the moment the use of \nthe gmm++ library.\n"
        <<
        "If you really feel that it is important to have the Matrix Market file reading \nas a built in feature, you can post it into our web page:\n"
        << "w3.mecanica.upm.es/lmx\n"
        << "Thank you very much for your interest and help in setting our priorities in the \ndevelopment schedule.\n"
        << "---------------------------------------------\n"
        << "gmm++ not defined.\nYou must set \"#define HAVE_GMM\" in your file in order to use this library." << endl;
        LMX_THROW(failure_error, message.str());
#else
        gmm::dense_matrix<T> gmm_dense_matrix( rows , cols );
        gmm::csc_matrix<T> cscmat( rows , cols );

        for (size_type i = 0; i<rows; ++i){
          for (size_type j = 0; j<cols; ++j){
            gmm_dense_matrix(i,j) = this->contents[i][j];
          }
        }

        gmm::copy( gmm_dense_matrix, cscmat );
        gmm::Harwell_Boeing_save(input_file, cscmat);

#endif
    }

    /** Traspose method.
      * Swaps elements with respect to the diagonal: A(i,j) = A(j,i) */
    void trn()
    {
        std::vector<std::vector<T> > temp;
        size_type old_rows;
        size_type i, j;

        temp.resize(cols);
        for (i = 0; i < cols; ++i) {
            temp[i].resize(rows);
        }
        for (i = 0; i < rows; ++i) {
            for (j = 0; j < cols; ++j) {
                temp[j][i] = contents[i][j];
            }
        }
        this->resize(cols, rows);
        for (i = 0; i < cols; ++i) {
            for (j = 0; j < rows; ++j) {
                contents[i][j] = temp[i][j];
            }
        }
        old_rows = rows;
        rows = cols;
        cols = old_rows;
    }

    /** Clean below method.
      * Makes equal to zero every element below given factor.
      * \param factor Reference value for cleaning. */
    void cleanBelow(const double factor)
    {
        for (size_type i = 0; i < rows; ++i) {
            for (size_type j = 0; j < cols; ++j) {
                if (std::abs(contents[i][j]) < std::abs(static_cast<T>(factor))) contents[i][j] = static_cast<T>(0);
            }
        }
    }

    /**
     * Clear method.
     * Wipes all data.
     */
    void clear()
    {
        contents.clear();
        resize(rows, cols);
    }

    //begin JCGO 18/03/09
    /**
     * Reset method.
     */
    void reset()
    {
        for (size_type i = 0; i < rows; ++i) {
            for (size_type j = 0; j < cols; ++j) {
                contents[i][j] = static_cast<T>(0);
            }
        }
    }
    //end JCGO

    /** Data pointer method
     * Gives the direction in memory of (pointer to) the object.
     * @return A pointer to the matrix's contents (Type_stdmatrix).
     */
    std::vector<std::vector<T> > * data_pointer() { return &(this->contents); }

    bool exists(size_type, size_type) { return 1; }

}; // class Type_stdmatrix definitions.

////////////////////////////////////////////////// Implementing the methods defined previously:

/** Standard constructor.
  * Creates a new object with parameter contents resized to (rows, columns) dimension.
  * \param rows_in Rows of dense matrix.
  * \param columns_in Columns of dense matrix. */
template<typename T>
Type_stdmatrix<T>::Type_stdmatrix(size_type rows_in, size_type columns_in) :
        Data_mat<T>(), rows(rows_in), cols(columns_in) { resize(rows_in, columns_in); }

/** Resize method.
  * Changes the size of the contents parameter.
  * \param mrows New value for rows of dense matrix.
  * \param ncolumns New value for columns of dense matrix. */
template<typename T>
void Type_stdmatrix<T>::resize(size_type mrows, size_type ncolumns)
{
    if (mrows != this->rows || ncolumns != this->cols) {

        for (size_type i = 0; i < this->rows; ++i) {
            contents[i].resize(ncolumns);
        }

        contents.resize(mrows);
        for (size_type i = rows; i < mrows; ++i) {
            contents[i] = std::vector<T>(ncolumns);
        }

        /*
        for (size_type i = this->rows; i < mrows; ++i) {
            contents[i].resize(ncolumns);
        }
        */

        rows = mrows;
        cols = ncolumns;
    }
}


}; // namespace lmx


#endif
