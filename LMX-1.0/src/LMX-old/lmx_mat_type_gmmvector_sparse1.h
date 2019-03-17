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
#ifndef LMXTYPE_GMMVECTOR_SPARSE_H
#define LMXTYPE_GMMVECTOR_SPARSE_H

// #include<fstream>
#include<iostream>
#include<cstdlib>
#include"lmx_mat_data_vec.h"

//////////////////////////////////////////// Doxygen file documentation entry:
    /*!
      \file lmx_mat_type_gmmvector_sparse1.h
      
      \brief This file contains both the declaration and implementation for Type_gmmVector_sparse (sparse vector) class member functions.
      
      \author Daniel Iglesias IbÔøΩez
      
    */
//////////////////////////////////////////// Doxygen file documentation (end)


namespace lmx {

    /**
    \class Type_gmmVector_sparse 
    \brief Template class Type_gmmVector_sparse
    
    This class implements the methods defined in virtual class data_mat so the Vector type "std::vector" can be used in lmx. For details about the caracteristics of this matrix type, see the STL library manual and code.
    
    @param contents Corresponds to a std::vector and it's the base of the methods implemented for this class.
    
    @author Daniel Iglesias .
    */
template <typename T> class Type_gmmVector_sparse : public Data_vec<T>
{
private:
  /** Matrix data contents */
    gmm::wsvector<T> contents;
//     T reference;

public:
    /// Empty constructor.
    Type_gmmVector_sparse(){}

    Type_gmmVector_sparse(size_type, size_type);

    Type_gmmVector_sparse(char* input_file);

    /// Destructor.
    ~Type_gmmVector_sparse(){}

  void resize(size_type, size_type);

    /** Read element method.
      * Implements a method for reading data of the sparse matrix.
      * \param mrows Row position in sparse matrix.
      * \param ncolumns Column position in sparse matrix.
      * \return Value of the element in the position given by the parameters. */
    const T& readElement(const size_type& mrows, const size_type& ncolumns) const
  { static T reference = 0;
    return reference = contents[mrows];
  }

    /** Write element method.
      * Implements a method for writing data on the sparse matrix.
      * \param mrows Row position in sparse matrix.
      * \param ncolumns Column position in sparse matrix.
      * \param value Numerical type value. */
  void writeElement(T value, const size_type& mrows, const size_type& ncolumns)
   { contents[mrows] = value; }

    /** Add element method.
      * Implements a method for adding data on the Harwell-Boeing matrix.
      * Copy-pasted from writeElement.
      * \param mrows Row position in sparse matrix.
      * \param ncolumns Column position in sparse matrix.
      * \param value Numerical type value. */
  inline void addElement(T value, const size_type& mrows, const size_type& ncolumns)
   { contents[mrows] = value; }

   /** Method for knowing the number of data rows. 
    * \returns Number of rows.
    */
   size_type getRows() const
   { return gmm::vect_size(contents); }

   /** Method for knowing the number of data columns. 
    * \returns Number of columns.
    */
   size_type getCols() const
   { return 1; }

    /** Copy method.
      * Equals the data in the object's contents to those given by the input vector's parameter.
      * \param vector_in pointer to an object that belongs to a class derived from Data. */
  void equals(const Data<T>* vector_in)
   { contents = static_cast<const Type_gmmVector_sparse*>(vector_in)->contents; }
  
    /** Add method.
      * Adds the the input vector parameter's elements to the object's contents.
      * Necessary for overloading the "+=" operator.
      * \param vector_in_1 pointer to an object that belongs to a class derived from Data. */
  void add(const Data<T>* vector_in_1)
   { gmm::add(static_cast<const Type_gmmVector_sparse*>(vector_in_1)->contents, contents); }

    /** Substract method.
      * Substracts the the input vector parameter's elements to the object's contents.
      * Necessary for overloading the "-=" operator.
      * \param vector_in_1 pointer to an object that belongs to a class derived from Data. */
  void substract(const Data<T>* vector_in_1)
   { gmm::scale(contents, -1.);
     gmm::add(static_cast<const Type_gmmVector_sparse*>(vector_in_1)->contents, contents);
     gmm::scale(contents, -1.);
   }

    /** Multiply scalar method.
      * Multiplies the object's vector (contents) with a scalar.
      * Necessary for overloading the "*" operator.
      * \param scalar A scalar factor of template's class. */
   void multiplyScalar(const T& scalar)
   { gmm::scale(contents, scalar);
   }
  
  /** Method multiplying element-by-element of two arrays. One would be the object's contents and the other the parameter's contents.
    * Necessary for implementing  Vector to Vector multElements.
    * \param vector_in pointer to an object that belongs to a class derived from Data.
   */
   void multiplyElements(const Data<T>* vector_in)
   {
     for (size_type i=0; i<contents.size(); ++i)
       contents[i] *= vector_in->readElement(i,1);
   }


//    /** Traspose method.
//      * Swaps elements with respect to the diagonal: A(i,j) = A(j,i) */
//    void trn()
//    { // Must throw an error? Two options:
//      //  1. Imposible to transpose a vector.
//      //  2. Treat it as a Nx1 matrix and make a 1xN matrix.
//      // Should agree with dimension check in matrix * vector operator (vector.h).
//    }

    /** Clean below method.
      * Makes equal to zero every element below given factor.
      * \param factor Reference value for cleaning. */
   void cleanBelow(const double factor)
   { gmm::clean(contents, factor);
   }

  /** Read data method.
   * Opens the file specified and reads the vector's data in it.
   * 
   * @param input_file Name of the file to be read.
   */
  void readDataFile(const char* input_file)
  { 
    std::ifstream file_input(input_file); /* New stream for reading data from file. */
    
    char file_type[100]; /* Array for storing file type. */
    char crap[100]; /* Array for temporary storing data not to be used. */
    T temp; /* Class T variable for temporary storing vector's values. */
    size_type rows; /* Variable for storing number of rows of file's vector. */

    file_input >> file_type;

    if (!strcmp(file_type,"%%MatrixMarket")) {
      file_input.getline(crap,99);
      file_input >> rows;
      file_input.getline(crap,99);
    }
    else (rows = atoi(file_type));

    this->resize(rows, 1);

    for(int i=0; i<rows ; ++i) {
      file_input >> temp;
      contents[i] = temp;
    }
  }

  /** Write data method.
   * Opens the file specified and writes the vector's data in it.
   *
   * @param output_file Name of the file to be read.
   */
  void writeDataFile(const char* output_file)
  {
    std::ofstream file_output(output_file);
    size_type rows = gmm::vect_size(contents);
    
    file_output << rows << endl;

    for(int i=0; i<rows ; ++i) {
      file_output << contents[i] << endl;
    }
  }

  /**
   * Clear method.
   * Wipes all data.
   */
  void clear()
  {
    gmm::clear(contents); /**< Vector data contents. */
  }

	//begin JCGO 18/03/09
  /**
   * Reset method.
   */
  void reset()
  {
  	gmm::clean(contents,0.0);
  }
  //end JCGO

  /** Data pointer method
   * Gives the direction in memory of (pointer to) the object.
   * @return A pointer to the vector's contents (Type_gmmVector_sparse).
   */
  gmm::wsvector<T>* data_pointer()
  {return &(this->contents);}


}; // class Type_gmmVector_sparse definitions.

  ////////////////////////////////////////////////
  // Implementing the methods defined previously:

    /** Standard constructor.
     * Creates a new object with parameter contents resized to (columns) dimension.
     * \param rows New size of dense stdVector.
     * \param columns Not used. */
  template <typename T> Type_gmmVector_sparse<T>::Type_gmmVector_sparse(size_type rows, size_type columns) : Data_vec<T>()
     { gmm::resize(contents, rows); }

    /** Constructor from a data file.
      * Creates a new object reading data from a file.
      * \param input_file Name of file to read from. */
  template <typename T> Type_gmmVector_sparse<T>::Type_gmmVector_sparse(char* input_file) : Data_vec<T>()
     { // read_file;
       // write data;
     }

    /** Resize method.
     * Changes the size of the contents parameter.
     * \param mrows New value for size of Vector.
     * \param ncolumns Not used. */
     template <typename T> void Type_gmmVector_sparse<T>::resize(size_type mrows, size_type ncolumns)
   {
     gmm::resize(contents, mrows); 
   }


}; // namespace lmx


#endif
