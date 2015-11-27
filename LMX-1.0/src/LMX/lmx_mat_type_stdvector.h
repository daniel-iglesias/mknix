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
#ifndef LMXTYPE_STDVECTOR_H
#define LMXTYPE_STDVECTOR_H

#include<fstream>
#include<iostream>
#include<cstdlib>
#include<vector>
#include<string.h>
#include"lmx_mat_data_vec.h"

//////////////////////////////////////////// Doxygen file documentation entry:
    /*!
      \file lmx_mat_type_stdvector.h

      \brief This file contains both the declaration and implementation for type_stdVector (dense vector) class member functions.

      \author Daniel Iglesias 

    */
//////////////////////////////////////////// Doxygen file documentation (end)


namespace lmx {

    /**
    \class Type_stdVector 
    \brief Template class Type_stdVector
    
    This class implements the methods defined in virtual class data_mat so the Vector type "std::vector" can be used in lmx. For details about the caracteristics of this matrix type, see the STL library manual and code.
    
    @param contents Corresponds to a std::vector and it's the base of the methods implemented for this class.
    
    @author Daniel Iglesias .
    */
template <typename T> class Type_stdVector : public Data_vec<T>
{
private:
  /** Matrix data contents */
    std::vector<T> contents;
    
public:
    /// Empty constructor.
    Type_stdVector(){}

    Type_stdVector(size_type, size_type);

    Type_stdVector(char*);

    /// Destructor.
    ~Type_stdVector(){}

    /** Resize method.
     * Changes the size of the contents parameter.
     * \param mrows New value for vector's size.
     * \param ncolumns Not used.
     */
    void resize(size_type mrows, size_type ncolumns)
    { contents.resize(mrows); }

    /** Read element method.
     * Implements a method for reading vector's data.
     * \param mrows Row position in dense matrix.
     * \param ncolumns Column position in dense matrix.
     * \return Value of the element in the position given by the parameter.
     */
    const T& readElement(const size_type& mrows, const size_type& ncolumns) const
   { return contents[mrows]; }

    /** Write element method.
      * Implements a method for writing data on the dense matrix.
      * \param mrows Row position in dense matrix.
      * \param ncolumns Column position in dense matrix.
      * \param value Numerical type value. */
  void writeElement(T value, const size_type& mrows, const size_type& ncolumns)
   { contents[mrows] = value; }
  
    /** Add element method.
      * Implements a method for adding data on the Harwell-Boeing matrix.
      * Copy-pasted from writeElement.
      * \param mrows Row position in dense matrix.
      * \param ncolumns Column position in dense matrix.
      * \param value Numerical type value. */
  inline void addElement(T value, const size_type& mrows, const size_type& ncolumns)
   { contents[mrows] += value; }
  
   /** Method for knowing the number of data rows. 
    * \returns Number of rows.
    */
   size_type getRows() const
   { return contents.size(); }

   /** Method for knowing the number of data columns. 
    * \returns Number of columns.
    */
   size_type getCols() const
   { return 1; }

    /** Copy method.
      * Equals the data in the object's contents to those given by the input vector's parameter.
      * \param vector_in pointer to an object that belongs to a class derived from Data. */
  void equals(const Data<T>* vector_in)
  {
    contents = static_cast<const Type_stdVector*>(vector_in)->contents;

   }

    /** Add method.
      * Adds the the input vector parameter's elements to the object's contents.
      * Necessary for overloading the "+=" operator.
      * \param vector_in_1 pointer to an object that belongs to a class derived from Data. */
   void add(const Data<T>* vector_in_1)
   {
     // This will be changed by a generalized STL algorithm... when I find the best way.
#ifdef HAVE_GMM
     gmm::add(static_cast<const Type_stdVector*>(vector_in_1)->contents, contents);
#else
     for(size_type i=0; i<this->contents.size(); ++i)
       this->contents[i] += static_cast<const Type_stdVector*>(vector_in_1)->contents[i];
#endif
   }

    /** Substract method.
      * Substracts the the input vector parameter's elements to the object's contents.
      * Necessary for overloading the "-=" operator.
      * \param vector_in_1 pointer to an object that belongs to a class derived from Data. */
   void substract(const Data<T>* vector_in_1)
   {
    // This will be changed by a generalized STL algorithm... when I find the best way.
     for (unsigned int i=0; i < this->contents.size(); ++i)
       this->contents[i] -= static_cast<const Type_stdVector*>(vector_in_1)->contents[i];
   }

    /** Multiply scalar method.
      * Multiplies the object's vector (contents) with a scalar.
      * Necessary for overloading the "*" operator.
      * \param scalar A scalar factor of template's class. */
  void multiplyScalar(const T& scalar)
  {
    // This will be changed by a generalized STL algorithm... when I find the best way.
    for (unsigned int i=0; i<contents.size(); ++i) contents[i] *= scalar;
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


   void cleanBelow(const double factor)
   {
     // This will be changed by a generalized STL algorithm... when I find the best way.
#ifdef HAVE_GMM
      gmm::clean(contents, factor);
#else
     for(size_type i=0; i<this->contents.size(); ++i){
        if (this->contents[i] < static_cast<T>(factor) )
          this->contents[i] = static_cast<T>(0);
     }
#endif

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

    for(unsigned int i=0; i<rows ; ++i) {
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
    
    file_output << contents.size() << endl;

    for(unsigned int i=0; i<contents.size() ; ++i) {
      file_output << contents[i] << endl;
    }
  }

  /**
   * Clear method.
   * Wipes all data.
   */
  void clear()
  { size_type rows = contents.size();
    contents.clear();
    contents.resize(rows);
  }

	//begin JCGO 18/03/09
  /**
   * Reset method.
   */
  void reset()
  {
  for(size_type i=0; i<this->contents.size(); ++i)
          this->contents[i] = static_cast<T>(0);
  }
  //end JCGO

  /** Data pointer method
   * Gives the direction in memory of (pointer to) the object.
   * @return A pointer to the vector's contents (Type_stdVector).
   */
  std::vector<T>* data_pointer()
  {
    return &(this->contents);
  }

  friend void mat_vec_mult<>( const Type_csc<T>*,
                              const Type_stdVector<T>*,
                                    Type_stdVector<T>*);

}; // End class Type_stdVector definitions.

  ////////////////////////////////////////////////
  // Implementing the methods defined previously:

/**
 * Standard constructor.
 *
 * Creates a new object with parameter contents resized to (columns) dimension.
 * @param rows Rows of dense stdVector.
 * @param columns Not used.
 * @return
 */
template <typename T>
    Type_stdVector<T>::Type_stdVector(size_type rows, size_type columns) : Data_vec<T>()
{
   contents.resize(rows);
}


/**
 * Constructor from a data file.
 * Creates a new object reading data from a file.
 * \param input_file Name of file to read from.
 */
template <typename T> Type_stdVector<T>::Type_stdVector(char* input_file) : Data_vec<T>()
{
  readDataFile(input_file);
}

}; // namespace lmx


#endif
