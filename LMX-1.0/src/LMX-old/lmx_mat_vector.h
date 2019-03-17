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

#ifndef LMXVECTOR_H
#define LMXVECTOR_H

#ifdef HAVE_GMM
#include"lmx_mat_type_gmmvector_sparse1.h"
#endif

// #include"lmx_base_selector.h"
#include"lmx_mat_elem_ref.h"
#include"lmx_mat_type_stdvector.h"

//////////////////////////////////////////// Doxygen file documentation entry:
    /*!
      \file lmx_mat_vector.h

      \brief This file contains both the declaration and implementation for Vector class member and friend functions.

      \author Daniel Iglesias

    */
//////////////////////////////////////////// Doxygen file documentation (end)


namespace lmx {

class Selector;
template <typename T> class Matrix;
template <typename T> class DenseMatrix;
template <typename T> class LinearSystem;

int setVectorType(int);
int getVectorType();

    /**
    \class Vector 
    \brief Template class Vector

    This class permits the creation of vector objects. The Vector class, which derives from Matrix class, has one parameter (nrows) that stores the dimension of the vector's container. The data is stored in an atribute (*type_vector) that points to some class which derives from the Data_vec class.

    \param elements The size of the Data Container (Data_vec).
    \param reference An Elem_ref object for r/w data access.
    \param *type_matrix The pointer to the vector data container.

    \author Daniel Iglesias.
    */
template <typename T> class Vector{
private:
  size_type elements;    /**< Number of rows in vector object. */
  Elem_ref<T>* reference; /**< Reference pointer to an element in type_matrix. */
  Data_vec<T>* type_vector; /**< Pointer to the container type. */
  static const size_type zero;

public:
  friend class LinearSystem<T>;

public:
  Vector();

  explicit Vector(size_type);

  Vector(const Vector&);

  ~Vector();

  /**
   * Read number of rows.
   * \return Dimension of vector.
   */
  inline size_type size() const { return this->elements; }

  void load(char*);

  void save(char*);

  void fillIdentity( T factor = static_cast<T>(1) );

  void fillRandom( T factor = static_cast<T>(1) );

  inline Elem_ref<T>& operator () (size_type);

  inline Vector& operator = (const Matrix<T>&);

  inline Vector& operator = (const Vector<T>&);

  template <typename C> inline Vector<T>& operator = (const Vector<C>&);

  inline Vector& operator += (const Vector&);

  inline Vector& operator -= (const Vector&);

  inline Vector& operator *= (const T&);

//   inline Vector operator * (const Matrix<T>&) const;

//   inline Vector operator * (const DenseMatrix<T>&) const;

  inline Vector& add(const Vector<T>&, const Vector<T>&);

  inline Vector& subs(const Vector<T>&, const Vector<T>&);

  inline Vector& mult(const Matrix<T>&, const Vector<T>&);

  inline Vector& mult(const DenseMatrix<T>&, const Vector<T>&);

  inline Vector& mult(const T&);

  inline Vector& mult(const Vector<T>&, const T&);

  inline Vector& mult(const T&, const Vector<T>&);

  inline Vector& mult(const Vector<T>&, const Vector<T>&);
  
  inline Vector& multElements(const Vector<T>&);

  inline Vector& multElements(const Vector<T>&, const Vector<T>&);

  inline T norm1 () const;

  inline T norm2 () const;

 /** Function for reading a Vector object's element from it's known position.
   * \param m vec's position (row) for lookup.
   * \return The element.
    *  */
  T readElement(size_type m) const
  { if (m >= elements){
      std::stringstream message;
      message << "Row index exceeds vector dimensions. Trying to access element " << m << " of vecor with dimension " << elements << "." << endl;
      LMX_THROW(dimension_error, message.str());
    }
    return type_vector->readElement(m,zero);
  }


  /** Function writes the value in specified position of Vector object.
   *  */
  inline
      void writeElement(T theValue, size_type m) const
  { this->type_vector->writeElement(theValue, m, zero); }


  /** Function writes the value in specified position of Vector object.
   *  */
  inline
      void addElement(T theValue, size_type m) const
  { this->type_vector->addElement(theValue, m, zero); }


  /** Cleans all numbers below given factor.
   *  */
  void clean(double factor)
  {
    this->type_vector->cleanBelow(factor);
  }

  //needs documentation
  /** Clears the object's contents.
   *  */
  inline
      void clear()
  { this->type_vector->clear(); }

	//begin JCGO 18/03/09
  //needs documentation
  /** Reset the object's contents.
   *  */
  inline
      void reset()
  { this->type_vector->reset(); }
	//end JCGO
	
  /** Resize the Vector with given size parameter.
   *  \param i New size.
   *  */
  void resize(size_type i)
  {
    this->elements = i;
    this->type_vector->resize(this->elements, 1);
  }

  /** Resize the Vector to the size of another Vector.
   *  \param vec Constant reference to a Vector.
   *  */
  inline
      void resize( const Vector& vec )
  {
    this->elements = vec.size();
    this->type_vector->resize(this->elements, 1);
  }

  /** Operator for adding two Vector objects.
   * Performs typicar vector operation: A + B = C, where A, B & C are vectors. A is Vector *this.
   * \param B Vector reference for second LHS term.
   * \return The numeric result of the operation: C = A + B.
      *  */
  Vector operator + (const Vector<T>& B) const
  {
    if ( elements != B.size() )
      LMX_THROW(dimension_error, "Vectors dimensions mismatch");
    Vector<T> sum = *this;
    sum += B;
    return sum;
  }

  /** Operator for substracting two Vector objects.
   * Performs typicar vector operation: A - B = C, where A, B & C are vectors. A is Vector *this.
   * \param B Vector reference for second LHS term.
   * \return The numeric result of the operation: C = A - B.
   *  */
  Vector operator - (const Vector<T>& B) const
  {
    if ( elements != B.size() )
      LMX_THROW(dimension_error, "Vectors dimensions mismatch");
    Vector<T> sum = *this;
    sum -= B;
    return sum;
  }

  /** Operator for scalar product of two Vector objects.
   * Performs typical vector operation: A * B = c, where A & B are vectors and c is a scalar.
   * \param B Vector reference for second LHS term.
   * \return The numeric result of the operation (type: \a double ).
      *  */
  T operator * (const Vector& B) const
  {
    if ( elements != B.size() )
      LMX_THROW(dimension_error, "Vectors dimensions mismatch");

    T scalar_product = 0;

    for (size_type i=0; i<elements; ++i)
      scalar_product += type_vector->readElement(i,zero) * B.readElement(i);
    return scalar_product;
  }

  /** Overload operator for multiplying a scalar and a Vector object.
   * \param a scalar reference.
   * \return A new Vector (mult = (*this) * a).
   *  */
  Vector operator * (const T& a) const
  {
    Vector<T> mult(*this);
    mult.type_vector->multiplyScalar(a);
    return mult;
  }



/////////////////////////////// Friend functions:


  /** Function to typeset a vector to a file stream.
   *  \param os File stream to write to.
   *  \param vec_name Name that will be given to vector in typesetting. It's not the vector variable name.
   *  \param vec The Vector variable to be typeset.
   *  \param prec Maximum numeric precision (max. number of digits to be used in output).
   *  */
  friend void latexPrint(std::ofstream& os, char* vec_name, Vector& vec, int prec)
  { os.precision(prec);
    os << "\\[" << vec_name << " = " << std::endl;
    os << "\\left\\lbrace  \\begin{array}{";
    for(int j=0; j<1; ++j) { os << "c"; }
    os << "}" << std::endl;
    for(int i=0; i<vec.size(); ++i){
      for(int j=0; j<1; ++j){
        if(j != 1-1) { os << vec(i) << " & "; }
        else if (i != vec.size()-1) { os << vec(i) << " \\\\ "; }
        else { os << vec(i) << " \\end{array} \\right\\rbrace \\]" << endl; }
      }
    }
  }

};

}; // namespace lmx

/////////////////////////////// Implementation of the methods defined previously

namespace lmx {

template <typename T>
   const size_type Vector<T>::zero = 0;

/**
 * Empty constructor.
 */
template <typename T>
    Vector<T>::Vector()
{
  elements = 0;

  switch (getVectorType()) {
    case 0 :
      type_vector = new Type_stdVector< T >;
      break;

#ifdef HAVE_GMM
    case 1 :
      type_vector = new Type_gmmVector_sparse< T >;
    break;
#else
    {
        std::stringstream message;
        message << "gmm++ not defined.\nYou must set \"#define HAVE_GMM\" in your file in order to use this library." << endl;
        LMX_THROW(failure_error, message.str() );
    }
#endif
  }

  reference = new Elem_ref<T>(this->type_vector);

}

  /**
   * Standard constructor.
   * \param rows Number of Vector's rows.
   */
template <typename T>
    Vector<T>::Vector(size_type rows)
{
//   if ( rows < 0 ) LMX_THROW(dimension_error, "Can't have a negative-sized vector");

  elements = rows;

  switch (getVectorType()) {
    case 0 :
      type_vector = new Type_stdVector< T >;
      break;

#ifdef HAVE_GMM
    case 1 :
      type_vector = new Type_gmmVector_sparse< T >;
      break;
#else
    {
        std::stringstream message;
        message << "gmm++ not defined.\nYou must set \"#define HAVE_GMM\" in your file in order to use this library." << endl;
        LMX_THROW(failure_error, message.str() );
    }
#endif
  }

  type_vector->resize(elements, 1);
  reference = new Elem_ref<T>(this->type_vector);

}


/**
 * Copy constructor.
 * \param A Vector to copy from.
 */
template <typename T>
    Vector<T>::Vector(const Vector& A) : elements(A.elements)
{
  switch (getVectorType()) {
    case 0 :
      type_vector = new Type_stdVector< T >;
      break;

#ifdef HAVE_GMM
    case 1 :
      type_vector = new Type_gmmVector_sparse< T >;
    break;
#else
    {
        std::stringstream message;
        message << "gmm++ not defined.\nYou must set \"#define HAVE_GMM\" in your file in order to use this library." << endl;
        LMX_THROW(failure_error, message.str() );
    }
#endif
  }
  this->type_vector->resize(this->elements, 1);
  this->type_vector->equals(A.type_vector);
  this->reference = new Elem_ref<T>(this->type_vector);
}


/**
 * Destructor.
 */
template <typename T>
    Vector<T>::~Vector()
{
      delete this->reference;
      this->reference = 0;

      delete this->type_vector;
      this->type_vector = 0;

}

/** Method for reading a Vector from a file.
 *  \param input_file Name of the file to read.
 *  */
template <typename T>
    void Vector<T>::load(char* input_file)
{
  this->type_vector->readDataFile(input_file);
  this->elements = this->type_vector->getRows();
}

/** Method for writing a Vector to a file.
 *  \param output_file Name of the file to write.
 *  */
template <typename T>
    void Vector<T>::save(char* output_file)
{
  this->type_vector->writeDataFile(output_file);
}


/**
 * Cleans vector and sets all terms to specified value.
 * \param factor Value of terms. Default is unit value.
 */
template <typename T>
    void Vector<T>::fillIdentity( T factor )
{
  for (size_type i=0; i<this->elements; ++i)
    this->type_vector->writeElement(static_cast<T>(factor),i,zero);
}

/**
 * Function for filling a vector with random numbers
 * \param factor Scales the random numbers (default value is unity).
 */
template <typename T>
    void Vector<T>::fillRandom( T factor )
{
  for (size_type i=0; i<this->elements; ++i)
      this->type_vector->writeElement( factor * static_cast<T>( std::rand() ) / static_cast<T>(RAND_MAX),i,1);
}


/**
 * Overload of element extraction method.
 * \param row Position in Vector for the element to extract.
 * \return element referenced.
 */
template <typename T>
    Elem_ref<T>& Vector<T>::operator () (size_type row)
{
  if ( row >= this->elements ){
    std::stringstream message;
    message << "Row index exceeds vector dimensions. \n"
        << "Index: (" << row << ")" << endl
        << "Vector dimension: (" << this->elements << ")" << endl;
    LMX_THROW(dimension_error, message.str() );
  }

  reference->write_pos(row,zero);
  return *reference;
}

/** Overload of equals operator.
 * Equals every element between the Vector and
 * the first column of a Matrix object.
 * \param A Matrix to equal to.
 * \return A reference to the Vector object.
 */
template <typename T>
     Vector<T>& Vector<T>::operator = (const Matrix<T>& A)
{
  this->elements = rows(A);
  for (size_type i=0; i<this->elements; ++i)
    this->operator ( )(i,zero) = readElement(A,i,zero);
  return *this;
}

/** Overload of equals operator.
 * Equals every element between two Vector objects of the same type.
 * \param v Vector to equal to.
 * \return A reference to the Vector object.
 */
template <typename T>
    Vector<T>& Vector<T>::operator = (const Vector<T>& v)
{
  if ( this->elements != v.size() ){
    std::stringstream message;
    message << "Vectors dimensions mismatch. \n"
        << "LHS vector dimension: (" << this->elements << ")" << endl
        << "RHS vector dimension: (" << v.size() << ")" << endl;
        LMX_THROW(dimension_error, message.str() );
  }
  type_vector->equals(v.type_vector);
  return *this;
}

/** Overload of equals operator.
 * Equals every element between two Vector objects of different types.
 * \param v Vector to equal to.
 * \return A reference to the Vector object.
 */
template <typename T>
    template <typename C>
    Vector<T>& Vector<T>::operator = (const Vector<C>& v)
{
  if ( this->elements != v.size() ){
    std::stringstream message;
    message << "Vectors dimensions mismatch. \n"
        << "LHS vector dimension: (" << this->elements << ")" << endl
        << "RHS vector dimension: (" << v.size() << ")" << endl;
    LMX_THROW(dimension_error, message.str() );
  }

  for (int i=0; i<(this->elements) ; ++i){
      this->operator ( )(i) = v.readElement(i);
  }
}

/** Overloaded operator for assigning the addition of elements
 *  between two Vector objects.
 * \param b Vector to equal to.
 * \return A reference to the Vector object.
 */
template <typename T>
    Vector<T>& Vector<T>::operator += (const Vector& b)
{
  if ( this->elements != b.size() ){
    std::stringstream message;
    message << "Vectors dimensions mismatch. \n"
        << "LHS vector dimension: (" << this->elements << ")" << endl
        << "RHS vector dimension: (" << b.size() << ")" << endl;
        LMX_THROW(dimension_error, message.str() );
  }
  type_vector->add(b.type_vector);
  return *this;
}

/** Overloaded operator for assigning the substraction of elements
 * between two Vector objects.
 * \param b Vector to equal to.
 * \return A reference to the Vector object.
 */
template <typename T>
    Vector<T>& Vector<T>::operator -= (const Vector& b)
{
  if ( this->elements != b.size() ){
    std::stringstream message;
    message << "Vectors dimensions mismatch. \n"
        << "LHS vector dimension: (" << this->elements << ")" << endl
        << "RHS vector dimension: (" << b.size() << ")" << endl;
        LMX_THROW(dimension_error, message.str() );
  }
  type_vector->substract(b.type_vector);
  return *this;
}

/**
 * Scaling operator.
 *
 * Scales each element in Vector with specified scalar factor.
 * \param scalar Scaling factor.
 * \return Reference to result.
 */
template <typename T> inline
    Vector<T>& Vector<T>::operator *= (const T& scalar)
{
  type_vector->multiplyScalar(scalar);
  return *this;
}

/**
 * (Cancelled function) Operator for multiplying Vector and Matrix objects.
 *
 * Performs typical vector operation: b * A = c, where b & c are vectors and A is a Matrix.
 * \param A Matrix reference for second LHS term.
 * \param b Vector reference for first LHS term.
 * \return A new Vector (c = b^T * A).
 */
// template <typename T> inline
//     Vector<T> Vector<T>::operator * (const Matrix<T>& A) const
// {
//   if ( A.cols() != elements ){
//     std::stringstream message;
//     message << "Can't multiply vector^T * matrix: Dimensions mismatch. \n"
//         << "Matrix dimension: (" << A.rows() << "," << A.cols() << ")" << endl
//         << "Vector dimension: (" << elements << ")" << endl;
//     LMX_THROW(dimension_error, message.str() );
//   }
//   Vector<T> c(A.rows());
//   c.mult(A,*this);
//   return c;
// }

/**
 * (Cancelled function) Operator for multiplying Vector and DenseMatrix objects.
 *
 * Performs typical vector operation: b * A = c, where b & c are vectors and A is a DenseMatrix.
 * \param A DenseMatrix reference for second LHS term.
 * \param b Vector reference for first LHS term.
 * \return A new Vector (c = b^T * A).
 */
// template <typename T> inline
//     Vector<T> Vector<T>::operator * (const DenseMatrix<T>& A) const
// {
//   if ( A.cols() != elements ){
//     std::stringstream message;
//     message << "Can't multiply Vector^T * DenseMatrix: Dimensions mismatch. \n"
//         << "Matrix dimension: (" << A.rows() << "," << A.cols() << ")" << endl
//         << "Vector dimension: (" << elements << ")" << endl;
//     LMX_THROW(dimension_error, message.str() );
//   }
//   Vector<T> c(A.rows());
//   c.mult(A,*this);
//   return c;
// }


/**
 * Element adding function.
 *
 * Adds two Vectors' elements and saves result in object.
 * 
 * \param A A vector.
 * \param B A vector.
 * \return Reference to result.
 */
template <typename T> inline
    Vector<T>& Vector<T>::add(const Vector<T>& A, const Vector<T>& B)
{
  if ( this->elements != A.size() || this->elements != B.size() ){
    std::stringstream message;
    message << "Vectors dimensions mismatch. \n"
        << "LHS vector dimension: (" << this->elements << ")" << endl
        << "First RHS vector dimension: (" << A.size() << ")" << endl
        << "Second RHS vector dimension: (" << B.size() << ")" << endl;
        LMX_THROW(dimension_error, message.str() );
  }
  this->type_vector->equals(A.type_vector);
  this->type_vector->add(B.type_vector);
  return *this;
}


/**
 * Element substracting function. Substract two Vectors' elements and saves result in object.
 * \param A A vector.
 * \param B A vector.
 * \return Reference to result.
 */
template <typename T> inline
    Vector<T>& Vector<T>::subs(const Vector<T>& A, const Vector<T>& B)
{
  if ( this->elements != A.size() || this->elements != B.size() ){
    std::stringstream message;
    message << "Vectors dimensions mismatch. \n"
        << "LHS vector dimension: (" << this->elements << ")" << endl
        << "First RHS vector dimension: (" << A.size() << ")" << endl
        << "Second RHS vector dimension: (" << B.size() << ")" << endl;
        LMX_THROW(dimension_error, message.str() );
  }
  this->type_vector->equals(A.type_vector);
  this->type_vector->substract(B.type_vector);
  return *this;
}

/**
 * Function for multiplying Matrix and Vector objects.
 * Performs typical vector operation: A * b = c, where b & c are vectors and A is a Matrix.
 *
 * c would be the object witch its invoked from.
 *
 * This function is specialized for some Matrix and Vector types.
 *
 * \param A Matrix reference for first LHS term.
 * \param b Const Vector reference for second LHS term.
 */
template <typename T>
    inline
    Vector<T>& Vector<T>::mult(const Matrix<T>& A, const Vector<T>& b)
{
  if ( A.cols() != b.size() ){
    std::stringstream message;
    message << "Can't multiply matrix * vector: Dimensions mismatch. \n"
        << "Matrix dimension: (" << A.rows() << "," << A.cols() << ")" << endl
        << "Vector dimension: (" << b.size() << ")" << endl;
    LMX_THROW(dimension_error, message.str() );
  }

  if( elements != A.rows() ) resize( A.rows() );

  if (getMatrixType()==1 && getVectorType()==0) {
    mat_vec_mult<T>( static_cast<const Type_csc<T>*>(A.type_matrix),
                     static_cast<const Type_stdVector<T>*>(b.type_vector),
                     static_cast<Type_stdVector<T>*>(this->type_vector) );
  }
//   else if (getMatrixType()==1 && getVectorType()==2) {
//     mat_vec_mult<T>( static_cast<const Type_csc<T>*>(A.type_matrix),
//                      static_cast<const Type_cVector<T>*>(b.type_vector),
//                      static_cast<Type_cVector<T>*>(this->type_vector) );
//   }

#ifdef HAVE_GMM
  else if (getMatrixType()==2 && getVectorType()==0) {
    mat_vec_mult<T>( static_cast<Type_gmm<T>*>(A.type_matrix),
          static_cast<Type_stdVector<T>*>(b.type_vector),
          static_cast<Type_stdVector<T>*>(this->type_vector) );
  }
  else if (getMatrixType()==3 && getVectorType()==0) {
    mat_vec_mult<T>( static_cast<Type_gmm_sparse<T>*>(A.type_matrix),
          static_cast<Type_stdVector<T>*>(b.type_vector),
          static_cast<Type_stdVector<T>*>(this->type_vector) );
  }
#endif
  else {
    for (size_type i=0; i<A.rows(); ++i){
      this->operator ( )(i) = static_cast<T>(0);
      for (size_type j=0; j<A.cols(); ++j)
        this->operator ( )(i) += A.readElement(i,j) * b.readElement(j);
    }
  }

  return *(this);
}

/**
 * Function for multiplying DenseMatrix and Vector objects.
 * Performs typical vector operation: A * b = c, where b & c are vectors and A is a DenseMatrix.
 *
 * c would be the object witch its invoked from.
 *
 * This function is specialized for some Matrix and Vector types.
 *
 * \param A DenseMatrix reference for first LHS term.
 * \param b Const Vector reference for second LHS term.
 */
template <typename T>
    inline
    Vector<T>& Vector<T>::mult(const DenseMatrix<T>& A, const Vector<T>& b)
{
  if ( A.cols() != b.size() ){
    std::stringstream message;
    message << "Can't multiply DenseMatrix * Vector: Dimensions mismatch. \n"
        << "DenseMatrix dimension: (" << A.rows() << "," << A.cols() << ")" << endl
        << "Vector dimension: (" << b.size() << ")" << endl;
    LMX_THROW(dimension_error, message.str() );
  }

  if( elements != A.rows() ) resize( A.rows() );

  for (size_type i=0; i<A.rows(); ++i){
    this->operator ( )(i) = static_cast<T>(0);
    for (size_type j=0; j<A.cols(); ++j)
      this->operator ( )(i) += A.readElement(i,j) * b.readElement(j);
  }

  return *(this);
}

/**
 * Scaling function.
 *
 * Scales each element in Vector with specified scalar factor.
 * \param scalar Scaling factor.
 * \return Reference to result.
 */
template <typename T> inline
    Vector<T>& Vector<T>::mult(const T& scalar)
{
  for(size_type i=0; i<this->elements; ++i) this->type_vector->multiplyScalar(scalar);
}

/**
 * Scaling function. (scalar in RHS)
 *
 * Scales each element in Vector B with specified scalar factor and saves the result in object.
 * \param B Vector to scale.
 * \param scalar Scaling factor.
 * \return Reference to result.
 */
template <typename T> inline
    Vector<T>& Vector<T>::mult(const Vector<T>& B, const T& scalar)
{
  if ( this->elements != B.size() ){
    std::stringstream message;
    message << "Vectors dimensions mismatch. \n"
        << "LHS vector dimension: (" << this->elements << ")" << endl
        << "RHS vector dimension: (" << B.size() << ")" << endl;
        LMX_THROW(dimension_error, message.str() );
  }
  this->type_vector->equals(B.type_vector) ;
  this->type_vector->multiplyScalar(scalar);
  return *(this);
}

/**
 * Scaling function (scalar in LHS)
 *
 * Scales each element in Vector B with specified scalar factor and saves the result in object.
 * \param scalar Scaling factor.
 * \param B Vector to scale.
 * \return Reference to result.
 */
template <typename T> inline
    Vector<T>& Vector<T>::mult(const T& scalar, const Vector<T>& B)
{
  if ( this->elements != B.size() ){
    std::stringstream message;
    message << "Vectors dimensions mismatch. \n"
        << "LHS vector dimension: (" << this->elements << ")" << endl
        << "RHS vector dimension: (" << B.size() << ")" << endl;
        LMX_THROW(dimension_error, message.str() );
  }
  this->type_vector->equals(B.type_vector) ;
  this->type_vector->multiplyScalar(scalar);
  return *(this);
}

/**
 * Cross product
 *
 * Computes the cross product A x B. Not optimized but useful
 * \param A Vector of size m.
 * \param B Vector of size m.
 * \return Reference to result.
 */
template <typename T> inline
    Vector<T>& Vector<T>::mult(const Vector<T>& A, const Vector<T>& B)
{
  if ( A.size() != 3 || B.size() != 3 || this->elements != 3 ){
    std::stringstream message;
    message << "Vectors dimensions mismatch for Cross product ( this->mult(A,B) ). \n"
        << "LHS vector dimension: (" << this->elements << ")" << endl
        << "RHS vector A dimension: (" << A.size() << ")" << endl
        << "RHS vector B dimension: (" << B.size() << ")" << endl;
        LMX_THROW(dimension_error, message.str() );
  }
    this->operator ( )(0) = A.readElement(1) * B.readElement(2) - A.readElement(2) * B.readElement(1);
    this->operator ( )(1) = A.readElement(2) * B.readElement(0) - A.readElement(0) * B.readElement(2);
    this->operator ( )(2) = A.readElement(0) * B.readElement(1) - A.readElement(1) * B.readElement(0);
  return *(this);
}


/**
 * Internal product
 *
 * Element-to-element multiplication between Vector objects.
 * \param B Vector to multiply to.
 * \return Reference to result.
 */
  template <typename T>
      inline
      Vector<T>& Vector<T>::multElements(const Vector<T>& B)
  {
    this->type_vector->multiplyElements(B.type_vector);
    return *(this);
  }

/**
   * Internal product of to Vectors
   *
   * Element-to-element multiplication between Vector objects.
   * \param A A Vector.
   * \param B A Vector.
   * \return Reference to result.
 */

  template <typename T>
      inline
      Vector<T>& Vector<T>::multElements(const Vector<T>& A, const Vector<T>& B)
  {
    if ( this->elements != A.size() || this->elements != B.size() ){
      std::stringstream message;
      message << "Vectors dimensions mismatch. \n"
          << "LHS vector dimension: (" << this->elements << ")" << endl
          << "First RHS vector dimension: (" << A.size() << ")" << endl
          << "Second RHS vector dimension: (" << B.size() << ")" << endl;
          LMX_THROW(dimension_error, message.str() );
    }
    this->type_vector->equals(A.type_vector);
    this->type_vector->multiplyElements(B.type_vector);
    return *(this);
  }

/**
 * Computes the first norm of the Vector.
 * \return The first norm of object.
 */
template <typename T>
    T Vector<T>::norm1 () const
{
  T norm = 0;
  for (size_type i=0; i<elements; ++i)
    norm += std::abs(type_vector->readElement(i,zero));

  return norm;
}

/** Computes de second (euclidean) norm of the Vector.
  * \return The second (euclidean) norm of object.*/
template <typename T>
    T Vector<T>::norm2 () const
{
  T norm = 0;
  for (size_type i=0; i<elements; ++i)
    norm += type_vector->readElement(i,zero) * type_vector->readElement(i,zero);
  return std::sqrt(norm);
}


 //////////////////////////////////////////////////////////////////
 // Overloaded operators with Vector as rvalue:  //////////////////
 //////////////////////////////////////////////////////////////////

/**
 * Operator for multiplying Matrix and Vector objects.
 * Performs typical vector operation: A * b = c, where b & c are vectors and A is a Matrix.
 * \param A Matrix reference for first LHS term.
 * \param b Vector reference for second LHS term.
 * \return A new Vector (c = A*b).
 */
template <typename T>
    Vector<T> operator * (const Matrix<T>& A, const Vector<T>& b)
{
  if ( A.cols() != b.size() ){
    std::stringstream message;
    message << "Can't multiply Matrix * Vector: Dimensions mismatch. \n"
        << "Matrix dimension: (" << A.cols() << "," << A.rows() << ")" << endl
        << "Vector dimension: (" << b.size() << ")" << endl;
    LMX_THROW(dimension_error, message.str() );
  }
  Vector<T> c(A.rows());
  c.mult(A,b);
  return c;
}

/**
 * Operator for multiplying DenseMatrix and Vector objects.
 * Performs typical vector operation: A * b = c, where b & c are vectors and A is a DenseMatrix.
 * \param A DenseMatrix reference for first LHS term.
 * \param b Vector reference for second LHS term.
 * \return A new Vector (c = A*b).
 */
template <typename T>
    Vector<T> operator * (const DenseMatrix<T>& A, const Vector<T>& b)
{
  if ( A.cols() != b.size() ){
    std::stringstream message;
    message << "Can't multiply DenseMatrix * vector: Dimensions mismatch. \n"
        << "DenseMatrix dimension: (" << A.cols() << "," << A.rows() << ")" << endl
        << "Vector dimension: (" << b.size() << ")" << endl;
    LMX_THROW(dimension_error, message.str() );
  }
  Vector<T> c(A.rows());
  c.mult(A,b);
  return c;
}

/**
 * Overload operator for multiplying a scalar and a Vector object.
 * \param a scalar reference.
 * \param B Vector reference for second LHS term.
 * \return A new Vector (mult = a*b).
 */
template <typename T>
    Vector<T> operator * (const T& a, const Vector<T>& B)
{
  Vector<T> mult(B);
  mult *= a;
  return mult;
}

/**
 * Overloaded operator that prints a Vector object into the output stream selected.
 * \param vec Vector reference to the object that is going to be printed.
 * \param os An output std stream.
 * \return (\a os ) The stream that was passed as a parameter.
 */
template <typename T>
    std::ostream& operator << (std::ostream& os, const Vector<T>& vec)
{
  os << "Vector (" << vec.size() << ") = " ;
  if ( !vec.size() ) os << "void";
  for (size_type i=0; i<vec.size(); ++i)
  {
    os << endl;
    os << vec.readElement(i) << " ";
  }
  os << endl;
  return os;
}


}; // namespace lmx


#endif
