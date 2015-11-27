/***************************************************************************
 *   Copyright (C) 2005 by Jaime Planas, Jose M Sancho                     *
 ***************************************************************************/
// TensorRank1.h
// COFE: continuum-oriented finite elements
// Copyright � 2003-2005 Jaime Planas and Jose M Sancho
// Universidad Politecnica de Madrid, Spain


#ifndef TensorRank1_H
#define TensorRank1_H

#include <ostream>
#include <sstream>


namespace cofe {

template <int Dim, class T> class TensorRank2;

template <int Dim, class T = double>
class TensorRank1 //: public virtual MathObject
{
public:
	friend class TensorRank2<Dim,T>;
public:
	typedef           T	                        scalar_type;
	typedef typename cofe::TensorRank1<Dim,T>	vector_type;
	typedef	 typename cofe::TensorRank2<Dim,T>	tensor_type;
	typedef	 typename cofe::TensorRank1<Dim,T>	own_type;

public:
	TensorRank1();
	TensorRank1(const T* ptr);
	TensorRank1(const T& aValue);
        TensorRank1(const std::string& stringOfComponents);
        TensorRank1(const char* c_s);
	TensorRank1(const own_type& inV);
	~TensorRank1();
        bool initializeFrom(const std::string & s);
        bool initializeFrom(std::istream & istre);
	int size();
	scalar_type& operator()(int i);
	const scalar_type& operator()(int i)const;
	void  operator *= (const scalar_type& aValue);
	own_type&  operator = (const scalar_type& aValue);
        own_type&  operator = (const std::string& stringOfComponents);
        own_type&  operator = (const char* c_s);
	own_type&  operator = (const own_type& inV);
	void zero();
	void copyFrom(const scalar_type* ptr);
	void  operator += (const vector_type aV);
	void  operator -= (const vector_type aV);
	scalar_type dot(const vector_type & aVector)const;
	void beProductOf(const tensor_type & aR2T,
					 const vector_type & aV);
	void beSolutionOf(const tensor_type & aR2T,
				   	  const vector_type & aV);
protected:
	scalar_type com[Dim];
};

//implementation
template <int Dim, class T>
inline
TensorRank1<Dim,T>::TensorRank1()
{
}

template <int Dim, class T>
inline
TensorRank1<Dim,T>::TensorRank1(const T* ptr)
{
	for(int i = 0; i<Dim;++i) com[i] = ptr[i];
}

template <int Dim, class T>
inline
TensorRank1<Dim,T>::TensorRank1(const T& aValue)
{
	for(int i = 0; i<Dim;++i) com[i] = aValue;
}

template <int Dim, class T>
inline
TensorRank1<Dim,T>::TensorRank1(const std::string& stringOfComponents)
{
    this->initializeFrom(stringOfComponents);
}

template <int Dim, class T>
inline
TensorRank1<Dim,T>::TensorRank1(const char* s)
{
    const std::string ss(s);
    this->initializeFrom(ss);
}

template <int Dim, class T>
inline
TensorRank1<Dim,T>::TensorRank1(const own_type& inV)
{
	for(int i = 0; i < Dim;++i) com[i] = inV.com[i];
}

template <int Dim, class T>
inline
TensorRank1<Dim,T>::~TensorRank1()
{
}

template <int Dim,class T>
inline
bool TensorRank1<Dim,T>::initializeFrom(const std::string & s)
{
    std::istringstream ist(s);
    return initializeFrom(ist);
}

template <int Dim,class T>
inline
bool TensorRank1<Dim,T>::initializeFrom(std::istream & istre)
{
    bool good = true;
    for(int i = 0; i < Dim ;++i) if ( ! (istre >> com[i]) ) good = false;
    return good;
}

template <int Dim, class T>
inline
int TensorRank1<Dim,T>::size()
{
	return Dim;
}

template <int Dim, class T>
inline
T & TensorRank1<Dim,T>::operator()(int i)
{
	return com[i];
}

template <int Dim, class T>
inline
const T& TensorRank1<Dim,T>::operator()(int i)const
{
	return com[i];
}

template <int Dim, class T>
inline
void  TensorRank1<Dim,T>::operator *= (const scalar_type& aValue)
{
	for(int i = 0; i<Dim;++i) com[i] *= aValue;
}

template <int Dim, class T>
inline
TensorRank1<Dim,T> &  TensorRank1<Dim,T>::operator = (const scalar_type& aValue)
{
    for(int i = 0; i<Dim;++i) com[i] = aValue;
    return *this;
}

template <int Dim, class T>
inline
TensorRank1<Dim,T> &  TensorRank1<Dim,T>::operator = (const std::string& s)
{
    this->initializeFrom(s);
    return *this;
}

template <int Dim, class T>
inline
TensorRank1<Dim,T> &  TensorRank1<Dim,T>::operator = (const char* c_s)
{
    std::string s(c_s);
    this->initializeFrom(s);
    return *this;
}

template <int Dim, class T>
inline
TensorRank1<Dim,T> &  TensorRank1<Dim,T>::operator = (const own_type& inV)
{
	for(int i = 0; i < Dim;++i) com[i] = inV.com[i];
	return *this;
}

template <int Dim, class T>
inline
void TensorRank1<Dim,T>::zero()
{
	for(int i = 0; i < Dim;++i) com[i] = 0.0;
}

template <int Dim, class T>
inline
void TensorRank1<Dim,T>::copyFrom(const scalar_type* ptr)
{
	for(int i = 0; i<Dim;++i) com[i] = ptr[i];
}

template <int Dim, class T>
inline
void  TensorRank1<Dim,T>::operator += (const vector_type aV)
{
	for(int i = 0; i < Dim;++i) com[i] += aV.com[i];
}

template <int Dim, class T>
inline
void  TensorRank1<Dim,T>::operator -= (const vector_type aV)
{
	for(int i = 0; i < Dim;++i) com[i] -= aV.com[i];
}


template <int Dim, class T>
inline
T TensorRank1<Dim,T>::dot(const vector_type & aVector)const
{
	scalar_type x = 0.0;
	for(int i = 0; i<Dim;++i) x += (com[i]*aVector.com[i]);
	return x;
}

template <int Dim, class T>
inline
void TensorRank1<Dim,T>::beProductOf
		(const tensor_type & aR2T,
		 const vector_type & aV)
{
	aR2T.operateOn(aV,*this);
}

template <int Dim, class T>
inline
void TensorRank1<Dim,T>::beSolutionOf
		(const tensor_type & aR2T,
		 const vector_type & aV)
{
	aR2T.solve(aV,*this);
}


template <int Dim,class T>
std::ostream & operator << (std::ostream & os, const cofe::TensorRank1<Dim,T> & a)
{
	os << "(";
	for(int i = 0; i< Dim-1; ++i) os << a(i)  << ",";
	os << a(Dim-1) << ")";
	return os;
}



}	// namespace cofe

#endif	// TensorRank1_H




