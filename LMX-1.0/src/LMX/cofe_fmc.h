/***************************************************************************
 *   Copyright (C) 2005 by Jaime Planas, Jose M Sancho                     *
 ***************************************************************************/
/*
 *  fmc.h
 *  Addapted from fmc
 *  Created by Jaime Planas on 2005/11/13.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef FMC_H
#define FMC_H

#include "cofe_TensorRank1.h"
#include "cofe_TensorRank2.h"
#include "cofe_TensorRank2Sym.h"

// namespace fmc
// {

/*    cofe::TensorRank1<DIM> makecofe::TensorRank1<DIM>(const std::string & s); //Pseudo constructor, uso:
                                          // cofe::TensorRank1<DIM> v = makecofe::TensorRank1<DIM>("1.0 1.1 1.2");
    cofe::TensorRank2<DIM> makeTensor(const std::string & s); //Pseudo constructor
    cofe::TensorRank2Sym<DIM> makecofe::TensorRank2Sym<DIM>(const std::string & s); //Pseudo constructor  */
// }//namespace fmc

//vector-valued operators & functions
/*template <int DIM> inline
cofe::TensorRank1<DIM>  fmc::makecofe::TensorRank1<DIM>(const std::string & s)
{
    cofe::TensorRank1<DIM> v;
    v.initializeFrom(s);
    return v;
}*/


template <int DIM> inline
typename cofe::TensorRank1<DIM>  operator + (const typename cofe::TensorRank1<DIM> & v1, const typename cofe::TensorRank1<DIM> & v2)
{
    cofe::TensorRank1<DIM> r(v1);
    r += v2;
    return r;
}

template <int DIM> inline
cofe::TensorRank1<DIM>  operator - (const cofe::TensorRank1<DIM> & v1, const cofe::TensorRank1<DIM> & v2)
{
    cofe::TensorRank1<DIM> r(v1);
    r -= v2;
    return r;
}

template <int DIM> inline
cofe::TensorRank1<DIM>  operator * (const double a, const cofe::TensorRank1<DIM> & v)
{
    cofe::TensorRank1<DIM> r(v);
    r *= a;
    return r;
}

//Este operador no añade nada porque no evito sentencias al usarlo
//ya que escribir v2 = A*v1 es lo mismo que escribir v2.beProductOf(A,v1);
//permite, sin embargo, escribir fórmulas: pero es caro
template <int DIM> inline
cofe::TensorRank1<DIM>  operator * (const cofe::TensorRank2<DIM> & T, const cofe::TensorRank1<DIM> & v)
{
    cofe::TensorRank1<DIM> r;
    r.beProductOf(T,v);
    return r;
}

//Tensor-valued operators & functions
/*template <int DIM> inline
cofe::TensorRank2<DIM>  fmc::makeTensor(const std::string & s)
{
    cofe::TensorRank2<DIM> v;
    v.initializeFrom(s);
    return v;
}*/

template <int DIM> inline
cofe::TensorRank2<DIM>  operator + (const cofe::TensorRank2<DIM> & T1, const cofe::TensorRank2<DIM> & T2)
{
    cofe::TensorRank2<DIM> r(T1);
    r += T2;
    return r;
}

template <int DIM> inline
cofe::TensorRank2<DIM>  operator - (const cofe::TensorRank2<DIM> & T1, const cofe::TensorRank2<DIM> & T2)
{
    cofe::TensorRank2<DIM> r(T1);
    r -= T2;
    return r;
}

template <int DIM, class C> inline
cofe::TensorRank2<DIM>  operator * (const C a, const cofe::TensorRank2<DIM> & T)
{
    cofe::TensorRank2<DIM> r(T);
    r *= a;
    return r;
}

template <int DIM, class C> inline
cofe::TensorRank2<DIM>  operator * (const cofe::TensorRank2<DIM> & T, const C a)
{
    cofe::TensorRank2<DIM> r(T);
    r *= a;
    return r;
}

template <int DIM, typename T> inline
cofe::TensorRank2<DIM, T>  operator * (const cofe::TensorRank2<DIM, T> & T1, const cofe::TensorRank2<DIM, T> & T2)
{
    cofe::TensorRank2<DIM, T> r;
    r.beProductOf(T1,T2);
    return r;
}

template <int DIM> inline
cofe::TensorRank2<DIM>  transposeOf(const cofe::TensorRank2<DIM> & T)
{
    cofe::TensorRank2<DIM> r;
    r.beTransposeOf(T);
    return r;
}

//Symmetric-Tensor-valued operators
/*template <int DIM> inline
cofe::TensorRank2Sym<DIM>  fmc::makecofe::TensorRank2Sym<DIM>(const std::string & s)
{
    cofe::TensorRank2Sym<DIM> v;
    v.initializeFrom(s);
    return v;
}*/

#include <limits>
template <int DIM> inline
cofe::TensorRank2<DIM> inverseOf(const cofe::TensorRank2<DIM> & T)
{
    cofe::TensorRank2<DIM> r;
    double det,invdet;
    det = T.determinant();
    if(std::abs(det) <= std::numeric_limits<double>::min()) cofe::CofeUtils::error("Near singular matrix");
    invdet= 1.0/det;
    
    r(0,0)=  T(1,1)*T(2,2) - T(1,2)*T(2,1);
    r(1,0)=  T(1,2)*T(2,0) - T(1,0)*T(2,2);
    r(2,0)=  T(1,0)*T(2,1) - T(1,1)*T(2,0);
    
    r(0,1)=  T(0,2)*T(2,1) - T(0,1)*T(2,2);
    r(1,1)=  T(0,0)*T(2,2) - T(0,2)*T(2,0);
    r(2,1)=  T(0,1)*T(2,0) - T(0,0)*T(2,1);
    
    r(0,2)=  T(0,1)*T(1,2) - T(0,2)*T(1,1);
    r(1,2)=  T(0,2)*T(1,0) - T(0,0)*T(1,2);
    r(2,2)=  T(0,0)*T(1,1) - T(0,1)*T(1,0);
    
    r *= invdet;
    return r;
}


template <int DIM> inline
cofe::TensorRank2Sym<DIM>  operator + (const cofe::TensorRank2Sym<DIM> & T1, const cofe::TensorRank2Sym<DIM> & T2)
{
    cofe::TensorRank2Sym<DIM> r(T1);
    r += T2;
    return r;
}

template <int DIM> inline
cofe::TensorRank2Sym<DIM>  operator - (const cofe::TensorRank2Sym<DIM> & T1, const cofe::TensorRank2Sym<DIM> & T2)
{
    cofe::TensorRank2Sym<DIM> r(T1);
    r -= T2;
    return r;
}

template <int DIM> inline
cofe::TensorRank2Sym<DIM>  operator * (const double a, const cofe::TensorRank2Sym<DIM> & T)
{
    cofe::TensorRank2Sym<DIM> r(T);
    r *= a;
    return r;
}

template <int DIM> inline
cofe::TensorRank2Sym<DIM>  transposeOf(const cofe::TensorRank2Sym<DIM> & T)
{
    cofe::TensorRank2Sym<DIM> r(T);
    return r;
}


#endif //FMC_H
