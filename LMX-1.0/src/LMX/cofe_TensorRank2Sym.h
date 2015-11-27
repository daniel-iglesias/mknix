/***************************************************************************
 *   Copyright (C) 2005 by Jaime Planas, Jose M Sancho                     *
 ***************************************************************************/
/* 
	TensorRank2Sym.h

	Author:			
	Description:	<describe the TensorRank2Sym class here>
*/

#ifndef TensorRank2Sym_H
#define TensorRank2Sym_H

//#include "MathObject.h"
#include "cofe_TensorRank1.h"
#include "cofe_TensorRank2.h"
#include <algorithm>
//#include "MathMatrix.h"
#include <cassert>
#include "cofe_jpjacobi.h"
#include <iostream>

namespace cofe {

template <int Dim, class T> class AbstractTensorRank4SS;


//tengo ++dudas sobre la implementacin + conveniente para el tensor simtrico
//he optado por la ms sencilla: guarda toda la matriz y fuerza la simetra

template <int Dim, class T = double>
class TensorRank2Sym:public TensorRank2<Dim,T>
{
public:
//	typedef typename cofe::MathMatrix<Dim,Dim,T> matrix;
//	typedef typename cofe::MathMatrix<(Dim*(Dim+1))/2,1,T> stress_likeVector;
	
	
	TensorRank2Sym(){}
	TensorRank2Sym(const T ptr[Dim][Dim]):TensorRank2<Dim,T>(ptr){assert(checkSymmetry());}
//	explicit TensorRank2Sym(const stress_likeVector & vec)
//		{for(int i = 0; i< Dim; ++i) for(int j = 0; j < Dim; ++j) this->com[i][j] = vec(index(i,j));}
	TensorRank2Sym(const T & aValue)
		:TensorRank2<Dim,T>(aValue){}
	TensorRank2Sym(const TensorRank2Sym<Dim,T> & aTen)
		{for(int i = 0; i < Dim; ++i) for(int j = 0; j < Dim; ++j) this->com[i][j] = aTen.com[i][j];}

	virtual ~TensorRank2Sym(){}
	void  operator = (const T& aValue){TensorRank2<Dim, T>::operator = (aValue);}
	TensorRank2Sym<Dim,T> &  operator = (const TensorRank2Sym<Dim,T> & aTen)
	{
		for(int i = 0; i < Dim; ++i) for(int j = 0; j < Dim; ++j) this->com[i][j] = aTen.com[i][j];
		return *this;
	}	
	virtual void beSymmetricPartOf(const cofe::TensorRank2<Dim, T>& r2T)
	{
		for(int i = 0; i< Dim; ++i)
		{
			for(int j = 0; j < Dim; ++j) this->com[i][j] = 0.5*(r2T(i,j)+r2T(j,i));
		}
	}

	virtual void beHydrostaticPartOf(const cofe::TensorRank2Sym<Dim, T>& r2T)
	{
	  this->beUnityTensor();
	  (*this) *= ((1.0/3.0)*r2T.trace());
	}

	virtual void beDeviatoricPartOf(const cofe::TensorRank2Sym<Dim, T>& r2T)
	{
	  cofe::TensorRank2Sym<Dim, T> aux;
	  aux.beHydrostaticPartOf(r2T);
	  (*this) = r2T;
	  (*this) -= aux;
	}

	void beProductOf(const cofe::AbstractTensorRank4SS<Dim, T>&  r4T,
						const cofe::TensorRank2Sym<Dim, T>& r2ST)
						{r4T.operateOn(r2ST,*this);}
	
	void beSymProductOf(const cofe::TensorRank1<Dim, T>&  lv,
						const cofe::TensorRank1<Dim, T>& rv)
	{
		for(int i = 0; i < Dim; ++i)
		{	
			this->com[i][i] = lv(i)*rv(i); //diagonal
			for(int j = i+1; j< Dim;++j) this->com[i][j]=this->com[j][i]= 0.5*(lv(i)*rv(j)+lv(j)*rv(i));
		}
	}
	
	bool checkSymmetry()
	{
		for(int i = 0; i< Dim; ++i)
			for(int j = 0; j<i;++j) 
				if(this->com[i][j] != this->com[j][i]) return false;
		return true;
	}
	
	double maxPrincipalValue()const;
	void computeMaxPrincipalDirection(TensorRank1<Dim,T>&)const;			 

protected:
	int index(int i, int j)const
	{if(j < i) std::swap(i,j); return i+((j-i)*(2*Dim+1-j+i))/2;}

};

template<int Dim, class T>
inline
double TensorRank2Sym<Dim,T>::maxPrincipalValue()const
{
	//implementado slo para dimension 2
	return 0.0;
}

template<>
inline
double TensorRank2Sym<3,double>::maxPrincipalValue()const
{
  double I1, I2, I3, aux, b, c, x = 0.0;
  const double locPI = 4.0*std::atan(1.0);//convendr� meterlo en un global
  double m, n, t;
   I1 = this->trace();
  I2 = this->secondInvariant();
  I3 = this->determinant();
  aux = I1*I1/3.0;
  b = I2-aux;
  c = I1*(-2.0*aux+3.0*I2)/9.0-I3;

  for (int i = 0; i < 1; ++i) 
    //ojo, s�o calcula el primer autovalor (supongo que es el mayor)
    {
      if ( std::abs(b) <= 1.0e-9)
	{
	  x = -std::pow(c,1.0/3.0);
	}
      else
	{
	  m = 2.0*std::sqrt(-b/3.0);
	  n = 3.0*c/(m*b);
	  t = std::atan2(std::sqrt(std::abs(1.0-n*n)),n)/3.0;
	  x = m*std::cos(t+2.0*i*locPI/3.0); 
	}
      x += I1/3.0;
    }
  return x;
}


template<>
inline
double TensorRank2Sym<2,double>::maxPrincipalValue()const
{
  double p = this->trace()*0.5;
  double d = 0.5*(this->com[0][0] - this->com[1][1]);
  double tau   = this->com[0][1];
  double q     = std::sqrt(d*d+tau*tau);
  return p + q;
}

template<int Dim, class T>
inline
void TensorRank2Sym<Dim,T>::computeMaxPrincipalDirection(TensorRank1<Dim,T>&)const
{
	//implementado slo para dimension 2
}

template<>
inline
void TensorRank2Sym<2,double>::computeMaxPrincipalDirection(TensorRank1<2,double>& prdir)const
{
  double d = 0.5*(this->com[0][0] - this->com[1][1]);
  double tau   = this->com[0][1];
  double angle = 0.5*std::atan2(tau,d);
  prdir(0) = std::cos(angle);
  prdir(1) = std::sin(angle);
}

template<>
inline
void TensorRank2Sym<3,double>::computeMaxPrincipalDirection(TensorRank1<3,double>& prdir)const
{
  TensorRank2Sym<3,double>  aux;
  aux = (*this);
  double eiva[3], eive[3][3];
  int nrot;

  cofe::jacobi<3,double>(aux.com,eiva,eive,nrot);

  double max=eiva[0];int i=0;
  if ( eiva[1] > eiva[0] ) {max = eiva[1]; i=1;}
  if ( eiva[2] > max     ) {i=2;}

  prdir(0) = eive[0][i];
  prdir(1) = eive[1][i];
  prdir(2) = eive[2][i];
}
			 


}	// namespace cofe

#endif	// TensorRank2Sym_H


