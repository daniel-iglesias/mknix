/***************************************************************************
 *   Copyright (C) 2005 by Jaime Planas, Jose M Sancho                     *
 ***************************************************************************/
// TensorRank2.h
// COFE: continuum-oriented finite elements
// Copyright  2003-2005 Jaime Planas and Jose M Sancho
// Universidad Politecnica de Madrid, Spain

#ifndef TensorRTwo_H
#define TensorRTwo_H

#include "cofe_TensorRank1.h"
#include "cofe_gausselim.h"
#include <limits>

namespace cofe
{

  template <int Dim, class T> class AbstractTensorRank4SS;

  template <int Dim, class T = double>
  class TensorRank2
  {
  public:
    friend class TensorRank1<Dim,T>;
  public:
	typedef	T	scalar_type;
	typedef	typename cofe::TensorRank1<Dim,T>	vector_type;
	typedef	typename cofe::TensorRank2<Dim,T>	tensor_type;
	typedef	typename cofe::TensorRank2<Dim,T>	own_type;
	typedef	typename cofe::AbstractTensorRank4SS<Dim, T> ar4tensor_type;
  public:
    TensorRank2();
    TensorRank2(const T ptr[Dim][Dim]);
    TensorRank2(const T & aValue);
    TensorRank2(const std::string& stringOfComponents);
    TensorRank2(const char* c_stringOfComponents);
    TensorRank2(const own_type & aTen);
    virtual ~TensorRank2();
    bool initializeFrom(const std::string & s);
    bool initializeFrom(std::istream & istre);
    int size();
    scalar_type & operator()(int i, int j);
    const scalar_type & operator()(int i, int j)const;
    scalar_type trace()const;
    tensor_type &  operator = (const T& aValue);
    tensor_type &  operator = (const std::string& stringOfComponents);
    tensor_type &  operator = (const char* c_stringOfComponents);
    tensor_type &  operator = (const tensor_type & aTen);
    void zero();
	void  operator *= (const scalar_type aValue);
	void  operator += (const TensorRank2<Dim, T> a);
	void  operator -= (const TensorRank2<Dim, T> a);
    scalar_type dot(const tensor_type & aTensor)const;
    void operateOn(const vector_type& operand, vector_type& answer)const;
    void beTransposeOf(const tensor_type& Ten);
    void beProductOf(const vector_type & leftVec,
                     const ar4tensor_type & r4T,
                     const vector_type & rightVec);
    void beProductOf(const vector_type& leftVec,
                     const vector_type& rightVec);
    void beProductOf(const tensor_type& leftTen,
                     const tensor_type& rightTen);
    void beTraProductOf(const tensor_type& leftTen, //Lt R
                     const tensor_type& rightTen);
    void beProductTraOf(const tensor_type& leftTen,//L Rt
                        const tensor_type& rightTen);
    void beUnityTensor();
    void solve(const vector_type& operand,
               vector_type & answer)const;  
    void beInverseOf(const tensor_type& atensor);
    scalar_type determinant() const;    
    scalar_type secondInvariant() const;
    
//////////////// Added by LMX-team:

    TensorRank2(const int & aValue);
    std::ostream & print_row(std::ostream &, const cofe::TensorRank2<Dim,T> &, int) const ;

    bool  operator == (const T& aValue) const {
      for(int i=0; i<Dim; ++i){
        for(int j=0; j<Dim; ++j)
          if(this->com[i][j] != aValue) return 0;
      }
      return 1;
    }

    bool  operator != (const T& aValue) const {
      for(int i=0; i<Dim; ++i){
        for(int j=0; j<Dim; ++j)
          if(this->com[i][j] != aValue) return 1;
      }
      return 0;
    }

    void  operator *= (const own_type aTensor){
      for(int i=0; i<Dim; ++i){
        for(int j=0; j<Dim; ++j)
          com[i][j] *= aTensor(i,j);
      }
    }

    operator T() {
      return this->com[0][0];
    }

/////////// Needed for gmm compatibility (overloaded operators...)
    friend own_type  operator * (const own_type & T1, const own_type & T2)
    {
        own_type r;
        r.beProductOf(T1,T2);
        return r;
    }

//////////////// Added by LMX-team (END)


  protected:
   	scalar_type com[Dim][Dim];
  };
  
  //implementation
  
   template <int Dim,class T>
   inline
   TensorRank2<Dim,T>::TensorRank2()
    {
    }
    
   template <int Dim,class T>
   inline
   TensorRank2<Dim,T>::TensorRank2(const T ptr[Dim][Dim])
    {
    	for(int i = 0; i < Dim; ++i) for(int j = 0; j < Dim; ++j) com[i][j] = ptr[i][j];
    }
    
    template <int Dim,class T>
	inline
    TensorRank2<Dim,T>::TensorRank2(const T & aValue)
    {
            for(int i = 0; i < Dim; ++i) for(int j = 0; j < Dim; ++j) com[i][j] =aValue;
    }
   
    template <int Dim,class T>
	inline
    TensorRank2<Dim,T>::TensorRank2(const int & aValue)
    {
            for(int i = 0; i < Dim; ++i) for(int j = 0; j < Dim; ++j) com[i][j] = T(aValue);
    }
   
   template <int Dim, class T>
       inline
   TensorRank2<Dim,T>::TensorRank2(const std::string& stringOfComponents)
   {
           this->initializeFrom(stringOfComponents);
   }
   
   template <int Dim, class T>
       inline
       TensorRank2<Dim,T>::TensorRank2(const char* c_s)
   {
           const std::string ss(c_s);
           this->initializeFrom(ss);
   }   
    
    template <int Dim,class T>
    inline
    TensorRank2<Dim,T>::TensorRank2(const own_type & aTen)
    {
            for(int i = 0; i < Dim; ++i) for(int j = 0; j < Dim; ++j) com[i][j] = aTen.com[i][j];
    }
    
    template <int Dim,class T>
    inline
    TensorRank2<Dim,T>::~TensorRank2()
    {
    }

   template <int Dim,class T>
    inline
    bool TensorRank2<Dim,T>::initializeFrom(const std::string & s)
   {
        std::istringstream ist(s);
        return initializeFrom(ist);
   }
   
   template <int Dim,class T>
    inline
    bool TensorRank2<Dim,T>::initializeFrom(std::istream & istre)
   {
        bool good = true;
        for(int i = 0; i < Dim; ++i) for(int j = 0; j < Dim; ++j) if (!(istre >> com[i][j])) good = false;
        return good;
   }
       
    template <int Dim,class T>
    inline
    int TensorRank2<Dim,T>::size()
    {
    	return Dim*Dim;
    }
	
  	template <int Dim,class T>
	inline
    T & TensorRank2<Dim,T>::operator()(int i, int j)
    {
    	return com[i][j];
    }
	
  	template <int Dim,class T>
	inline
    const T & TensorRank2<Dim,T>::operator()(int i, int j)const
    {
    	return com[i][j];
    }
	
  	template <int Dim,class T>
	inline
    T TensorRank2<Dim,T>::trace()const
    {
    	scalar_type x = 0.0;
    	for(int i = 0; i< Dim; ++i) x += com[i][i];
    	return x;
    }
	
    template <int Dim,class T>
    inline
    TensorRank2<Dim,T> &  TensorRank2<Dim,T>::operator = (const T& aValue)
    {
        for(int i = 0; i < Dim; ++i)
                for(int j = 0; j < Dim; ++j) 
                        com[i][j] =aValue;
        return *this;
    }
   
   template <int Dim,class T>
   inline
   TensorRank2<Dim,T> &  TensorRank2<Dim,T>::operator = (const std::string& s)
   {
       this->initializeFrom(s);
       return *this;
   }
   
   template <int Dim,class T>
       inline
       TensorRank2<Dim,T> &  TensorRank2<Dim,T>::operator = (const char* c_s)
   {
           std::string s(c_s);
           this->initializeFrom(s);
           return *this;
   }
	
  	template <int Dim,class T>
	inline
    TensorRank2<Dim,T> &  TensorRank2<Dim,T>::operator = (const tensor_type & aTen)
    {
      for(int i = 0; i < Dim; ++i)
      	for(int j = 0; j < Dim; ++j)
      		com[i][j] = aTen.com[i][j];
      return *this;
    }
	
  	template <int Dim,class T>
	inline
    void TensorRank2<Dim,T>::zero()
	{
		for(int i = 0; i < Dim; ++i)
			for(int j = 0; j < Dim; ++j)
				com[i][j] =0.0;
	}
	
  	template <int Dim,class T>
	inline
	void  TensorRank2<Dim,T>::operator *= (const scalar_type aValue)
	{
		for(int i = 0; i < Dim; ++i)
			for(int j = 0; j < Dim; ++j)
				com[i][j] *=aValue;
	}
	
  	template <int Dim,class T>
	inline
	void  TensorRank2<Dim,T>::operator += (const TensorRank2<Dim, T> a)
	{
		for(int i = 0; i < Dim; ++i)
			for(int j = 0; j < Dim; ++j)
				com[i][j] += a.com[i][j];
	}
	
  	template <int Dim,class T>
	inline
	void  TensorRank2<Dim,T>::operator -= (const TensorRank2<Dim, T> a)
	{
		for(int i = 0; i < Dim; ++i)
			for(int j = 0; j < Dim; ++j)
				com[i][j] -= a.com[i][j];
	}
	
  	template <int Dim,class T>
	inline
    T TensorRank2<Dim,T>::dot(const tensor_type & aTensor)const
    {
      scalar_type x = 0;
      for(int i = 0; i < Dim; ++i)
      	for(int j = 0; j < Dim; ++j)
      		x += (com[i][j] * aTensor.com[i][j]);
      return x;
    }
	
  	template <int Dim,class T>
	inline
    void TensorRank2<Dim,T>::operateOn(const vector_type& operand, vector_type& answer)const
    {
      for(int i = 0; i < Dim; ++i)
      {
        answer.com[i] = 0;
        for(int j = 0; j < Dim; ++j) answer.com[i] += (this->com[i][j])*(operand.com[j]);
      }
    }
   
   template <int Dim,class T>
       inline
   void TensorRank2<Dim,T>::beTransposeOf(const tensor_type& aTensor)
   {
       for(int i = 0; i < Dim; ++i) for(int j = 0; j < Dim; ++j) com[i][j] = aTensor.com[j][i]; 
   }
  	template <int Dim,class T>
	inline
    void TensorRank2<Dim,T>::beProductOf(const vector_type & leftVec,
                     const ar4tensor_type & r4T,
                     const vector_type & rightVec)
    {
      r4T.extremeContract(leftVec,rightVec, *this);
    }
	
  	template <int Dim,class T>
	inline
    void TensorRank2<Dim,T>::beProductOf(const vector_type& leftVec,
                     const vector_type& rightVec)
    {
      for(int i = 0; i < Dim; ++i)
      {
        for(int j = 0; j < Dim; ++j) com[i][j] = (leftVec.com[i])*(rightVec.com[j]);
      }
    }
	
  	template <int Dim,class T>
	inline
    void TensorRank2<Dim,T>::beProductOf(const tensor_type& leftTen,
                     const tensor_type& rightTen)
    {
      for(int i = 0; i < Dim; ++i)
      {
        for(int j = 0; j < Dim; ++j)
        {
            com[i][j] = 0.0;
            for(int k = 0; k < Dim; ++k)  com[i][j] += (leftTen.com[i][k])*(rightTen.com[k][j]);
        }
      }
    }
   
   template <int Dim,class T>
       inline
       void TensorRank2<Dim,T>::beTraProductOf(const tensor_type& leftTen,
                                            const tensor_type& rightTen)
   {
           for(int i = 0; i < Dim; ++i)
           {
               for(int j = 0; j < Dim; ++j)
               {
                   com[i][j] = 0.0;
                   for(int k = 0; k < Dim; ++k)  com[i][j] += (leftTen.com[k][i])*(rightTen.com[k][j]);
               }
           }
   }
   
   template <int Dim,class T>
       inline
       void TensorRank2<Dim,T>::beProductTraOf(const tensor_type& leftTen,
                                               const tensor_type& rightTen)
   {
           for(int i = 0; i < Dim; ++i)
           {
               for(int j = 0; j < Dim; ++j)
               {
                   com[i][j] = 0.0;
                   for(int k = 0; k < Dim; ++k)  com[i][j] += (leftTen.com[i][k])*(rightTen.com[j][k]);
               }
           }
   }
   
   
	
  	template <int Dim,class T>
	inline
	void TensorRank2<Dim,T>::beUnityTensor()
	{
		this->zero();
		for(int i = 0; i < Dim;++i) com[i][i] = 1.0;
	}
	
	template <int Dim,class T>
	inline
	T TensorRank2<Dim,T>::determinant() const
	{
		//implemented only for Dim = 1,2,3
		return com[0][0];
	}

  template<>
  inline
  double TensorRank2<2,double>::determinant() const
  {
    return com[0][0]*com[1][1]-com[0][1]*com[1][0];
  }

  template<>
  inline
  double TensorRank2<3,double>::determinant() const
  {

    return com[0][0]*(com[1][1]*com[2][2]-com[2][1]*com[1][2])+
           com[0][1]*(com[1][2]*com[2][0]-com[1][0]*com[2][2])+
           com[0][2]*(com[1][0]*com[2][1]-com[1][1]*com[2][0]);
  }

  template <int Dim,class T>
  inline
  void TensorRank2<Dim, T>::solve(const vector_type& rh,
                                  vector_type & lh)const
  {
    TensorRank2<Dim,double>  auxm;
    auxm = (*this);
    TensorRank1<Dim,double>  auxv;
    auxv = rh;
    cofe::gausselim(auxm.com,auxv.com);
    lh = auxv;
  }

  template<>
  inline
  void TensorRank2<2,double>::solve(const cofe::TensorRank1<2,double>& rh,
                                    cofe::TensorRank1<2,double> & lh)const
  {
    double detInv= 1.0/(com[0][0]*com[1][1]-com[0][1]*com[1][0]);
    lh(0) = (com[1][1]*rh(0)-com[0][1]*rh(1))*detInv;
    lh(1) = (com[0][0]*rh(1)-com[1][0]*rh(0))*detInv;
  }

  template <int Dim,class T>
  inline
  void TensorRank2<Dim, T>::beInverseOf(const tensor_type& )
  {}

  template <>
  inline
  void TensorRank2<2, double>::beInverseOf(const cofe::TensorRank2<2, double>& aT)
  {
    double detInv= 1.0/(aT(0,0)*aT(1,1)-aT(0,1)*aT(1,0));
    com[0][0] = detInv*aT(1,1);
    com[0][1] = -detInv*aT(0,1);
    com[1][0] = -detInv*aT(1,0);
    com[1][1] = detInv*aT(0,0);
  }
  

  template <>
  inline
  void TensorRank2<3, double>::beInverseOf(const cofe::TensorRank2<3, double>& T)
  {
    double det,invdet;
    det = T.determinant();    
    if(std::abs(det) <= std::numeric_limits<double>::min()) cofe::CofeUtils::error("Near singular matrix");
    invdet= 1.0/det;
    
    (*this)(0,0)=  T(1,1)*T(2,2) - T(1,2)*T(2,1);
    (*this)(1,0)=  T(1,2)*T(2,0) - T(1,0)*T(2,2);
    (*this)(2,0)=  T(1,0)*T(2,1) - T(1,1)*T(2,0);
    
    (*this)(0,1)=  T(0,2)*T(2,1) - T(0,1)*T(2,2);
    (*this)(1,1)=  T(0,0)*T(2,2) - T(0,2)*T(2,0);
    (*this)(2,1)=  T(0,1)*T(2,0) - T(0,0)*T(2,1);
    
    (*this)(0,2)=  T(0,1)*T(1,2) - T(0,2)*T(1,1);
    (*this)(1,2)=  T(0,2)*T(1,0) - T(0,0)*T(1,2);
    (*this)(2,2)=  T(0,0)*T(1,1) - T(0,1)*T(1,0);

    (*this) *= invdet;
}

  template <int Dim, typename T>
  inline
  T TensorRank2<Dim,T>::secondInvariant() const
  {
  	T tr = this->trace();
  	return 0.5*(tr*tr - this->dot(*this));
  }

  template <int Dim,class T>
  std::ostream & operator << (std::ostream & os, const cofe::TensorRank2<Dim,T> & a)
  {
    for(int i = 0; i< Dim; ++i)
    {
      os <<"|";
      for(int j = 0; j< Dim-1;++j)os << a(i,j) << " ";
      os << a(i,Dim-1) << "|\n";
    }
    return os;
  }

//////////////// Added by LMX-team:
  template <int Dim,class T>
  std::ostream & TensorRank2<Dim,T>::print_row(std::ostream & os, const cofe::TensorRank2<Dim,T> & a, int row) const
  {
      os <<"|";
      for(int j = 0; j< Dim-1;++j) os << a(row,j) << " ";
      os << a(row,Dim-1) ;
      os <<"|";
    return os;
  }


}	// namespace cofe

#endif	// TensorRTwo_H

