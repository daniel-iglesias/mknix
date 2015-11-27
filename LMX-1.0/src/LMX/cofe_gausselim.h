/***************************************************************************
 *   Copyright (C) 2005 by Jaime Planas, Jose M Sancho                     *
 ***************************************************************************/

#ifndef GAUSSELIM_H
#define GAUSSELIM_H

#include <iostream>
#include "cofe_CofeUtils.h"

namespace cofe
{
template<int N, class T>
void gausselim(T a[N][N], T b[N] );
}

template<int N, class T>
void cofe::gausselim(T a[N][N], T b[N])
{
  int n = N;
  for(int k = 0; k < n - 1 ;++k)
    {
      if (a[k][k] == 0)
	{cofe::CofeUtils::error("pivote nulo");}//	error;
      for(int i = k + 1; i < n; ++i)
	{

	  if(a[i][k] != 0)
	    {
	      T mult = a[i][k]/a[k][k];
	      a[i][k]=mult;
	      for(int j = k + 1; j < n; ++j)
		a[i][j]-= mult*a[k][j];
	    }
	}
    }

  for(int i = 1; i < n; ++i)
    for(int j = 0; j < i; ++j) b[i] -= a[i][j]*b[j];
  for(int i = n - 1; i >= 0; i--){
    for(int j = i+1; j<n; ++j) b[i] -= a[i][j]*b[j];
    b[i] /= a[i][i];
  }
}

#endif
