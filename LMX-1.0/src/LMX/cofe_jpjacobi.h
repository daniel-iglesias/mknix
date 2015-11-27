/***************************************************************************
 *   Copyright (C) 2005 by Jaime Planas, Jose M Sancho                     *
 ***************************************************************************/

#ifndef JACOBI_H
#define JACOBI_H

#include "cofe_CofeUtils.h"

namespace cofe
{


template<int N, class T>
void jacobi(T a[N][N], T d[N], T v[N][N], int & nrot);
//a[N][N]= input matrix IS DESTROYED -> pass a copy
//d[N]= eigenvalue vector
//v[N][N]= v[][i] eigenvector i (column)
//nrot = number of Jacobi rotations

template<int N, class T>
void rot(T a[N][N],  T s,  T tau,  int i,  int j,  int k,  int l)
	{
		T g,h;

		g=a[i][j];
		h=a[k][l];
		a[i][j]=g-s*(h+g*tau);
		a[k][l]=h+s*(g-h*tau);
	}





template<int N, class T>
void jacobi(T a[N][N], T d[N], T v[N][N], int & nrot)
{
	using namespace std;
	int i,j,ip,iq;
	T tresh,theta,tau,t,sm,s,h,g,c;

	int n= N;
	T b[N],z[N];
	for (ip=0;ip<n;ip++) {
		for (iq=0;iq<n;iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	for (ip=0;ip<n;ip++) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
	}
	nrot=0;
	for (i=1;i<=50;++i) {
		sm=0.0;
		for (ip=0;ip<n-1;ip++) {
			for (iq=ip+1;iq<n;iq++)
				sm += fabs(a[ip][iq]);
		}
		if (sm == 0.0)
			return;
		if (i < 4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=0;ip<n-1;ip++) {
			for (iq=ip+1;iq<n;iq++) {
				g=100.0*fabs(a[ip][iq]);
				if (i > 4 && (fabs(d[ip])+g) == fabs(d[ip])
					&& (fabs(d[iq])+g) == fabs(d[iq]))
						a[ip][iq]=0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if ((fabs(h)+g) == fabs(h))
						t=(a[ip][iq])/h;
					else {
						theta=0.5*h/(a[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;
					for (j=0;j<ip;++j)
						cofe::rot<N,T>(a,s,tau,j,ip,j,iq);
					for (j=ip+1;j<iq;++j)
						cofe::rot<N,T>(a,s,tau,ip,j,j,iq);
					for (j=iq+1;j<n;++j)
						cofe::rot<N,T>(a,s,tau,ip,j,iq,j);
					for (j=0;j<n;++j)
						cofe::rot<N,T>(v,s,tau,j,ip,j,iq);
					++nrot;
				}
			}
		}
		for (ip=0;ip<n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
	cofe::CofeUtils::error("Too many iterations in routine jacobi");
}

}//namespace cofe

#endif// JACOBI_H


