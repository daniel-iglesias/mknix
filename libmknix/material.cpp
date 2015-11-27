/***************************************************************************
 *   Copyright (C) 2013 by Daniel Iglesias                                 *
 *   http://code.google.com/p/mknix                                        *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "material.h"

namespace mknix {

Material::Material()
    : dim(0)
    , capacity(0)
    , kappa(0)
    , young(0)
    , poisson(0)
    , beta(0)
    , density(0)
{
}

Material::~Material()
{
}

void Material::setThermalProps( double capacity_in, double kappa_in, double beta_in, double density_in )
{
    capacity = capacity_in;
    kappa = kappa_in;
    beta = beta_in;
    density = density_in;
}

// Mechanical needs dimension
void Material::setMechanicalProps( int dim_in, double young_in, double poisson_in, double density_in )
{
    dim = dim_in;
    young = young_in;
    poisson = poisson_in;
    density = density_in;
    lambda = (poisson*young) / ( (1+poisson)*(1-2*poisson) );
    mu = young / ( 2*(1+poisson) );
    computeD();
    computeC();
}

void Material::computeD( )
{
    double comFacD; // Common factor for matrix D.

    if(dim==2) {
        D.resize(3,3);
        // Plain stress case:
        comFacD = young / (1 - std::pow(poisson, 2) ); // = E/(1-mu^2)
        D(0,0) = 1;
        D(0,1) = poisson;
        D.writeElement(D.readElement(0,1),1,0);
        D(1,1) = 1;
        D(2,2) = (1-poisson) / 2;
        D *= comFacD; //Apply common factor.
    }
    else if(dim==3) {
        D.resize(6,6);
        comFacD = young / ( (1+poisson)*(1-2*poisson) ); // = E/((1+mu)(1-2*mu))
        D.writeElement(1.-poisson,0,0);
        D.writeElement(1.-poisson,1,1);
        D.writeElement(1.-poisson,2,2);
        D.writeElement(   poisson,0,1);
        D.writeElement(   poisson,0,2);
        D.writeElement(   poisson,1,2);
        D.writeElement(D.readElement(0,1),1,0);
        D.writeElement(D.readElement(0,2),2,0);
        D.writeElement(D.readElement(2,1),1,2);
        D.writeElement( 0.5-poisson,3,3);
        D.writeElement( 0.5-poisson,4,4);
        D.writeElement( 0.5-poisson,5,5);
        D *= comFacD; //Apply common factor.
    }
}

void Material::computeC( )
{
    if(dim==2) {
        C.resize(3,3);
        // Plain stress case:
        C(0,0) = lambda+2*mu;
        C(0,1) = lambda;
        C.writeElement(C.readElement(0,1),1,0);
        C(1,1) = lambda+2*mu;
        C(2,2) = mu;
    }
    else if(dim==3) {
        C.resize(6,6);
        // Plain stress case:
        C(0,0) = lambda+2*mu;
        C(0,1) = lambda;
        C(0,2) = lambda;
        C(1,2) = lambda;
        C.writeElement(C.readElement(0,1),1,0);
        C.writeElement(C.readElement(0,2),2,0);
        C.writeElement(C.readElement(1,2),2,1);
        C(1,1) = lambda+2*mu;
        C(2,2) = lambda+2*mu;
        C(3,3) = mu;
        C(4,4) = mu;
        C(5,5) = mu;
    }
    //   cout << lambda << ", " << mu << endl;
//   cout << C << endl;
//   int kk=1; cout << this->getCsym(kk,kk,kk,kk) << endl;
}

double Material::getCsym( int& i, int& j, int& k, int& l )
{
    return (0.25*( Cijkl(i,j,k,l) + Cijkl(i,j,l,k) + Cijkl(j,i,k,l) + Cijkl(j,i,l,k) ) );
}

double Material::Cijkl( int& i, int& j, int& k, int& l )
{
    double res=0.;
    if(i==j && k==l) res += lambda;
    if(i==k && j==l) res += 2*mu;
    return res;
}


void Material::computeS(cofe::TensorRank2Sym<2,double> & S,
                        const cofe::TensorRank2<2,double> & F,
                        double temperature_in
                       )
{
    cofe::TensorRank2Sym<2,double> one;
    one.beUnityTensor();
    // E = 1/2 (F'F - 1)
    E.beTraProductOf(F,F);
    E -= one;
    E *= 0.5;
    // Adding thermal expansion
    one *= beta*temperature_in;
    E -= one;
    one.beUnityTensor();
    // St. Venant Kirchoff:
    // S = lambda*tr(E)*1 + 2*mu*E
    S = E;
    S *= 2*mu;
    one *= lambda*E.trace();
    S += one;
//     // Adding thermal expansion
//     one.beUnityTensor();
//     one *= beta*temperature_in;
//     S -= one;
}

void Material::computeS(cofe::TensorRank2Sym<3,double> & S,
                        const cofe::TensorRank2<3,double> & F)
{
    cofe::TensorRank2Sym<3,double> one;
    one.beUnityTensor();
    // E = 1/2 (F'F - 1)
    E3.beTraProductOf(F,F);
    E3 -= one;
    E3 *= 0.5;
    // St. Venant Kirchoff:
    // S = lambda*tr(E)*1 + 2*mu*E
    S = E3;
    S *= 2*mu;
    one *= lambda*E3.trace();
    S += one;

}

double Material::computeEnergy( const cofe::TensorRank2<2,double> & F )
{
    cofe::TensorRank2Sym<2,double> one;
    one.beUnityTensor();
    // E = 1/2 (F'F - 1)
    E.beTraProductOf(F,F);
    E -= one;
    E *= 0.5;

    // W = 1/2 lambda * tr(E)^2 + mu * E:E
    double energy = 0.5 * lambda * pow( E.trace(), 2 ) + mu * E.dot( E );
    return energy;
}

double Material::computeEnergy( const cofe::TensorRank2<3,double> & F )
{
    cofe::TensorRank2Sym<3,double> one;
    one.beUnityTensor();
    // E = 1/2 (F'F - 1)
    E3.beTraProductOf(F,F);
    E3 -= one;
    E3 *= 0.5;

    // W = 1/2 lambda * tr(E)^2 + mu * E:E
    double energy = 0.5 * lambda * pow( E.trace(), 2 ) + mu * E.dot( E );
    return energy;
}

}
