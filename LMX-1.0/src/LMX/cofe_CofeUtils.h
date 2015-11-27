/***************************************************************************
 *   Copyright (C) 2005 by Jaime Planas, Jose M Sancho                     *
 ***************************************************************************/
/* 
	CofeUtils.h
  Created by Jaime Planas.
	Description:
*/

#ifndef CofeUtils_H
#define CofeUtils_H

#include <cmath>
#include <string>
#include <iostream>
#include <cstdlib>

/** \namespace cofe{}
 *
 *  \brief Continuum Oriented Finite Elements.
 *
 *  This namespace contains all the Tensor related facilities.
 *
 * \author Jaime Planas & J. M. Sancho.
 *
 */
namespace cofe {

class CofeUtils
{
public:
    static void error(const std::string & stri)
    {
        std::cout << stri << std::endl;
        std::cin.get();
        std::exit(0);
    }
    
    static void error(const std::string & stri, int code)
    {
        std::cout << stri << " " << code << std::endl;
        std::cin.get();
        std::exit(0);
    }
};

template<class T>
double norm(const T& a){ return std::sqrt(a.dot(a));}

}	// namespace cofe

#endif	// CofeUtils_H


