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
#include "common.h"

namespace mknix
{

double interpolate1D(double key, const std::map<double, double>& theMap)
{
    typedef std::map<double, double>::const_iterator i_t;

    i_t i = theMap.upper_bound(key);
    if (i == theMap.end())
    {
        return (--i)->second;
    }
    if (i == theMap.begin())
    {
        return i->second;
    }
    i_t l = i;
    --l;

    const double delta = (key - l->first) / (i->first - l->first);
    return (delta * i->second + (1 - delta) * l->second);
}


boxFIR::boxFIR(int _numCoeffs) :
    numCoeffs(_numCoeffs * 2)
{
    if (numCoeffs < 1)   //Must be > 0 or bad stuff happens
    {
        numCoeffs = 2;
    }

    double val = 1. / numCoeffs;
    for (std::size_t ii = 0; ii < numCoeffs; ++ii)
    {
        b.push_back(val);
        m.push_back(0.);
    }
}

void boxFIR::filter(vector<double>& a)
{
    double output;

    // init with all memories equal to first values:
    for (std::size_t ii = 0; ii < numCoeffs; ++ii)
    {
        m[ii] = a[ii];
    }

    for (std::size_t nn = 0; nn < a.size(); ++nn)
    {
        output = 0;
        if (nn < numCoeffs / 2)
        {
            //Apply smoothing filter to signal
            //     m[0] = a[nn+numCoeffs/2];
            for (std::size_t ii = 0; ii < numCoeffs; ++ii)
            {
                output += b[ii] * m[ii];
            }
        }
        else if ((a.size() - nn) < numCoeffs)
        {
            //Apply smoothing filter to signal
            m[0] = a[nn];
            for (std::size_t ii = 0; ii < numCoeffs; ++ii)
            {
                output += b[ii] * m[ii];
            }
        }
        else
        {
            //Apply smoothing filter to signal
            m[0] = a[nn + numCoeffs / 2];
            for (std::size_t ii = 0; ii < numCoeffs; ++ii)
            {
                output += b[ii] * m[ii];
            }
        }

        //Reshuffle memories
        if (nn > numCoeffs / 2)
        {
            for (std::size_t ii = numCoeffs - 1; ii != 0; --ii)
            {
                m[ii] = m[ii - 1];
            }
        }
        a[nn] = output;
    }
}

std::vector<double> doubles_in_vector(const std::string& str)
{
    std::istringstream stm(str); // input stringstream to read from the line

    // create a vector containing doubles in the line (left to right)
    using iterator = std::istream_iterator<double>;
    std::vector<double> seq{ iterator(stm), iterator() };

    return seq; // and return it
}

std::vector<std::vector<double>> read_lines(std::istream& stm)
{
    std::vector<std::vector<double>> result;
    std::string line;

    while (std::getline(stm, line)) result.push_back(doubles_in_vector(line));

    return result;
}


}
