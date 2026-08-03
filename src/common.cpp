/******************************************************************************
 *  Copyright (C) 2015 by Daniel Iglesias                                     *
 *                                                                            *
 *  This file is part of MkniX.                                             *
 *                                                                            *
 *  MkniX is free software: you can redistribute it and/or modify           *
 *  it under the terms of the GNU Lesser General Public License as            *
 *  published by the Free Software Foundation, either version 3 of the        *
 *  License, or (at your option) any later version.                           *
 *                                                                            *
 *  MkniX is distributed in the hope that it will be useful,                *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU Lesser General Public License for more details.                       *
 *                                                                            *
 *  You should have received a copy of the GNU Lesser General Public          *
 *  License along with MkniX.  If not, see <http://www.gnu.org/licenses/>.  *
 *****************************************************************************/

#include "common.h"

#ifdef HAVE_VTK
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellLocator.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkGenericCell.h>
#include <vtkIdList.h>
#include <vtkCell.h>
#endif

namespace mknix
{

/**
 * @brief Performs linear interpolation on a 1D map of (key, value) pairs.
 *
 * Given a lookup key and a sorted map of data points, returns the linearly
 * interpolated value. If the key is beyond the upper bound of the map the
 * last value is returned; if it is below the lower bound the first value
 * is returned.
 *
 * @param key       The lookup key to interpolate at.
 * @param theMap    A sorted map of (key, value) data points.
 * @return The interpolated value at the given key.
 */
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

/**
 * @brief Performs bilinear interpolation on a 2D nested map of data points.
 *
 * Interpolates along the first key axis by finding the two bracketing rows
 * and calling interpolate1D on each row for the second key, then linearly
 * combines the two results. Returns 0.0 if the map is empty.
 *
 * @param key1      The first-axis lookup key.
 * @param key2      The second-axis lookup key.
 * @param the2DMap  A sorted nested map: outer key → inner (key, value) map.
 * @return The bilinearly interpolated value at (key1, key2).
 */
double interpolate2D(
    double key1,
    double key2,
    const std::map<double, std::map<double, double>>& the2DMap)
{
    typedef std::map<double, std::map<double, double>>::const_iterator i2_t;

    if (the2DMap.empty())
    {
        return 0.0;
    }

    i2_t i = the2DMap.upper_bound(key1);
    if (i == the2DMap.end())
    {
        --i;
        if (i->second.empty())
        {
            return 0.0;
        }
        return interpolate1D(key2, i->second);
    }
    if (i == the2DMap.begin())
    {
        if (i->second.empty())
        {
            return 0.0;
        }
        return interpolate1D(key2, i->second);
    }

    i2_t l = i;
    --l;

    const double vL = l->second.empty() ? 0.0 : interpolate1D(key2, l->second);
    const double vU = i->second.empty() ? vL : interpolate1D(key2, i->second);

    if (i->first == l->first)
    {
        return vL;
    }

    const double delta = (key1 - l->first) / (i->first - l->first);
    return (delta * vU + (1.0 - delta) * vL);
}

/**
 * @brief Performs trilinear interpolation on a 3D nested map of data points.
 *
 * Interpolates along the first key axis by finding the two bracketing slices
 * and calling interpolate2D on each slice for (key2, key3), then linearly
 * combines the two results. Returns 0.0 if the map is empty.
 *
 * @param key1      The first-axis lookup key (x coordinate).
 * @param key2      The second-axis lookup key (y coordinate).
 * @param key3      The third-axis lookup key (z coordinate).
 * @param the3DMap  A sorted nested map: outer key → 2D nested map of (key, value) pairs.
 * @return The trilinearly interpolated value at (key1, key2, key3).
 */
double interpolate3D(
    double key1,
    double key2,
    double key3,
    const std::map<double, std::map<double, std::map<double, double>>>& the3DMap)
{
    typedef std::map<double, std::map<double, std::map<double, double>>>::const_iterator i3_t;

    if (the3DMap.empty())
    {
        return 0.0;
    }

    i3_t i = the3DMap.upper_bound(key1);
    if (i == the3DMap.end())
    {
        --i;
        if (i->second.empty())
        {
            return 0.0;
        }
        return interpolate2D(key2, key3, i->second);
    }
    if (i == the3DMap.begin())
    {
        if (i->second.empty())
        {
            return 0.0;
        }
        return interpolate2D(key2, key3, i->second);
    }

    i3_t l = i;
    --l;

    const double vL = l->second.empty() ? 0.0 : interpolate2D(key2, key3, l->second);
    const double vU = i->second.empty() ? vL : interpolate2D(key2, key3, i->second);

    if (i->first == l->first)
    {
        return vL;
    }

    const double delta = (key1 - l->first) / (i->first - l->first);
    return (delta * vU + (1.0 - delta) * vL);
}


/**
 * @brief Constructs a box (moving-average) FIR filter with the given number of coefficients.
 *
 * The effective filter length is doubled internally. If the resulting length
 * is less than 1 it is clamped to 2. All coefficients are set to 1/numCoeffs
 * so that the filter computes a uniform moving average.
 *
 * @param _numCoeffs Desired half-length of the filter; the internal length
 *                   will be 2 * _numCoeffs.
 */
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

/**
 * @brief Applies the box FIR filter in-place to a signal vector.
 *
 * Smooths the input signal using a centred moving-average strategy:
 * the beginning and end regions of the signal are handled by holding
 * the boundary values, while the interior is filtered using a look-ahead
 * window of numCoeffs/2 samples. The result overwrites the input vector.
 *
 * @param a Signal vector to be filtered in-place.
 */
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

/**
 * @brief Parses all whitespace-separated double values from a string.
 *
 * Uses a string stream and an istream_iterator to extract every floating-point
 * number present in the string, from left to right.
 *
 * @param str Input string containing whitespace-separated numeric values.
 * @return A vector of doubles parsed from the string, in order.
 */
std::vector<double> doubles_in_vector(const std::string& str)
{
    std::istringstream stm(str); // input stringstream to read from the line

    // create a vector containing doubles in the line (left to right)
    using iterator = std::istream_iterator<double>;
    std::vector<double> seq{ iterator(stm), iterator() };

    return seq; // and return it
}

/**
 * @brief Reads all lines from an input stream and parses each as a row of doubles.
 *
 * Iterates over lines until EOF, calling doubles_in_vector on each line and
 * collecting the resulting rows into a 2D vector.
 *
 * @param stm Input stream to read lines from.
 * @return A 2D vector where each inner vector contains the doubles from one line.
 */
std::vector<std::vector<double>> read_lines(std::istream& stm)
{
    std::vector<std::vector<double>> result;
    std::string line;

    while (std::getline(stm, line)) result.push_back(doubles_in_vector(line));

    return result;
}

/// @brief Checks if a string represents a number. 
/// @param str The string to check. Can be empty.
/// @return True if the string represents a number, False otherwise. 
bool isNumber(const std::string& str) {
    if (str.empty()) return false; // Handle empty strings

    for (char c : str) {
        if (!std::isdigit(c)) return false;
    }
    return true;
}

/**
 * @brief Reads a 2D source distribution file into a nested map.
 *
 * Reads the file line by line. Lines that cannot be parsed as exactly three
 * whitespace-separated floating-point numbers (x, y, value) are silently
 * skipped, allowing the file to contain an arbitrary header with units or
 * other human-readable information.
 *
 * @param fileName Path to the file containing the 2D source data.
 * @param dest     Output nested map to populate: dest[x][y] = value.
 */
void readFile2D(const std::string& fileName,
                std::map<double, std::map<double, double>>& dest)
{
    std::ifstream file(fileName);
    if (!file.is_open())
    {
        cerr << "ERROR: FILE2D NOT FOUND: " << fileName << endl;
        return;
    }

    std::string line;
    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        double x, y, val;
        if (iss >> x >> y >> val)
        {
            dest[x][y] = val;
        }
        // Lines that do not yield 3 doubles (headers, comments, etc.) are skipped.
    }
}

/**
 * @brief Reads a 3D source distribution file into a nested map.
 *
 * Reads the file line by line. Lines that cannot be parsed as exactly four
 * whitespace-separated floating-point numbers (x, y, z, value) are silently
 * skipped, allowing the file to contain an arbitrary header with units or
 * other human-readable information.
 *
 * @param fileName Path to the file containing the 3D source data.
 * @param dest     Output nested map to populate: dest[x][y][z] = value.
 */
void readFile3D(const std::string& fileName,
                std::map<double, std::map<double, std::map<double, double>>>& dest)
{
    std::ifstream file(fileName);
    if (!file.is_open())
    {
        cerr << "ERROR: FILE3D NOT FOUND: " << fileName << endl;
        return;
    }

    std::string line;
    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        double x, y, z, val;
        if (iss >> x >> y >> z >> val)
        {
            dest[x][y][z] = val;
        }
        // Lines that do not yield 4 doubles (headers, comments, etc.) are skipped.
    }
}

#ifdef HAVE_VTK
/**
 * @brief Reads a VTU file and builds the VTK objects needed for spatial interpolation.
 *
 * Reads the file using vtkXMLUnstructuredGridReader, builds a vtkCellLocator
 * for fast spatial queries, and resolves the first scalar point-data array.
 * Returns true on success, false if the file could not be read or contains
 * no scalar point data.
 *
 * @param fileName Path to the .vtu file to load.
 * @param grid     Output: the loaded unstructured grid.
 * @param locator  Output: cell locator built on the grid.
 * @param cell     Output: reusable generic-cell object for queries.
 * @param scalars  Output: pointer to the scalar data array (non-owning).
 * @return True if all objects were successfully initialised, false otherwise.
 */
bool readFileVTK(const std::string& fileName,
                 vtkSmartPointer<vtkUnstructuredGrid>& grid,
                 vtkSmartPointer<vtkCellLocator>& locator,
                 vtkSmartPointer<vtkGenericCell>& cell,
                 vtkDataArray*& scalars)
{
    auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(fileName.c_str());
    reader->Update();

    grid = reader->GetOutput();
    if (!grid || grid->GetNumberOfPoints() == 0)
    {
        cerr << "ERROR: FILE_VTK failed to read or empty grid: " << fileName << endl;
        return false;
    }

    locator = vtkSmartPointer<vtkCellLocator>::New();
    locator->SetDataSet(grid);
    locator->BuildLocator();

    cell = vtkSmartPointer<vtkGenericCell>::New();

    scalars = grid->GetPointData()->GetScalars();
    if (!scalars && grid->GetPointData()->GetNumberOfArrays() > 0)
    {
        scalars = grid->GetPointData()->GetArray(0);
    }

    if (!scalars)
    {
        cerr << "ERROR: FILE_VTK no scalar point data found in: " << fileName << endl;
        return false;
    }

    return true;
}

/**
 * @brief Interpolates a scalar value from a VTK unstructured grid at a given point.
 *
 * If the point lies inside the grid, uses VTK shape-function weights from
 * vtkCellLocator::FindCell. If the point is outside, projects it onto the
 * nearest cell face via FindClosestPoint + EvaluatePosition and interpolates
 * there. Returns 0.0 if no cell can be found at all.
 *
 * @param x       X coordinate of the query point.
 * @param y       Y coordinate of the query point.
 * @param z       Z coordinate of the query point.
 * @param locator Cell locator built on the source grid.
 * @param cell    Reusable generic-cell object (updated by the locator).
 * @param scalars Scalar point-data array to interpolate from.
 * @return The interpolated scalar value at (x, y, z).
 */
double interpolateVTK(double x, double y, double z,
                      vtkCellLocator* locator,
                      vtkGenericCell* cell,
                      vtkDataArray* scalars)
{
    double pos[3] = { x, y, z };
    double pcoords[3];
    std::vector<double> weights(VTK_CELL_SIZE);
    vtkIdType cellId = locator->FindCell(pos, 0.0, cell, pcoords, weights.data());

    if (cellId >= 0)
    {
        // Point is inside a cell: interpolate using VTK shape functions.
        vtkIdList* ptIds = cell->GetPointIds();
        double load = 0.0;
        for (vtkIdType i = 0; i < ptIds->GetNumberOfIds(); ++i)
        {
            load += weights[i] * scalars->GetComponent(ptIds->GetId(i), 0);
        }
        return load;
    }

    // Point is outside the grid: project onto the nearest cell face and
    // interpolate there.
    double closestPoint[3];
    double dist2;
    int subId;
    locator->FindClosestPoint(pos, closestPoint, cell, cellId, subId, dist2);
    if (cellId >= 0)
    {
        double distSq;
        std::vector<double> bweights(VTK_CELL_SIZE);
        cell->EvaluatePosition(closestPoint, nullptr, subId, pcoords, distSq, bweights.data());
        vtkIdList* ptIds = cell->GetPointIds();
        double load = 0.0;
        for (vtkIdType i = 0; i < ptIds->GetNumberOfIds(); ++i)
        {
            load += bweights[i] * scalars->GetComponent(ptIds->GetId(i), 0);
        }
        return load;
    }

    return 0.0;
}
#endif

}
