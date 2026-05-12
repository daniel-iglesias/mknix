/***************************************************************************
 *   Copyright (C) 2013 by Daniel Iglesias                                 *
 *   https://github.com/daniel-iglesias/mknix                                        *
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <simulation/mknixWrapper.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "LMX/lmx.h"

using namespace std;

std::vector<std::vector<double> > transpose(const std::vector<std::vector<double> >& data) {
    // this assumes that all inner vectors have the same size and
    // allocates space for the complete result in advance
    std::vector<std::vector<double> > result(data[0].size(),
                                          std::vector<double>(data.size()));
    for (std::vector<double>::size_type i = 0; i < data[0].size(); i++) 
        for (std::vector<double>::size_type j = 0; j < data.size(); j++) {
            result[i][j] = data[j][i];
        }
    return result;
}

// Function to split a string by a delimiter and convert to double
std::vector<double> split_to_double(const std::string& str, char delimiter) {
    std::vector<double> tokens;
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, delimiter)) {
        try {
            tokens.push_back(std::stod(token));
        } catch (const std::invalid_argument& e) {
            tokens.push_back(0.0);  // Handle non-numeric values as 0.0
        }
    }
    return tokens;
}


int readInputFile( const std::string & input_file_name, std::vector<std::vector<double>> & data, int & size ) 
{
    std::ifstream file(input_file_name);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << input_file_name << std::endl;
        return 1;
    }

    std::string line;
    std::vector<std::string> headers;

    // Read the header line
    if (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string header;
        while (std::getline(ss, header, ',')) {
            headers.push_back(header);
        }
        size = headers.size();
    }

    // Read the data lines
    while (std::getline(file, line)) {
        data.push_back(split_to_double(line, ','));
    }

    file.close();

    // Display the headers
    std::cout << "Headers:" << std::endl;
    for (const auto& header : headers) {
        std::cout << header << " ";
    }
    std::cout << std::endl << std::endl;

    // Display the data
    std::cout << "Data:" << std::endl;
    for (const auto& row : data) {
        for (const auto& value : row) {
            std::cout << value << " ";
        }
        std::cout << std::endl;
    }
    return EXIT_SUCCESS;
}

void printLastVector(const std::vector<double*> &signals, int last_vector_size) {
    if (signals.empty()) {
        std::cerr << "Error: signals vector is empty." << std::endl;
        return;
    }

    // Get the pointer to the last vector
    double* last_vector = signals.back();

    // Print each element
    std::cout << "Last vector contents: ";
    for (int i = 0; i < last_vector_size; ++i) {
        std::cout << last_vector[i] << " ";
    }
    std::cout << std::endl;
}

int main(int argc, char * argv[])
{
        double initialTemperature = 100;
        int output_files_detail = 2;
                // 0 none
                // 1 only times and output.reader
                // 2 all output files
//   try{
    if (argc >= 3 && argc < 5){
        // cout << "EXECUTING: $ mknixwrapper " << argv[1] 
        //      << " " << argv[2] << endl;
        // system("echo -n '1. Current Directory is '; pwd");
    if (argc > 3){
        lmx::setMatrixType(atoi(argv[3]));
        if (argc == 5){ lmx::setLinSolverType( atoi(argv[4]) ); }
    }
        std::vector<std::vector<double>> data;
        int column_size ; 
        std::string filename(argv[2]);
        readInputFile(filename, data, column_size);
        int steps = (column_size-1)/2; // assuming each signal represents a step
        cout << "Number of signals read = " << data.size() << endl;
        cout << "Number of steps read = " << steps << endl;

        std::vector<double*> temperatures;
        std::vector<double*> signals;
        std::vector<std::vector<double>> transposed_data = transpose(data);

        mknix::MknixWrapper simulationWrapper(argv[1], output_files_detail);
        simulationWrapper.init( initialTemperature );
        for(int i = 0; i < steps; ++i){
            temperatures.push_back( new double[4] );
            signals.push_back( transposed_data[2*i+2].data() );
            cout << "Temperature: " << endl;
            cout << "Heat signal: " << endl;
            printLastVector( signals, data.size() );
            simulationWrapper.run( signals.back(), temperatures.back() );
            printLastVector( temperatures, 4 );
        }
        // Cleaning up: delete allocated memory for output temperatures
        for (auto& ptr : temperatures) {
            delete[] ptr;
        } 
    }
    else if (argc < 3){
        cout << "Need at least two parameters: " << endl
             << "1. Name of the mknix input file" << endl
             << "2. Signals (input heat) file name" << endl
             << "3. (Optional) Matrix type (0 for dense, 1 for sparse)" << endl
             << "4. (Optional) Linear solver type (0 for direct, 1 for iterative)" << endl
             << "Example: $./mknixwrapper \"input_2triangles.fem.mknix\"  0 0" << endl;
    }
//   } catch (const std::exception& e) {
//     cerr << "Error: " << e.what() << endl;
//     return EXIT_FAILURE;
//   }
    return EXIT_SUCCESS;
}
