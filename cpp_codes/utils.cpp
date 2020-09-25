#include "utils.hpp"
#include "solver.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <armadillo>
#include <time.h>
#include <stdio.h>  // remove(filename)

using namespace std;
using namespace arma;

// This file contains various (large) functions for running different parts of
// the project in order to prevent main.cpp from becoming messy. 

void run_buckling_beam(){
    ofstream ofile;

    string filename = "../results/buckling_beam/iterations.txt"; 
    string filename_timing = "../results/buckling_beam/timing.txt";

    // The "remove(...)" function below requires the input file name to be
    // of the type const char*.
    const char * filename_c_str = filename.c_str(); 

    // Remove the file if it already exists.
    cout << "remove(filename_c_str) = " << remove(filename_c_str) << endl;   
    
    // Do the same for the timing file.
    const char * filename_t_c_str = filename_timing.c_str();
    cout << "remove(filename_t_c_str) = " << remove(filename_t_c_str) << endl; 

    ofile.open(filename);
    // Choose the data column titles:
    ofile << "n, " << "number_of_transformations" << endl;
    // The file will be re-opened in the function write_to_file.
    ofile.close();  

    // Do the same for the timing file.
    ofile.open(filename_timing);
    ofile << "n,time_used,time_used_arma" << endl;
    ofile.close();  
    
    clock_t start, end;
    double timeused;
    double timeused_arma;

    // Initialise the dimensionality schedule, the number of repeats and
    // the tolerance to be reached before terminating the jacobi 
    // rotations. 
    double tol = 1e-8;
    int repetition = 1;
    int max_n = 150;

    for(int j=0; j<repetition; j++){
        for(int n=50; n< max_n; n+=50){
            double h = 1/double(n);
            double hh = h*h;
            double d = 2/hh;
            double a = -1/hh;

            mat A = zeros<mat>(n, n);
            A.diag().fill(d);
            A.diag(-1).fill(a);
            A.diag(1).fill(a);

            Solver my_solver;
            my_solver.init(n, A, tol);

            // Run the algorithm and time it.
            start = clock();
            my_solver.run();
            end = clock();
            timeused = ((double)end-(double)start)/CLOCKS_PER_SEC;

            // Run armadillo solver and time it.
            mat D = zeros<mat>(n, n);
            D = diagmat(A);
            vec eigval;
            mat eigvec;
            start = clock();
            eig_sym(eigval, eigvec, A);
            end = clock();
            timeused_arma = ((double)end-(double)start)/CLOCKS_PER_SEC;

            // Write timings to file.
            cout <<scientific<< "time used by solver class" << timeused << endl;
            ofile.open(filename_timing, ios::app);
            ofile << n << "," << timeused << "," << timeused_arma << endl;
            ofile.close();

            // Sort the eigenvectors and eigenvalues from solver class.
            my_solver.sort_eigvec_and_eigval();
            string filename_R = "../results/buckling_beam/eig_vec_"+ to_string(n) + ".txt";
            my_solver.write_to_file(filename, filename_R);  // Write 

            // Create files for saving analytic eigenvalues and vectors to.
            string filename_analytic_eigvec = "../results/buckling_beam/analytic_eig_vec" 
                                              + to_string(n) + ".txt";
            string filename_analytic_eigval = "../results/buckling_beam/analytic_eig_val" 
                                              + to_string(n) + ".txt";

            // Run analytic eigenvector solver. 
            my_solver.analytic_eigvec(filename_analytic_eigvec, filename_analytic_eigval);
        }
    }
}

