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


void run_q_dots_one_electron(){
    ofstream ofile;
    
    string filename_N = "2d_eigenvalues_vs_N.txt";
    string filename_rhoMax = "2d_eigenvalues_vs_rhoMax.txt";

    // Choose some different values for rho_max to run the code for:
    vec rhoMaxList = vec("100 1000 10000 100000");
    vec NList = vec("10 20 30 40 50 60 80 100 125 150 175 200");
    rhoMaxList.print();

    // Remove the files if they already exists:
    const char * filename_N_c_str = filename_N.c_str(); // The "remove(...)" function below requires the input file name to be
    // of the type const char*.
    cout << "remove(filename_N_c_str) = " << remove(filename_N_c_str) << endl; 

    const char * filename_rhoMax_c_str = filename_rhoMax.c_str();
    cout << "remove(filename_rhoMax_c_str) = " << remove(filename_rhoMax_c_str) << endl;


    ofile.open(filename_N);
    // Choose the data column titles:
    ofile << "n, " << "eigenvalues" << endl;
    ofile.close();  // The file will be re-opened in the function write_to_file.

    ofile.open(filename_rhoMax);
    ofile << "rho_max, " << "eigenvalues" << endl;
    ofile.close();  
    
    clock_t start, end;
    double timeused;

    //int maxPower = 2;
    double tol = 1e-8;
    int repetition = 1;

    // Run algorithm as function of n.
    for(int i=1; i <= maxPower; i++){
        int n = pow(10, i);
        mat A = zeros<mat>(n, n);
        A.diag().fill(2);
        A.diag(-1).fill(-1);
        A.diag(1).fill(-1);
        
        Solver my_solver;
        my_solver.init(n, A, tol);

        // Run the algorithm and time it:
        start = clock();
        my_solver.run();
        end = clock();
        timeused = ((double)end-(double)start)/CLOCKS_PER_SEC;
        cout <<scientific<< "time used" << timeused << endl;
        ofile.open(filename_timing, ios::app);
        ofile << "\n" << n << "," << timeused << endl;
        ofile.close();

        my_solver.sort_eigvec_and_eigval();
        string filename_R = "eig_vec_"+ to_string(n) + ".txt";
        my_solver.write_to_file(filename, filename_R);  // Write 
    }

    // Run algorithm as function of rhoMax.
}