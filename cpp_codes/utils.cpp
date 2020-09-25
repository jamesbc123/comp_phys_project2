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
    
    string filename = "iterations.txt"; // This file saves number of iterations for
    string filename_timing = "timing.txt";

    // different matrix sizes n. It is a .csv file.
    const char * filename_c_str = filename.c_str(); // The "remove(...)" function below requires the input file name to be
    // of the type const char*.
    cout << "remove(filename_c_str) = " << remove(filename_c_str) << endl;   // Remove the file if it already exists.
    
    // Do the same for the timing file.
    const char * filename_t_c_str = filename_timing.c_str();
    cout << "remove(filename_t_c_str) = " << remove(filename_t_c_str) << endl; 

    ofile.open(filename);
    // Choose the data column titles:
    ofile << "n, " << "number_of_transformations" << endl;
    ofile.close();  // The file will be re-opened in the function write_to_file.

    ofile.open(filename_timing);
    ofile << "n, " << "time used by solver class (seconds)" << endl;
    ofile.close();  
    
    clock_t start, end;
    double timeused;

    int maxPower = 2;
    double tol = 1e-8;
    int repetition = 1;

    for(int j=0; j<repetition; j++){
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