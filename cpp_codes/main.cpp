#include "solver.hpp"
#include <armadillo>
#include <time.h>
#include <fstream>
#include <stdio.h>  // remove(filename)

using namespace std;
using namespace arma;

ofstream ofile;

int main(){
    string filename = "iterations.txt"; // This file saves number of iterations for
    string filename_timing = "timing.txt";

    // different matrix sizes n. It is a .csv file.
    const char * filename_c_str = filename.c_str(); // The "remove(...)" function below requires the input file name to be
    // of the type const char*.
    cout << "remove(filename_c_str) = " << remove(filename_c_str) << endl;   // Remove the file if it already exists.
    
    ofile.open(filename);
    // Choose the data column titles:
    ofile << "n, " << "number_of_transformations" << endl;
    ofile.close();  // The file will be re-opened in the function write_to_file.
    
    clock_t start, end;
    double timeused;

    int pow = 5;
    double tol = 1e-8;

    for(int i=1; i < pow; i++){
        int n = 10^pow;
        mat A = zeros<mat>(n, n);
        A.diag().fill(2);
        A.diag(-1).fill(-1);
        A.diag(1).fill(-1);
        
        Solver my_solver;
        cout << "n (before .init()): " << n << endl;
        my_solver.init(n, A, tol);

        // Run the algorithm and time it:
        start = clock();
        my_solver.run();
        end = clock();
        timeused = (end-start)/CLOCKS_PER_SEC;

        ofile.open(filename_timing);
        ofile << "n = " << n << ", time used by my_solver.run() (seconds): " << timeused << endl;
        ofile.close();

        my_solver.write_to_file(filename);  // Write 
    return 0;
}