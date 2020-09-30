#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <armadillo>
#include <time.h>
#include <stdio.h>  // remove(filename)
#include <sstream>  // stringstream objects
#include "utils.hpp"
#include "solver.hpp"

using namespace std;
using namespace arma;

// This file contains various (large) functions for running different parts of
// the project in order to prevent main.cpp from becoming messy. 

void run_buckling_beam(){
    ofstream ofile;

    string filename = "../results/buckling_beam/iterations.txt"; 
    string filename_timing = "../results/buckling_beam/timing.txt";

    // Remove the files if they already exists:
    // The "remove(...)" function below requires the input file name to be
    // of the type const char*.
    const char * filename_c_str = filename.c_str(); 
    remove(filename_c_str);   
    const char * filename_t_c_str = filename_timing.c_str();
    remove(filename_t_c_str); 

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

    // Initialise the dimensionality schedule, the number of 
    // repeats(to get the average time) and the tolerance 
    // to be reached before terminating the jacobi rotations. 
    double tol = 1e-8;
    int repetition = 1;
    int max_n = 150;

    for(int j=0; j<repetition; j++){
        for(int n=50; n<=max_n; n+=50){
            cout << endl;
            double h = 1/double(n);
            double hh = h*h;
            double d = 2/hh;
            double a = -1/hh;

            mat A = zeros<mat>(n-1, n-1);
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
            cout << scientific << "Time used by solver class: " << timeused << endl;
            ofile.open(filename_timing, ios::app);
            ofile << n << "," << timeused << "," << timeused_arma << endl;
            ofile.close();

            // Sort the eigenvectors and eigenvalues from solver class.
            my_solver.sort_eigvec_and_eigval();

            // Create files for the numerical eigenvectors and eigenvalues.
            string filename_num_eigvec = "../results/buckling_beam/num_eigvec_"
                                         + to_string(n) + ".txt";
            string filename_num_eigval = "../results/buckling_beam/num_eigval_"
                                         + to_string(n) + ".txt";
            my_solver.write_to_file(filename, filename_num_eigvec, filename_num_eigval);  // Write            

            // Create files for saving analytic eigenvalues and vectors to.
            string filename_analytic_eigvec = "../results/buckling_beam/analytic_eigvec_" 
                                              + to_string(n) + ".txt";
            string filename_analytic_eigval = "../results/buckling_beam/analytic_eigval_" 
                                              + to_string(n) + ".txt";

            // Run analytic eigenvector solver. 
            my_solver.analytic_eigvec(filename_analytic_eigvec, filename_analytic_eigval);
        }
    }
}

void run_q_dots_one_electron(){
    // Choose some different values for n and rho_max to run the code for:
    //vec nList = vec("10 20 30 40 50 60 80 100 125 150 175 200");
    vec nList = vec("10 20 30 50 100 200 300 400 500");
    vec rhoMaxList = vec("1 5 10");
    
    nList.print("nList:");
    rhoMaxList.print("rhoMaxList:");
    
    ofstream ofile;
    string filename_n = "2d_eigenvalues_vs_n.txt";
    string filename_rhoMax = "2d_eigenvalues_vs_rhoMax.txt";

    // Remove the files if they already exists:
    const char * filename_n_c_str = filename_n.c_str(); // The "remove(...)" function below requires the input file name to be
    // of the type const char*.
    cout << "remove(filename_N_c_str) = " << remove(filename_n_c_str) << endl; 

    const char * filename_rhoMax_c_str = filename_rhoMax.c_str();
    cout << "remove(filename_rhoMax_c_str) = " << remove(filename_rhoMax_c_str) << endl;

    // Create the text files to contain the results:
    ofile.open(filename_n);
    // Choose the data column titles:
    ofile << "n, " << "eigenvalues" << endl;
    ofile.close();  // The file will be re-opened in the function write_to_file.

    ofile.open(filename_rhoMax);
    ofile << "rho_max, " << "eigenvalues" << endl;
    ofile.close();  
    
    clock_t start, end;
    double timeused;

    double tol = 1e-8;

    // Run algorithm for different n:
    cout << "Calculating eigenvalues for different values of n...\n";
    // Choose a max value for rho (rho_max, approximation for rho=infinity):
    double rho_max = 7;
    double rho_0 = 0;
    double rho_n = rho_max; // rho_n = rho_max
    cout << "rho_max: " << rho_max << endl;
    for (auto n : nList){ // For all chosen values of n.
        cout << "n: " << n << endl;
        mat A = zeros<mat>(n-1, n-1);
        // In the quantum case we must add the potential V(rho). It contains different
        // values along the diagonal.

        //double h = 1/double(n);
        double h = (rho_n - rho_0)/double(n);
        cout << "h = " << h << endl;
        double hh = h*h;
        double d = 2/hh;
        double a = -1/hh; // Same as buckling beam case.

        vec rhoList(n-1); // Position list of length n-1 (end points rho_0 and rho_n are
        // excluded).
        for (int i=0; i<=n-2; i++){rhoList(i) = rho_0 + (i+1)*h;} // Fill rhoList, and 
        // exclude end points.
        //rhoList.print("rhoList: ");

        // Off-diagonal elements:
        A.diag(-1).fill(a);
        A.diag(1).fill(a);

        double V_i;
        double rho_i;
        for (int i=0; i<=n-2; i++){ // Fill the diagonal elements of the matrix.
        // In the quantum case, these elements vary because of the potential V_i.
            V_i = rhoList(i)*rhoList(i); // V_i = rho_i^2
            A(i,i) = d + V_i; // The diagonal elements (d_i in the project text).
        }
        // Now the matrix A represents the QM equation, and all that is needed is to 
        // solve that system of equations.

        Solver my_solver;
        my_solver.init(n, A, tol);

        // Run the algorithm and time it:
        start = clock();
        my_solver.run();
        end = clock();
        timeused = ((double)end-(double)start)/CLOCKS_PER_SEC;
        // Sort the eigenvalues in ascending order:
        my_solver.sort_eigvec_and_eigval();

        // Print the eigenvalues:
        vec eigenVals = my_solver.get_sorted_eigenvalues();
        cout << "eigenVals.n_elem: " << eigenVals.n_elem;
        eigenVals(span(0,8)).print("\neigenVals (first 9 elements):"); // Print the first 9 eigenvalues.
        cout << "\t.\n\t.\n\t.\n\n";
        
        /*
        cout <<scientific<< "time used" << timeused << endl;
        ofile.open(filename_timing, ios::app);
        ofile << "\n" << n << "," << timeused << endl;
        ofile.close();

        my_solver.sort_eigvec_and_eigval();
        string filename_R = "eig_vec_"+ to_string(n) + ".txt";
        my_solver.write_to_file(filename, filename_R);  // Write 
        */
    }
    cout << "Calculations of eigenvalues for different values of n complete.\n\n";

    /*
    // Run algorithm as function of rhoMax:
    cout << "Calculating eigenvalues for different values of rho_max...\n";
    // Choose one value of n for all of these rho_max values:
    int n = 100;
    cout << "n: " << n << endl;
    for (auto rhoMax : rhoMaxList){ // For all chosen values of rho_max.
        cout << endl << "rhoMax: " << rhoMax << endl;
        mat A = zeros<mat>(n-1, n-1);
        // In the quantum case we must add the potential V(rho). It contains different
        // values along the diagonal.

        double rho_0 = 0;
        double rho_n = rhoMax; // rho_n = rho_max

        //double h = 1/double(n);
        double h = (rho_n - rho_0)/double(n);
        cout << "h = " << h << endl;
        double hh = h*h;
        double d = 2/hh;
        double a = -1/hh; // Same as buckling beam case.
        
        vec rhoList(n-1); // Position list of length n-1 (end points rho_0 and rho_n are
        // excluded).
        for (int i=0; i<=n-2; i++){rhoList(i) = rho_0 + (i+1)*h;} // Fill rhoList, and 
        // exclude end points.
        //rhoList.print("rhoList: ");

        // Off-diagonal elements:
        A.diag(-1).fill(a);
        A.diag(1).fill(a);
        // Diagonal elements:
        double V_i;
        double rho_i;
        for (int i=0; i<=n-2; i++){ // Fill the diagonal elements of the matrix.
            V_i = rhoList(i)*rhoList(i); // V_i = rho_i^2.
            A(i,i) = d + V_i; // The diagonal elements (d_i in the project text).
        }

        Solver my_solver;
        my_solver.init(n, A, tol);

        // Run the algorithm and time it:
        start = clock();
        my_solver.run();
        end = clock();
        timeused = ((double)end-(double)start)/CLOCKS_PER_SEC;
        // Sort the eigenvalues in ascending order:
        my_solver.sort_eigvec_and_eigval();

        // Print the eigenvalues:
        vec eigenVals = my_solver.get_sorted_eigenvalues();
        cout << "eigenVals.n_elem: " << eigenVals.n_elem;
        //eigenVals.print("\neigenVals:");
        eigenVals(span(0,8)).print("\neigenVals (first 9 elements):"); // Print the first 9 eigenvalues.
        //cout << ".\n.\n.\n\n";
        cout << "\t.\n\t.\n\t.\n\n";
    }
    cout << "Calculations of eigenvalues for different values of rho_max complete.\n";
    */
}

void run_q_dots_two_electrons(){  
    ofstream ofile;
    
    clock_t start, end;
    double timeused;
    double tol = 1e-8;
    vec omegaList = vec("0.01 0.5 1 5"); // The values of omega_r specified in
    // the project text.

    // Run algorithm for different n:
    cout << "Calculating eigenvalues for different values of n...\n";
    // Choose a max value for rho (rho_max, approximation for rho=infinity):
    double rho_max = 7;
    double rho_0 = 0;
    double rho_n = rho_max; // rho_n = rho_max
    double n = 200;
    for (auto omega : omegaList){ // For all chosen values of omega_r
        cout << "omega_r: " << omega << endl;

        double omegaSq = omega*omega;  // omega_r^2
        mat A = zeros<mat>(n-1, n-1);
        // In the quantum case we must add the potential V(rho). It contains different
        // values along the diagonal.

        //double h = 1/double(n);
        double h = (rho_n - rho_0)/double(n);
        cout << "h = " << h << endl;
        double hh = h*h;
        double d = 2/hh;
        double a = -1/hh; // Same as buckling beam case.

        vec rhoList(n-1); // Position list of length n-1 (end points rho_0 and rho_n are
        // excluded).
        for (int i=0; i<=n-2; i++){rhoList(i) = rho_0 + (i+1)*h;} // Fill rhoList, and 
        // exclude end points.

        // Off-diagonal elements:
        A.diag(-1).fill(a);
        A.diag(1).fill(a);

        double rho_i;
        double rhoSq_i; // (rho_i)^2
        for (int i=0; i<=n-2; i++){ // Fill the diagonal elements of the matrix.
        // In the quantum case, these elements vary because of the potential V_i.
            rho_i = rhoList(i);
            rhoSq_i = rho_i*rho_i;
            A(i,i) = d + (omegaSq*rhoSq_i + 1/rho_i); // The diagonal elements.
        }
        // Now the matrix A represents the QM equation for two electrons, and all
        // that needs to be done now is to solve that system of equations.

        Solver my_solver;
        my_solver.init(n, A, tol);

        // Run the algorithm and time it:
        start = clock();
        my_solver.run();
        end = clock();
        timeused = ((double)end-(double)start)/CLOCKS_PER_SEC;
        // Sort the eigenvalues in ascending order:
        my_solver.sort_eigvec_and_eigval();

        // Print the eigenvalues:
        vec eigenVals = my_solver.get_sorted_eigenvalues();
        cout << "eigenVals.n_elem: " << eigenVals.n_elem;
        eigenVals(span(0,8)).print("\neigenVals (first 9 elements):"); // Print the first 9 eigenvalues.
        cout << "\t.\n\t.\n\t.\n\n";
        
        // Store the results in vectors:
        double eigenvalueGS = my_solver.get_sorted_eigenvalues()(0);  // Get the ground state eigenvalue.
        vec eigenvectorGS = my_solver.get_sorted_R().col(0); // Get the ground state eigenvector.
        

        // ******** Write the results to file: ********
        // Create one file for each value of omega_r.
        // Each file name includes the omega_r value:
        stringstream omegaStringstream; omegaStringstream << setprecision(2)
        << fixed << omega;
        string omegaStr = omegaStringstream.str();
        string directory = "../results/problem_2e/"; // Save the files in the results folder.
        string filename = directory + "2e_eigvecGS_omega" + omegaStr + ".txt";  // "GS" stands for ground state.
        cout << "filename: " << filename << endl;

        // Remove the files if they already exists:
        const char * filename_c_str = filename.c_str(); // The "remove(...)" function below requires the input file name to be
        // of the type const char*.
        cout << "remove(filename_c_str) = " << remove(filename_c_str) << endl;
        ofile.open(filename);//, ios::app);
        // Choose the data column titles:
        ofile << "rho, " << "eigenvector (ground state)";
        ofile.close();  // The file will be re-opened in the function write_to_file.

        // Write the data to file:
        ofile.open(filename, ios::app);
        for(int i = 0; i <= n-2; i++){
            ofile << endl << rhoList(i) << ", " << eigenvectorGS(i);
        }
        ofile.close();
    }
}
