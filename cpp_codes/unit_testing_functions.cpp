#include "unit_testing_functions.hpp"
#include "solver.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <armadillo>

using namespace arma;
using namespace std;

void test_Armadillo_vs_analytical(int n) {
    // This function finds the eigenvalues and eigenvectors of an nxn Toeplitz matrix
    // first numerically using Armadillo and then from the analytical expression.
    // The results are then compared, and the eigenvectors and eigenvalues are tested
    // through matrix multiplication.
    double x0 = 0;
    double xn = 1;
    double h = (xn-x0) / (double)n; // Step width
    double hh = h*h;

    // Diagonal elements:
    double d = 2/hh;
    double a = -1/hh;

    cout << setprecision(3); // Print formatting 
    mat A = zeros<mat>(n-1, n-1);

    A.diag().fill(d);
    A.diag(-1).fill(a);
    A.diag(1).fill(a);
    A.raw_print(cout, "\nA: ");

    // Get the eigenvalues using armadillo:
    vec eigenValues;
    mat eigenVectors;
    eig_sym(eigenValues, eigenVectors, A); // A is a symmetric matrix, so we can use eig_sym to get the eigenvalues.

    //cout << fixed << setprecision(3);
    eigenValues.raw_print(cout, "\neigenValues (from Armadillo solver): ");
    eigenVectors.raw_print(cout, "\neigenVectors (from Armadillo solver): ");
    
    /*
    cout << setprecision(2) << "eigenValues: " << endl << eigenValues << endl;
    //cout << fixed << setprecision(2);
    //eigenVectors.raw_print(cout, "eigenVectors:");
    cout << "eigenVectors: " << endl << eigenVectors << endl << endl;
    */
    
    cout << "\nTest Armadillo eigenvalue and eigenvector: " << endl;
    vec eigVec1 = eigenVectors.col(0);
    cout << "eigVec1: " << endl << eigVec1 << endl;

    double eigVal1 = eigenValues[0];
    cout << "eigVal1: " << eigVal1 << endl;
    cout << "eigVal1*eigVec1: " << endl << eigVal1*eigVec1 << endl;
    cout << "A*eigVec1: " << endl << A*eigVec1 << endl;
    
    

    // Now, print the analytical solutions for the eigenvalues and eigenvectors:
    // Analytical eigenvalues:
    vec lambdaVec = vec(n-1);
    for(int i=0; i<=n-2; i++){
        lambdaVec(i) = d + 2*a*cos((i+1)*M_PI / (double)n);
    }

    // Analytical eigenvectors:
    mat u = mat(n-1, n-1);   
    for(int i=0; i<=n-2; i++){
        for(int j=0; j<=n-2; j++){
            u(i,j) = sin((i+1)*(j+1)*M_PI / (double)n);  // Adding 1 to i and j because of C++'s 0-indexing (in order to match the analytical equation).
        }
    }

    lambdaVec.raw_print(cout, "\nlambdaVec (eigenvalues from analytical expression): ");
    u.raw_print(cout, "\nu (eigenvectors from analytical expression): ");

    // Test one analytical eigenvector and its corresponding eigenvalue:
    cout << "\nTest analytical eigenvalue and eigenvector: " << endl;
    vec eigVec1_an = u.col(0);
    cout << "eigVec1_an: " << endl << eigVec1_an << endl;

    double eigVal1_an = lambdaVec[0];
    cout << "eigVal1_an: " << eigVal1_an << endl;
    cout << "eigVal1_an*eigVec1_an: " << endl << eigVal1_an*eigVec1_an << endl;
    cout << "A*eigVec1_an: " << endl << A*eigVec1_an << endl;

    
    cout << "eigVec1 = " << endl << eigVec1 << endl;
    cout << "2.239130435*eigVec1 = " << endl << 2.235890014*eigVec1 << endl;
}

void test_finding_largest_off_diagonal(int n){
    // This function tests the function Solver::max_off_diag().

    // Create a matrix with random off-diagonal elements:
    //mat A = randi<mat>(n, n, distr_param(0, 100));
    arma_rng::set_seed_random();  // Set the seed to a random value
    mat A = randu<mat>(n, n)*100;
    //A.raw_print(cout, "\nA: ");
    A = trimatu(A, 1);  // Make A an upper triangular matrix, with the
    // diagonal as zero as well.
    A.print();

    Solver test_solver;
    double tol = 1e-8;  // Just set some tolerance limit. Not used here, though.
    test_solver.init(n, A, tol);

    // Find the largest element:
    double max = test_solver.max_off_diag();
    cout << "\ntest_solve.max_off_diag() output: " << max << "\n\n";
    // Visually inspect the matrix A to check if max_off_diag() was correct.
    // Do this multiple times in order to check that the function works.
}

void test_orthogonal_eigenvectors(int n){
    // This function checks if the eigenvectors output from Solver::run()
    // are orthogonal, as they should be for a symmetric matrix.
    //int n = 10;
    mat A = zeros<mat>(n, n);
    A.diag().fill(2);
    A.diag(-1).fill(-1);
    A.diag(1).fill(-1);
    
    Solver test_solver;
    double tol = 1e-8;
    test_solver.init(n, A, tol);
    test_solver.run();
    mat R = test_solver.get_R();  // Get the eigenvector solutions (along the rows of R).
    mat vi, vj; // (Different) eigenvectors from Jacobi's algorithm.

    cout << "\ntest_orthogonal_eigenvectors(): \n";
    mat dotProducts(n,n);   // An nxn matrix containing a table of dot products.
    for (int i=0; i<=n-1; i++){
        //for (int j=i; j<=n-1; j++){
        for (int j=0; j<=n-1; j++){
            //vi = R.row(i).t();  // .t(): Transposed.
            //vj = R.row(j).t();
            vi = R.row(i);
            vj = R.row(j);
            double dotProd = dot(vi,vj);
            dotProducts(i,j) = dotProd; // Add the dot product to the table.
            cout << "\nv" << i << " dot v" << j << " = " << dotProd << endl;
        }
    }
    cout << "dotProducts.print(): \n";
    cout << setprecision(3);
    dotProducts.raw_print(cout);

}