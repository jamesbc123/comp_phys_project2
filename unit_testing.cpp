// This program is for performing small tests comparing the numerical results vs. analytical/mathematical results.

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <armadillo>

using namespace arma;
using namespace std;

int main()
{
    double x0 = 0; double xN = 1; // End points in x
    int N = 10;  // Number of points = N+1
    double h = (xN-x0) / (double)N; // Step width
    double hh = h*h;

    // Diagonal elements:
    double d = 2/hh;
    double a = -1/hh;

    cout << setprecision(3); // Print formatting 
    mat A = zeros<mat>(N-1, N-1);
    //Mat<int> A = zeros<mat>(N-1, N-1);

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
    vec lambdaVec = vec(N-1);
    for(int i=0; i<=N-2; i++){
        lambdaVec(i) = d + 2*a*cos((i+1)*M_PI / (double)N);
    }

    // Analytical eigenvectors:
    mat u = mat(N-1, N-1);   
    for(int i=0; i<=N-2; i++){
        for(int j=0; j<=N-2; j++){
            u(i,j) = sin((i+1)*(j+1)*M_PI / (double)N);  // Adding 1 to i and j because of C++'s 0-indexing (in order to match the analytical equation).
        }
    }

    lambdaVec.raw_print(cout, "\nlambdaVec (eigenvalues from analytical expression): ");
    u.raw_print(cout, "\nu (eigenvectors from analytical expression): ");

   //   ### (Hvorfor er det analytiske uttrykket annerledes enn Armadillo-uttrykket? Ovenfor er det allerede vist at egenvektorene
   // og egenvektorene stemmer for Armadillo-løseren, så da må det være det analytiske uttrykket som er feil. Hvor er feilen?)

    return 0;
}

