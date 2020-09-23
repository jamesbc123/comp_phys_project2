// This program is for performing small tests comparing the numerical results vs. analytical/mathematical results.

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <armadillo>
#include "unit_testing_functions.hpp"

using namespace arma;
using namespace std;

int main()
{
    //double x0 = 0; double xN = 1; // End points in x
    /*
    int n = 10;  // Number of points = N+1
    double h = (xN-x0) / (double)n; // Step width
    double hh = h*h;

    // Diagonal elements:
    double d = 2/hh;
    double a = -1/hh;
    */


    // Test the numerical (Armadillo) vs analytical eigelvalues/eigenvectors:
    /*
    int n = 10;
    test_Armadillo_vs_analytical(n);
    */

    int n = 5;
    test_finding_largest_off_diagonal(n);

    return 0;
}

