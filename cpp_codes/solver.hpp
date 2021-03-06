#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <iostream>
#include <cmath>
#include <fstream>
#include <armadillo>

class Solver{
private:
    int m_n;        // Size of matrix A
    arma::mat m_A;  // Matrix A, the system we want to solve for the eigenvalues.
    arma::mat m_R;  // Rotation matrix
    arma::mat m_sorted_R; // Sorted rotation matrix
    double m_tol;   // Tolerance of the sum of the off-diagonal elements
    int m_max_iter; // Maximum accepted number of iterations
    int m_k;        // First matrix index k
    int m_p;        // Second matrix index p
    double m_tau;   // tau = cot(2*theta)
    bool m_tol_reached;
    int m_i;
    arma::vec m_eigval;
    arma::vec m_sorted_eigval;
    std::ofstream m_ofile;

public:
    void init(int n, arma::mat A, double tol);
    double max_off_diag();
    void calc_tau();
    void rotate();
    void run();
    void sort_eigvec_and_eigval();
    void write_to_file(std::string filename_iter, std::string filename_num_eigvec, std::string filename_num_eigval);
    void analytic_eigvec(std:: string filename_eigvec, std::string filename_eigval);
    arma::mat get_R(); // Solver::run() should be ran first in order
    // to get the solution.
    arma::mat get_sorted_R();
    arma::vec get_eigenvalues();
    arma::vec get_sorted_eigenvalues();
};
#endif
