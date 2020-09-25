#include "solver.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
#include <armadillo>

using namespace arma;
using namespace std;


void Solver::init(int n, mat A, double tol){
    /* Initialise the member variables in the class */
    m_n = n;
    m_A = A;
    m_R = zeros<mat>(m_n-1, m_n-1);
    m_R.diag().fill(1.0);
    m_eigval;
    m_tol = tol;
    m_max_iter = (m_n-1)*(m_n-1)*(m_n-1);
    m_tol_reached = false;
}

double Solver::max_off_diag(){
    /* Find the max element in A and change the 
    indices of p and k to those found. */
    double max = 0.0;
    for (int i = 0; i < m_n-1; ++i){
        for (int j = i+1; j < m_n-1; ++j){
            double aij = fabs(m_A(i,j));
            if ( aij > max){
                max = aij; m_k = i; m_p = j;
            }
        }
    }
    if(max <= m_tol){
        m_tol_reached = true;
    }
    return max;
}

void Solver::calc_tau(){
    /* Calculate tau */
    m_tau = (m_A(m_p,m_p) - m_A(m_k,m_k))/(2*m_A(m_k,m_p));
}

void Solver::rotate(){
    /* Perform Jacobi method by making similarity transformations. */
    double s, c;
    if(m_A(m_k, m_p) != 0.0){
        double t;
        Solver::calc_tau(); 
        if(m_tau >= 0){
            t = 1.0/(m_tau + sqrt(1.0 + m_tau*m_tau));
        }else{
            t = -1.0/(-m_tau +sqrt(1.0 + m_tau*m_tau));
        }
        c = 1/sqrt(1+t*t);
        s = c*t;
    }
    else{
        c = 1.0;
        s = 0.0;
    }

    double a_kk, a_pp, a_ik, a_ip, r_ik, r_ip;
    a_kk = m_A(m_k,m_k);
    a_pp = m_A(m_p,m_p);
    m_A(m_k,m_k) = c*c*a_kk - 2.0*c*s*m_A(m_k,m_p) + s*s*a_pp;
    m_A(m_p,m_p) = s*s*a_kk + 2.0*c*s*m_A(m_k,m_p) + c*c*a_pp;
    m_A(m_k,m_p) = 0.0; // Hard-coding non-diagonal elements by hand.
    m_A(m_p,m_k) = 0.0;
    for(int i = 0; i < m_n-1; i++ ){
        if ( i != m_k && i != m_p ) {
            a_ik = m_A(i,m_k);
            a_ip = m_A(i,m_p);
            m_A(i,m_k) = c*a_ik - s*a_ip;
            m_A(m_k,i) = m_A(i,m_k);
            m_A(i,m_p) = c*a_ip + s*a_ik;
            m_A(m_p,i) = m_A(i,m_p);
        }
        // And finally the new eigenvectors
        r_ik = m_R(i,m_k);
        r_ip = m_R(i,m_p);
        m_R(i,m_k) = c*r_ik - s*r_ip;
        m_R(i,m_p) = c*r_ip + s*r_ik;
    }
    return;
}

void Solver::run(){
    // This function has no output, but effectively it does since
    // it changes m_A to the solved (diagonalized) version.
    // And the rows of m_R contains the eigenvectors.
    for(m_i = 0; m_i < m_max_iter; m_i++){
            Solver::max_off_diag();
            if(m_tol_reached == true){
                cout << "Tolerance limit reached; sufficiently good results obtained.\n";
                return;
            }
            Solver::calc_tau();
            Solver::rotate();
            
            // Update m_eigval to contain the found eigenvalues:
            m_eigval = m_A.diag();
    }

    // If we have reached here then we are at max iter.
    cout << "Maximum number of accepted iterations reached; aborting run.\n";
}

/*
void Solver::write_to_file(string filename, string filename_R){
    // Write the information to file.
    
    // Columns in the text file: n, number_of_transformations
    ofstream m_ofile;
    
    // Append the data to the file.
    m_ofile.open(filename, ios::app);
    m_ofile << m_n << "," << m_i <<endl;  
    m_ofile.close();

    m_ofile.open(filename_R);
    m_ofile << m_R << endl;
    m_ofile.close();
}
*/
void Solver::write_to_file(string filename_iter, string filename_num_eigvec, string filename_num_eigval){
    /* Write the information to file */
    
    // Columns in the text file: n, number_of_transformations
    ofstream m_ofile;

    cout << "m_n (inside write_to_file): " << m_n << endl;
    cout << "m_i (inside write_to_file): " << m_i << endl;
    
    // Append the data to the file.
    m_ofile.open(filename_iter, ios::app);
    m_ofile << m_n << "," << m_i <<endl;  
    m_ofile.close();

    m_ofile.open(filename_num_eigvec);
    m_ofile << m_R << endl;
    m_ofile.close();

    m_ofile.open(filename_num_eigval);
    m_ofile << m_eigval << endl;
    m_ofile.close();
}


void Solver::sort_eigvec_and_eigval(){
    // Sort the eigenvaules by value and the eigenvectors by this ordering.
    vec eigval = m_A.diag();

    uvec indices = sort_index(eigval, "ascend");
    sort (eigval.begin(), eigval.begin()+m_n-1);
    mat sorted_R = zeros<mat>(m_n-1,m_n-1);

    for (int i=0; i<m_n-1; i++){
        sorted_R.row(i) = m_R.row(indices(i));
    }

    // Change m_R to the new sorted R.
    m_R = sorted_R;
}

/*
void Solver::sort_eigvec_and_eigval(){
    // Sort the eigenvaules by value and the eigenvectors by this ordering.
    m_eigval = m_A.diag();

    uvec indices = sort_index(m_eigval, "ascend");
    sort (m_eigval.begin(), m_eigval.begin()+m_n);
    mat sorted_R = zeros<mat>(m_n,m_n);

    for (int i=0; i<m_n; i++){
        sorted_R.row(i) = m_R.row(indices(i));
    }

    // Change m_R to the new sorted R.
    m_R = sorted_R;
}
*/

void Solver::analytic_eigvec(string filename_eigvec, string filename_eigval){
    /* Calculate the analytic eigenvectors and write them to file.
    
    */
    double h = 1/double(m_n);
    double hh = h*h;
    double d = 2/hh;
    double a = -1/hh;

    vec eigval_a = vec(m_n-1);
    for(int i=0; i<=m_n-2; i++){
        eigval_a(i) = d + 2*a*cos((i+1)*M_PI / (double)m_n);
    }

    // Analytical eigenvectors:
    mat u = mat(m_n-1, m_n-1);   
    for(int i=0; i<=m_n-2; i++){
        for(int j=0; j<=m_n-2; j++){
            // Adding 1 to i and j because of C++'s 0-indexing (in order to match the analytical equation).
            u(i,j) = sin((i+1)*(j+1)*M_PI / (double)m_n);  
            
        }
    }

    // Write eigenvalues and eigenvectors to file.
    m_ofile.open(filename_eigvec);
    m_ofile << u << endl;
    m_ofile.close();

    m_ofile.open(filename_eigval);
    m_ofile << eigval_a << endl;
    m_ofile.close();
}

arma::mat Solver::get_R(){return m_R;} // Solver::run() should be ran first. This
// returns the final rotation matrix.

arma::vec Solver::get_eigenvalues(){return m_A.diag();} // Solver::run() should be ran first.
// The diagonal elements are the eigenvalues.
