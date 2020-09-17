#include <iostream>
#include <cmath>
#include <fstream>
#include <armadillo>
#include "solver.hpp"

using namespace arma;
using namespace std;

void Solver::init(int n, mat A, double tol){
    /* Initialise the member variables needed in the class
    */
    int m_n = n;
    mat m_A = A;
    mat m_R;
    double m_tol = tol;
    int m_k;
    int m_k;
    double m_tau;
    int m_max_iter = n*n*n;
    bool m_tol_reached = False;
    ofstream m_ofile;
    int m_i;
}

void Solver::max_off_diag(){
    /* Find the max element in A and change the 
    indices of p and k to those found.
    */
    double max = 0.0;
    for (int i = 0; i < n; ++i){
        for (int j = i+1; j < n; ++j){
            double aij = fabs(m_A(i,j));
            if ( aij > max){
                max = aij; m_k = i; m_p = j;
            }
        }
    }
    if(max <= m_tol){
        m_tol_reached = True;
    }
}

void Solver::calc_tau(){
    /* Calculate tau
    */
    m_tau = (m_A(m_p,m_p) - m_A(m_k,m_k))/(2*m_A(m_k,m_p));
}

void Solver::rotate(){
    /* Perform Jacobi method by making similarity transformations.

    */
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
    }else{
        c = 1.0;
        s = 0.0;
        }

    double a_kk, a_pp, a_ik, a_ip, r_ik, r_ip;
    a_kk = m_A(m_k,m_k);
    a_pp = m_A(m_p,m_p);
    m_A(m_k,m_k) = c*c*a_kk - 2.0*c*s*m_A(k,l) + s*s*a_pp;
    m_A(m_p,m_p) = s*s*a_kk + 2.0*c*s*m_A(k,l) + c*c*a_pp;
    m_A(m_k,m_p) = 0.0; // hard-coding non-diagonal elements by hand
    m_A(m_p,m_k) = 0.0; // same here
    for(int i = 0; i < n; i++ ){
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
        r_il = m_R(i,m_p);
        m_R(i,m_k) = c*r_ik - s*r_ip;
        m_R(i,m_p) = c*r_ip + s*r_ik;
    }
    return;
}
void Solver::run(){
    
    for(m_i = 0; m_i < m_max_iter; m_i++){
            Solver.max_off_diag();
            if(m_rol_reached == True){
                //print some stuff using a class function.
                 
                Solver::print_out();
                return;
            }

            Solver::calc_tau();
            Solver::rotate();
    }
    
    // If we have reached here then we are at max iter. 
    // Hence print some stuff using a class function.
    // Or write to file.
    Solver::print_out();
    return;
}

void Solver::print_out(){

    // print to file the number of iterations or something. 
    return;
}

void Solver::write_to_file(string filename){
    /* Write the information to file
    */
    m_ofile.open(filename);
    m_ofile << m_n << " " << m_i << "iterations" << end1;
    m_ofile.close();

}