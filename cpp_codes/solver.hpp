#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <iostream>
#include <cmath>
#include <fstream>
#include <armadillo>

using namespace arma;
using namespace std;

class Solver{

private:
    int m_n;
    mat m_A;
    mat m_R;
    double m_tol;
    int m_max_iter;
    int m_k;
    int m_p;
    double m_tau;
    bool m_tol_reached;
    int m_i;
    ofstream m_ofile;

public:
    void init(int n, mat A, double tol);
    void max_off_diag();
    void calc_tau();
    void rotate();
    void write_to_file(string filename);
    void print_out();
protected:


};
#endif
