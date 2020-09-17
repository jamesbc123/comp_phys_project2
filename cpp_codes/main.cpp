#include "solver.hpp"
#include <armadillo>
#include <time.h>
#include <fstream>

using namespace std;
using namespace arma;

int main(){
    clock_t start, end;
    double timeused;

    string filename = "./iterations.txt.";
    ofsteam ofile;

    int n = 100;
    double tol = 1e-8;
    mat A = zeros<mat>(n, n);
    A.diag().fill(2);
    A.diag(-1).fill(-1);
    A.diag(1).fill(-1);
    
    Solver my_solver;
    my_solver.init(n, A, tol);
    start = clock();
    my_solver.run();
    end = clock();
    timeused = (end-start)/CLOCKS_PER_SEC; 
    my_solver.write_to_file(filename);

    // Open here again and add timing
    // This could be done inside the class too. What do 
    // you think amund?
    ofile.open(filename);
    m_ofile << timused << "time_used" << end1;
    m_ofile.close();
    return 0;
}