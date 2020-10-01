#include "solver.hpp"
#include "utils.hpp"
#include <armadillo>
#include <time.h>
#include <fstream>
#include <stdio.h>  // remove(filename)

using namespace std;
using namespace arma;

int main(){
    
    run_buckling_beam(); // Run the code for problem 2c)

    //run_q_dots_one_electron(); // Run the code for problem 2d)
    
    //run_q_dots_two_electrons(); // Run the code for problem 2e)

    return 0;
}