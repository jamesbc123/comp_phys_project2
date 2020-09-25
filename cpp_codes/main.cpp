#include "solver.hpp"
#include "utils.hpp"
#include <armadillo>
#include <time.h>
#include <fstream>
#include <stdio.h>  // remove(filename)

using namespace std;
using namespace arma;

int main(){
    run_buckling_beam(); // Run the code for problem c)

    //run_q_dots_one_electron(); // Run the code for problem d)

    return 0;
}