#include "Solver/BVP.h"
using namespace std;

int main(int argc, char* argv[]){
    BVP<2> bvp;
    bvp.read(argv[1]);
    bvp.printProblem();
    bvp.solve();
    bvp.output("result.txt");
    vector<double>err = bvp.checkError();
    return 0;
}