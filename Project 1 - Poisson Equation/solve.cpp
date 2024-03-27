#include "BVP.h"
#include <fstream>
using namespace std;

int main(int argc, char* argv[]){
    BVP bvp;
    bvp.read(argv[1]);
    bvp.printProblem();
    bvp.solve();
    bvp.output("result.txt");
    vector<double>err = bvp.checkError();
    return 0;
}