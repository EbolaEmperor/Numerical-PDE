#include "BVP.h"

int main(int argc, char* argv[]){
    BVP bvp;
    bvp.read(argv[1]);
    bvp.printProblem();
    bvp.solve();
    bvp.output("result.txt");
    bvp.checkError();
    return 0;
}