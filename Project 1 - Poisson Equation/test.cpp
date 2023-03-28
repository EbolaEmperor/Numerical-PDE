#include "BVP.h"
#include <fstream>
using namespace std;

int main(int argc, char* argv[]){
    BVP bvp[4];
    Json::Reader reader;
    Json::Value problem;
    ifstream ifs;
    ifs.open(argv[1]);
    if(!ifs.is_open()){
        cerr << "cannot read file " << argv[1] << endl;
        exit(-1);
    }
    if(!reader.parse(ifs, problem)){
        cerr << "parse error" << endl;
        exit(-1);
    }

    vector<double> err[4];
    int grid = 16;
    for(int i = 0; i < 4; i++){
        problem["Grid Size"] = grid;
        bvp[i].read(problem);
        if(i==0) bvp[i].printProblem();
        bvp[i].solve();
        cout << "Grid Size: " << grid << endl;
        //bvp[i].output("result.txt");
        err[i] = bvp[i].checkError();
        grid *= 2;
    }
    cout << "----------------------------------------------------------" << endl;
    cout << "Converge Rate (1-norm): " << log2(err[2][0]/err[3][0]) << endl;
    cout << "Converge Rate (2-norm): " << log2(err[2][1]/err[3][1]) << endl;
    cout << "Converge Rate (max-norm): " << log2(err[2][2]/err[3][2]) << endl;
    return 0;
}