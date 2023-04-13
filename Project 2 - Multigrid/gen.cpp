#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include "json.h"
#include <fstream>
using namespace std;
using filesystem::directory_iterator;

bool issub(const string &s, const string &t){
    for(int i = 0; i < s.size()-t.size(); i++){
        if(s.substr(i,t.size()) == t)
            return true;
    }
    return false;
}

int main() {
    string path = "examples/";
    for (const auto & file : directory_iterator(path)){
        string s = file.path();
        if(!issub(s, "test3-5")) continue;
        cout << s << endl;

        Json::Reader reader;
        Json::Value problem;
        ifstream ifs;
        ifs.open(s);
        if(!ifs.is_open()){
            cerr << "cannot read file " << s << endl;
            exit(-1);
        }
        if(!reader.parse(ifs, problem)){
            cerr << "parse error" << endl;
            exit(-1);
        }
        ifs.close();

        problem["Cycle Type"] = "FMG";
        problem["eps"] = 1e-8;
        problem["Max Iter"] = 20;

        Json::StyledWriter writer;
        ofstream os;
        os.open(s);
        os << writer.write(problem);
        os.close();
    }
    return 0;
}