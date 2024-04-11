#pragma once

#include <iostream>
#include <string>
#include <map>

class BType{
private:
    std::map<std::string, std::string> types;
public:
    void setBondary(const std::string& bon, const std::string& typ){
        if(typ != "Dirichlet" && typ != "Neumann"){
            std::cerr << "Unrecignized bondary '" << typ << "'" << std::endl;
            exit(-1);
        }
        types[bon] = typ;
    }

    std::string operator () (const std::string& bon) const{
        if(!types.count(bon)){
            std::cerr << "Undefined bondary '" << bon << "'" << std::endl;
            exit(-1);
        }
        return types.find(bon)->second;
    }
};