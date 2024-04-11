#include "mathExpr.h"
#include <stack>
#include <cmath>
#include <iostream>

const char LOG = 1;
const char LOG2 = 2;
const char LOG10 = 3;
const char EXP = 4;
const char SIN = 5;
const char COS = 6;
const char TAN = 7;
const char ASIN = 8;
const char ACOS = 9;
const char ATAN = 10;
const char SINH = 11;
const char COSH = 12;
const char TANH = 13;
const char SQRT = 14;
const char FLOOR = 15;
const char CEIL = 16;
const char ROUND = 17;
const char ABS = 18;

const double E = exp(1);
const double PI = acos(-1);

int priority(const char &op){
    switch(op){
        case '^': return 2;
        case '*': case '/': return 3;
        case '+': case '-': return 4;
        default: return 1;
    }
}

char recognize(std::string &s){
    if(s=="log") return LOG;
    if(s=="log2") return LOG2;
    if(s=="log10") return LOG10;
    if(s=="exp") return EXP;
    if(s=="sin") return SIN;
    if(s=="cos") return COS;
    if(s=="tan") return TAN;
    if(s=="asin") return ASIN;
    if(s=="acos") return ACOS;
    if(s=="atan") return ATAN;
    if(s=="sinh") return SINH;
    if(s=="cosh") return COSH;
    if(s=="tanh") return TANH;
    if(s=="sqrt") return SQRT;
    if(s=="floor") return FLOOR;
    if(s=="ceil") return CEIL;
    if(s=="round") return ROUND;
    if(s=="abs") return ABS;
    std::cerr << "[Error] Unrecognized expression '" << s << "'" << std::endl;
    exit(-1);
    return 0;
}

double calc(const char &op, const double &a, const double &b = 0){
    switch(op){
        case '*': return a * b;
        case '/': return a / b;
        case '+': return a + b;
        case '-': return a - b;
        case '^': return pow(a, b);
        case LOG: return log(a);
        case LOG2: return log2(a);
        case LOG10: return log10(a);
        case EXP: return exp(a);
        case SIN: return sin(a);
        case COS: return cos(a);
        case TAN: return tan(a);
        case ASIN: return asin(a);
        case ACOS: return acos(a);
        case ATAN: return atan(a);
        case SINH: return sinh(a);
        case COSH: return cosh(a);
        case TANH: return tanh(a);
        case SQRT: return sqrt(a);
        case FLOOR: return floor(a);
        case CEIL: return ceil(a);
        case ROUND: return round(a);
        case ABS: return fabs(a);
    }
    return 0;
}

void deal(std::stack<char> &opt, std::stack<double> &num, const int &pri = 5){
    while(!opt.empty() && priority(opt.top()) <= pri){
        char p = opt.top(); opt.pop();
        if(priority(p) == 1){
            if(num.empty()){
                std::cerr << "[Error] Syntax Error!" << std::endl;
                exit(-1);
            }
            double a = num.top(); num.pop();
            num.push(calc(p, a));
        } else {
            if(num.size() < 2){
                std::cerr << "[Error] Syntax Error!" << std::endl;
                exit(-1);
            }
            double b = num.top(); num.pop();
            double a = num.top(); num.pop();
            num.push(calc(p, a, b));
        }
    }
}

void Expression::addConstant(const std::string &conName, const double &val){
    if(constant.count(conName) || conName=="e" || conName=="pi"){
        std::cerr << "[Error] Constant '" << conName << "' has already existed" << std::endl;
        exit(-1);
    }
    constant[conName] = val;
}

void Expression::addVariable(const std::string &varName){
    if(variable.count(varName)){
        std::cerr << "[Error] Variable '" << varName << "' has already existed" << std::endl;
        exit(-1);
    }
    variable[varName] = varcnt++;
}

void Expression::setExpr(const std::string &s){
    expr = s;
    expr += '\0';
}

double Expression::eval(const std::vector<double> &parm) const{
    if(parm.size() != varcnt){
        std::cerr << "[Error] The number of parameters does't match." << std::endl;
        exit(-1);
    }
    int curPos = -1;
    return evaluate(curPos, parm);
}

double Expression::eval() const{
    std::vector<double> parm;
    return eval(parm);
}

double Expression::operator() (const std::vector<double> &parm) const{
    return eval(parm);
}

double Expression::evaluate(int &curPos, const std::vector<double> &parm) const{
    std::stack<char> op;
    std::stack<double> num;
    curPos++;
    while(expr[curPos]){
        if(expr[curPos] == '-' && (!curPos || expr[curPos-1] == '(')){
            num.push(0);
            op.push('-');
            curPos++;
        } else if(isdigit(expr[curPos])){
            std::string tmp;
            while(isdigit(expr[curPos])) tmp += expr[curPos++];
            if(expr[curPos] == '.'){
                tmp += expr[curPos++];
                while(isdigit(expr[curPos])) tmp += expr[curPos++];
            }
            num.push(stod(tmp));
        } else {
            if(expr[curPos] == '('){
                num.push(evaluate(curPos, parm));
                curPos++;
            }
            else if(expr[curPos] == ')') break;
            else{
                char curop = expr[curPos];
                if(isalpha(curop)){
                    std::string tmp;
                    while(isalpha(expr[curPos]) || isdigit(expr[curPos]))
                        tmp += expr[curPos++];
                    if(variable.count(tmp)){
                        num.push(parm[variable.find(tmp)->second]);
                        continue;
                    } else if(constant.count(tmp)){
                        num.push(constant.find(tmp)->second);
                        continue;
                    } else if(tmp=="e"){
                        num.push(E);
                        continue;
                    } else if(tmp == "pi"){
                        num.push(PI);
                        continue;
                    }
                    curop = recognize(tmp);
                } else curPos++;
                deal(op, num, priority(curop));
                op.push(curop);
            }
        }
    }
    deal(op, num);
    if(num.size() != 1){
        std::cerr << "[Error] Syntax Error!" << std::endl;
        exit(-1);
    }
    return num.top();
}