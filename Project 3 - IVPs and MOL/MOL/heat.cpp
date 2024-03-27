#include <bits/stdc++.h>
#include "IVP.h"
#include <jsoncpp/json/json.h>
using namespace std;
typedef complex<double> Complex;

const double PI = acos(-1);

class HeatEquation : public TimeFunction{
private:
    double h;
    // 在这里设置边值条件
    double b0(const double &t) const{
        return 0.0;
    }
    double b1(const double &t) const{
        return 0.0;
    }
public:
    HeatEquation(const int &n){
        __isLinear = true;
        h = 1.0 / n;
    }
    ColVector operator () (const ColVector &U, const double &t) const{
        ColVector U1(U.size());
        for(int i = 0; i < U.size(); i++){
            U1(i) = -2*U(i)/(h*h);
            if(i-1 >= 0) U1(i) += U(i-1)/(h*h);
            if(i+1 < U.size()) U1(i) += U(i+1)/(h*h);
            if(i==0) U1(i) += b0(t) / (h*h);
            if(i+1==U.size()) U1(i) += b1(t) / (h*h);
        }
        return U1;
    }
    ColVector solve(const Matrix &a, const ColVector &c, const ColVector &U0, const double &t0, const double &k) const {
        int s = c.size(), m = U0.size();
        ColVector rhs(s*m);
        for(int i = 0; i < s; i++)
            rhs.setSubmatrix(i*m, 0, (*this)(U0, t0+c(i)*k));
        Matrix coef(s*m, s*m);
        for(int i = 0; i < s; i++)
            for(int j = 0; j < s; j++)
                for(int l = 0; l < m; l++){
                    coef(i*m+l, j*m+l) = a(i,j)*k*2.0/(h*h);
                    if(i==j) coef(i*m+l, j*m+l) += 1;
                    if(l) coef(i*m+l, j*m+l-1) = -a(i,j)*k*1.0/(h*h);
                    if(l<m-1) coef(i*m+l, j*m+l+1) = -a(i,j)*k*1.0/(h*h);
                }
        return coef.solve(rhs);
    }
};

// 在这里设置初值条件
double g0(const double &x){
    if(0.45 <= x && x < 0.5) return 20*(x-0.45);
    else if(0.5 <= x && x < 0.55) return -20*(x-0.55);
    else return 0;
}

// Fourier Solver 使用分离变量后的Fourier级数展开，是解析方法，但是只支持齐次边值条件
class FourierSolver{
private:
    int maxn;
    vector<double> A;

    double F(const double &x) const{
        return 2 * initial(x) * sin(maxn*PI*x);
    }

    double squareIntegrate(const double &l, const double &r) const{
        double mid = (l + r) / 2;
        return (F(l) + 4*F(mid) + F(r)) * (r - l) / 6;
    }

    // Simpson自适应积分
    double integrate(const double &l, const double &r, const double &V, const double &eps) const{
        double mid = (l + r) / 2;
        double L = squareIntegrate(l,mid);
        double R = squareIntegrate(mid,r);
        if(fabs(L+R-V) <= 15*eps) return L + R + (L+R-V) / 15.0;
        return integrate(l,mid,L,eps) + integrate(mid,r,R,eps);
    }

    double integrate(const double &l, const double &r, const double &eps) const{
        // 将积分区间随机划分为9个小区间，防止出现因为5个等分点为0导致整个积分返回0的情况
        int points[10], m;
        for(int i = 1; i < 9; i++) points[i] = rand();
        points[0] = 0; points[9] = RAND_MAX;
        sort(points, points+10);
        m = unique(points, points+10) - points;
        double sum = 0;
        for(int i = 0; i < m-1; i++){
            double hl = l + (double)points[i]/RAND_MAX * (r-l);
            double hr = l + (double)points[i+1]/RAND_MAX * (r-l);
            sum += integrate(hl, hr, squareIntegrate(hl,hr), max(1e-4*eps, 1e-15));
        }
        return sum;
    }

public:
    FourierSolver(const int &N){
        cout << "--------------------------------------------------------------------------------" << endl;
        cout << "Solving with Fourier solver..." << endl;
        maxn = 0;
        do{
            ++maxn;
            A.push_back(integrate(0, 1, 1e-15));
        } while(maxn < N);
        cout << maxn << " terms has been retained in the Fourier series." << endl;
    }

    double initial(const double &x) const{
        return g0(x);
    }

    double operator () (const double &x, const double &t) const{
        double sum = 0;
        for(int k = 1; k <= maxn; k++){
            sum += A[k-1] * exp(-k*k*PI*PI*t) * sin(k*PI*x);
        }
        return sum;
    }

    void output(const string &fname, const double &T, const double &h){
        cout << "--------------------------------------------------------------------------------" << endl;
        ofstream fout(fname);
        fout << setprecision(12) << T << " ";
        for(double x = h; x <= 1-h+1e-14; x += h){
            fout << (*this)(x, T) << " ";
        }
        fout << endl;
        fout.close();
        cout << "Output: Results has been saved to " << fname << endl;
    }

    void denseDiscreteOutput(const string &fname, const double &step, const double &T, const double &h){
        cout << "--------------------------------------------------------------------------------" << endl;
        ofstream fout(fname);
        fout << setprecision(12);
        for(double x = h/2; x <= 1; x += h){
            fout << initial(x) << " ";
        }
        fout << endl;
        for(double t = step; t <= T+1e-14; t += step){
            for(double x = h/2; x <= 1; x += h){
                fout << (*this)(x, t) << " ";
            }
            fout << endl;
        }
        fout.close();
        cout << "Dense-Discrete Output: Results has been saved to " << fname << endl;
    }
};

// DFT Solver 使用FFT将问题转化为频域上若干个独立的线性ODE，只支持周期边界条件，且空间划分数必须为2的幂次
class DFTSolver{
private:
    int m;
    vector<Complex> initFcoef;

    vector<Complex> fft(const vector<Complex> &in_a, const int &v) const{
        const int n = in_a.size();
        static const double pi = acos(-1);
        int l = 0, *r = new int[n];
        memset(r, 0, sizeof(int) * n);
        for(int i = 1; i < n; i <<= 1) l++;
        for(int i = 0; i < n; i++)
            r[i] = ( r[i/2] / 2 ) | ( (i&1) << (l-1) );
        vector<Complex> a = in_a;
        for(int i = 0; i < n; i++)
            if(i < r[i]) swap(a[i], a[r[i]]);
        for(int i = 1; i < n; i <<= 1){
            Complex wn(cos(pi/i), v*sin(pi/i));
            int p = i << 1;
            for(int j = 0; j < n; j += p)
            {
                Complex w(1, 0);
                for(int k = 0; k < i; k++)
                {
                    Complex x = a[j+k], y = w*a[i+j+k];
                    a[j+k] = x + y;
                    a[i+j+k] = x - y;
                    w = w * wn;
                }
            }
        }
        if(v==-1) for(auto& z : a) z/=n;
        return a;
    }

    double initial(const double &x) const{
        return g0(x);
    }

public:
    DFTSolver(const int &N){
        cout << "--------------------------------------------------------------------------------" << endl;
        cout << "Solving diffusion equation with DFT solver (with periodic boundary)..." << endl;
        m = N;
        double h = 1.0/m;
        for(int i = 0; i < m; i++){
            initFcoef.push_back(initial(h/2 + i*h));
        }
        initFcoef = fft(initFcoef, 1);
    }

    RowVector operator () (const double &t) const{
        double sum = 0;
        auto tmpCoef = initFcoef;
        for(int k = 0; k < m; k++){
            //cout << tmpCoef[k] << " ";
            int n = (k+m/2)%m - m/2;
            tmpCoef[k] *= exp(-4*PI*PI*n*n*t);
            //cout << tmpCoef[k] << endl;
        }
        tmpCoef = fft(tmpCoef, -1);
        RowVector res(m);
        for(int i = 0; i < m; i++){
            res(i) = tmpCoef[i].real();
        }
        return res;
    }

    void output(const string &fname, const double &T, const double &h){
        cout << "--------------------------------------------------------------------------------" << endl;
        ofstream fout(fname);
        auto res = (*this)(T);
        for(int i = 0; i < res.size(); i++)
            fout << res(i) << " ";
        fout.close();
        cout << "Output: Results has been saved to " << fname << endl;
    }

    void denseDiscreteOutput(const string &fname, const double &step, const double &T, const double &h){
        cout << "--------------------------------------------------------------------------------" << endl;
        ofstream fout(fname);
        fout << setprecision(12);
        for(double x = h/2; x <= 1; x += h){
            fout << initial(x) << " ";
        }
        fout << endl;
        for(double t = step; t <= T+1e-14; t += step){
            auto res = (*this)(t);
            for(int i = 0; i < res.size(); i++)
                fout << res(i) << " ";
            fout << endl;
        }
        fout.close();
        cout << "Dense-Discrete Output: Results has been saved to " << fname << endl;
    }
};

int main(int argc, char* argv[]){
    Json::Reader reader;
    Json::Value problem;
    ifstream ifs;
    if(argc < 2){
        cerr << "Please provide input file name!" << endl;
        exit(-1);
    }
    ifs.open(argv[1]);
    if(!ifs.is_open()){
        cerr << "cannot read file " << argv[1] << endl;
        exit(-1);
    }
    if(!reader.parse(ifs, problem)){
        cerr << "parse error" << endl;
        exit(-1);
    }

    int n = problem["Space Section"].asInt();
    double h = 1.0 / n;
    double T = problem["End Time"].asDouble();

    if(problem["Method"].asString() == "Fourier"){
        FourierSolver solver(problem["Fourier Terms"].asInt());
        if(problem["Output"].asBool())
            solver.output("result.txt", T, h);
        if(problem["Dense-Discrete Output"].asBool())
            solver.denseDiscreteOutput("result-dense.txt", problem["Dense-Discrete Output Step"].asDouble(), T, h);
        return 0;
    }

    if(problem["Method"].asString() == "DFT"){
        DFTSolver solver(problem["Space Section"].asInt());
        if(problem["Output"].asBool())
            solver.output("result.txt", T, h);
        if(problem["Dense-Discrete Output"].asBool())
            solver.denseDiscreteOutput("result-dense.txt", problem["Dense-Discrete Output Step"].asDouble(), T, h);
        return 0;
    }

    ColVector U0(n-1);
    for(int i = 1; i < n; i++)
        U0(i-1) = g0(i*h);
    HeatEquation f(n);

    auto& factory = TimeIntegratorFactory::Instance();
    auto solver = factory.createTimeIntegrator(problem["Method"].asString(), problem["Order"].asInt());
    if(problem.isMember("Time Section")){
        solver->solveWithInfo(f, U0, T, problem["Time Section"].asInt());
    } else {
        if(problem.isMember("Max Step"))
            solver->setMaxStep(problem["Max Step"].asDouble());
        solver->solveWithInfo(f, U0, T, problem["Tolerance"].asDouble());
    }
    
    if(problem["Output"].asBool())
        solver->output("result.txt");
    if(problem["Dense-Discrete Output"].asBool())
        solver->denseDiscreteOutput("result-dense.txt", problem["Dense-Discrete Output Step"].asDouble());
    return 0;
}