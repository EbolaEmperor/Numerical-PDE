#include <bits/stdc++.h>
using namespace std;

const double pi = acos(-1);

class FourierSolver{
private:
    int maxn;
    vector<double> A;

    double F(const double &x) const{
        return 2 * initial(x) * sin(maxn*pi*x);
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
    FourierSolver(const double &eps){
        maxn = 0;
        do{
            ++maxn;
            A.push_back(integrate(0, 1, eps));
        } while(maxn < 2 || fabs(A[maxn-1]) > eps || fabs(A[maxn-2]) > eps);
        cout << maxn << " items has been retained in the Fourier series." << endl;
    }

    double initial(const double &x) const{
        if(0.45 <= x && x < 0.5) return 20*(x-0.45);
        else if(0.5 <= x && x < 0.55) return -20*(x-0.55);
        else return 0;
    }

    double operator () (const double &x, const double &t) const{
        double sum = 0;
        for(int k = 1; k <= maxn; k++){
            sum += A[k-1] * exp(-k*k*pi*pi*t) * sin(k*pi*x);
        }
        return sum;
    }
};

int main(int argc, char* argv[]){
    double t = stod(argv[1]), eps = stod(argv[2]);
    ofstream fout("result.txt");
    FourierSolver truesol(eps);
    for(int i = 0; i <= 1000; i++){
        double val = t ? truesol(0.001*i, t) : truesol.initial(0.001*i);
        fout << 0.001*i << " " << val << "\n";
    }
    fout.close();
    cout << "The result at t=" << t << " has been saved to result.txt." << endl;
    return 0;
}