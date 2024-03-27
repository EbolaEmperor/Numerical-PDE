#ifndef _PYLYNOMIAL_H_
#define _POLYNOMIAL_H_

#include <vector>
using namespace std;

template<class Type> class Polynomial{
public:
    int n;
    vector<Type> a;
    static Polynomial one;
    void removeLeadingZero(){
        while(n && a[n]==0) n--;
        if(n>=0) a.resize(n+1);
        else{
            n = 0;
            a.resize(1);
            a[0] = Type(0);
        }
    }
    Polynomial(): n(0){a.resize(1);}
    Polynomial(const int &_n){
        n = _n;
        a.resize(n+1);
    }
    Polynomial(const int &_n, const vector<Type> &p){
        n = _n;
        a = p;
    }
    Polynomial(const Polynomial &rhs) = default;
    Polynomial& operator = (const Polynomial &rhs) = default;
    Polynomial operator * (const Polynomial &rhs) const{
        Polynomial res(n+rhs.n);
        for(int i = 0; i <= n; i++)
            for(int j = 0; j <= rhs.n; j++)
                res.a[i+j] += a[i] * rhs.a[j];
        res.removeLeadingZero();
        return res;
    }
    Polynomial operator + (const Polynomial &rhs) const{
        int m = max(n, rhs.n);
        Polynomial res(m);
        for(int i = 0; i <= n; i++)
            res.a[i] += a[i];
        for(int i = 0; i <= rhs.n; i++)
            res.a[i] += rhs.a[i];
        res.removeLeadingZero();
        return res;
    }
    Polynomial operator - (const Polynomial &rhs) const{
        int m = max(n, rhs.n);
        Polynomial res(m);
        for(int i = 0; i <= n; i++)
            res.a[i] -= a[i];
        for(int i = 0; i <= rhs.n; i++)
            res.a[i] -= rhs.a[i];
        res.removeLeadingZero();
        return res;
    }
    Polynomial operator - () const{
        Polynomial res(*this);
        for(int i = 0; i <= n; i++)
            res.a[i] = -res.a[i];
        return res;
    }
    Polynomial integral () const{
        Polynomial res(n+1);
        for(int i = 0; i <= n; i++)
            res.a[i+1] = a[i] / Type(i+1);
        return res;
    }
    friend std::ostream& operator << (std::ostream &out, const Polynomial &p){
        int sti = p.n;
        while(sti >= 0 && p.a[sti] == 0) sti--;
        for(int i = sti; i >= 0; i--){
            if(p.a[i]==0) continue;
            if(i<sti && p.a[i]>=0) out << "+";
            out << p.a[i];
            if(i>=2) out << "*z^" << i;
            if(i==1) out << "z";
        }
        if(sti < 0) out << '0';
        return out;
    }
    bool operator == (const Polynomial & rhs) const{
        if(n!=rhs.n) return false;
        for(int i = 0; i <= n; i++)
            if(a[i] != rhs.a[i])
                return false;
        return true;
    }
};

template<class Type>Polynomial<Type> constPolynomial(const Type &a){
    Polynomial<Type> res(0);
    res.a[0] = a;
    return res;
}

template<class Type>Polynomial<Type> constPolynomial(const int &a){
    Polynomial<Type> res(0);
    res.a[0] = Type(a);
    return res;
}

template<class Type>Polynomial<Type> identityPolynomial(){
    Polynomial<Type> res(1);
    res.a[1] = Type(1);
    return res;
}

template<class Type> Polynomial<Type> Polynomial<Type>::one = constPolynomial<Type>(1);

template<class Type> class fracPolynomial{
private:
    Polynomial<Type> a, b;
public:
	fracPolynomial(){
        a = constPolynomial<Type>(0);
        b = constPolynomial<Type>(1);
    }
	fracPolynomial(const int &_a){
        a = constPolynomial<Type>(_a);
        b = constPolynomial<Type>(1);
    }
    fracPolynomial(const Polynomial<Type> &_a){
        a = _a;
        b = constPolynomial<Type>(1);
    }
	fracPolynomial(const Polynomial<Type> &_a, const Polynomial<Type> &_b){
        a = _a;
        b = _b;
    }
	fracPolynomial(const fracPolynomial &rhs) = default;
	friend std::ostream & operator << (std::ostream &out, const fracPolynomial &x) {
		if(x.b == 1) out << x.a;
		else out << '(' << x.a << ")/(" << x.b << ")";
		return out;
	}
	bool operator == (const fracPolynomial &y) const{
		return a == y.a && b == y.b;
	}
	bool operator != (const fracPolynomial &y) const{
		return !(*this==y);
	}
	fracPolynomial operator * (const fracPolynomial &y) const{
		return fracPolynomial(a * y.a, b * y.b);
	}
	fracPolynomial operator / (const fracPolynomial &y) const{
		return fracPolynomial(a * y.b , b * y.a);
	}
	fracPolynomial operator + (const fracPolynomial &y) const{
		return fracPolynomial(a * y.b + y.a * b, b * y.b);
	}
	fracPolynomial operator - (const fracPolynomial &y) const{
		return fracPolynomial(a * y.b - y.a * b, b * y.b);
	}
	fracPolynomial operator - () const{
		return fracPolynomial(-a, b);
	}
	fracPolynomial operator *= (const fracPolynomial &y){
		return (*this)=(*this)*y;
	}
	fracPolynomial operator /= (const fracPolynomial &y){
		return (*this)=(*this)/y;
	}
	fracPolynomial operator += (const fracPolynomial &y){
		return (*this)=(*this)+y;
	}
	fracPolynomial operator -= (const fracPolynomial &y){
		return (*this)=(*this)-y;
	}
};

#endif