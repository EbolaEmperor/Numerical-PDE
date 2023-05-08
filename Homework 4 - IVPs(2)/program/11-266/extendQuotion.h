#include "fraction.h"

// a + b * sq(6)

const fraction six(6,1);

class extendQuotion6{
private:
    fraction a, b;
public:
    extendQuotion6(): a(0), b(0){}
    extendQuotion6(const int &n): a(n), b(0){}
    extendQuotion6(const fraction &a): a(a), b(0){}
    extendQuotion6(const fraction &a, const fraction &b): a(a), b(b){}
    extendQuotion6(const extendQuotion6 & rhs) = default;
    extendQuotion6 operator + (const extendQuotion6 &rhs) const{
        extendQuotion6 res;
        res.a = a + rhs.a;
        res.b = b + rhs.b;
        return res;
    }
    extendQuotion6 operator - (const extendQuotion6 &rhs) const{
        extendQuotion6 res;
        res.a = a - rhs.a;
        res.b = b - rhs.b;
        return res;
    }
    extendQuotion6 operator - () const{
        extendQuotion6 res;
        res.a = -a;
        res.b = -b;
        return res;
    }
    extendQuotion6 operator * (const extendQuotion6 &rhs) const{
        extendQuotion6 res;
        res.a = a*rhs.a + six*b*rhs.b;
        res.b = a*rhs.b + b*rhs.a;
        return res;
    }
    extendQuotion6 operator ^ (const int &rhs) const{
        extendQuotion6 res(1), a(*this);
        for(int b = rhs; b; b >>= 1, a *= a)
            if(b & 1) res *= a;
        return res;
    }
    extendQuotion6 operator / (const fraction &rhs) const{
        extendQuotion6 res;
        res.a = a / rhs;
        res.b = b / rhs;
        return res;
    }
    extendQuotion6 operator / (const extendQuotion6 &rhs) const{
        extendQuotion6 tmp;
        tmp.a = rhs.a; tmp.b = -rhs.b;
        return (*this) * tmp / (rhs.a*rhs.a - six*rhs.b*rhs.b);
    }
    extendQuotion6& operator += (const extendQuotion6 &rhs){
        return (*this) = (*this) + rhs;
    }
    extendQuotion6& operator -= (const extendQuotion6 &rhs){
        return (*this) = (*this) - rhs;
    }
    extendQuotion6& operator *= (const extendQuotion6 &rhs){
        return (*this) = (*this) * rhs;
    }
    extendQuotion6& operator /= (const extendQuotion6 &rhs){
        return (*this) = (*this) / rhs;
    }
    bool operator == (const extendQuotion6 &rhs) const{
        return a==rhs.a && b==rhs.b;
    }
    bool operator != (const extendQuotion6 &rhs) const{
        return a!=rhs.a || b!=rhs.b;
    }
    bool operator >= (const extendQuotion6 &rhs) const{
        extendQuotion6 res = (*this) - rhs;
        if(res.a>=0 && res.b>=0) return true;
        else if(res.a<0 && res.b<0) return false;
        else if(res.a>=0 && res.b<0) return res.a*res.a >= res.b*res.b*6;
        else return res.a*res.a <= res.b*res.b*6;
    }
    friend std::ostream & operator << (std::ostream &out, const extendQuotion6 &x) {
		if(x.b==0) out<<x.a;
		else{
            out << "(" << x.a;
            if(x.b > 0) out << '+';
            out << x.b << "*sq[6])";
        }
		return out;
	}
};