#ifndef PRIME_FIELD_INCLUDED
#define PRIME_FIELD_INCLUDED
#include<bits/stdc++.h>
using namespace std;
template< int p >
class PrimeField {
    public:
    int value;

    static int from_int(int v) {
        return ( v < 0 ? p - (-v % p) : (v % p) );
    }

    PrimeField<p> binpow(int a, int n){
        PrimeField<p> res(a);
        return binpow(res, n);
    }
    //binary exponentation
    PrimeField<p> binpow(PrimeField<p> a, int n){
        if(n==0)return 1;
        if(n==1)return a;
        if(n & 1)return a * binpow(a, n-1);
        auto res = binpow(a, n/2);
        return res * res;
    }

    PrimeField(int v = 0) : value( from_int(v) ) {}
    PrimeField<p> operator = (int n){
        value = n;
    }
    //reverse value in module P
    PrimeField<p>& operator~ (){
        return *(new PrimeField<p>(binpow(value, p - 2)));
    }

    //https://cp-algorithms.com/algebra/primitive-root.html
    PrimeField<p>& get_generator(){
        vector<int> fact;
        int phi = p-1,  n = phi;
        for (int i=2; i*i<=n; ++i)
            if (n % i == 0) {
                fact.push_back (i);
                while (n % i == 0)
                    n /= i;
            }
        if (n > 1)
            fact.push_back (n);

        for (int res=2; res<=p; ++res) {
            bool ok = true;
            for (size_t i=0; i<fact.size() && ok; ++i)
                ok &= binpow (res, phi / fact[i]) != 1;
            if (ok)  return *(new PrimeField<p>(res));
        }
        return *(new PrimeField<p>(-1));
    }
    PrimeField<p> & operator += (PrimeField<p> v) {

        value += v.value;
        if (value >= p) {
            value -= p;
        }

        return *this;
    }

    PrimeField<p> & operator -= (PrimeField<p> v) {

        value -= v.value;
        if (value < 0) {
            value += p;
        }

        return *this;
    }

    PrimeField<p> & operator *= (PrimeField<p> v) {

        value = static_cast< long long >( value ) * v.value % p;

        return *this;
    }

    PrimeField<p> & operator /= (PrimeField<p> a) {
        return *this *= binpow(a, p - 2);
    }
    PrimeField<p>& operator /= (int a) {
        return *this *= binpow(a, p - 2);
    }

    PrimeField<p> & operator - () {

        value = p - value;

        return *this;
    }

    int to_int() const {
        return value;
    }
    operator int() const {
        return value;
    }

    operator bool() const {
        return value;
    }

    PrimeField<p> & operator ++ () {
        value = value + 1;
        return *this;
    }

    PrimeField<p> & operator -- () {
        value = value - 1;
        return *this;
    }

    PrimeField<p> & operator ++ (int) {
        PrimeField<p> o(*this);
        *this += 1;
        return o;
    }

    PrimeField<p> & operator -- (int) {
        PrimeField<p> o(*this);
        *this -= 1;
        return o;
    }

};
//exponentiation
template< int p >
PrimeField< p > operator ^ (PrimeField<p> a, int n) {
    return a.binpow(a, n);
}
template< int p >
PrimeField< p > operator + (PrimeField<p> a, PrimeField<p> b) {
    a += b;
    return a;
}

template< int p >
PrimeField< p > operator - (PrimeField<p> a, PrimeField<p> b) {
    a -= b;
    return a;
}

template< int p >
PrimeField< p > operator * (PrimeField<p> a, PrimeField<p> b) {
    a *= b;
    return a;
}

template< int p >
PrimeField< p > operator / (PrimeField<p> a, PrimeField<p> b) {
    a /= b;
    return a;
}
template< int p >
PrimeField< p > operator / (PrimeField<p> a, int b) {
    a /= b;
    return a;
}
template< int p >
PrimeField< p > operator / (int a, PrimeField<p> b) {
    PrimeField<p> fa(a);
    fa /= b;
    return fa;
}
template< int p >
bool operator == (PrimeField<p> a, PrimeField<p> b) {
    return a.to_int() == b.to_int();
}

template< int p >
bool operator != (PrimeField<p> a, PrimeField<p> b) {
    return a.to_int() != b.to_int();
}

template< int p >
std::ostream& operator<<( std::ostream& out, PrimeField< p > a ) {
    out << a.to_int();

    return out;
}



#endif
