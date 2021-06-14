#include<bits/stdc++.h>
#include "PrimeField.hpp"
#define ll long long
using namespace std;

size_t power2(size_t n){
    size_t res = 1;
	while(n > res)
		res <<= 1;
	res <<= 1;
	return res;
}
template<int p>
void fft (vector<PrimeField<p> > &a, bool invert) {
    PrimeField<p> root=root.get_generator();
	int n = (int) a.size();
	if (n == 1)  return;

	vector<PrimeField<p> > a0 (n/2),  a1 (n/2);
	for (int i=0, j=0; i<n; i+=2, ++j) {
		a0[j] = a[i];
		a1[j] = a[i+1];
	}
	fft <p>(a0, invert);
	fft <p>(a1, invert);

	PrimeField<p> w = 1;
	PrimeField<p> wlen = root;
	wlen = wlen^((p-1)/(int)a.size());
	if(invert)
		wlen = ~wlen;
	for (int i=0; i<n/2; ++i) {
		a[i] = a0[i] + w * a1[i];
		a[i+n/2] = a0[i] - w * a1[i];
		if(invert){
            a[i] /= 2;
            a[i+n/2] /= 2;
		}
        w = w * wlen;
	}
}
template<int p>
vector<PrimeField<p> >& multiply(vector<PrimeField<p> > &a, vector<PrimeField<p> > &b){
	vector<PrimeField<p> > fa(a.begin(), a.end()), fb(b.begin(), b.end());
    size_t n = power2(max(a.size(), b.size()));
	fa.resize(n);
	fb.resize(n);
	fft<p>(fa, 0);
	fft<p>(fb, 0);
	for(size_t i=0;i<n;i++)
		fa[i] = fa[i] * fb[i];
	fft<p>(fa, 1);
	return  *(new vector<PrimeField<p> > (fa.begin(), fa.end()));
}
vector<unsigned> pollard3(vector<unsigned> &a, vector<unsigned> &b){
    unsigned n = power2(a.size() + b.size());
    const unsigned
        p1 = 27 * (1 << 26) + 1,
        p2 =  7 * (1 << 26) + 1,
        p3 = 33 * (1 << 25) + 1;

    PrimeField<p1> w1 = w1.get_generator();
    PrimeField<p2> w2 = w2.get_generator();
    PrimeField<p3> w3 = w3.get_generator();

    vector<PrimeField<p1>> newa(a.begin(), a.end());
    vector<PrimeField<p1>> newb(b.begin(), b.end());

    vector<PrimeField<p2>> newa2(a.begin(), a.end());
    vector<PrimeField<p2>> newb2(b.begin(), b.end());

    vector<PrimeField<p3>> newa3(a.begin(), a.end());
    vector<PrimeField<p3>> newb3(b.begin(), b.end());

    vector<PrimeField<p1>> v1 = multiply<p1>(newa, newb);
    vector<PrimeField<p2>> v2 = multiply<p2>(newa2, newb2);
    vector<PrimeField<p3>> v3 = multiply<p3>(newa3, newb3);

    auto N2 = p1;
    auto N3 = (unsigned long long)p1 * p2;

    auto C2 = PrimeField<p2>(1) / PrimeField<p2>(N2);
    auto C3 = PrimeField<p3>(1) / PrimeField<p3>(N3);

    vector<unsigned> res(n);
    unsigned d = 0;
    for(int i = 0; i < n; i++){
        unsigned  w = v1[i].to_int();
        w += (unsigned long long)N2 * ((v2[i].to_int() - w)*C2);
        w += (unsigned          )N3 * ((v3[i].to_int() - w)*C3);
        res[i] = d + w;
        d = (d + w)>>32;
    }

    return res;
}
int main(){
    const int p = 998244353;//998244353
    PrimeField<p> root=root.get_generator();

    vector<PrimeField<p> > a = {1,1,1};
    vector<PrimeField<p> > b = {1,1,1};
	for(auto i:multiply<p>(a,b)){
		cout<<i<<" ";
	}
}

