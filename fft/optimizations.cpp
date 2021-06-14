#include<bits/stdc++.h>
using namespace std;
const double PI = acos(-1);
typedef complex<double> base;

const int maxn = 1<<15;
//const int n = 1 << 13;
const int cs = 1 << 16;
const int numThreads = 4;

int rev[maxn];
base a[maxn], b[maxn];
thread threads[numThreads];
base fa[maxn], fb[maxn];
base w_pow[maxn];

void fft(base a[], int n, bool invert){
//    int n = (int)a.size();

    for (int i=1, j=0; i<n; ++i) {
		int bit = n >> 1;
		for (; j>=bit; bit>>=1)
			j -= bit;
		j += bit;
		if (i < j)
			swap (a[i], a[j]);
	}


    for(int len=2;len<=n;len<<=1){
        double ang = 2*PI/len*(invert?-1:1);
        base wn(cos(ang), sin(ang));

        w_pow[0] = base(1, 0);
        for(int i=1;i<len/2;i++)
            w_pow[i] = w_pow[i-1]*wn;

        for(int i=0;i<n;i+=len){
            for(int j=0;j<len/2;j++){
                base a0 = a[i+j], a1 = w_pow[j]*a[i+j+len/2];
                a[i+j]       = a0 + a1;
                a[i+j+len/2] = a0 - a1;
            }
        }
    }
    if(invert)
        for(int i=0;i<n;i++)
            a[i] /= n;
}
void mul(base a[], base b[],const int lena,const int lenb, vector<int> &res){
    int n=1;
    while(n<max(lena, lenb))
        n<<=1;
    n<<=1;
    for(int i=0;i<n;i++){
        fa[i] = a[i];
        fb[i] = b[i];
    }

    res.resize(n);
    fft(fa, n, 0);
    fft(fb, n, 0);
    for(int i=0;i<n;i++)
        fa[i] *= fb[i];
    fft(fa, n, 1);
    for(int i=0;i<n;i++)
        res[i] = int(fa[i].real() + 0.5);
}
void calc_rev();
int main(){

    calc_rev();

    vector<int> res;
    int n = 1<<15;
    for(int i=0;i<n;i++)a[i] = rand(), b[i] = rand();

    auto start = std::chrono::steady_clock::now();

    mul(a, b, n, n,res);

    auto finish = std::chrono::steady_clock::now();
    std::chrono::duration<double> time = finish - start;
    cout<<std::chrono::duration_cast<std::chrono::milliseconds>(time).count();
    cout<<endl;
//    for(int i:res){
//        cout<<i<<" ";
//    }
    return 0;
}
void calc_rev() {
    for(int i=1,j=0;i<maxn;i++){
        int bit = maxn>>1;
        for(;j>=bit;bit>>=1)
            j -= bit;
        j += bit;
        rev[i] = j;
    }
}
