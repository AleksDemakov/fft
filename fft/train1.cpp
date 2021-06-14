#include<bits/stdc++.h>
using namespace std;
typedef complex<double> base;
vector<int> rev;
base w    [14][(1<<14)];
base w_inv[14][(1<<14)];

int logw(int n){
    return sizeof(n)*8 - __builtin_clz(n)-1;
}

void bit(int n){
    rev.resize(n);
    for(int i=1, j=0;i<n;i++){
        int bit = n>>1;
        for(;j>=bit;bit>>=1)
            j -= bit;
        j += bit;
        rev[i] = j;
    }
}

void precalc_w(int n){
    int log_n = log2(n);
    const double PI = acos(-1);
    double ang = 2*PI/n;
    double ang_inv = -2*PI/n;
    base wn(cos(ang), sin(ang)), wn_inv(cos(ang_inv), sin(ang_inv));
    w    [log_n][0] = 1;
    w_inv[log_n][0] = 1;
    for(int i=1;i<n;i++){
        w    [log_n][i] = w    [log_n][i-1] * wn;
        w_inv[log_n][i] = w_inv[log_n][i-1] * wn_inv;
    }

    for(int i=log_n-1;i>0;i--){
        for(int j=0;j<(1<<i);j++){
            w    [i][j] = w    [i+1][j<<1];
            w_inv[i][j] = w_inv[i+1][j<<1];
        }
    }

}
void fft_thread(vector<base> &a, int starti, int startj, int cnt,  int len, bool invert){

}
void fft(vector<base> &a, bool invert){
    int n = (int)a.size();
    int but_cnt = n/2;
    int but_in_threads = max(1, but_cnt/numThreads);
    vector<thread> threads(numThreads);

    for(int i=0;i<n;i++)
        if(i < rev[i])
            swap(a[i], a[rev[i]]);

    for(int len=2;len<n;len<<=1){
        threads.clear();
        if(len/2 > but_in_threads){
            for(int i=0;i<n;i+=len){
                for(int j=0;j<len/2;j+=but_in_thread){

                }
            }
        }else{
            for(int i=0;i<n;i+=but_in_thread*2){

            }
        }
    }

    if(invert)
        for(base &i:a)
            i /= n;

}

vector<int> mul(const vector<int> &a, const vector<int> &b){
    vector<base> ta(a.begin(), a.end()),
                 tb(b.begin(), b.end());
    int n = 1;
    while(n < max(a.size(), b.size()))
        n<<=1;
    n<<=1;
    bit(n);
    ta.resize(n);
    tb.resize(n);

    fft(ta, 0);
    fft(tb, 0);
    for(int i=0;i<n;i++)
        ta[i] *= tb[i];
    fft(ta, 1);

    vector<int> res(n);
    for(int i=0;i<n;i++)
        res[i] = (int)(ta[i].real() + 0.5);
    return res;
}








int main(){

}
//TODO: fft one thread
//TODO: fft function

