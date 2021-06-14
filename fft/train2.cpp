#include<bits/stdc++.h>
using namespace std;
const double PI = acos(-1);
typedef complex<double> base;

const int maxn = 1<<20;
//const int n = 1 << 13;
const int cs = 4;
const int numThreads = 2;

int rev[maxn];
base a[maxn], b[maxn];
base fa[maxn], fb[maxn];
base w_pow[maxn];
base w_pow2[maxn];

void calc_rev(int n) {
    for(int i=1, j=0;i<n;i++){
        int bit = n>>1;
        for(;j>=bit;bit>>=1)
            j -= bit;
        j += bit;
        rev[i] = j;
    }
}
void syncronize(int i, int other, atomic<bool> *done){
    done[i] = 1;
    if(i == other)
        return;

    while(done[other] == 0)
        this_thread::yield();
}
void one_thread_sequence(base *a, int i, int len, int thread_id){
    for(int j=0;j<len/2;j++){
        string s = to_string(thread_id)+"<-id "+to_string(len)+" "+to_string(i)+" "+to_string(j)+"\n";cout<<s;
        base a0 = a[i+j], a1 = w_pow[j]*a[i+j+len/2];
        a[i+j]       = a0 + a1;
        a[i+j+len/2] = a0 - a1;
    }
}
void one_thread(base *a, int i, int len, int thread_id, int thread_size, atomic<bool> *done){
    int but_in_thread = cs/(2*numThreads);//сколько бабочек в одном потоке
    //делим j на порции по размеру кэша(счетчик k)
    for(int k=0;k < (len/2);k+=but_in_thread){
        done[thread_id] = 0;
        for(int j=k;j<min(len/2, k+but_in_thread);j++){
            string s = to_string(thread_id)+"<-id "+to_string(len)+" "+to_string(i)+"\n";cout<<s;
            base a0 = a[i+j], a1 = w_pow[j]*a[i+j+len/2];
            a[i+j]       = a0 + a1;
            a[i+j+len/2] = a0 - a1;
        }
        if(thread_id == thread_size)
            syncronize(thread_id, 0, done);
        syncronize(thread_id, thread_size, done);


        int i0 = (thread_id%2==0)? i : i-len;
        for(int j=k + (((i/len)%2==0)?0:len/2) ;j<min(len, k + (((i/len)%2==0)?0:len/2) + but_in_thread);j++){
            string s = to_string(thread_id)+"<-id "+to_string(len*2)+" "+to_string(i0)+" "+to_string(j)+"\n";cout<<s;
            base a0 = a[i0+j], a1 = w_pow2[j]*a[i0+j+len];
            a[i0+j]     = a0 + a1;
            a[i0+j+len] = a0 - a1;
        }
        syncronize(thread_id, min(thread_size, thread_id+1), done);
    }
}
void fft(base a[], int n, bool invert){
    calc_rev(n);
    for (int i=1, j=0; i<n; ++i)
        if(i < rev[i])
            swap (a[i], a[rev[i]]);

    vector<thread> threads;
    atomic<bool> done[numThreads] = {};

    for(int len=2;len<=n;len<<=2){
            double ang = 2*PI/len*(invert?-1:1);
            base wn(cos(ang), sin(ang));
            w_pow[0] = base(1, 0);
            for(int i=1;i<len/2;i++)
                w_pow[i] = w_pow[i-1]*wn;
            double ang2 = 2*PI/(2*len)*(invert?-1:1);
            base wn2(cos(ang2), sin(ang2));
            w_pow2[0] = w_pow[0];
            for(int i=1;i<len;i++)
                w_pow2[i] = w_pow2[i-1]*wn2;
            threads.clear();

            //TODO: костыль
//            if(len <= cs/2){
//                for(int i=0;i<n;i+=len)
//                    threads.emplace_back( one_thread_sequence, a, i, len, 0);
//                for(int i=0;i<n;i+=len*2)
//                    threads.emplace_back( one_thread_sequence, a, i, len*2, 0);
//            }else{

                for(int i=0, cur_th=0;i<n;i+=len){
                    threads.emplace_back( one_thread, a, i, len, cur_th, min(n/len, numThreads-1), done);

                    if(cur_th ==  min(n/len, numThreads-1)){
                        for(int th=0;th<numThreads;th++)
                            threads[th].join();
                        cur_th = 0;
                        threads.clear();
                    }else
                        cur_th++;
                }

//            }
            for(thread &th:threads)
                th.join();

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

int main(){


    vector<int> res;
    int n = 8;
    for(int i=0;i<n;i++)a[i] = 1, b[i] = 1;

    auto start = std::chrono::steady_clock::now();

    mul(a, b, n, n,res);

    auto finish = std::chrono::steady_clock::now();
    std::chrono::duration<double> time = finish - start;
    cout<<std::chrono::duration_cast<std::chrono::milliseconds>(time).count();
    cout<<endl;
    for(int i:res)    cout<<i<<" ";

    return 0;
}

