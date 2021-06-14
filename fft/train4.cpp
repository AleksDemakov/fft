#include<bits/stdc++.h>
using namespace std;
const double PI = acos(-1);
typedef complex<double> base;

const int maxn = 1<<20;
//const int n = 1 << 13;
const int cs = 8;
const int numThreads = 4;

int rev[maxn];
base a[maxn], b[maxn];
base fa[maxn], fb[maxn];
base w_pow[maxn];
base w_pow2[maxn];
condition_variable done;
mutex mut;

atomic<int> done_cnt;
atomic<int> out_cnt;


void calc_rev(int n) {
    for(int i=1, j=0;i<n;i++){
        int bit = n>>1;
        for(;j>=bit;bit>>=1)
            j -= bit;
        j += bit;
        rev[i] = j;
    }
}
void syncronize(int id, int thread_size){
//    cout<<to_string(id)+" "+to_string(thread_size)<<endl;
//    done_cnt++;
//    out_cnt--;
//    unique_lock<mutex> guard(mut);
//    done.wait(guard, [&thread_size, &id](){
//
//    cout<<"|"+to_string(done_cnt)+" "+to_string(id)+" "+to_string(thread_size)+"|predicate"<<endl;
//              return (done_cnt == thread_size);});
//cout<<to_string(id)+"done"<<endl;

}
void one_thread_sequence(base *a, int i, int len, int thread_id, base *w_pow){
    for(int j=0;j<len/2;j++){
//        string s = to_string(thread_id)+"<-id "+to_string(len)+" "+to_string(i)+" "+to_string(j)+"\n";cout<<s;
        base a0 = a[i+j], a1 = w_pow[j]*a[i+j+len/2];
        a[i+j]       = a0 + a1;
        a[i+j+len/2] = a0 - a1;
    }
}

void first(base *a, int i, int k, int len, int thread_id, int thread_size){
    int but_in_thread = cs/(2*thread_size);//сколько бабочек в одном потоке
    for(int j=k;j<min(k+len/2, k+but_in_thread);j++){
//        string s = to_string(thread_id)+"<-id "+to_string(len)+" "+to_string(i)+"\n";cout<<s;
        base a0 = a[i+j], a1 = w_pow[j]*a[i+j+len/2];
        a[i+j]       = a0 + a1;
        a[i+j+len/2] = a0 - a1;
    }
}
void second(base *a, int i, int k, int len, int thread_id, int thread_size){
    int but_in_thread = cs/(2*thread_size);//сколько бабочек в одном потоке
//    cout<<to_string(i)+" "+to_string(len)<<endl;
    int i0 = (thread_id%2==0)? i : i-len;
    for(int j=k + (((i/len)%2==0)?0:len/2) ;j<min(k+len, k + (((i/len)%2==0)?0:len/2) + but_in_thread);j++){
//        string s = to_string(thread_id)+"<-id "+to_string(len*2)+" "+to_string(i0)+" "+to_string(j)+"\n";cout<<s;
        base a0 = a[i0+j], a1 = w_pow2[j]*a[i0+j+len];
        a[i0+j]     = a0 + a1;
        a[i0+j+len] = a0 - a1;
    }
}
void one_thread(base *a, int i, int len, int thread_id, int thread_size){

}
void fft(base a[], int n, bool invert){
    calc_rev(n);
    for (int i=1, j=0; i<n; ++i)
        if(i < rev[i])
            swap (a[i], a[rev[i]]);

    vector<thread> threads;


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
            if(len == n){
                for(int i=0;i<n;i+=len)
                    one_thread_sequence(a, i, len, 0, w_pow);
                break;
            }


            if(len <= cs/2){
                 for(int i=0;i<n;i+=len)
                    one_thread_sequence(a, i, len, 0, w_pow);
                for(int i=0;i<n;i+=len*2)
                    one_thread_sequence(a, i, 2*len, 0, w_pow2);
            }else{
                done_cnt = 0;
                int but_in_thread = cs/(2*min(n/len, numThreads));//сколько бабочек в одном потоке

                for(int k=0;k < (len/2);k+=but_in_thread){
                        threads.clear();
                        for(int i=0,j=0, cur_th=0;i<n;i+=len){
//                            for(int cur_th=0;cur_th<min(n/len, numThreads);cur_th++)
                            threads.emplace_back( first, a, i, k, len, cur_th, min(n/len, numThreads));


                            if(cur_th ==  min(n/len, numThreads)-1){
                                for(thread &th:threads)
                                    th.join();
                                threads.clear();
                                cur_th = 0;
                                for(;j<=i;j+=len){
                                    threads.emplace_back( second, a, j, k, len, cur_th, min(n/len, numThreads));
                                    cur_th++;
                                }
                                for(thread &th:threads)
                                    th.join();
                                threads.clear();
                                cur_th = 0;
                            }else
                                cur_th++;

//                            for(int cur_th=0;cur_th<min(n/len, numThreads);cur_th++)
//                                threads.emplace_back( second, a, i, k, len, cur_th, min(n/len, numThreads));
                        }
//                        threads.clear();
//                        for(int i=0;i<n;i+=len){
//                            for(int cut_th=0;cur_th<min(n/len, numThreads)-1;cur_th++)
//                                threads.emplace_back( second, a, i, k, len, cur_th, min(n/len, numThreads));
//
//                            if(cur_th ==  but_in_thread){
//                                for(thread &th:threads)
//                                    th.join();
//                                cur_th = 0;
//                                threads.clear();
//                            }else
//                                cur_th++;
//                        }




//                    threads.emplace_back( one_thread, a, i, len, cur_th, min(n/len, numThreads));
                }

            }
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
    int n = 1<<13;
    for(int i=0;i<n;i++)a[i] = rand(), b[i] = rand();

    auto start = std::chrono::steady_clock::now();

    mul(a, b, n, n,res);

    auto finish = std::chrono::steady_clock::now();
    std::chrono::duration<double> time = finish - start;
    cout<<std::chrono::duration_cast<std::chrono::milliseconds>(time).count();
    cout<<endl;
//    for(int i:res)    cout<<i<<" ";

    return 0;
}

