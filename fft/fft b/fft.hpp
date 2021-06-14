 #ifndef FFT_HPP_INCLUDED
 #define FFT_HPP_INCLUDED

 #include <mutex>
 #include <atomic>
 #include <thread>

int log2ceil(unsigned n) {
    assert(n);
    int t = 0;
    for ( n = n-1; n; n >>= 1 )
        t++;
    return t;
 }

 template < typename Field, typename ConstPtr = const Field* >
 void fft(Field* y, ConstPtr a, const Field* w, const int N, const int i0 = 0, const int step = 1) {

    if ( N == 1 ) {
        y[0] = a[i0];
        return;
    }

    fft(y      , a, w, N/2, i0                            , step * 2);

    fft(y + N/2, a, w, N/2, i0 + (step < 0 ? -step : step), step * 2);

    auto t = y[N/2];
    y[N/2] = y[0] - t;
    y[0  ] = y[0] + t;

    for ( int i = 1, j = step < 0 ? N*(-step) : 0; i <  N/2; i++ ) {
        auto t = w[j + i*step] * y[i+N/2];
        y[i+N/2] = y[i] - t;
        y[i    ] += t;
    }
 }

 template < typename Field >
 void fft_inv(Field* y, const Field* a, const Field* w, int N ) {
    fft(y, a, w, N, 0, -1);
    Field N_inv = Field(1) / N;
    for ( int i = 0; i < N; i++ )
        y[i] *= N_inv;
 }

 template < typename T >
 void fft_parallel_stage2(T* y, const T* omegas, unsigned t, int p) {
    const unsigned n = 1U << t;
    // di - размер "бабочки".
    for ( unsigned di = 1, s = 0; s < t - 1; s++, di <<= 1 ) {
        auto omega = omegas[s];
        for ( unsigned i0 = 0; i0 < n; i0 += 2*di ) {
            auto tmp = y[i0+di];
            y[i0+di] = y[i0] - tmp;
            y[i0   ] = y[i0] + tmp;
            auto w = omega;
            for ( unsigned i = 1; i < di; i++, w *= omega ) {
                auto tmp = w * y[i0+i+di];
                y[i0+i+di] = y[i0+i] - tmp;
                y[i0+i   ] = y[i0+i] + tmp;
            }
        }

    }
 }

 void syncronize(unsigned p, unsigned q, unsigned i, std::atomic<unsigned>* current_iteration) {
    current_iteration[p] = i;
    if ( p == q )
        return;
    for ( int j = 1000; j--; )
        if ( current_iteration[q] >= i )
            return;
    while ( current_iteration[q] < i ) {
        std::this_thread::yield();
    }
 }

 template < typename T >
 void fft_parallel_stage3(T* y, const T& omega, const T& base_omega, unsigned s, unsigned k, unsigned p, unsigned threads_size, unsigned sb, unsigned sync, std::atomic<unsigned>* current_iteration, bool force_sync) {
    if ( force_sync ) {
        syncronize(p, 2 * p % threads_size, sync + k, current_iteration);
        syncronize(p, (2 * p + 1) % threads_size, sync + k, current_iteration);
    }
    else if ( k )
        syncronize(p, p ^ (1U << (k-1)), sync + k, current_iteration);

    // sb - размер блока.
    const unsigned di = 1U << (s + k); // размер "бабочки"
    const unsigned mask = (2U << k) - 1;
    const unsigned i0 = (p << 1 & ~mask ^ p & mask >> 1) << s;
    auto w = base_omega;
    for ( unsigned i = 0; i < sb; i++, w *= omega ) {
        auto tmp = w * y[i0+i+di];
        y[i0+i+di] = y[i0+i] - tmp;
        y[i0+i   ] = y[i0+i] + tmp;
    }
 }

 template < typename T >
 void fft_one_thread(T* y, T* omegas, unsigned log_n, unsigned p, unsigned threads_log, unsigned cache_log, std::atomic<unsigned>* current_iteration) {
    unsigned sync = 0;
    syncronize(p, p, sync, current_iteration);

    const unsigned n = 1U << log_n;
    const unsigned cache_size = 1U << cache_log;
    const unsigned threads_size = 1U << threads_log;

    assert( cache_size % (2 * threads_size) == 0 );
    assert( cache_log <= log_n );

    auto t_reduced = cache_log - threads_log - 1;
    unsigned su = 2U << (t_reduced + threads_log); // размер участка
    const unsigned sb = 1U << t_reduced;
    for ( unsigned iu = 0; iu < n; iu += su, sync += threads_log + 1 ) {
        auto y_shifted = y + iu;
        assert(y_shifted + p * (2U << t_reduced) + (1U << t_reduced + 1) <= y + n);
        fft_parallel_stage2(y_shifted + p * (2U << t_reduced), omegas, t_reduced + 1, p);
        std::vector<T> base_omega(threads_log + 1, 1);
        for ( unsigned k = 1; k <= threads_log; k++ )
            base_omega[k] = dicho(omegas[k], p % (1 << k));
        for ( unsigned k = 0; k <= threads_log; k++ ) {
            fft_parallel_stage3(y_shifted, omegas[t_reduced + k], base_omega[k], t_reduced, k, p, threads_size, sb, sync, current_iteration, k == 0);
        }
    }
    for ( unsigned s = t_reduced + threads_log + 1; s < log_n; s += threads_log + 1) {
        const unsigned k_min = s + threads_log + 1 <= log_n ? 0 : s + threads_log + 1 - log_n;
        const unsigned pow_s = 1 << s - k_min;
        const unsigned su = 2U << s + threads_log - k_min;
        for ( unsigned iu = 0; iu < n; iu += su ) {
            std::vector<T> base_omega(threads_log + 1, 1);
            for ( unsigned k = 1; k <= threads_log; k++ )
                base_omega[k] = dicho(omegas[k], p % (1 << k));
            for ( unsigned ib = 0; ib < pow_s; ib += sb, sync += threads_log + 1 ) {
                auto y_shifted = y + iu + ib;
                for ( unsigned k = k_min; k <= threads_log; k++ ) {
                    assert(y_shifted + sb <= y + n);
                    fft_parallel_stage3(y_shifted, omegas[s + k - k_min], base_omega[k], s - k_min, k, p, threads_size, sb, sync, current_iteration, k == k_min);
                    if ( ib + sb < pow_s )
                        base_omega[k] *= omegas[s + k - k_min - t_reduced];
                }
            }
        }
    }
 }

 unsigned* get_inv_indexes(unsigned n);

 unsigned next_inv_index(unsigned j, unsigned i, unsigned b) {
    return j ^ (i ^ i+1) << b - __builtin_popcount(i ^ i+1);
 }

 template < typename Field, typename ConstPtr = const Field* >
 void fft_parallel(Field* y, ConstPtr a, const Field* w, const unsigned N, const int i0 = 0, int step = 1) {

    const unsigned log_N = log2ceil(N);

    for ( unsigned i = 0, j = 0; i < N; j = next_inv_index(j, i, log_N), i++ ) {
        if ( i == j )
            y[i] = a[j];
        else if ( i < j ) {
            auto t = a[i];
            y[i] = a[j];
            y[j] = t;
        }
    }

    const unsigned log_cs = 21; // логарифм от размера кеша.
    const unsigned log_P = 4; log2ceil(std::thread::hardware_concurrency());
    const unsigned log_q = log_cs - log_P - 1;

    std::atomic<unsigned> current_iteration[1 << log_P] = {};

    assert(log_q + log_P < log_N);

    std::vector<Field> omegas(log_N);
    auto omega = step < 0 ? w[N-1] : w[1];
    for ( unsigned i = log_N; i--; omega *= omega )
        omegas[i] = omega;

    std::vector<std::thread> threads;

    for ( int p = 1; p < (1 << log_P); p++ )
        threads.emplace_back(fft_one_thread<Field>, y, &omegas[0], log_N, p, log_P, log_cs, current_iteration);

    fft_one_thread(y, &omegas[0], log_N, 0, log_P, log_cs, current_iteration);

    for ( auto& t : threads )
        t.join();
 }

 template < typename Field, typename ConstPtr = const Field* >
 void fft_iterative(Field* y, ConstPtr a, const Field* w, const unsigned N, const int i0 = 0, int step = 1) {

    unsigned* inv = get_inv_indexes(N);
    for ( unsigned i = 0; i < N; i++ ) {
        if ( i == inv[i] )
            y[i] = a[i];
        else if ( i < inv[i] ) {
            auto t = a[i];
            y[i] = a[inv[i]];
            y[inv[i]] = t;
        }
    }

    unsigned j = step < 0 ? N : 0;
    for ( unsigned bat_size = 1, st = N; st >>= 1; bat_size <<= 1) {
        step = step < 0 ? -st : st;
        for ( unsigned first_index = 0; first_index < N; first_index += 2 * bat_size) {
            auto i0 = first_index;
            auto di = bat_size;
            auto t = y[i0+0+di];
            y[i0+0+di] = y[i0+0] - t;
            y[i0+0   ] = y[i0+0] + t;
            for ( int i = 1; i < bat_size; i++ ) {
                auto t = w[j + i*step] * y[i0+i+di];
                y[i0+i+di] = y[i0+i] - t;
                y[i0+i   ] += t;
            }
        }
    }
 }

 template < typename Field >
 void fft_inv_iterative(Field* y, const Field* a, const Field* w, int N ) {
    fft_iterative(y, a, w, N, 0, -1);
    Field N_inv = Field(1) / N;
    for ( int i = 0; i < N; i++ )
        y[i] *= N_inv;
 }

 template < typename Field >
 void fft_inv_parallel(Field* y, const Field* a, const Field* w, int N ) {
    fft_parallel(y, a, w, N, 0, -1);
    Field N_inv = Field(1) / N;
    for ( int i = 0; i < N; i++ )
        y[i] *= N_inv;
 }

 template < int t, typename T >
 void gen_inv_indexes(T arr_inv[]) {

    assert( t < sizeof(unsigned)*8 );

    const unsigned N = 1 << t;

    unsigned* inv = new unsigned[N];

    for ( unsigned i = 0; i < N; i++ ) {
        inv[i] = 0;
        for ( unsigned bit = 0, j = i; bit < t; bit++ ) {
            inv[i] <<= 1;
            inv[i] |= j & 1;
            j >>= 1;
        }
    }

    arr_inv[t] = inv;
 }

 template < int t, typename T >
 T& get_or_gen_inv_indexes(T inv[]) {
    static std::once_flag f;
    if ( inv[t] )
        return inv[t];
    std::call_once(f, gen_inv_indexes<t, T>, inv);
    return inv[t];
 }

 unsigned* get_inv_indexes(unsigned n) {

    static std::atomic<unsigned*> w[32] = {0};

    switch ( n ) {
        case 1U <<  0: return get_or_gen_inv_indexes< 0>(w);
        case 1U <<  1: return get_or_gen_inv_indexes< 1>(w);
        case 1U <<  2: return get_or_gen_inv_indexes< 2>(w);
        case 1U <<  3: return get_or_gen_inv_indexes< 3>(w);
        case 1U <<  4: return get_or_gen_inv_indexes< 4>(w);
        case 1U <<  5: return get_or_gen_inv_indexes< 5>(w);
        case 1U <<  6: return get_or_gen_inv_indexes< 6>(w);
        case 1U <<  7: return get_or_gen_inv_indexes< 7>(w);
        case 1U <<  8: return get_or_gen_inv_indexes< 8>(w);
        case 1U <<  9: return get_or_gen_inv_indexes< 9>(w);
        case 1U << 10: return get_or_gen_inv_indexes<10>(w);
        case 1U << 11: return get_or_gen_inv_indexes<11>(w);
        case 1U << 12: return get_or_gen_inv_indexes<12>(w);
        case 1U << 13: return get_or_gen_inv_indexes<13>(w);
        case 1U << 14: return get_or_gen_inv_indexes<14>(w);
        case 1U << 15: return get_or_gen_inv_indexes<15>(w);
        case 1U << 16: return get_or_gen_inv_indexes<16>(w);
        case 1U << 17: return get_or_gen_inv_indexes<17>(w);
        case 1U << 18: return get_or_gen_inv_indexes<18>(w);
        case 1U << 19: return get_or_gen_inv_indexes<19>(w);
        case 1U << 20: return get_or_gen_inv_indexes<20>(w);
        case 1U << 21: return get_or_gen_inv_indexes<21>(w);
        case 1U << 22: return get_or_gen_inv_indexes<22>(w);
        case 1U << 23: return get_or_gen_inv_indexes<23>(w);
        case 1U << 24: return get_or_gen_inv_indexes<24>(w);
        case 1U << 25: return get_or_gen_inv_indexes<25>(w);
        case 1U << 26: return get_or_gen_inv_indexes<26>(w);
        case 1U << 27: return get_or_gen_inv_indexes<26>(w);
        case 1U << 28: return get_or_gen_inv_indexes<26>(w);
        case 1U << 29: return get_or_gen_inv_indexes<27>(w);
        case 1U << 30: return get_or_gen_inv_indexes<28>(w);
        case 1U << 31: return get_or_gen_inv_indexes<29>(w);
    }
    assert ( false );
    return 0;
 }


 #endif
