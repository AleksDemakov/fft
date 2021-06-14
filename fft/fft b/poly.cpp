 #include <iostream>
 #include <vector>
 #include <cassert>

 #include "polynomial.hpp"
 #include "prime-field.hpp"
 #include "../mm/time_measure.hpp"

 template < typename T >
 std::vector<T> conv_mul(const std::vector<T>& a, const std::vector<T>& b) {
    const int n = a.size();
    const int m = b.size();

    if ( n < m )
        return std::move(conv_mul(b, a));

    std::vector<T> c(n+m-1, 0);

    for ( int i = 0; i < m-1; i++ )
        for ( int j = 0; j <= i; j++ )
            c[i] += a[j] * b[i-j];

    for ( int i = m-1; i < n-1; i++ )
        for ( int j = 0; j < m; j++ )
            c[i] += a[i-j] * b[j];

    for ( int i = n-1; i < n+m-1; i++ )
        for ( int j = 1; j < m+n-i; j++ )
            c[i] += a[n-j] * b[i+j-n];

    return std::move(c);
 }

 int main() {

    const int N = 1 << 24;
    const int p = 469762049;

    const int n = N / 2;
    assert((p-1) % N == 0);

    const PrimeField<p> xi = PrimeField<p>::primitive_element();
    std::cout
        << " w = " << xi
        << "; w^((p-1)/2) = " << dicho(xi, (p-1)/2)
        << "; w^(p-1) = "     << dicho(xi, p-1)
        << std::endl
    ;

    PrimeField<p> omega = dicho(xi, (p-1)/N);

    PrimeField<p>* w = new PrimeField<p>[N];
    w[0] = 1;
    for ( int i = 1; i < N; i++ )
        w[i] = w[i-1] * omega;

    std::vector<PrimeField<p>> a(n), b(n);
    for ( auto& t : a )
        t = rand();
    for ( auto& t : b )
        t = rand();

    time_measure tm1;
    auto c = poly_mul(a, b, w, N);
    tm1.stop();

    time_measure tm2;
    auto d = a; // conv_mul(a, b);
    tm2.stop();

    delete [] w;
    w = 0;

    std::cout << " N = " << N << " om = " << omega << "; fft: " << tm1 << "; conv: " << tm2 << std::endl;

    const_poly_stub<PrimeField<p>> cc(&c[0], c.size());
    const_poly_stub<PrimeField<p>> dd(&d[0], d.size());
    bool flag = true;
    for ( int i = 0; i < N; i++ )
        flag = flag && cc[i] == dd[i];
    assert(flag);

    return 0;
 }
