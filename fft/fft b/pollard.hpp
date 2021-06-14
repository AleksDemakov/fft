 #ifndef POLLARD_HPP_INCLUDED
 #define POLLARD_HPP_INCLUDED

 #include <vector>
 #include <cassert>
 
 #include "fft.hpp"
 #include "polynomial.hpp"
 #include "prime-field.hpp"

 template <
    unsigned p1 = 27 * (1 << 26) + 1,
    unsigned p2 =  7 * (1 << 26) + 1,
    unsigned p3 = 33 * (1 << 25) + 1
 >
 std::vector<unsigned> pollard3(const std::vector<unsigned>& a, const std::vector<unsigned>& b) {

    const unsigned n = 1 << log2ceil(a.size() + b.size());

    assert(a.size());
    assert(b.size());
    assert( (p1-1) % n == 0 );
    assert( (p2-1) % n == 0 );
    assert( (p3-1) % n == 0 );

    auto w1 = get_roots_of_unity<p1>(n);
    auto w2 = get_roots_of_unity<p2>(n);
    auto w3 = get_roots_of_unity<p3>(n);

    // std::clog << " w<" << p1 << ">:"; for ( int i = 0; i < n; i++ ) std::clog << ' ' << w1[i]; std::clog << std::endl;
    // std::clog << " w<" << p2 << ">:"; for ( int i = 0; i < n; i++ ) std::clog << ' ' << w2[i]; std::clog << std::endl;
    // std::clog << " w<" << p3 << ">:"; for ( int i = 0; i < n; i++ ) std::clog << ' ' << w3[i]; std::clog << std::endl;

    auto v1 = poly_mul(a, b, w1, n);
    auto v2 = poly_mul(a, b, w2, n);
    auto v3 = poly_mul(a, b, w3, n);

    // std::clog << " v<" << p1 << ">:"; for ( auto t: v1 ) std::clog << ' ' << t; std::clog << std::endl;
    // std::clog << " v<" << p2 << ">:"; for ( auto t: v2 ) std::clog << ' ' << t; std::clog << std::endl;
    // std::clog << " v<" << p3 << ">:"; for ( auto t: v3 ) std::clog << ' ' << t; std::clog << std::endl;

    auto N2 = p1;
    auto N3 = (unsigned long long)p1 * p2;

    auto C2 = PrimeField<p2>(1) / PrimeField<p2>(N2);
    auto C3 = PrimeField<p3>(1) / PrimeField<p3>(N3);

    std::vector<unsigned> u(n);
    unsigned __int128 d = 0;
    for ( int i = 0; i < n; i++ ) {
        unsigned __int128 w = v1[i].to_int();
        w += (unsigned long long)N2 * ((v2[i] - w) * C2).to_int();
        w += (unsigned __int128 )N3 * ((v3[i] - w) * C3).to_int();
        u[i] = d + w; // Само округлится, так как unsigned 32-битный.
        d = ( d + w ) >> 32;
    }

    return std::move(u);
 }

 template <
    unsigned p1 = 27 * (1 << 26) + 1,
    unsigned p2 =  7 * (1 << 26) + 1,
    unsigned p3 = 33 * (1 << 25) + 1
 >
 std::vector<unsigned> pollard3_iterative(const std::vector<unsigned>& a, const std::vector<unsigned>& b) {

    const unsigned n = 1 << log2ceil(a.size() + b.size());

    assert(a.size());
    assert(b.size());
    assert( (p1-1) % n == 0 );
    assert( (p2-1) % n == 0 );
    assert( (p3-1) % n == 0 );

    auto w1 = get_roots_of_unity<p1>(n);
    auto w2 = get_roots_of_unity<p2>(n);
    auto w3 = get_roots_of_unity<p3>(n);

    // std::clog << " w<" << p1 << ">:"; for ( int i = 0; i < n; i++ ) std::clog << ' ' << w1[i]; std::clog << std::endl;
    // std::clog << " w<" << p2 << ">:"; for ( int i = 0; i < n; i++ ) std::clog << ' ' << w2[i]; std::clog << std::endl;
    // std::clog << " w<" << p3 << ">:"; for ( int i = 0; i < n; i++ ) std::clog << ' ' << w3[i]; std::clog << std::endl;

    auto v1 = poly_mul_iterative(a, b, w1, n);
    auto v2 = poly_mul_iterative(a, b, w2, n);
    auto v3 = poly_mul_iterative(a, b, w3, n);

    // std::clog << " v<" << p1 << ">:"; for ( auto t: v1 ) std::clog << ' ' << t; std::clog << std::endl;
    // std::clog << " v<" << p2 << ">:"; for ( auto t: v2 ) std::clog << ' ' << t; std::clog << std::endl;
    // std::clog << " v<" << p3 << ">:"; for ( auto t: v3 ) std::clog << ' ' << t; std::clog << std::endl;

    auto N2 = p1;
    auto N3 = (unsigned long long)p1 * p2;

    auto C2 = PrimeField<p2>(1) / PrimeField<p2>(N2);
    auto C3 = PrimeField<p3>(1) / PrimeField<p3>(N3);

    std::vector<unsigned> u(n);
    unsigned __int128 d = 0;
    for ( int i = 0; i < n; i++ ) {
        unsigned __int128 w = v1[i].to_int();
        w += (unsigned long long)N2 * ((v2[i] - w) * C2).to_int();
        w += (unsigned __int128 )N3 * ((v3[i] - w) * C3).to_int();
        u[i] = d + w; // Само округлится, так как unsigned 32-битный.
        d = ( d + w ) >> 32;
    }

    return std::move(u);
 }

 template <
    unsigned p1 = 27 * (1 << 26) + 1,
    unsigned p2 =  7 * (1 << 26) + 1,
    unsigned p3 = 33 * (1 << 25) + 1
 >
 std::vector<unsigned> pollard3_parallel(const std::vector<unsigned>& a, const std::vector<unsigned>& b) {

    const unsigned n = 1 << log2ceil(a.size() + b.size());

    assert(a.size());
    assert(b.size());
    assert( (p1-1) % n == 0 );
    assert( (p2-1) % n == 0 );
    assert( (p3-1) % n == 0 );

    auto w1 = get_roots_of_unity<p1>(n);
    auto w2 = get_roots_of_unity<p2>(n);
    auto w3 = get_roots_of_unity<p3>(n);

    // std::clog << " w<" << p1 << ">:"; for ( int i = 0; i < n; i++ ) std::clog << ' ' << w1[i]; std::clog << std::endl;
    // std::clog << " w<" << p2 << ">:"; for ( int i = 0; i < n; i++ ) std::clog << ' ' << w2[i]; std::clog << std::endl;
    // std::clog << " w<" << p3 << ">:"; for ( int i = 0; i < n; i++ ) std::clog << ' ' << w3[i]; std::clog << std::endl;

    auto v1 = poly_mul_parallel(a, b, w1, n);
    auto v2 = poly_mul_parallel(a, b, w2, n);
    auto v3 = poly_mul_parallel(a, b, w3, n);

    // std::clog << " v<" << p1 << ">:"; for ( auto t: v1 ) std::clog << ' ' << t; std::clog << std::endl;
    // std::clog << " v<" << p2 << ">:"; for ( auto t: v2 ) std::clog << ' ' << t; std::clog << std::endl;
    // std::clog << " v<" << p3 << ">:"; for ( auto t: v3 ) std::clog << ' ' << t; std::clog << std::endl;

    auto N2 = p1;
    auto N3 = (unsigned long long)p1 * p2;

    auto C2 = PrimeField<p2>(1) / PrimeField<p2>(N2);
    auto C3 = PrimeField<p3>(1) / PrimeField<p3>(N3);

    std::vector<unsigned> u(n);
    unsigned __int128 d = 0;
    for ( int i = 0; i < n; i++ ) {
        unsigned __int128 w = v1[i].to_int();
        w += (unsigned long long)N2 * ((v2[i] - w) * C2).to_int();
        w += (unsigned __int128 )N3 * ((v3[i] - w) * C3).to_int();
        u[i] = d + w; // Само округлится, так как unsigned 32-битный.
        d = ( d + w ) >> 32;
    }

    return std::move(u);
 }

 #endif
