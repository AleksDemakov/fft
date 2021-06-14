 #include <iostream>
 #include <vector>
 #include <string>
 #include <cstdlib>

 #include "pollard.hpp"
 #include "../mm/time_measure.hpp"

 template < unsigned long long b, class T  >
 std::string to_dec(T a) {
    std::string s;
    for ( int i = a.size(); i; ) {
        unsigned d = 0;
        for ( int j = i; j--; ) {
            auto t = b * d + a[j];
            d = t % 10;
            a[j] = t / 10;
            if ( a[j] == 0 && i == j+1 )
                i--;
        }
        s.push_back('0' + d);
    }
    for ( int j = s.length()-1, i = (j+1) / 2; i--; )
        std::swap(s[i], s[j-i]);
    return s;
 }

 int main() {

    const unsigned long long B = 1ULL << 32;
    const int n = 2000000;

    std::vector<std::vector<unsigned>> a(10), b(10);
    for ( int i = 0; i < a.size() / 2; i++ ) {
        a[2*i+0].resize(n);
        a[2*i+1].resize(n);
        for ( int j = 0; j < n; j++ )
            a[2*i+0][j] = a[2*i+1][j] = rand() % B;
    }
    for ( int i = 0; i < b.size() / 2; i++ ) {
        b[2*i+0].resize(n);
        b[2*i+1].resize(n);
        for ( int j = 0; j < n; j++ )
            b[2*i+0][j] = b[2*i+1][j] = rand() % B;
    }

    for ( int i = 0; i < a.size() / 2; i++ ) {

        time_measure tm1;
        auto u = pollard3_parallel(a[2*i+0], b[2*i+0]);
        tm1.stop();
        std::cout << " Умножение двух " << n << "-разрядных "<< B << "-ичных чисел параллельным БПФ заняло " << tm1.get() << " мс" << std::endl;

        time_measure tm2;
        auto v = pollard3_iterative(a[2*i+1], b[2*i+1]);
        tm2.stop();
        std::cout << " Умножение двух " << n << "-разрядных "<< B << "-ичных чисел итеративным БПФ заняло " << tm2.get() << " мс" << std::endl;

        if ( n == 2 ) {
            std::clog << " a[" << 2 * i + 0 << "]:"; for ( auto& t : a[2*i+0] ) std::clog << ' ' << t; std::clog << std::endl;
            std::clog << " a[" << 2 * i + 1 << "]:"; for ( auto& t : a[2*i+1] ) std::clog << ' ' << t; std::clog << std::endl;
            std::clog << " b[" << 2 * i + 0 << "]:"; for ( auto& t : b[2*i+0] ) std::clog << ' ' << t; std::clog << std::endl;
            std::clog << " b[" << 2 * i + 1 << "]:"; for ( auto& t : b[2*i+1] ) std::clog << ' ' << t; std::clog << std::endl;
            std::clog << " u:"; for ( auto& t : u ) std::clog << ' ' << t; std::clog << std::endl;
            std::clog << " v:"; for ( auto& t : v ) std::clog << ' ' << t; std::clog << std::endl;
        }

        assert( u.size() == v.size() );
        for ( int j = 0; j < u.size(); j++ )
            if ( u[j] != v[j] ) {
                std::cerr << " " << j << " ";
                assert( u[j] == v[j] );
            }
    }

    return 0;
 }
