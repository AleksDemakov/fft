 #include <iostream>
 #include <vector>

 #include "prime-field.hpp"
 #include "dicho.hpp"
 #include "fft.hpp"

 int main() {
    const int p = 5;

    {

    int N = 4;
    std::vector<PrimeField<p>> a = { 0, 1, 2, 3 };
    std::vector<PrimeField<p>> y(N), z(N);
    std::vector<PrimeField<p>> w = { 1, 2, 2*2, 2*2*2 };

    for ( auto t : a )
        std::cout << ' ' << t;
    std::cout << std::endl;

    fft(&y[0], &a[0], &w[0], N);

    for ( auto t : y )
        std::cout << ' ' << t;
    std::cout << std::endl;

    fft_inv(&z[0], &y[0], &w[0], N);

    for ( auto t : z )
        std::cout << ' ' << t;
    std::cout << std::endl;

    }

    std::cout << std::endl;

    PrimeField<p> a(0), b(0);

    std::cout << " * |";
    for ( int j = 0; j < p; ++j, ++a )
        std::cout << ' ' << a;
    std::cout << std::endl;

    std::cout << "---+";
    for ( int j = 0; j < p; ++j, ++a )
        std::cout << "--";
    std::cout << std::endl;

    for ( int i = 0; i < p; ++i, ++b ) {
        std::cout << ' ' << b << " |";
        a = 0;
        for ( int j = 0; j < p; ++j, ++a )
            std::cout << ' ' << a*b;
        std::cout << std::endl;
    }

    std::cout << std::endl;

    b = 1;

    std::cout << " / |";
    for ( int j = 0; j < p; ++j, ++a )
        std::cout << ' ' << a;
    std::cout << std::endl;

    std::cout << "---+";
    for ( int j = 0; j < p; ++j, ++a )
        std::cout << "--";
    std::cout << std::endl;

    for ( int i = 1; i < p; ++i, ++b ) {
        std::cout << ' ' << b << " |";
        a = 0;
        for ( int j = 0; j < p; ++j, ++a )
            std::cout << ' ' << a/b;
        std::cout << std::endl;
    }

    std::cout << std::endl;

    a = 0;
    for ( int j = 0; j < p; ++j, ++a ) {
        for ( int n = 0; n < 20; ++n )
            std::cout << ' ' << dicho(a, n);
        std::cout << std::endl;
    }

    return 0;
 }
