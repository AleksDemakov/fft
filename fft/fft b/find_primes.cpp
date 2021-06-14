 #include <iostream>
 #include <vector>
 #include <cassert>

 bool is_prime(unsigned n) {
    static std::vector<unsigned> primes;
    if ( primes.empty() ) {
        primes.push_back(3);
        for ( unsigned p = 5; p < 50'000; p += 2 )
            if ( is_prime(p) )
                primes.push_back(p);
    }
    for ( auto p : primes ) {
        if ( n == p )
            return true;
        if ( p * p > n )
            return true;
        if ( n % p == 0 )
            return false;
    }
    assert(false);
    return true;
 }

 int main() {

    int counter = 10;

    for ( int t = 30; t--; )
        for ( unsigned long long k = 1; k * (1 << t) < 2'000'000'000; k += 2 )
            if ( is_prime(k * (1 << t) + 1) ) {
                std::cout << ' ' << k * (1 << t) + 1 << " = " << k << " x 2^" << t << " + 1" << std::endl;
                if ( !--counter )
                    return 0;
            }

    return 0;
 }
