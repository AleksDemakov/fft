 #ifndef POLYNOMIAL_HPP_INCLUDED
 #define POLYNOMIAL_HPP_INCLUDED

 #include <vector>

 #include "fft.hpp"

 template < class Field >
 class const_poly_stub {
    const Field* arr;
    int size;
 public:
    const_poly_stub(const Field* a, int n) : arr(a), size(n) {}
    Field operator[](int i) const {
        if ( i >= size )
            return 0;
        return arr[i];
    }
 };

 template < class Poly, class Field >
 std::vector<Field> poly_mul(const Poly& a, const Poly& b, const Field* w, int N) {
    Field* x = new Field[N];
    Field* y = new Field[N];
    std::vector<Field> c(N);

    fft(x, const_poly_stub(&a[0], a.size()), w, N);
    fft(y, const_poly_stub(&b[0], b.size()), w, N);

    for ( int i = 0; i < N; i++ )
        x[i] *= y[i];

    fft_inv(&c[0], x, w, N);

    delete [] x;
    delete [] y;

    return std::move(c);
 }

 template < class Poly, class Field >
 std::vector<Field> poly_mul_iterative(const Poly& a, const Poly& b, const Field* w, int N) {
    Field* x = new Field[N];
    Field* y = new Field[N];
    std::vector<Field> c(N);

    fft_iterative(x, const_poly_stub(&a[0], a.size()), w, N);
    fft_iterative(y, const_poly_stub(&b[0], b.size()), w, N);

    for ( int i = 0; i < N; i++ )
        x[i] *= y[i];

    fft_inv_iterative(&c[0], x, w, N);

    delete [] x;
    delete [] y;

    return std::move(c);
 }

 template < class Poly, class Field >
 std::vector<Field> poly_mul_parallel(const Poly& a, const Poly& b, const Field* w, int N) {
    Field* x = new Field[N];
    Field* y = new Field[N];
    std::vector<Field> c(N);

    fft_parallel(x, const_poly_stub(&a[0], a.size()), w, N);
    fft_parallel(y, const_poly_stub(&b[0], b.size()), w, N);

    for ( int i = 0; i < N; i++ )
        x[i] *= y[i];

    fft_inv_parallel(&c[0], x, w, N);

    delete [] x;
    delete [] y;

    return std::move(c);
 }

 template < typename T >
 T gorner(const std::vector<T>& a, const T& x) {
    const int n = a.size();
    T v = a[n-1];
    for ( int i = n-1; i--; )
        (v *= x) += a[i];
    return v;
 }

 #endif
