#ifndef QPROP_LINALG_HH
#define QPROP_LINALG_HH
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
using qprop_lapack_complex = __CLPK_doublecomplex;
#define qprop_zgttrf zgttrf_
#define qprop_zgttrs zgttrs_
#define qprop_sincos vvsincos
#else
#include <cblas.h>
#include <lapacke.h>
#include <cmath>
#include <type_traits>
static_assert(std::is_same<lapack_int, int>::value, "Use LP64 BLAS/LAPACK (32-bit integers)");
using qprop_lapack_complex = lapack_complex_double;
inline void qprop_zgttrf(int* n, qprop_lapack_complex* dl,
                        qprop_lapack_complex* d, qprop_lapack_complex* du,
                        qprop_lapack_complex* du2, int* piv, int* info) {
    *info = LAPACKE_zgttrf_work(*n, dl, d, du, du2, piv);
}
inline void qprop_zgttrs(char* trans, int* n, int* nrhs,
                        qprop_lapack_complex* dl, qprop_lapack_complex* d,
                        qprop_lapack_complex* du, qprop_lapack_complex* du2,
                        int* piv, qprop_lapack_complex* b, int* ldb, int* info) {
    *info = LAPACKE_zgttrs_work(LAPACK_COL_MAJOR, *trans, *n, *nrhs,
                               dl, d, du, du2, piv, b, *ldb);
}
inline void qprop_sincos(double* s, double* c, const double* x, const int* n) {
    for (int i = 0; i < *n; ++i) { s[i] = std::sin(x[i]); c[i] = std::cos(x[i]); }
}
#endif
#endif
