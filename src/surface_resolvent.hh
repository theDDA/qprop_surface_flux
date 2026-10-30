#ifndef QPROP_SURFACE_RESOLVENT_HH
#define QPROP_SURFACE_RESOLVENT_HH
#include "qprop_linalg.hh"
#include <algorithm>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

// Responses use the ordinary transpose, without complex conjugation.
class SurfaceResolvent {
    using z = std::complex<double>;
    int N, S, L, C, j, dimension;
    double dr, charge, corner_d, corner_m;
    std::vector<z> radial, rhs;
    std::vector<double> derivative, centrifugal;
    std::vector<z> lo, di, up, potential_values, core_lo, core_di, core_up;
    std::vector<z> sources, dl, d, du, du2, solution, residual;
    std::vector<int> piv;
    long fallbacks = 0, refinements = 0;
    long direct_solves = 0, condensed_solves = 0, response_solves = 0;
    double largest_residual = 0;
    bool full_grid = std::getenv("QPROP_SURFACE_FULL_GRID") != nullptr;

    static qprop_lapack_complex *lapack(std::vector<z>& a) {
        return reinterpret_cast<qprop_lapack_complex*>(a.data());
    }

    const std::vector<z>& solve_system(const std::vector<z>& lower,
                                      const std::vector<z>& diagonal,
                                      const std::vector<z>& upper,
                                      const std::vector<z>& right_hand_sides,
                                      int nrhs, char transpose) {
        int n = int(diagonal.size()), info = 0;
        dl = lower; d = diagonal; du = upper;
        du2.resize(std::max(1,n-2)); piv.resize(n);
        qprop_zgttrf(&n, lapack(dl), lapack(d), lapack(du), lapack(du2), piv.data(), &info);
        if (info) throw std::runtime_error("Singular radial resolvent factorization");
        solution = right_hand_sides;
        residual.resize(nrhs*n);
        qprop_zgttrs(&transpose, &n, &nrhs, lapack(dl), lapack(d), lapack(du),
                lapack(du2), piv.data(), lapack(solution), &n, &info);
        if (info) throw std::runtime_error("Radial solve failed");
        const auto& sub = transpose == 'T' ? upper : lower;
        const auto& super = transpose == 'T' ? lower : upper;
        for (int pass = 0; pass < 2; ++pass) {
            double error = 0, scale = 0;
            for (int q = 0; q < nrhs; ++q)
                for (int i = 0; i < n; ++i) {
                    const int k = q*n+i;
                    z az = diagonal[i]*solution[k];
                    double row_scale = std::abs(diagonal[i])*std::abs(solution[k]);
                    if (i) { az += sub[i-1]*solution[k-1]; row_scale += std::abs(sub[i-1])*std::abs(solution[k-1]); }
                    if (i+1<n) { az += super[i]*solution[k+1]; row_scale += std::abs(super[i])*std::abs(solution[k+1]); }
                    residual[k] = right_hand_sides[k] - az;
                    error = std::max(error, std::abs(residual[k]));
                    scale = std::max(scale, row_scale+std::abs(right_hand_sides[k]));
                }
            const double eta = error / std::max(scale, std::numeric_limits<double>::min());
            if (pass || eta <= 8*std::numeric_limits<double>::epsilon()) {
                largest_residual = std::max(largest_residual, eta);
                if (eta > 1e-11) throw std::runtime_error("Radial solve residual exceeds tolerance");
                break;
            }
            qprop_zgttrs(&transpose, &n, &nrhs, lapack(dl), lapack(d), lapack(du),
                    lapack(du2), piv.data(), lapack(residual), &n, &info);
            if (info) throw std::runtime_error("Radial query refinement failed");
            for (int i = 0; i < nrhs*n; ++i) solution[i] += residual[i];
            ++refinements;
        }
        return solution;
    }

public:
    template<class Scalar, class Absorber>
    SurfaceResolvent(grid g, grid small, const wavefunction& wf,
                     double Z, Scalar potential, Absorber absorber, double radius)
        : N(int(g.ngps_x())), L(int(g.ngps_y())), j(int(g.rindex(radius))),
          dimension(g.dimens()), dr(g.delt_x()), charge(Z) {
        C = dimension == 44 ? L*L : L;
        corner_d = -2/(dr*dr)*(1-charge*dr/(12-10*charge*dr));
        corner_m = -2*(1+dr*dr/12*corner_d);
        int extent = 0;
        for (int l = 0; l < L; ++l) {
            const int ml = dimension == 44 ? l : 0;
            for (int m = -ml; m <= ml; ++m)
                for (int r = 0; r < small.ngps_x(); ++r)
                    if (wf[small.index(r,l,m,0)] != z(0)) extent = std::max(extent,r+1);
        }
        S = std::min(N,std::max(j+1,extent+1)); // One extra Numerov source row.
        rhs.assign(C*S,z(0));
        for (int l = 0; l < L; ++l) {
            const int ml = dimension == 44 ? l : 0;
            for (int m = -ml; m <= ml; ++m) {
                const int c = dimension == 44 ? l*(l+1)+m : l;
                auto value = [&](int r) { return r >= 0 && r < small.ngps_x() ? wf[small.index(r,l,m,0)] : z(0); };
                for (int r = 0; r < S; ++r)
                    rhs[c*S+r] = (r == 0 && l == 0 ? corner_m : -10./6.)*value(r)
                               - (value(r-1)+value(r+1))/6.;
            }
        }
        radial.resize(N); centrifugal.resize(N);
        for (int r = 0; r < N; ++r) {
            const double x = g.r(r);
            radial[r] = z(potential(x,0,0,0,0),-absorber(r,0,0,0,g));
            centrifugal[r] = 0.5/(x*x);
        }
        // Row j of the compact derivative D = C_1^{-1} Delta_1.
        const double y = sqrt(3.)-2.;
        std::vector<double> diag(N,4./6.), v(N,0), factor(N,0);
        diag.front() = diag.back() = (4+y)/6.;
        v[j] = 1;
        for (int i = 1; i < N; ++i) {
            factor[i] = (1./6.)/diag[i-1];
            diag[i] -= factor[i]/6.;
            v[i] -= factor[i]*v[i-1];
        }
        v[N-1] /= diag[N-1];
        for (int i = N-2; i >= 0; --i) v[i] = (v[i]-v[i+1]/6.)/diag[i];
        derivative.resize(N);
        derivative[0] = (y*v[0]-v[1])*0.5/dr;
        for (int i = 1; i < N-1; ++i) derivative[i] = (v[i-1]-v[i+1])*0.5/dr;
        derivative[N-1] = (v[N-2]-y*v[N-1])*0.5/dr;
        std::cout << "Surface resolvent: N=" << N << ", support=" << S
                  << ", channels=" << C << ", exact source extent=" << extent << '\n';
    }

    void apply(double energy, wavefunction& value, wavefunction& deriv) {
        lo.resize(N-1); di.resize(N); up.resize(N-1); potential_values.resize(N);
        for (int l = 0; l < L; ++l) {
            for (int r = 0; r < N; ++r) potential_values[r] = radial[r] + double(l*(l+1))*centrifugal[r] - energy;
            for (int r = 0; r < N; ++r) di[r] = -2/(dr*dr) - (10./6.)*potential_values[r];
            if (l == 0) di[0] = corner_d + corner_m*potential_values[0];
            for (int r = 0; r < N-1; ++r) {
                up[r] = 1/(dr*dr) - potential_values[r+1]/6.;
                lo[r] = 1/(dr*dr) - potential_values[r]/6.;
            }
            int n = full_grid ? N : S;
            z p(0), q(0);
            bool safe = true;
            if (S < N && !full_grid) {
                p = di[N-1]; q = derivative[N-1];
                for (int r = N-2; r >= S-1; --r) {
                    if (std::abs(p) <= 1e-14*(std::abs(di[r+1])+std::abs(lo[r]))) { safe = false; break; }
                    if (r == S-1) break;
                    q = derivative[r] - lo[r]*q/p;
                    p = di[r] - up[r]*lo[r]/p;
                }
                if (!safe) { n = N; ++fallbacks; }
            }
            core_di.assign(di.begin(),di.begin()+n);
            core_lo.assign(lo.begin(),lo.begin()+n-1);
            core_up.assign(up.begin(),up.begin()+n-1);
            z derivative_correction(0);
            if (n < N) {
                core_di[n-1] -= up[n-1]*lo[n-1]/p;
                derivative_correction = lo[n-1]*q/p;
            }
            const int first = dimension == 44 ? l*l : l;
            const int last = dimension == 44 ? (l+1)*(l+1) : l+1;
            if (dimension == 34) {
                sources.assign(n,z(0));
                std::copy(rhs.begin()+first*S,rhs.begin()+(first+1)*S,sources.begin());
                const auto& state = solve_system(core_lo,core_di,core_up,sources,1,'N');
                value[first] = state[j];
                z slope(0);
                for (int r = 0; r < n; ++r) slope += derivative[r]*state[r];
                deriv[first] = slope-derivative_correction*state[n-1];
                if (n < N) ++condensed_solves;
                else ++direct_solves;
            } else {
                sources.assign(2*n,z(0));
                sources[j] = 1;
                for (int r = 0; r < n; ++r) sources[n+r] = derivative[r];
                sources[2*n-1] -= derivative_correction;
                const auto& rows = solve_system(core_lo,core_di,core_up,sources,2,'T');
                for (int c = first; c < last; ++c) {
                    z a(0), b(0);
                    for (int r = 0; r < S; ++r) {
                        a += rows[r]*rhs[c*S+r];
                        b += rows[n+r]*rhs[c*S+r];
                    }
                    value[c] = a; deriv[c] = b;
                }
                ++response_solves;
            }
        }
    }

    ~SurfaceResolvent() {
        std::cout << "Surface query residual max=" << largest_residual
                  << ", full-grid fallbacks=" << fallbacks
                  << ", refinements=" << refinements << '\n';
        std::cout << "Surface branches: direct=" << direct_solves
                  << ", condensed=" << condensed_solves
                  << ", shared=" << response_solves << '\n';
    }
};
#endif
