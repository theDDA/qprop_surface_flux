vector<cplxd> fft_amplitude;
#ifdef QPROP_DIRECT_ANGULAR
const bool fft_use_direct = true;
#else
const bool fft_use_direct = false;
#endif
struct SurfaceCoupling {
    int index[4] = {-1, -1, -1, -1};
    double coefficient[4] = {0, 0, 0, 0};
};
vector<SurfaceCoupling> fft_couplings;
vector<double> fft_bessel, fft_angular, fft_cos_phi, fft_sin_phi;
vector<double> fft_cos_theta, fft_sin_theta, fft_energy;
vector<cplxd> fft_weights, fft_synthesis, fft_orbital_phase;
vector<size_t> fft_weight_offsets;
vector<cplxd> fft_u, fft_h, fft_product, fft_rings, fft_magnetic;

struct SourceSlice {
    cplxd *pointer;
    cplxd *begin() { return pointer; }
    cplxd *data() { return pointer; }
};
vector<SourceSlice> fft_packed;
vector<cplxd> fft_shared_serial;
cplxd *fft_shared_base = nullptr;
size_t fft_shared_bytes = 0;

void allocate_shared_sources(size_t count) {
    fft_shared_bytes = count*sizeof(cplxd);
    fft_shared_serial.resize(count);
    fft_shared_base = fft_shared_serial.data();
    cout << "SHARED_SOURCE_BYTES=" << fft_shared_bytes
         << " SOURCE_COPIES=1 SOURCE_OWNER=" << (i_proc == 0) << endl;
}

vector<double> fft_phase, fft_sine, fft_cosine;
vector<double> fft_times, fft_ax, fft_ay, fft_az, fft_window;
vector<cplxd> fft_field;
int fft_block_k = 1;
std::unique_ptr<gsl_fft_complex_wavetable, decltype(&gsl_fft_complex_wavetable_free)>
    fft_table{nullptr, gsl_fft_complex_wavetable_free};
std::unique_ptr<gsl_fft_complex_workspace, decltype(&gsl_fft_complex_workspace_free)>
    fft_work{nullptr, gsl_fft_complex_workspace_free};

void prep_fft_constants() {
    const bool full_m = qprop_dim == 44;
    const int L = int(ell_grid_size), Lout = full_m ? L : L+1;
    const int channels = full_m ? L*L : Lout;
    const int np = full_m ? int(num_phi_surff-1) : 1;
    const int mlo = full_m ? 1-L : int(initial_m), mhi = full_m ? L-1 : int(initial_m);
    const int nm = mhi-mlo+1;
    const size_t bytes_per_k = size_t(cache_size_t) *
        (16 + 40*size_t(np) + (fft_use_direct ? 16*size_t(nm) : 0));
    const size_t capacity = std::max(size_t(1), (size_t(1)<<30)/bytes_per_k);
    fft_block_k = 1;
    while (size_t(fft_block_k)*2 <= capacity) fft_block_k *= 2;
    if (size_t(num_k_proc) <= capacity) fft_block_k = std::max(1, int(num_k_proc));
    const int block_k = fft_block_k;
    const int blocks = (num_k_proc+block_k-1)/block_k;
    fft_bessel.resize(size_t(num_k_proc)*(L+1));
    fft_energy.resize(num_k_proc);
    for (int k = 0; k < num_k_proc; ++k) {
        fft_energy[k] = 0.5*k_values[k]*k_values[k];
        fill_bessel_array(k_values[k], R_tsurff, fft_bessel.data()+size_t(k)*(L+1));
    }
    fft_cos_phi.resize(np); fft_sin_phi.resize(np);
    for (int p = 0; p < np; ++p) {
        const double phi = full_m ? phis_surff[p] : 0;
        fft_cos_phi[p] = cos(phi); fft_sin_phi[p] = sin(phi);
    }
    fft_cos_theta.resize(num_theta_surff); fft_sin_theta.resize(num_theta_surff);
    fft_angular.assign(size_t(num_theta_surff)*channels, 0);
    fft_orbital_phase.resize(Lout);
    for (int l = 0; l < Lout; ++l) fft_orbital_phase[l] = pow_neg_i(l+1);
    for (int th = 0; th < num_theta_surff; ++th) {
        fft_cos_theta[th] = cos(thetas_surff[th]);
        fft_sin_theta[th] = sin(thetas_surff[th]);
        for (int l = 0; l < Lout; ++l) {
            const int lo = full_m ? -l : int(initial_m), hi = full_m ? l : int(initial_m);
            for (int m = lo; m <= hi; ++m) {
                if (abs(m) > l) continue;
                const int c = full_m ? l*(l+1)+m : l;
                fft_angular[size_t(th)*channels+c] = Ylm_prefactor(m) *
                    gsl_sf_legendre_sphPlm(l, abs(m), fft_cos_theta[th]);
            }
        }
    }
    const size_t active_channels = full_m ? size_t(L)*L : Lout-abs(int(initial_m));
    fft_weights.resize(2*size_t(num_k_proc)*num_theta_surff*active_channels);
    fft_weight_offsets.reserve(size_t(num_theta_surff)*blocks*nm);
    size_t offset = 0;
    for (int th = 0; th < num_theta_surff; ++th)
        for (int kb = 0; kb < num_k_proc; kb += block_k) {
            const int nk = std::min(block_k, int(num_k_proc)-kb);
            for (int m = mlo; m <= mhi; ++m) {
                const int first = abs(m), nl = Lout-first, inner = 2*nl;
                fft_weight_offsets.push_back(offset);
                for (int a = 0; a < nl; ++a) {
                    const int l = first+a, c = full_m ? l*(l+1)+m : l;
                    const cplxd angular = fft_orbital_phase[l]*fft_angular[size_t(th)*channels+c];
                    for (int k = 0; k < nk; ++k) {
                        const double *j = fft_bessel.data()+size_t(kb+k)*(L+1);
                        fft_weights[offset+k*inner+a] = angular*j[l];
                        fft_weights[offset+k*inner+nl+a] = l < L ?
                            angular*(k_values[kb+k]*j[l+1]-l/R_tsurff*j[l]) : cplxd(0);
                    }
                }
                offset += size_t(nk)*inner;
            }
        }
    if (i_proc == 0) fft_couplings.resize(channels);
    if (i_proc == 0) for (int l = 0; l < Lout; ++l) {
        const int lo = full_m ? -l : int(initial_m), hi = full_m ? l : int(initial_m);
        for (int m = lo; m <= hi; ++m) {
            if (abs(m) > l) continue;
            auto &entry = fft_couplings[full_m ? l*(l+1)+m : l];
            if (full_m) {
                if (l > 0 && m >= -l+2) {
                    entry.index[0]=(l-1)*l+m-1;
                    entry.coefficient[0]=clebsch_mm1_ellm1(l,m);
                }
                if (l > 0 && m <= l-2) {
                    entry.index[1]=(l-1)*l+m+1;
                    entry.coefficient[1]=clebsch_mp1_ellm1(l,m);
                }
                if (l+1 < L) {
                    entry.index[2]=(l+1)*(l+2)+m-1;
                    entry.coefficient[2]=clebsch_mm1_ellp1(l,m);
                    entry.index[3]=(l+1)*(l+2)+m+1;
                    entry.coefficient[3]=clebsch_mp1_ellp1(l,m);
                }
            } else {
                if (l > abs(m)) {
                    entry.index[0]=l-1;
                    entry.coefficient[0]=sqrt(double(l*l-m*m)/((2.*l-1)*(2.*l+1)));
                }
                if (l+1 < L) {
                    entry.index[1]=l+1;
                    entry.coefficient[1]=sqrt(double((l+1)*(l+1)-m*m)/((2.*l+1)*(2.*l+3)));
                }
            }
        }
    }
    if (fft_use_direct) {
        fft_synthesis.resize(size_t(nm)*np);
        for (int m = mlo; m <= mhi; ++m)
            for (int p = 0; p < np; ++p)
                fft_synthesis[size_t(m-mlo)*np+p] =
                    exp(cplxd(0,full_m ? m*phis_surff[p] : 0));
    } else if (np > 1) {
        fft_table.reset(gsl_fft_complex_wavetable_alloc(np));
        fft_work.reset(gsl_fft_complex_workspace_alloc(np));
    }

    if (i_proc == 0) {
        fft_u.resize(size_t(channels)*cache_size_t);
        fft_h.resize(fft_u.size());
    }
    allocate_shared_sources(2*active_channels*cache_size_t);
    fft_packed.resize(nm);
    size_t source_offset = 0;
    for (int m = mlo; m <= mhi; ++m) {
        fft_packed[m-mlo].pointer = fft_shared_base+source_offset;
        source_offset += size_t(2*(Lout-abs(m)))*cache_size_t;
    }
    fft_product.resize(size_t(block_k)*cache_size_t);
    fft_rings.resize(size_t(block_k)*cache_size_t*np);
    fft_phase.resize(fft_rings.size());
    fft_sine.resize(fft_rings.size()); fft_cosine.resize(fft_rings.size());
    if (fft_use_direct) fft_magnetic.resize(size_t(block_k)*cache_size_t*nm);
    fft_times.resize(cache_size_t); fft_ax.resize(cache_size_t);
    fft_ay.resize(cache_size_t); fft_az.resize(cache_size_t);
    fft_window.resize(cache_size_t);
    fft_field.resize(cache_size_t);
}

void prep_fft_space() {
    const long angles = num_theta_surff * (qprop_dim == 44 ? num_phi_surff : 1);
    fft_amplitude.assign(num_k_proc * angles, cplxd(0));
    if (i_proc == 0) {
        psi_R.reset(new cplxd[ell_m_grid_size * cache_size_t]);
        psi_R_deriv.reset(new cplxd[ell_m_grid_size * cache_size_t]);
    }
    const auto setup_tick = std::chrono::steady_clock::now();
    prep_fft_constants();
    cout.precision(17);
    cout << "QPROP_SETUP_WALL="
         << std::chrono::duration<double>(std::chrono::steady_clock::now()-setup_tick).count()
         << " QPROP_WEIGHT_CACHE_BYTES=" << fft_weights.size()*sizeof(cplxd)
         << " QPROP_MOMENTUM_TILE=" << fft_block_k << endl;
    if (i_proc == 0)
        cout << "FFT scheme 1: " << fft_amplitude.size()
             << " accumulated complex amplitudes; time block " << cache_size_t << endl;
}

long read_fft_psi_data() {
    long count = 0;
    if (i_proc == 0) {
        const std::streamsize bytes = size_t(cache_size_t)*ell_m_grid_size*sizeof(cplxd);
        psi.read(reinterpret_cast<char*>(psi_R.get()), bytes);
        const std::streamsize actual = psi.gcount();
        dpsi_dr.read(reinterpret_cast<char*>(psi_R_deriv.get()), actual);
        count = actual/sizeof(cplxd);
        if (actual % (ell_m_grid_size*sizeof(cplxd)) || dpsi_dr.gcount() != actual)
            count = -1;
    }
    if (count < 0) throw std::runtime_error("Mismatched surface records");

    return count;
}

void time_integration_fft(long start_t, long stop_t) {
    const int nt = int(stop_t - start_t);
    if (!nt) return;
    const bool full_m = qprop_dim == 44;
    const int np = full_m ? int(num_phi_surff - 1) : 1;
    const int np_out = full_m ? int(num_phi_surff) : 1;
    const int L = int(ell_grid_size);
    const int Lout = full_m ? L : L + 1; // dim34 retains the top A_z-generated harmonic.
    const int channels = full_m ? L * L : Lout;
    auto &u = fft_u; auto &h = fft_h;
    auto &times = fft_times; auto &ax = fft_ax; auto &ay = fft_ay;
    auto &az = fft_az; auto &win = fft_window;
    for (int t = 0; t < nt; ++t) {
        const double time = delta_t * double(start_t + t);
        times[t] = time;
        ax[t] = vecpot_x.integral(time);
        ay[t] = vecpot_y.integral(time);
        az[t] = vecpot_z.integral(time);
        win[t] = hanning_win(time);
        fft_field[t] = cplxd(vecpot_x(time, 0), vecpot_y(time, 0));
        const double Az = vecpot_z(time, 0);
        if (i_proc != 0) continue;
        const cplxd *p = psi_R.get() + t * ell_m_grid_size;
        const cplxd *d = psi_R_deriv.get() + t * ell_m_grid_size;
        if (!full_m) for (int l = 0; l < Lout; ++l) {
            if (l >= labs(initial_m)) {
                cplxd value = l < L ? d[l] - p[l] / R_tsurff : cplxd(0);
                const auto &entry = fft_couplings[l];
                if (entry.index[0] >= 0) {
                    const double c = entry.coefficient[0];
                    value += cplxd(0, 2 * Az * c) * p[l-1];
                }
                if (entry.index[1] >= 0) {
                    const double c = entry.coefficient[1];
                    value += cplxd(0, 2 * Az * c) * p[l+1];
                }
                u[l * nt + t] = l < L ? p[l] : cplxd(0);
                h[l * nt + t] = value;
            }
        }
    }
    if (full_m && i_proc == 0)
        for (int tb = 0; tb < nt; tb += 4)
            for (int cb = 0; cb < channels; cb += 32)
                for (int t = tb; t < std::min(tb+4,nt); ++t) {
                    const cplxd A = fft_field[t];
                    const cplxd *p = psi_R.get()+t*ell_m_grid_size;
                    const cplxd *d = psi_R_deriv.get()+t*ell_m_grid_size;
                    for (int c = cb; c < std::min(cb+32,channels); ++c) {
                        const auto &entry = fft_couplings[c];
                        cplxd coupling(0);
                        if (entry.index[0] >= 0)
                            coupling -= conj(A)*entry.coefficient[0]*p[entry.index[0]];
                        if (entry.index[1] >= 0)
                            coupling += A*entry.coefficient[1]*p[entry.index[1]];
                        if (entry.index[2] >= 0)
                            coupling += conj(A)*entry.coefficient[2]*p[entry.index[2]];
                        if (entry.index[3] >= 0)
                            coupling -= A*entry.coefficient[3]*p[entry.index[3]];
                        u[c*nt+t] = p[c];
                        h[c*nt+t] = d[c]-p[c]/R_tsurff+cplxd(0,sqrt(2.0))*coupling;
                    }
                }
    const auto &cp = fft_cos_phi; const auto &sp = fft_sin_phi;
    const int block_k = fft_block_k;
    const cplxd one(1), zero(0);
    const int mlo = full_m ? 1-L : int(initial_m);
    const int mhi = full_m ? L-1 : int(initial_m);
    const int nm = mhi-mlo+1, blocks = (num_k_proc+block_k-1)/block_k;
    auto &packed = fft_packed;
    if (i_proc == 0) for (int m = mlo; m <= mhi; ++m) {
        const int first = abs(m), nl = Lout-first;
        for (int a = 0; a < nl; ++a) {
            const int l = first+a, c = full_m ? l*(l+1)+m : l;
            std::copy(h.begin()+c*nt, h.begin()+(c+1)*nt, packed[m-mlo].begin()+a*nt);
            std::copy(u.begin()+c*nt, u.begin()+(c+1)*nt, packed[m-mlo].begin()+(nl+a)*nt);
        }
    }
    auto &product = fft_product; auto &rings = fft_rings;
    auto &phase = fft_phase; auto &sine = fft_sine; auto &cosine = fft_cosine;
    for (int th = 0; th < num_theta_surff; ++th) {
        const double ct = fft_cos_theta[th], st = fft_sin_theta[th];
        for (int kb = 0; kb < num_k_proc; kb += block_k) {
            const int nk = std::min(block_k, int(num_k_proc) - kb);
            std::fill(rings.begin(), rings.begin() + nk * nt * np, zero);
            for (int m = mlo; m <= mhi; ++m) {
                const int first = abs(m), nl = Lout - first, inner = 2 * nl;
                const cplxd *weights = fft_weights.data()+
                    fft_weight_offsets[(size_t(th)*blocks+kb/block_k)*nm+m-mlo];
                cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                            nk, nt, inner, &one, weights, inner,
                            packed[m-mlo].data(), nt, &zero, product.data(), nt);
                const int bin = full_m ? (m % np + np) % np : 0;
                for (int kt = 0; kt < nk * nt; ++kt) {
                    if (fft_use_direct) fft_magnetic[kt*nm+m-mlo] = product[kt];
                    else rings[kt*np+bin] += product[kt]; // Add aliased signed modes.
                }
            }
            if (fft_use_direct)
                cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                            nk*nt, np, nm, &one, fft_magnetic.data(), nm,
                            fft_synthesis.data(), np, &zero, rings.data(), np);
            for (int k = 0; k < nk; ++k) {
                const double kval = k_values[kb+k], energy = fft_energy[kb+k];
                const double k_transverse = kval*st, k_axial = kval*ct;
                for (int t = 0; t < nt; ++t) {
                    const int offset = (k*nt+t)*np;
                    if (!fft_use_direct && np > 1)
                        gsl_fft_complex_backward(reinterpret_cast<double*>(rings.data()+offset),
                                                 1, np, fft_table.get(), fft_work.get());
                    const double free_phase = energy*times[t];
                    if (full_m) {
                        const double phase_x = k_transverse*ax[t];
                        const double phase_y = k_transverse*ay[t];
                        for (int p = 0; p < np; ++p)
                            phase[offset+p] = free_phase+cp[p]*phase_x+sp[p]*phase_y;
                    } else {
                        phase[offset] = free_phase+k_axial*az[t];
                    }
                }
            }
            const int nphase = nk*nt*np;
            qprop_sincos(sine.data(), cosine.data(), phase.data(), &nphase);
            for (int k = 0; k < nk; ++k)
                for (int t = 0; t < nt; ++t)
                    for (int p = 0; p < np; ++p) {
                        const int offset = (k*nt+t)*np+p;
                        fft_amplitude[((kb+k)*num_theta_surff+th)*np_out+p] +=
                            win[t] * cplxd(cosine[offset], sine[offset]) * rings[offset];
                    }
        }
    }
}

void eval_fft_spectrum() {
    const bool full_m = qprop_dim == 44;
    const int np_out = full_m ? int(num_phi_surff) : 1;
    const int np = full_m ? np_out - 1 : 1;
    std::ifstream tail("./dat/isurfv.dat");
    vector<cplxd> value(ell_m_grid_size), deriv(ell_m_grid_size);
    double energy = 0, kval = 0, re, im, dre, dim;
    auto read_tail = [&]() {
        tail >> energy >> kval;
        for (int c = 0; c < ell_m_grid_size; ++c) {
            tail >> re >> im >> dre >> dim;
            value[c] = cplxd(re, im); deriv[c] = cplxd(dre, dim);
        }
        if (!tail) throw std::runtime_error("Incomplete isurfv.dat");
    };
    if (tsurff_method == 2)
        for (int k = 0; k < i_proc*(num_k_surff/num_proc); ++k) read_tail();
    ofstream finite("./dat/tsurff_polar_"+fixed_length_number(num_proc,i_proc)+".dat");
    ofstream total("./dat/isurfv_polar_"+fixed_length_number(num_proc,i_proc)+".dat");
    finite.precision(17); total.precision(17);
    const int channels = full_m ? ell_m_grid_size : ell_grid_size+1;
    const int tail_mlo = full_m ? 1-int(ell_grid_size) : int(initial_m);
    const int tail_nm = full_m ? 2*int(ell_grid_size)-1 : 1;
    vector<cplxd> tail_ring(np);
    vector<cplxd> tail_modes(fft_use_direct ? tail_nm : 0);
    vector<cplxd> tail_radial(ell_m_grid_size);
    const double normalization = R_tsurff/sqrt(2*M_PI);
    for (int k = 0; k < num_k_proc; ++k) {
        kval = k_values[k]; energy = fft_energy[k];
        if (tsurff_method == 2) read_tail();
        const double *bessel = fft_bessel.data()+size_t(k)*(ell_grid_size+1);
        cplxd tail_prefactor(0);
        if (tsurff_method == 2) {
            tail_prefactor = normalization*cplxd(0,1)*exp(cplxd(0,energy*duration));
            for (int l = 0; l < ell_grid_size; ++l) {
                const int mlo = full_m ? -l : int(initial_m), mhi = full_m ? l : int(initial_m);
                const double dj = -kval*bessel[l+1]+l/R_tsurff*bessel[l];
                for (int m = mlo; m <= mhi; ++m) {
                    if (abs(m) > l) continue;
                    const int c = full_m ? l*(l+1)+m : l;
                    tail_radial[c] = fft_orbital_phase[l]*
                        (bessel[l]*(deriv[c]-value[c]/R_tsurff)-dj*value[c]);
                }
            }
        }
        for (int th = 0; th < num_theta_surff; ++th) {
            const double theta = thetas_surff[th];
            if (tsurff_method == 2) {
                std::fill(tail_ring.begin(), tail_ring.end(), cplxd(0));
                if (fft_use_direct) std::fill(tail_modes.begin(), tail_modes.end(), cplxd(0));
                for (int l = 0; l < ell_grid_size; ++l) {
                    const int mlo = full_m ? -l : int(initial_m), mhi = full_m ? l : int(initial_m);
                    for (int m = mlo; m <= mhi; ++m) {
                        if (abs(m) > l) continue;
                        const int c = full_m ? l*(l+1)+m : l;
                        const int bin = full_m ? (m % np + np) % np : 0;
                        cplxd &target = fft_use_direct ? tail_modes[m-tail_mlo] : tail_ring[bin];
                        target += fft_angular[size_t(th)*channels+c]*tail_radial[c];
                    }
                }
                if (fft_use_direct) {
                    const cplxd one(1), zero(0);
                    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                                1, np, tail_nm, &one, tail_modes.data(), tail_nm,
                                fft_synthesis.data(), np, &zero, tail_ring.data(), np);
                } else if (np > 1)
                    gsl_fft_complex_backward(reinterpret_cast<double*>(tail_ring.data()),
                                             1, np, fft_table.get(), fft_work.get());
            }
            for (int p = 0; p < np_out; ++p) {
                const int source_p = full_m && p == np_out-1 ? 0 : p;
                const double phi = full_m ? phis_surff[p] : 0;
                cplxd b = normalization*delta_t*fft_amplitude[(k*num_theta_surff+th)*np_out+source_p];
                finite << energy << ' ' << kval << ' ' << theta << ' ' << phi << ' '
                       << norm(b)*kval << ' ' << real(b) << ' ' << imag(b) << '\n';
                if (tsurff_method == 2) {
                    b += tail_prefactor*tail_ring[source_p];
                    total << energy << ' ' << kval << ' ' << theta << ' ' << phi << ' '
                          << norm(b)*kval << ' ' << real(b) << ' ' << imag(b) << '\n';
                }
            }
            if (full_m) { finite << '\n'; total << '\n'; }
        }
        finite << '\n'; total << '\n';
    }
}
