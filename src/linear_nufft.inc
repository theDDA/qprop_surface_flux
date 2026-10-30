bool use_linear_nufft() const {
#ifdef HAVE_FINUFFT
    return expansion_scheme == 1 && qprop_dim == 34 && initial_m == 0 && !fft_use_direct;
#else
    return false;
#endif
}

void prep_nufft_space() {
    const int L = int(ell_grid_size), nc = L+1;
    fft_amplitude.assign(size_t(num_k_proc)*num_theta_surff, cplxd(0));
    fft_energy.resize(num_k_proc);
    fft_bessel.resize(size_t(num_k_proc)*nc);
    fft_angular.resize(size_t(num_theta_surff)*nc);
    fft_orbital_phase.resize(nc);
    for (int l=0; l<nc; ++l) fft_orbital_phase[l]=pow_neg_i(l+1);
    for (int k=0; k<num_k_proc; ++k) {
        fft_energy[k]=0.5*k_values[k]*k_values[k];
        fill_bessel_array(k_values[k], R_tsurff, fft_bessel.data()+size_t(k)*nc);
    }
    for (int th=0; th<num_theta_surff; ++th)
        for (int l=0; l<nc; ++l)
            fft_angular[size_t(th)*nc+l]=gsl_sf_legendre_sphPlm(l,0,cos(thetas_surff[th]));
}

void time_integration_nufft() {
#ifdef HAVE_FINUFFT
    const auto start=std::chrono::steady_clock::now();
    const int L=int(ell_grid_size);
    const size_t nq=size_t(num_k_proc)*num_theta_surff;
    psi.seekg(0,std::ios::end);
    const std::streamoff bytes=psi.tellg();
    dpsi_dr.seekg(0,std::ios::end);
    if (bytes<=0 || bytes%(L*sizeof(cplxd)) || dpsi_dr.tellg()!=bytes)
        throw std::runtime_error("A1: empty or mismatched surface histories");
    const size_t nt=bytes/(L*sizeof(cplxd));
    vector<cplxd> u(nt*L), d(nt*L);
    psi.seekg(0); dpsi_dr.seekg(0);
    psi.read(reinterpret_cast<char*>(u.data()),bytes);
    dpsi_dr.read(reinterpret_cast<char*>(d.data()),bytes);
    if (!psi || !dpsi_dr) throw std::runtime_error("A1: incomplete history read");
    vector<double> times(nt), alpha(nt), field(nt), weights(nt), energy(nq), kz(nq);
    double amin=0,amax=0;
    for (size_t t=0; t<nt; ++t) {
        times[t]=delta_t*t;
        if (vecpot_x(times[t],0)!=0 || vecpot_y(times[t],0)!=0 ||
            vecpot_x.integral(times[t])!=0 || vecpot_y.integral(times[t])!=0)
            throw std::runtime_error("A1 requires a z-polarized field");
        alpha[t]=vecpot_z.integral(times[t]);
        field[t]=vecpot_z(times[t],0);
        weights[t]=hanning_win(times[t]);
        if (!t) amin=amax=alpha[t];
        amin=std::min(amin,alpha[t]); amax=std::max(amax,alpha[t]);
    }
    const double center=(amin+amax)/2;
    for (size_t t=0; t<nt; ++t) alpha[t]-=center;
    for (int k=0; k<num_k_proc; ++k)
        for (int th=0; th<num_theta_surff; ++th) {
            const size_t q=size_t(k)*num_theta_surff+th;
            energy[q]=fft_energy[k]; kz[q]=k_values[k]*cos(thetas_surff[th]);
        }
    struct Plan {
        finufft_plan p=nullptr;
        ~Plan() { if (p) finufft_destroy(p); }
    } plan;
    auto check=[](int status) {
        if (status) throw std::runtime_error("A1: FINUFFT status "+std::to_string(status));
    };
    finufft_opts opts; finufft_default_opts(&opts);
    opts.nthreads=1; opts.spread_thread=2;
    check(finufft_makeplan(3,2,nullptr,1,2,1e-12,&plan.p,&opts));
    check(finufft_setpts(plan.p,nt,times.data(),alpha.data(),nullptr,nq,
                        energy.data(),kz.data(),nullptr));
    vector<cplxd> coeff(2*nt), transformed(2*nq);
    double transform_seconds=0, sum_seconds=0;
    for (int l=0; l<=L; ++l) {
        for (size_t t=0; t<nt; ++t) {
            const cplxd *p=u.data()+t*L;
            cplxd h=l<L ? d[t*L+l]-p[l]/R_tsurff : cplxd(0);
            if (l>0) h+=cplxd(0,2*field[t]*l/sqrt((2.*l-1)*(2.*l+1)))*p[l-1];
            if (l+1<L) h+=cplxd(0,2*field[t]*(l+1)/sqrt((2.*l+1)*(2.*l+3)))*p[l+1];
            // The output routine applies delta_t and R/sqrt(2*pi).
            coeff[t]=weights[t]*h;
            coeff[nt+t]=l<L ? weights[t]*p[l] : cplxd(0);
        }
        auto tick=std::chrono::steady_clock::now();
        check(finufft_execute(plan.p,coeff.data(),transformed.data()));
        transform_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-tick).count();
        tick=std::chrono::steady_clock::now();
        for (int k=0; k<num_k_proc; ++k) {
            const double *j=fft_bessel.data()+size_t(k)*(L+1);
            const double minus_dj=l<L ? k_values[k]*j[l+1]-l/R_tsurff*j[l] : 0;
            for (int th=0; th<num_theta_surff; ++th) {
                const size_t q=size_t(k)*num_theta_surff+th;
                fft_amplitude[q]+=fft_orbital_phase[l]*fft_angular[size_t(th)*(L+1)+l]*
                    (j[l]*transformed[q]+minus_dj*transformed[nq+q]);
            }
        }
        sum_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-tick).count();
    }
    for (size_t q=0; q<nq; ++q) fft_amplitude[q]*=exp(cplxd(0,kz[q]*center));
    cout.precision(17);
    cout << "A1_NUFFT_WALL=" << transform_seconds << " A1_SUM_WALL=" << sum_seconds
         << " A1_TOTAL_WALL=" << std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()
         << " A1_SAMPLES=" << nt << " A1_CHANNELS=" << L+1 << endl;
#else
    throw std::runtime_error("A1 requires a build with FINUFFT");
#endif
}
