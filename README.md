# Qprop surface-flux extraction

Code accompanying *Reducing repeated work in surface-flux photoelectron-spectrum
extraction*, by Damir Agačević. The repository contains four added source files
and an integration patch for Qprop 3.2.

| Algorithm | Implementation |
| --- | --- |
| A1: temporal type-3 NUFFT, linear polarization with $m=0$ | `src/linear_nufft.inc` |
| B1: angular contraction with direct azimuthal synthesis | `src/fft_scheme1.inc` |
| B2: angular contraction with azimuthal FFT | `src/fft_scheme1.inc` |
| C1: exact elimination of the source-free radial exterior | `src/surface_resolvent.hh` |
| C2: direct radial solve for conserved $m$ | `src/surface_resolvent.hh` |
| C3: shared surface responses for the full magnetic basis | `src/surface_resolvent.hh` |

## Build

Obtain [Qprop 3.2](https://www.qtmps.physik.uni-rostock.de/forschung/qprop/qprop/)
separately. Requirements: a C++11 compiler, Make, `patch`, GSL, FINUFFT 2.5.1,
and LP64 BLAS/LAPACK with LAPACKE. On macOS, the build uses Accelerate for
BLAS/LAPACK.

```sh
make \
  -j 4 \
  QPROP_SRC=/path/to/qprop_v3.2/src \
  FINUFFT_CFLAGS=-I/path/to/finufft/include \
  FINUFFT_LIBS="-L/path/to/finufft/lib -Wl,-rpath,/path/to/finufft/lib -lfinufft" \
  all \
  b1
```

The patch is applied to copies under `build/overlay`; the supplied Qprop tree is
unchanged. Executables are written to `build/bin/{attoclock,long-wavelength,vortex,lissajous}/`:

- `tsurff`: A1 for dimension 34 with $m=0$, B1 for other conserved-$m$
  cases, and B2 for dimension 44.
- `tsurff_b1`: B1, for the direct angular comparisons in the paper.
- `isurfv`: C2 for dimension 34, C3 for dimension 44, with C1 wherever a
  source-free exterior is present.

The modified t-SURFF algorithms use `expansion-scheme = 1`. Dimension 34 uses a
field along $z$; dimension 44 uses a field in the $xy$ plane. These are serial
executables. A1 uses one FINUFFT thread and tolerance $10^{-12}$.

For a build without A1, set `FINUFFT_LIBS=`; `tsurff` then uses B1/B2. Use a separate
`BUILD` directory when changing compiler flags or dependencies.

## Run

Run from the Qprop calculation directory containing `initial.param`,
`propagate.param`, `tsurff.param`, and the saved propagation data in `dat/`.
Use the executable corresponding to that example's field.

```sh
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
/path/to/qprop_surface_flux/build/bin/attoclock/isurfv
/path/to/qprop_surface_flux/build/bin/attoclock/tsurff
```

The programs write Qprop's native `dat/isurfv.dat`, `dat/tsurff_polar_0.dat`, and
`dat/isurfv_polar_0.dat` files. Run `isurfv` first when using `tsurff-version = 2`.
Set `QPROP_SURFACE_FULL_GRID=1` to disable C1 for the full-grid comparison.

Qprop remains subject to its [upstream terms of use](https://www.qtmps.physik.uni-rostock.de/forschung/qprop/qprop/).
Please cite the relevant Qprop papers, including
[Tulsky and Bauer, *Comput. Phys. Commun.* **251**, 107098 (2020)](https://doi.org/10.1016/j.cpc.2019.107098).
