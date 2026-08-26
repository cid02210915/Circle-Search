# Circle-Search

`Circle-Search` is a research pipeline for topology-independent searches for matched circles in cosmic microwave background (CMB) temperature maps. It contains the original Fourier-space $S$-statistic, two additional real-space statistics, a MASTER-like mask-correction method, map-inpainting utilities, and the notebooks used for validation and null tests.

## Scientific motivation

A Galactic mask removes part of each circle. Iterative inpainting fills the missing region, but the reconstructed temperatures depend on the random inpainting realisation. Reliable results therefore require many inpainted realisations for every CMB map, which makes the full search computationally expensive.

This repository investigates four connected approaches:

| Method | Main purpose | Geometry search | Dependence on inpainting |
|---|---|---|---|
| Original $S$ | Baseline matched-circle detection | Maximised over $(i,j,\beta)$ | Uses the inpainted map directly |
| $D_{ij}(\beta_S)$ | Test a candidate selected by $S$ using filtered temperatures and jointly valid samples | Uses $(i_S,j_S,\beta_S)$ from $S$ | Inherits the $S$-selected geometry and the filtered input map |
| $D_{ij}^{\Delta}$ | Improve phase localisation and small-circle detectability using neighbouring temperature differences | Independently maximised over $(i,j,\beta)$ | Only jointly valid differences enter the statistic |
| MASTER-$S$ | Correct mask-induced mode coupling before calculating $S$ | Independently maximised over $(i,j,\beta)$ | The unfiltered branch uses the masked map directly and avoids stochastic inpainting |

The $D^{\Delta}$ and MASTER-$S$ developments address different problems: $D^{\Delta}$ sharpens the matched-circle response, while unfiltered MASTER-$S$ provides an alternative to repeating a stochastic and expensive inpainting procedure.

## Search methods

### 1. Original $S$-statistic

For two circles of opening angle $\alpha$, the code evaluates

```math
S_{ij}(\alpha,\beta)
=
\frac{
2\sum_m m\,\operatorname{Re}\!\left[C_{ij}(m)e^{-im\beta}\right]
}{
\sum_m m\left[C_{ii}(m)+C_{jj}(m)\right]
},
```

and maximises over the allowed circle centres $(i,j)$ and relative phase $\beta$. The lowest ring modes can be removed using `m_min`.

### 2. Filtered $D$ at the $S$-selected geometry

The filtered statistic is calculated only at the geometry selected by $S$:

```math
(i_S,j_S,\beta_S)
\longrightarrow
D_{ij}(\beta_S).
```

Using only jointly valid samples,

```math
D_{ij}(\beta_S)
=
\frac{
2\sum_{k\in\mathrm{joint}}T_i(k)T_j(k_{\beta_S})
}{
\sum_{k\in\mathrm{joint}}
\left[T_i^2(k)+T_j^2(k_{\beta_S})\right]
}.
```

This is a candidate-validation statistic, not an independent $D$-geometry search.

### 3. Independent temperature-difference search $D^{\Delta}$

Neighbouring temperature differences are formed at separations $s=1,2$:

```math
\Delta_sT(k)
=
T(k+s)-T(k).
```

Taking differences gives greater weight to small-scale changes along a circle and produces a sharper response around the matched phase. The code independently maximises $D_{ij}^{\Delta}$ over $(i,j,\beta)$ at every $\alpha$.

A difference is used only when both of its endpoints are valid. The output therefore includes $N_{\mathrm{joint}}^{\Delta}$, and a minimum support can be imposed with `d_min_joint_samples`.

Optional tools provide:

- a bias correction depending on $\alpha$ and $N_{\mathrm{joint}}^{\Delta}$;
- null histograms indexed by $\alpha$, $N_{\mathrm{joint}}^{\Delta}$, and $D^{\Delta}$.

### 4. Mask-reconstructed MASTER-$S$

`MASTER_circles.F90` works directly with a temperature map and its mask. For every circle it transforms

```math
T,
\qquad
wT,
\qquad
w,
```

where $w$ is the binary circle window. It then constructs the masked pseudo-mode products

```math
\widehat C_{ij},
\qquad
\widehat C_{ii},
\qquad
\widehat C_{jj},
```

and the mask-coupling eigenvalues

```math
\mu_q
\equiv
\mathcal{F}\!\left[M(:,0)\right]_q.
```

Each product is reconstructed using Tikhonov regularisation:

```math
\widetilde C_q
=
\frac{
\mu_q^*
}{
|\mu_q|^2+\lambda^2
}
\widehat C_q,
\qquad
\lambda
=
f_{\mathrm{reg}}\max_q|\mu_q|.
```

Ordinary generalized cross-validation (GCV) selects $f_{\mathrm{reg}}$. The initial grid contains 21 points from $10^{-9}$ to $10^{-7}$ and is extended when the minimum lies on a boundary.

For efficiency, the $ii$ and $jj$ reconstructions and denominator contributions are cached once per circle. Only the geometry-dependent $ij$ reconstruction is repeated for each candidate pair.

The corrected products are inserted into $S(\beta)$. Phase bins satisfying

```math
f_{\mathrm{lost}}(\beta)>0.5
```

are excluded before maximisation.

Filtered and unfiltered MASTER-$S$ searches can be run separately using `use_sw_filter`. Because the harmonic filter is applied before the mask, the filtered branch requires a defined full-sky input map, whereas the unfiltered branch can use the mask directly without stochastic filling.

## Repository structure

```text
Circle-Search/
├── README.md
├── src/
│   ├── circles.F90
│   ├── MASTER_circles.F90
│   ├── Makefile.am
│   └── Makefile.in
├── scripts/
│   └── inpaint/
│       ├── down_to_512.py
│       ├── apply_mask.py
│       ├── iterative_inpaint.py
│       └── conditional_inpaint.py
├── notebooks/
│   ├── E1_filter.ipynb
│   ├── MASTER_test_1.ipynb
│   ├── MASTER_test_2.ipynb
│   ├── results.ipynb
│   ├── S_null_thresholds_and_injection.ipynb
│   ├── f_lost.ipynb
│   └── test_CMB_inpaint.ipynb
├── config/
│   ├── s_plus_d_at_s_geometry.ini
│   └── d_diff.ini
└── data/
    ├── filters/
    ├── masks/
    └── spectra/
```

## Source files

### `src/circles.F90`

Contains the original $S$ search and the added direct-search functions:

- configurable `m_min`;
- return of the best circle centres and $\beta$;
- filtered $D_{ij}(\beta_S)$;
- independent $D_{ij}^{\Delta}$;
- joint-validity counting and thresholds;
- optional $D^{\Delta}$ bias correction and null histograms;
- independent switches for $S$, $D$, and $D^{\Delta}$;
- symmetry and phase-convention diagnostics.

### `src/MASTER_circles.F90`

Contains the complete mask-reconstruction pipeline:

- Fourier transforms of $T$, $wT$, and $w$;
- pseudo-mode products and circulant coupling matrices;
- GCV-selected Tikhonov inversion;
- cached auto reconstructions for $ii$ and $jj$;
- pair-dependent reconstruction of $ij$;
- phase-dependent $f_{\mathrm{lost}}$ rejection;
- optional filtered and unfiltered MASTER-$S$ searches.

## Inpainting scripts

| Script | Purpose |
|---|---|
| `down_to_512.py` | Downgrade input maps to the working resolution |
| `apply_mask.py` | Apply the chosen Galactic mask |
| `iterative_inpaint.py` | Generate iterative harmonic inpainting realisations |
| `conditional_inpaint.py` | Generate conditional inpainting realisations |

Because the inpainted region is stochastic, the inpainting scripts should be run with multiple seeds for each underlying CMB realisation.

## Notebooks

| Notebook | Main role |
|---|---|
| `E1_filter.ipynb` | Test the filtered matched-circle response |
| `MASTER_test_1.ipynb` | Construct and validate the MASTER-like inverse problem |
| `MASTER_test_2.ipynb` | Test regularisation and reconstructed $S$ |
| `results.ipynb` | Compare search outputs and detection performance |
| `S_null_thresholds_and_injection.ipynb` | Calculate null critical levels and test injected circles |
| `f_lost.ipynb` | Study signal recovery as a function of the lost-circle fraction |
| `test_CMB_inpaint.ipynb` | Validate the CMB inpainting pipeline and its realisation dependence |

## Requirements

### Fortran search code

- GNU Fortran or another modern Fortran compiler
- MPI
- OpenMP
- HEALPix Fortran libraries
- FFTW3
- CFITSIO
- Autotools when regenerating the build system

The repository tracks `Makefile.am` and `Makefile.in`. A standalone build also requires the top-level Autotools configuration from the parent `circles-in-the-sky` project. After configuring the compiler and dependency paths for the target system, build the executables with `make`.

### Python analysis

The scripts and notebooks primarily require:

```text
numpy
pandas
matplotlib
healpy
```

## Running the searches

### Original $S$ with $D$ at the selected geometry

```bash
mpirun -np <MPI_RANKS> ./circles config/s_plus_d_at_s_geometry.ini
```

### Independent $D^{\Delta}$ search

```bash
mpirun -np <MPI_RANKS> ./circles config/d_diff.ini
```

### MASTER-$S$ search

```bash
mpirun -np <MPI_RANKS> ./MASTER_circles <master_configuration.ini>
```

The MASTER configuration specifies the map, mask, search geometry, and optional harmonic filter. Exact executable paths and MPI launch commands may differ between systems or cluster schedulers.

## Important configuration options

### Common geometry settings

| Option | Purpose |
|---|---|
| `opening_angle_min`, `opening_angle_max` | Range of opening angles $\alpha$ |
| `n_opening_angles` | Number of opening-angle grid points |
| `nphi_ring` | Number of samples around each circle |
| `nside_search` | HEALPix resolution of the circle-centre grid |
| `min_separation_angle` | Minimum separation between circle centres |
| `orientable` | Select the circle-orientation convention |
| `m_min` | Remove the lowest ring Fourier modes |
| `mapfilename` | Input CMB map |

### Original $S$ and $D(\beta_S)$

| Option | Purpose |
|---|---|
| `use_s_statistics` | Enable the original $S$-statistic search |
| `use_d_statistics` | Evaluate filtered $D$ at the $S$-selected geometry |
| `d_filter_filename` | Harmonic filter applied before calculating $D$ |
| `d_filter_lmax` | Maximum multipole of the filter |
| `maskfilename` | Mask used to identify jointly valid samples |
| `mask_threshold` | Threshold used to binarise the mask |
| `d_min_joint_samples` | Minimum number of jointly valid samples |
| `write_position_smax` | Write the best result for every source centre |

`use_d_statistics = T` requires `use_s_statistics = T` because $D$ uses the geometry selected by $S$.

### Independent $D^{\Delta}$

| Option | Purpose |
|---|---|
| `use_d_diff_statistics` | Enable the independent temperature-difference search |
| `d_filter_filename` | Harmonic filter applied before taking differences |
| `d_filter_lmax` | Maximum multipole of the filter |
| `maskfilename` | Mask used to validate both endpoints of each difference |
| `mask_threshold` | Threshold used to binarise the mask |
| `d_min_joint_samples` | Minimum $N_{\mathrm{joint}}^{\Delta}$ |
| `use_d_diff_bias_correction` | Subtract the calibrated null bias |
| `write_d_diff_null_histogram` | Write conditional null-histogram counts |
| `d_diff_hist_nbins` | Number of $D^{\Delta}$ histogram bins |

The $S$ and $D^{\Delta}$ searches can be enabled independently.

### MASTER-$S$

| Option | Purpose |
|---|---|
| `maskfilename` | Mask used to construct the circle windows and coupling matrices |
| `mask_threshold` | Threshold used to binarise each circle window |
| `use_sw_filter` | Select filtered or unfiltered MASTER-$S$ input |
| `sw_filter_filename` | Harmonic filter used when `use_sw_filter = T` |
| `sw_filter_lmax` | Maximum multipole of the MASTER filter |
| `write_position_smax` | Write the best MASTER-$S$ result for every source centre |

The GCV grid, adaptive extensions, Tikhonov inversion, and $f_{\mathrm{lost}}\leq0.5$ requirement are implemented internally in `MASTER_circles.F90`.

## Outputs

### Original $S$

The main `.dat` table contains:

```text
alpha_deg, Smax,
pixel_i, theta_i_deg, phi_i_deg,
pixel_j, theta_j_deg, phi_j_deg,
separation_deg, beta_max_bin, beta_max_deg
```

When `write_position_smax = T`, the position-specific files contain the best $S$, $\beta_S$, and—when enabled—$D_{ij}(\beta_S)$ and `n_joint` for every source centre.

### Filtered $D(\beta_S)$

The `_d.dat` table contains:

```text
alpha_deg, D_at_betaS,
pixel_i, theta_i_deg, phi_i_deg,
pixel_j, theta_j_deg, phi_j_deg,
separation_deg, beta_S_bin, beta_S_deg, n_joint
```

The circle pair and phase in this file are inherited from the maximum-$S$ geometry.

### Independent $D^{\Delta}$

The `_d_diff.dat` table contains:

```text
alpha_deg, D_diff_max,
pixel_i, theta_i_deg, phi_i_deg,
pixel_j, theta_j_deg, phi_j_deg,
separation_deg, beta_max_bin, beta_max_deg, n_joint_diff
```

### MASTER-$S$

The final MASTER table contains the reconstructed maximum statistic and geometry:

```text
alpha_deg, MASTER_Smax,
pixel_i, theta_i_deg, phi_i_deg,
pixel_j, theta_j_deg, phi_j_deg,
separation_deg, beta_max_bin, beta_max_deg
```

Position-specific files contain the best reconstructed $S$ and $\beta$ for each source centre. Per-MPI diagnostic files additionally record the selected regularisation fractions:

```text
f_reg_ij, f_reg_ii, f_reg_jj
```

## Included data

- `SW_Doppler_ISW_filter_lmax300.dat`: harmonic filter used by $D$, $D^{\Delta}$, and optional filtered MASTER-$S$ runs.
- `COM_Mask_CMB-common-Mask-Int_128_R3.00.fits`: Planck common CMB mask at $N_{\mathrm{side}}=128$.
- `COM_PowerSpect_CMB-TT-full_R3.01.txt`: Planck temperature power spectrum used for simulations and inpainting.

## Run records

For every run, record the common map and search settings:

- CMB realisation and random seed;
- $N_{\mathrm{side}}$, $\ell_{\max}$, and $N_\phi$;
- opening-angle and centre-separation ranges;
- mask file and threshold;
- filter file and $\ell_{\max}$;
- MPI ranks and OpenMP threads.

For the inpainted $S$ pipeline, also record the number of inpainting realisations, every inpainting seed, and the number of iterations.

For $D^{\Delta}$, record the joint-sample threshold and whether bias correction was enabled.

For MASTER-$S$, record whether filtering was enabled, the selected $f_{\mathrm{reg}}$ values, and the applied $f_{\mathrm{lost}}$ limit.

All detection claims should be calibrated against null simulations processed through the same pipeline as the injected or observed maps.
