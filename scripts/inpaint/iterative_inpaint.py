import os
import healpy as hp
import numpy as np

# --------------------------------------------------
# Configuration
# --------------------------------------------------
full_map_file = (
    "/home/agnes/Skyclean/docs/circle/"
    "ilc_lmax800_nside512_masked.fits"
)

power_spectrum_file = "COM_PowerSpect_CMB-TT-full_R3.01.txt"

output_dir = (
    "/home/agnes/Skyclean/docs/circle/"
    "iterative_inpainted_ilc_lmax800_nside512"
)

nside = 512
lmax = 800
n_iter = 50
n_realizations = 1
seed_base = 12345

os.makedirs(output_dir, exist_ok=True)

# --------------------------------------------------
# Load already-masked ILC map
# --------------------------------------------------
full_map = hp.read_map(full_map_file, field=0, verbose=False)

assert hp.get_nside(full_map) == nside, "Input map is not NSIDE=512"

# The input map is already masked:
# observed pixels have data; masked pixels are hp.UNSEEN.
valid_mask = full_map != hp.UNSEEN

# Use zero only internally in the unknown region.
observed_map = full_map.copy()
observed_map[~valid_mask] = 0.0

print("Loaded already-masked ILC map")
print("NSIDE:", hp.get_nside(full_map))
print("Valid sky fraction:", np.mean(valid_mask))
print("Observed mean:", np.mean(full_map[valid_mask]))
print("Observed std:", np.std(full_map[valid_mask]))

# -------------------------------------------------------------
# Load Planck theory TT spectrum: D_l [microK^2] -> C_l [K^2]
# -------------------------------------------------------------
data = np.loadtxt(power_spectrum_file)

ell_data = data[:, 0].astype(int)
Dl_microK2 = data[:, 1]

cl_theory = np.zeros(lmax + 1)

for ell, Dl in zip(ell_data, Dl_microK2):
    if 2 <= ell <= lmax:
        Dl_K2 = Dl * 1e-12
        cl_theory[ell] = Dl_K2 * 2.0 * np.pi / (ell * (ell + 1.0))

# No monopole/dipole in random prior maps
cl_theory[0] = 0.0
cl_theory[1] = 0.0

# --------------------------------------------------
# Inpainting function
# --------------------------------------------------
def generate_inpainted_realization(
    observed_map,
    valid_mask,
    cl,
    nside,
    lmax,
    n_iter,
    seed,
):
    np.random.seed(seed)

    alm_rand = hp.synalm(cl, lmax=lmax, new=True)
    work_map = hp.alm2map(
        alm_rand,
        nside=nside,
        lmax=lmax,
        verbose=False,
    )

    # Keep observed ILC sky fixed.
    work_map[valid_mask] = observed_map[valid_mask]

    for _ in range(n_iter):
        alm = hp.map2alm(work_map, lmax=lmax, iter=3)

        work_map = hp.alm2map(
            alm,
            nside=nside,
            lmax=lmax,
            verbose=False,
        )

        # Re-impose observed ILC sky.
        work_map[valid_mask] = observed_map[valid_mask]

    return work_map

# --------------------------------------------------
# Generate realizations
# --------------------------------------------------
for i in range(n_realizations):
    seed = seed_base + i

    inpainted = generate_inpainted_realization(
        observed_map=observed_map,
        valid_mask=valid_mask,
        cl=cl_theory,
        nside=nside,
        lmax=lmax,
        n_iter=n_iter,
        seed=seed,
    )

    outname = os.path.join(
        output_dir,
        f"iterative_inpainted_ilc_{i + 1:03d}.fits",
    )

    # Default output ordering is RING.
    hp.write_map(outname, inpainted, overwrite=True)

    print(f"[{i + 1:03d}/{n_realizations}] wrote {outname}")
    print("Full output mean:", np.mean(inpainted))
    print("Observed output std:", np.std(inpainted[valid_mask]))
    print("Filled-region std:", np.std(inpainted[~valid_mask]))
    print("Full output std:", np.std(inpainted))

print("Done.")