import os
import sys
import numpy as np
import healpy as hp

# --------------------------------------------------
# Configuration
# --------------------------------------------------
cmb64_file = "ilc_null_64_nd.fit"         # full finite SMICA map / null simulation 
mask64_file = "mask_64.fits"          # binary 0/1 mask
power_spectrum_file = "COM_PowerSpect_CMB-TT-full_R3.01.txt"

lmax = 192
n_realizations = 20

work_dir = "cmbpaint_run"
input_dir = os.path.join(work_dir, "input")
matrix_dir = os.path.join(work_dir, "matrix")
output_dir = "inpainted_maps"

In_Name = "smica_for_paint.fits"
I_Mask = "mask_for_paint.fits"

repo_path = "/home/agnes/Skyclean/docs/circle/CMB-PAInT"
cmbpaint_global_path = os.path.join(repo_path, "CMBPAInT") + "/"

os.makedirs(input_dir, exist_ok=True)
os.makedirs(matrix_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)

# --------------------------------------------------
# Import CMB-PAInT
# --------------------------------------------------
sys.path.insert(0, repo_path)
os.environ["PYTHONPATH"] = repo_path + ":" + os.environ.get("PYTHONPATH", "")

from CMBPAInT import CMB_PAInT as CMBP

print("CMB-PAInT path exists:", os.path.exists(cmbpaint_global_path))

# --------------------------------------------------
# Load map and mask directly
# --------------------------------------------------
map_for_paint = hp.read_map(cmb64_file)
mask_for_paint = hp.read_map(mask64_file)

nside = hp.get_nside(map_for_paint)

print("Loaded input map and mask")
print("Nside:", nside)
print("Mask unique values:", np.unique(mask_for_paint))
print("Valid sky fraction:", np.mean(mask_for_paint > 0.5))
print("Observed map std:", np.std(map_for_paint[mask_for_paint > 0.5]))
print("Number of hp.UNSEEN in input map:", np.sum(map_for_paint == hp.UNSEEN))

# --------------------------------------------------
# Load Planck TT spectrum: Dl -> Cl
# --------------------------------------------------
spec = np.loadtxt(power_spectrum_file, comments="#")

ell = spec[:, 0].astype(int)
Dl = spec[:, 1]

cl_fid = np.zeros(lmax + 1)

for l, d_l in zip(ell, Dl):
    if 2 <= l <= lmax:
        cl_fid[l] = 2.0 * np.pi * d_l / (l * (l + 1))

# Convert microK^2 -> K^2
cl_fid *= 1e-12

print("Loaded fiducial Planck TT Cl")
print("Cl shape:", cl_fid.shape)
print("Cl[2:10]:", cl_fid[2:10])

# --------------------------------------------------
# Save CMB-PAInT inputs
# --------------------------------------------------
map_path = os.path.join(input_dir, In_Name)
mask_path = os.path.join(input_dir, I_Mask)
cl_path = os.path.join(input_dir, "cl_tt.npy")

hp.write_map(map_path, map_for_paint, overwrite=True, dtype=np.float64)
hp.write_map(mask_path, mask_for_paint, overwrite=True, dtype=np.float64)
np.save(cl_path, np.array([cl_fid]))

print("Saved CMB-PAInT input files")
print("Map:", map_path)
print("Mask:", mask_path)
print("Cl:", cl_path)

# --------------------------------------------------
# Configure CMB-PAInT
# --------------------------------------------------
Inpainting_Data = CMBP(
    cl_path,             # Cls
    nside,               # Nside
    lmax,                # Lmax
    matrix_dir + "/",    # Out_PATH
    input_dir + "/",     # Inp_In_PATH
    output_dir + "/",    # Inp_Out_PATH
    In_Name,             # input map file name
    n_realizations,      # Num_Sims
    True,                # Single
    False,               # No_Z
    False,               # Cons_Uncons
    1,                   # Fields: temperature only
    None,                # Zbar_PATH
    mask_path,           # I_Mask
    1536,                  # Chunks
    None,                # job_name
    None,                # env
    local=True,
    P_Mask=None,
    DP=False,
    Noise_Level=0.011
)

# --------------------------------------------------
# Run CMB-PAInT
# --------------------------------------------------
Inpainting_Data.generate_config_file()

print("Computing covariance matrix...")
Inpainting_Data.calculate_covariance_matrix(global_path=cmbpaint_global_path)

print("Computing Cholesky decomposition...")
Inpainting_Data.calculate_chol(global_path=cmbpaint_global_path)

print("Generating inpainted maps...")
Inpainting_Data.calculate_inp(global_path=cmbpaint_global_path)

print("Done.")
print("Inpainted maps written to:", output_dir)