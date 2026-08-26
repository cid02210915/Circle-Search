import healpy as hp
import numpy as np

# -------------------------
# Input files
# -------------------------
cmb_file = "COM_CMB_IQU-smica_2048_R3.00_full.fits"
mask_file = "COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits"

# -------------------------
# Output files
# -------------------------
cmb64_file = "smica_64.fits"
mask64_file = "mask_64.fits"
masked64_file = "smica_masked_64.fits"

# -------------------------

# Read input map and mask
# -------------------------
cmb2048 = hp.read_map(cmb_file, field=0)
mask2048 = hp.read_map(mask_file, field=0)

# -------------------------
# Downgrade CMB map
# -------------------------
# Smooth before downgrading to NSIDE=64.
# NSIDE=64 has pixel size ~1.8 deg, so use stronger smoothing than 5 arcmin.
fwhm_rad = np.radians(2.0)   # 2 degrees
cmb2048_smooth = hp.smoothing(cmb2048, fwhm=fwhm_rad)

cmb64 = hp.ud_grade(cmb2048_smooth, nside_out=64)

# -------------------------
# Downgrade mask
# -------------------------
mask64_float = hp.ud_grade(mask2048, nside_out=64)

# Conservative binary mask
mask64 = (mask64_float > 0.9).astype(np.float64)

# -------------------------
# Apply mask
# -------------------------
masked64 = cmb64.copy()
masked64[mask64 < 0.5] = hp.UNSEEN

# -------------------------
# Save outputs
# -------------------------
hp.write_map(cmb64_file, cmb64, overwrite=True)
hp.write_map(mask64_file, mask64, overwrite=True)
hp.write_map(masked64_file, masked64, overwrite=True)

# -------------------------
# Quick checks
# -------------------------
print("Original nside:", hp.get_nside(cmb2048))
print("Downgraded nside:", hp.get_nside(cmb64))
print("Mask sky fraction at 64:", np.mean(mask64 > 0.5))
print("Masked pixels at 64:", np.sum(masked64 == hp.UNSEEN))
print("Saved:", cmb64_file, mask64_file, masked64_file)