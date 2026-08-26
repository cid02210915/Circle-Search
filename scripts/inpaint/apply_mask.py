import healpy as hp
import numpy as np

cmb_file = "/home/agnes/Skyclean/docs/circle/ilc_lmax1023.fits"
mask_file = "COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits"

output_file = (
    "/home/agnes/Skyclean/docs/circle/"
    "ilc_lmax800_nside512_masked.fits"
)

nside_out = 512

# Read already processed ILC map: RING, NSIDE=512, lmax=800
cmb = hp.read_map(cmb_file, field=0, nest=False)

# Read high-resolution common mask
mask = hp.read_map(mask_file, field=0, nest=False)

# Only downgrade the mask
mask = hp.ud_grade(
    mask,
    nside_out=nside_out,
    order_in="RING",
    order_out="RING",
    power=0,
)

print("CMB nside:", hp.get_nside(cmb))
print("Mask nside:", hp.get_nside(mask))
print("CMB npix:", len(cmb))
print("Mask npix:", len(mask))
print("Mask min/max:", np.min(mask), np.max(mask))

if hp.get_nside(cmb) != nside_out:
    raise ValueError("Input ILC map is not NSIDE=512.")

# Apply mask
valid = mask > 0.5

masked = cmb.copy()
masked[~valid] = hp.UNSEEN

# Visualize
hp.mollview(mask, title="Common intensity mask NSIDE=512")
hp.mollview(masked, title="Masked ILC map, NSIDE=512")

print("Total number of pixels:", len(mask))
print("Number of masked pixels:", np.sum(~valid))
print("Valid sky fraction:", np.mean(valid))

# Save masked RING map
hp.write_map(
    output_file,
    masked,
    overwrite=True,
)

print("Saved:", output_file)