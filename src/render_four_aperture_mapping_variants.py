from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio


root = Path(__file__).resolve().parents[1]
data_path = root / "ImageForShow" / "four_aperture_preview_data.mat"
out_path = root / "ImageForShow" / "four_aperture_mapping_variants.png"

mat = sio.loadmat(data_path)
intensity = mat["I_norm"]
x_mm = mat["X"][0, :] * 1e3
y_mm = mat["Y"][:, 0] * 1e3

variants = [
    (99.95, 18, 0.90),
    (99.90, 35, 0.80),
    (99.85, 60, 0.72),
    (99.80, 100, 0.65),
]

fig, axes = plt.subplots(2, 2, figsize=(10, 10), constrained_layout=True)
last_image = None

for ax, (percentile, alpha, gamma) in zip(axes.ravel(), variants):
    clip = np.percentile(intensity, percentile)
    visible = np.minimum(intensity, clip) / clip
    visible = np.log1p(alpha * visible) / np.log1p(alpha)
    visible = visible**gamma

    last_image = ax.imshow(
        visible,
        extent=[x_mm.min(), x_mm.max(), y_mm.min(), y_mm.max()],
        origin="lower",
        cmap="gray",
        vmin=0,
        vmax=1,
    )
    ax.set_aspect("equal")
    ax.set_title(f"p{percentile} alpha={alpha} gamma={gamma}")
    ax.set_xlabel("x / mm")
    ax.set_ylabel("y / mm")

fig.colorbar(last_image, ax=axes.ravel().tolist(), shrink=0.8)
fig.savefig(out_path, dpi=180)
print(out_path)
