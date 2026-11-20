from pathlib import Path
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
from scipy.io import loadmat

module = Path(__file__).resolve().parents[2]
folder = module / "Output" / "case01" / "array_pupil_comparison"
trials = [loadmat(folder / f"trial{k}.mat", simplify_cells=True) for k in (1, 2, 3)]
labels = ("4×4 · 入瞳 0.20 mm", "8×8 · 入瞳 0.20 mm", "8×8 · 入瞳 0.40 mm")
colors = ("#5E6B82", "#2475A5", "#C5652D")
channel_colors = ("#B53B45", "#33814F", "#3464B1")
coordinate = trials[0]["x"] * 1e6
frequency = trials[0]["frequency"]
center = np.abs(coordinate).argmin()
freq_positive = frequency[center:]

for trial in trials[1:]:
    np.testing.assert_array_equal(trial["frequency"], frequency)
    np.testing.assert_array_equal(trial["x"] * 1e6, coordinate)
for key in ("psf", "mtf", "fixedExposure", "captured"):
    difference = np.max(np.abs(trials[0][key] - trials[1][key]))
    print(f"{key}: maximum difference with fixed pupil = {difference:.4g}")
    np.testing.assert_allclose(trials[0][key], trials[1][key], rtol=0, atol=1e-10)

def crossing(x, y, level):
    j = np.flatnonzero(y <= level)[0]
    return np.interp(level, y[j-1:j+1][::-1], x[j-1:j+1][::-1])

captured = np.stack([trial["captured"] for trial in trials])
area = np.array([trial["incidentArea"] for trial in trials])
relative_signal = captured * area[:, None] / (captured[0] * area[0])
mtf_green = np.stack([trial["mtf"][center, center:, 1] for trial in trials])
mtf50 = [crossing(freq_positive, mtf, 0.5) for mtf in mtf_green]
green_psf = [trial["psf"][:, :, 1] for trial in trials]
fwhm = [
    2 * crossing(coordinate[center:], psf[center, center:], psf[center, center] * 0.5)
    for psf in green_psf
]
mask_fraction = [trial["openSamples"] / trial["pupilSamples"] for trial in trials]
with (folder / "summary.csv").open("w", newline="", encoding="utf-8") as stream:
    writer = csv.writer(stream)
    writer.writerow(["PitchesPerAxis", "PupilDiameterMm", "PupilAreaMm2", "OpenFraction",
                     "GreenPSFFWHMUm", "GreenMTF50CyclesPerMm", "CapturedR", "CapturedG", "CapturedB",
                     "RelativeSignalR", "RelativeSignalG", "RelativeSignalB"])
    for i, (periods, diameter) in enumerate(((4, 0.2), (8, 0.2), (8, 0.4))):
        writer.writerow([periods, diameter, area[i] * 1e6, mask_fraction[i], fwhm[i], mtf50[i],
                         *captured[i], *relative_signal[i]])

font_path = "C:/Windows/Fonts/msyh.ttc"
font_manager.fontManager.addfont(font_path)
plt.rcParams.update({
    "font.family": font_manager.FontProperties(fname=font_path).get_name(),
    "font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
    "axes.unicode_minus": False, "text.color": "#14273C",
    "axes.labelcolor": "#35475B", "xtick.color": "#35475B", "ytick.color": "#35475B",
})
fig = plt.figure(figsize=(13.2, 11.5), facecolor="white")
grid = fig.add_gridspec(3, 4, left=0.065, right=0.938, bottom=0.044, top=0.868,
                       width_ratios=(1, 1, 1, 0.045), height_ratios=(1.03, 0.84, 1.18),
                       hspace=0.57, wspace=0.28)
fig.text(0.065, 0.973, "case01 · 阵列周期与入瞳对成像的影响", weight="bold", size=23, va="top")
fig.text(0.065, 0.925, "相同 Mask 与 nk 膜层 · 3 mm 焦距 · 41 波长 D65 · 1 μm 像元",
         size=11, color="#537087", va="top")
peak = max(psf.max() for psf in green_psf)
for i, (psf, label) in enumerate(zip(green_psf, labels)):
    ax = fig.add_subplot(grid[0, i])
    log_psf = 10 * np.log10(np.maximum(psf[center-60:center+61, center-60:center+61] / peak, 1e-5))
    heatmap = ax.imshow(log_psf, extent=[-60.5, 60.5] * 2, origin="lower",
                       cmap="magma", vmin=-50, vmax=0, interpolation="nearest")
    ax.set_title(label, loc="left", weight="bold", pad=8)
    ax.set(xlabel="传感器 x (μm)", ylabel="传感器 y (μm)", xticks=[-60,-30,0,30,60],
           yticks=[-60,-30,0,30,60])
    ax.text(0.05, 0.94, f"G 通道 PSF · 半高全宽 {fwhm[i]:.2f} μm",
            va="top", transform=ax.transAxes, color="white", fontsize=9)
fig.colorbar(heatmap, cax=fig.add_subplot(grid[0, 3]), label="相对共同峰值 (dB)")

ax = fig.add_subplot(grid[1, :2])
for i, (label, color) in enumerate(zip(labels, colors)):
    ax.plot(freq_positive, mtf_green[i], color=color, linestyle=("--" if i == 1 else "-"),
            lw=2.0, label=f"{label}   MTF50 {mtf50[i]:.1f}")
ax.set_title("G 通道水平 MTF", loc="left", weight="bold", pad=8)
ax.set(xlabel="空间频率 (cycles/mm)", ylabel="MTF", xlim=(0,330), ylim=(0,1.02))
ax.grid(alpha=0.2)
ax.legend(loc="upper right", frameon=False, fontsize=8.5)

ax = fig.add_subplot(grid[1, 2])
pos = np.arange(3)
for i, (label, color) in enumerate(zip("RGB", channel_colors)):
    ax.bar(pos+(i-1)*0.21, relative_signal[:, i], width=0.19, color=color, label=label)
ax.set_title("入瞳面积修正后的相对信号", loc="left", weight="bold", pad=8, fontsize=11)
ax.set(ylabel="相对第 1 组 (×)", ylim=(0,4.65), xticks=pos, xticklabels=["第 1 组","第 2 组","第 3 组"])
ax.grid(axis="y", alpha=0.2)
ax.legend(loc="upper left", frameon=False, ncol=3, fontsize=8)

for i, trial in enumerate(trials):
    ax = fig.add_subplot(grid[2, i])
    linear = np.clip(trial["fixedExposure"] / captured[i][None,None,:], 0, 1)
    image = np.where(linear <= 0.0031308, 12.92 * linear, 1.055 * linear**(1/2.4) - 0.055)
    ax.imshow(image, interpolation="nearest")
    ax.set_title(f"亮度匹配成像 · 第 {i+1} 组", loc="left", weight="bold", pad=8)
    ax.set_axis_off()
fig.text(0.065, 0.014,
         "图像按通道匹配亮度以比较清晰度。相对信号按相同入射辐照度和曝光估算。",
         size=9.5, color="#537087")
fig.savefig(folder / "comparison.png", dpi=155, facecolor="white")
print(f"Green PSF FWHM (μm): {fwhm}")
print(f"Green MTF50 (cycles/mm): {mtf50}")
print(f"Pupil open fractions: {mask_fraction}")
print(f"Sampled pupil area ratios: {area/area[0]}")
print(f"Relative collected RGB signals: {relative_signal}")
