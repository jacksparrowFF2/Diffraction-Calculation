"""Prepare the Case5–8 inputs using the periodic centers in Case1."""

import json
from pathlib import Path
import shutil

import numpy as np
from PIL import Image


inputs = Path(__file__).resolve().parents[2] / "Input"


def settings(case):
    return json.loads((inputs / case / "config.json").read_text(encoding="utf-8"))


def smaller_openings(original, pitch_um):
    size = original.shape[0]
    half = size // 2
    axis = (np.arange(size) + 0.5 - half) * pitch_um / size
    y, x = np.meshgrid(axis, axis, indexing="ij")
    nominal_radius = np.hypot(x, y)
    reduced = np.zeros_like(original)
    # Roll the four corners to the center before fitting their shared circle.
    for shift in (0, half):
        opening = np.roll(original, (shift, shift), axis=(0, 1))
        interior = opening.copy()
        for offset, dimension in ((1, 0), (-1, 0), (1, 1), (-1, 1)):
            interior &= np.roll(opening, offset, axis=dimension)
        # Limit the fit to the existing circular boundary near the nominal radius.
        edge = opening & ~interior & (abs(nominal_radius - 59) < 1)
        col, row = x[edge], y[edge]
        cx, cy, intercept = np.linalg.lstsq(
            np.column_stack((2 * col, 2 * row, np.ones_like(col))),
            col**2 + row**2, rcond=None,
        )[0]
        radius = np.sqrt(intercept + cx**2 + cy**2)
        assert abs(radius - 59) < 0.3 and np.hypot(cx, cy) < 0.5
        resized = np.hypot(x - cx, y - cy) <= 54
        assert np.all(opening[resized])
        print(f"Roll {shift}: center offset ({cx:.3f}, {cy:.3f}) um, "
              f"source inner edge {radius:.3f} um, generated radius 54 um")
        reduced |= np.roll(resized, (-shift, -shift), axis=(0, 1))
    return reduced


def write_config(case, config):
    folder = inputs / case
    folder.mkdir(exist_ok=True)
    (folder / "config.json").write_text(
        json.dumps(config, indent=2) + "\n", encoding="utf-8")
    return folder


def main():
    case1 = settings("case1")
    case4 = settings("case4")
    with Image.open(inputs / "case1/BM.png") as image:
        pixels = np.asarray(image.convert("L"))
    size = case1["stack"]["geometry"]["inputSamplesPerPeriod"]
    assert pixels.shape == (size, size) and set(np.unique(pixels)) == {0, 255}
    openings = smaller_openings(pixels == 255, case1["cfg"]["pitch"] * 1e6)
    mask = Image.fromarray(openings.astype(np.uint8) * 255)
    for case, config in (("case5", settings("case2")), ("case6", case1)):
        config["layerImageDir"] = f"topic_udc_imagingV2/Input/{case}"
        folder = write_config(case, config)
        materials = folder / "Materials"
        materials.mkdir(exist_ok=True)
        shutil.copyfile(inputs / "case1/Materials/BM_nk.csv", materials / "BM_nk.csv")
        if case == "case5":
            shutil.copyfile(inputs / "case1/BM.png", folder / "BM.png")
            shutil.copyfile(inputs / "case2/Materials/Cath_nk.csv", materials / "Cath_nk.csv")
            mask.save(folder / "Cath.png")
        else:
            mask.save(folder / "BM.png")
    write_config("case7", {**case4, "reuseProject": "case4"})
    write_config("case8", {**case4, "reuseProject": "case7", "neutralTransmittance": 0.5})
    pitch_um = case1["cfg"]["pitch"] * 1e6
    for samples, raster in ((1000, openings), (500, openings[1::2, 1::2])):
        fraction = raster.mean()
        equivalent_radius = pitch_um * np.sqrt(fraction / (2 * np.pi))
        assert abs(equivalent_radius - 54) < 0.05
        print(f"{samples} points/Pitch: {fraction:.6%} open, "
              f"area-equivalent radius {equivalent_radius:.3f} um")


if __name__ == "__main__":
    main()
