"""Shared conventions for the static and interactive comparisons."""

import csv
from io import BytesIO

import numpy as np
from PIL import Image


DIRECTIONS = ('Horizontal', 'Vertical', 'Diagonal45', 'Diagonal135')
PSF_FLOOR_DB = -70
PSF_COLORMAP = 'magma'
RGB_TONE_STRENGTH = 60
RGB_LUMINANCE_GAMMA = 0.72
RGB_LUMINANCE_WEIGHTS = np.array([0.2126, 0.7152, 0.0722])


def rgb_peak(kernels):
    """Find the shared exposure reference, including available clear apertures."""
    return max((float(kernel.max()) for kernel in kernels if kernel is not None), default=0.0)


def rgb_psf_png(kernel, reference):
    """Apply the COE display mapping to UDC's existing linear RGB channels."""
    rgb = np.asarray(kernel).transpose(2, 1, 0)[::-1]
    if reference > 0:
        rgb = (rgb / reference) * 0.9
    luminance = rgb @ RGB_LUMINANCE_WEIGHTS
    target = (np.log1p(RGB_TONE_STRENGTH * luminance) /
              np.log1p(RGB_TONE_STRENGTH)) ** RGB_LUMINANCE_GAMMA
    scale = np.divide(target, luminance, out=np.ones_like(luminance), where=luminance > 0)
    rgb = rgb * scale[..., None]
    rgb /= np.maximum(1, rgb.max(axis=2, keepdims=True))
    mapped = np.where(rgb <= 0.0031308, 12.92 * rgb, 1.055 * rgb**(1 / 2.4) - 0.055)
    png = BytesIO()
    Image.fromarray(np.rint(255 * np.clip(mapped, 0, 1)).astype(np.uint8)).save(png, format='PNG', optimize=True)
    return png.getvalue()


def read_csv(folder, name):
    with (folder / name).open(encoding='utf-8', newline='') as stream:
        return list(csv.DictReader(stream))


def green_mtf(profiles):
    paths = [[r for r in profiles if r['Channel'] == 'G' and r['Direction'] == direction]
             for direction in DIRECTIONS]
    frequencies = np.array([[float(r['CyclesPerMm']) for r in path] for path in paths])
    for frequency in frequencies[1:]:
        np.testing.assert_array_equal(frequencies[0], frequency)
    mtf = np.array([[float(r['MTF']) for r in path] for path in paths])
    clear = np.array([[float(r['ClearMTF']) for r in path] for path in paths])
    return frequencies[0], mtf.min(axis=0), clear.min(axis=0)
