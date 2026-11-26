"""Recreate the effective amplitude pupil panel from MATLAB result data."""

from dataclasses import dataclass
from io import BytesIO

import h5py
import numpy as np
from PIL import Image


MAX_PUPIL_SIDE = 1500
WAVELENGTH_NM = 550


@dataclass
class PupilPreview:
    extent_um: tuple[float, float, float, float]
    wavelength_nm: float
    amplitude: bytes
    clear: bytes


def gray_png(values):
    # HDF5 stores MATLAB image arrays as x,y. Display positive y at the top.
    pixels = np.rint(255 * np.clip(values.T[::-1], 0, 1)).astype(np.uint8)
    png = BytesIO()
    Image.fromarray(pixels).save(png, format='PNG')
    return png.getvalue()


def sample_grid(dataset, stride):
    if stride == 1 or dataset.chunks is None or dataset.chunks[1] != 1:
        return dataset[::stride, ::stride]
    # MATLAB's column chunks can be read sequentially with much less overhead
    # than an HDF5 selection strided across both axes. Sample x in memory.
    shape = tuple((n + stride - 1) // stride for n in dataset.shape)
    sampled = np.empty(shape, dtype=dataset.dtype)
    for j, col in enumerate(range(0, dataset.shape[1], stride)):
        sampled[:, j] = dataset[:, col][::stride]
    return sampled


def read_pupil_preview(result):
    fields = ('optics/pupil', 'optics/regionMap', 'optics/transmissionLookup',
              'optics/x', 'optics/y', 'cfg/wavelengths')
    datasets = [result.get(key) for key in fields]
    if not all(isinstance(dataset, h5py.Dataset) for dataset in datasets):
        return None
    pupil, region, lookup, axis_x, axis_y, wavelength = datasets
    x = axis_x[()].ravel() * 1e6
    y = axis_y[()].ravel() * 1e6
    wavelengths = wavelength[()].ravel() * 1e9
    if (len(x) < 2 or len(y) < 2 or pupil.shape != (len(x), len(y)) or region.shape != pupil.shape
            or lookup.ndim != 2 or lookup.shape[0] != len(wavelengths) or not len(wavelengths)
            or not all(np.isfinite(v).all() for v in (x, y, wavelengths))):
        return None

    wave = int(np.argmin(np.abs(wavelengths - WAVELENGTH_NM)))
    transmission = lookup[wave]
    if transmission.dtype.names is not None:
        if set(transmission.dtype.names) != {'real', 'imag'}:
            return None
        amplitude = np.hypot(transmission['real'], transmission['imag'])
    else:
        amplitude = np.abs(transmission)
    if not np.isfinite(amplitude).all():
        return None

    # Match plotDiffractionDiagnostics.m, including its preview-only stride.
    stride = (max(pupil.shape) + MAX_PUPIL_SIDE - 1) // MAX_PUPIL_SIDE
    clear = sample_grid(pupil, stride)
    regions = sample_grid(region, stride)
    if not np.issubdtype(regions.dtype, np.integer) or regions.min() < 0 or regions.max() >= len(amplitude):
        return None
    extent = []
    for axis in (x, y):
        step = stride * (axis[1] - axis[0])
        if step <= 0:
            return None
        extent.extend((float(axis[0] - step / 2), float(axis[::stride][-1] + step / 2)))
    return PupilPreview(tuple(extent), float(wavelengths[wave]), gray_png(clear * amplitude[regions]), gray_png(clear))
