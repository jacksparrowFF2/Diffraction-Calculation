"""Read display previews from the material coverage saved by the UDC solver."""

from dataclasses import dataclass
from io import BytesIO

import h5py
import numpy as np
from PIL import Image


MASK_LAYERS = ('BM', 'Cath')
MAX_SIDE = 750


@dataclass
class MaskPreviews:
    width_um: float
    height_um: float
    tiled: bool
    subsampled: bool
    layers: dict[str, bytes]


def read_mask_previews(result):
    geometry = result.get('optics/layerGeometry')
    if not isinstance(geometry, h5py.Group):
        return None
    masks = geometry.get('layerMasks')
    if not isinstance(masks, h5py.Group) or 'gridStep' not in geometry:
        return None
    datasets = {name: masks[name] for name in MASK_LAYERS if name in masks}
    grid = next(iter(datasets.values()), geometry.get('unionMask'))
    step = float(geometry['gridStep'][()].item()) * 1e6
    if grid is None or grid.ndim != 2 or not np.isfinite(step) or step <= 0:
        return None

    # The solver requires all layer images to start on the same cell. A full
    # coverage layer also repeats on that cell, regardless of how it was input.
    nx, ny = grid.shape
    tiled = False
    if 'repeatCounts' in geometry:
        repeats = geometry['repeatCounts'][()].ravel().tolist()
        if (len(repeats) == 2 and all(v >= 1 and float(v).is_integer() for v in repeats)
                and nx % int(repeats[0]) == ny % int(repeats[1]) == 0):
            nx //= int(repeats[0])
            ny //= int(repeats[1])
            tiled = max(repeats) > 1
    stride = (max(nx, ny) + MAX_SIDE - 1) // MAX_SIDE
    pngs = {}
    for name, dataset in datasets.items():
        # MATLAB's HDF5 order is x,y. Flip y for a Cartesian preview, like PSF.
        covered = dataset[:nx:stride, :ny:stride].T[::-1].astype(bool)
        pixels = np.where(covered, np.uint8(32), np.uint8(255))
        png = BytesIO()
        Image.fromarray(pixels).save(png, format='PNG')
        pngs[name] = png.getvalue()
    return MaskPreviews(nx * step, ny * step, tiled, stride > 1, pngs)
