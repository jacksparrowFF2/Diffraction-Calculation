"""Render the shared self-contained HTML comparison page."""

from base64 import b64encode
from io import BytesIO
import json
from pathlib import Path

from matplotlib import colormaps
import numpy as np
from PIL import Image

from .viewer_data import PSF_COLORMAP, PSF_FLOOR_DB, RGB_LUMINANCE_GAMMA, RGB_TONE_STRENGTH


HERE = Path(__file__).resolve().parent


def png_url(data):
    return 'data:image/png;base64,' + b64encode(data).decode('ascii')


def mask_payload(previews):
    if previews is None:
        return None
    return dict(widthUm=previews.width_um, heightUm=previews.height_um, tiled=previews.tiled,
                subsampled=previews.subsampled,
                layers={name: png_url(png) for name, png in previews.layers.items()})


def pupil_payload(preview, clear):
    if preview is None:
        return None
    return dict(extent=preview.extent_um, wavelengthNm=preview.wavelength_nm,
                png=png_url(preview.clear if clear else preview.amplitude))


def psf_url(psf, reference):
    scale = 10**(PSF_FLOOR_DB / 10)
    db = 10 * np.log10(np.maximum(psf / reference, scale))
    index = np.rint(255 * np.clip(1 - db / PSF_FLOOR_DB, 0, 1)).astype(np.uint8)
    png = BytesIO()
    # Store dB indices as grayscale so the browser can recolor without rescaling.
    # Positive y is at the top, as in the static plot.
    Image.fromarray(index[::-1]).save(png, format='PNG', optimize=True)
    return png_url(png.getvalue())


def psf_palettes():
    levels = np.linspace(0, 1, 256)
    palettes = {name: colormaps[name](levels, bytes=True)[:, :3].tolist()
                for name in (PSF_COLORMAP, 'viridis', 'gray')}
    matlab = json.loads((HERE / 'matlab_colormaps.json').read_text(encoding='utf-8'))['palettes']
    palettes['jet'] = matlab['jet']
    palettes.update(matlab)
    return palettes


def render_page(data):
    rgb_display = {'toneStrength': RGB_TONE_STRENGTH, 'gamma': RGB_LUMINANCE_GAMMA}
    data = json.dumps(dict(data, rgbDisplay=rgb_display), ensure_ascii=False, allow_nan=False, separators=(',', ':'))
    template = (HERE / 'interactive_comparison.html').read_text(encoding='utf-8')
    script = (HERE / 'interactive_comparison.js').read_text(encoding='utf-8')
    html = template.replace('/* UDC_SCRIPT */', script, 1)
    return html.replace('<!-- UDC_DATA -->', data.replace('<', '\\u003c'), 1)
