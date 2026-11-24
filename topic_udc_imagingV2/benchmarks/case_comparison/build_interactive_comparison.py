"""Build a self-contained HTML comparison from existing UDC outputs."""

import argparse
from base64 import b64encode
from io import BytesIO
import json
from pathlib import Path

import h5py
from matplotlib import colormaps
import numpy as np
from PIL import Image

from comparison_data import (PSF_COLORMAP, PSF_FLOOR_DB, RGB_LUMINANCE_GAMMA, RGB_TONE_STRENGTH,
                             green_mtf, read_csv, rgb_peak, rgb_psf_png)


HERE = Path(__file__).resolve().parent
OUTPUT = HERE.parents[1] / 'Output'
CASES = (
    ('case1', 'Case1 · BM 59 μm', '#176b91'),
    ('case2', 'Case2 · BM + 图案阴极', '#be512f'),
    ('case3', 'Case3 · BM', '#427252'),
    ('case4', 'Case4 · BM + 连续阴极', '#8261bd'),
    ('case5', 'Case5 · BM 59 / 阴极 54 μm', '#ad4b72'),
    ('case6', 'Case6 · BM 54 μm', '#168478'),
    ('case7', 'Case7 · BM 59 μm + 连续阴极', '#3348a0'),
    ('case8', 'Case8 · Case7 + 50% 衰减', '#a26a12'),
)
DEFAULT = {'case1', 'case5', 'case6', 'case7', 'case8'}


def png_url(data):
    return 'data:image/png;base64,' + b64encode(data).decode('ascii')


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


def scene_pixels(png):
    with Image.open(BytesIO(png)) as image:
        return np.asarray(image.convert('RGB'))


def read_case(folder, case):
    with h5py.File(folder / 'udc_results.mat') as f:
        expected = {'pitch': 236.2e-6, 'samplesPerPeriod': 500, 'pupilRadius': 0.002,
                    'focalLength': 0.003, 'sensorPitch': 1e-6, 'outputSize': 201}
        for key, value in expected.items():
            np.testing.assert_allclose(f[f'cfg/{key}'][()].item(), value, atol=0, rtol=1e-12)
        np.testing.assert_array_equal(f['cfg/periods'][()].ravel(), [17, 17])
        np.testing.assert_allclose(f['cfg/wavelengths'][()].ravel(), np.arange(380, 781, 10) * 1e-9,
                                   atol=0, rtol=1e-12)
        np.testing.assert_allclose(f['optics/coordinate'][()].ravel()*1e6, np.arange(-100, 101),
                                   atol=1e-10, rtol=0)
        psf = f['optics/psf'][1].T
        kernel_rgb = f['optics/kernel'][()]
        absolute_psf = kernel_rgb[1].T
        clear_psf = f['optics/quality/clearPsf'][1].T
        inputs = {key: f[key][()] for key in ('cfg/spectralWeights', 'cfg/channelResponse',
                                             'images/inputLinear', 'optics/clearKernel')}
        means = f['images/fixedExposure'][()].mean(axis=(1, 2))
        clear_means = f['images/clear'][()].mean(axis=(1, 2))
    for array in (psf, clear_psf):
        assert array.shape == (201, 201) and np.all(array >= 0), case
        np.testing.assert_allclose(array.sum(), 1, atol=1e-10, rtol=0)
    channel = {r['Channel']: r for r in read_csv(folder, 'udc_channel_metrics.csv')}
    radii = [10, 20, 50]
    tail = {float(r['RadiusUm']): r for r in read_csv(folder, 'udc_psf_tail.csv')
            if r['Channel'] == 'G' and r['RadiusFitsWindow'] == '1'}
    mtf_summary = [r for r in read_csv(folder, 'udc_mtf_directional_summary.csv') if r['Channel'] == 'G']
    freq, mtf, clear_mtf = green_mtf(read_csv(folder, 'udc_mtf_directional.csv'))

    def metrics(reference):
        capture_key = 'ClearCapturedThroughput' if reference else 'CapturedThroughput'
        fraction_key = 'ClearTailFractionWithinWindow' if reference else 'TailFractionWithinWindow'
        tail_key = 'ClearTailThroughput' if reference else 'TailThroughput'
        band_key = 'ClearBandMinMTF' if reference else 'BandMinMTF'
        curve = clear_mtf if reference else mtf
        heatmap = clear_psf if reference else psf
        kernel = inputs['optics/clearKernel'][1].T if reference else absolute_psf
        captured = [float(channel[c][capture_key]) for c in 'RGB']
        np.testing.assert_allclose(kernel, heatmap * captured[1], rtol=1e-12, atol=0,
                                   err_msg=f'{case}: inconsistent G-channel PSF and throughput')
        return {
            'freq': freq.tolist(), 'mtf': curve.tolist(), 'psf': heatmap,
            'absolutePsf': kernel,
            'rgbKernel': inputs['optics/clearKernel'] if reference else kernel_rgb,
            'relative': [1 if reference else float(channel[c]['RelativeToClear']) for c in 'RGB'],
            'captured': captured,
            'total': [1 if reference else float(channel[c]['Throughput']) for c in 'RGB'],
            'tailFraction': [float(tail[r][fraction_key]) for r in radii],
            'tailAbsolute': [float(tail[r][tail_key]) for r in radii],
            'bandMin': min(float(r[band_key]) for r in mtf_summary),
            'mtf300': float(np.interp(300, freq, curve)),
            'center': float(heatmap[100, 100]),
            'mean': ((clear_means if reference else means) / clear_means).tolist(),
            'scene': (folder / ('udc_clear.png' if reference else 'udc_fixed_exposure.png')).read_bytes(),
        }

    return metrics(False), metrics(True), inputs


def payload(results):
    records = []
    baseline = None
    for case, name, color in CASES:
        data, clear, inputs = read_case(results / case / 'nk', case)
        if baseline is None:
            baseline = clear
            baseline_inputs = inputs
            baseline_scene = scene_pixels(clear['scene'])
        else:
            for key, array in inputs.items():
                np.testing.assert_array_equal(array, baseline_inputs[key], err_msg=f'{case}: {key}')
            np.testing.assert_allclose(clear['psf'], baseline['psf'], atol=1e-12, rtol=0)
            np.testing.assert_array_equal(clear['freq'], baseline['freq'])
            np.testing.assert_allclose(clear['mtf'], baseline['mtf'], atol=1e-12, rtol=0)
            np.testing.assert_array_equal(scene_pixels(clear['scene']), baseline_scene, err_msg=case)
            assert all(clear[key] == baseline[key] for key in
                       ('captured', 'tailFraction', 'tailAbsolute', 'mean')), case
        np.testing.assert_array_equal(data['freq'], baseline['freq'])
        records.append(dict(data, id=case, name=name, color=color, selected=case in DEFAULT))
        print(f'Loaded {case}')
    records.insert(0, dict(baseline, id='clear', name='无遮挡', color='#596678', selected=True))
    peak = max(record['psf'].max() for record in records)
    rgb_reference = rgb_peak(record['rgbKernel'] for record in records)
    for record in records:
        record['psf'] = psf_url(record['psf'], peak)
        # The kernel stores energy per sensor bin as a fraction of incident energy.
        record['absolutePsf'] = psf_url(record['absolutePsf'], 1)
        record['rgbPsf'] = png_url(rgb_psf_png(record.pop('rgbKernel'), rgb_reference))
        record['scene'] = png_url(record['scene'])
    return {'records': records, 'palettes': psf_palettes(), 'floorDb': PSF_FLOOR_DB}


def render_page(data):
    rgb_display = {'toneStrength': RGB_TONE_STRENGTH, 'gamma': RGB_LUMINANCE_GAMMA}
    data = json.dumps(dict(data, rgbDisplay=rgb_display), ensure_ascii=False, allow_nan=False, separators=(',', ':'))
    template = (HERE / 'interactive_comparison.html').read_text(encoding='utf-8')
    script = (HERE / 'interactive_comparison.js').read_text(encoding='utf-8')
    html = template.replace('/* UDC_SCRIPT */', script, 1)
    return html.replace('<!-- UDC_DATA -->', data.replace('<', '\\u003c'), 1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, default=OUTPUT, help='Directory containing case1–case8.')
    parser.add_argument('--output', type=Path, help='Destination HTML file.')
    args = parser.parse_args()
    output = args.output or args.results / 'interactive_comparison' / 'udc_comparison.html'
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(render_page(payload(args.results)), encoding='utf-8')
    print(output.resolve())


if __name__ == '__main__':
    main()
