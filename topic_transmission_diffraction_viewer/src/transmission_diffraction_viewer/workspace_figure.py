"""Adapt selected V2 cases to the shared offline comparison page."""

import numpy as np

from .case_compatibility import RULES, common_band, common_tail_radii, evaluate
from .case_workspace import RGB
from .comparison_page import mask_payload, png_url, psf_palettes, psf_url, pupil_payload, render_page
from .viewer_data import PSF_FLOOR_DB, rgb_peak, rgb_psf_png


CASE_COLORS = ('#2379ab', '#be552b', '#228475', '#8756b2', '#b34d72', '#996811')
REFERENCE_COLORS = ('#596678', '#776f88', '#567b79')


def reference_groups(cases):
    groups = {}
    for index, case in enumerate(cases):
        key = tuple(case.conditions.get(field) for field in ('clear', 'grid', 'scene', 'clear_preview'))
        groups.setdefault(key if all(key) else (index,), []).append(case)
    for group in groups.values():
        name = ('共同无遮挡参考' if len(group) == len(cases) and len(cases) > 1 else
                '无遮挡 · ' + ' / '.join(case.name for case in group))
        representative = max(group, key=lambda case: (case.clear_rgb_kernel is not None,
                                                      case.pupil_preview is not None,
                                                      case.clear_rgb_psf is not None))
        yield name, representative


def workspace_page(cases, show_clear):
    full_radius = round(min(min(-c.axis_um[0], c.axis_um[-1]) for c in cases), 9)
    if full_radius <= 0:
        raise ValueError('所选 PSF 没有以零点为中心的共同接收范围。')
    band = common_band(cases)
    radii = common_tail_radii(cases)
    rgb_reference = rgb_peak(kernel for case in cases for kernel in (case.rgb_kernel, case.clear_rgb_kernel))

    def make_record(case, name, key, color, clear=False, dash=()):
        axis = case.axis_um
        step = float(axis[1] - axis[0])
        if not np.allclose(np.diff(axis), step, rtol=1e-9, atol=1e-9):
            raise ValueError(f'{name} 的像面采样间距不均匀。')
        keep = np.flatnonzero(np.abs(axis) <= full_radius + step / 2 + 1e-9)
        kernel = case.clear_kernel if clear else case.kernel
        crop = kernel[np.ix_(keep, keep)]
        energy = kernel.sum()
        normalized = crop / energy if energy > 0 else crop
        mtf = case.clear_mtf if clear else case.mtf
        freq = []
        curve = []
        band_min = mtf_300 = None
        if band is not None:
            lo, hi = band
            freq = np.r_[lo, case.freq[(case.freq > lo) & (case.freq < hi)], hi].tolist()
            curve = np.interp(freq, case.freq, mtf).tolist()
            if lo <= 10 and hi >= 100:
                nodes = np.r_[10, case.freq[(case.freq > 10) & (case.freq < 100)], 100]
                band_min = float(np.min(np.interp(nodes, case.freq, mtf)))
            if lo <= 300 <= hi:
                mtf_300 = float(np.interp(300, case.freq, mtf))
        captured = [case.channels[ch][2 if clear else 1] for ch in RGB]
        denominators = [case.channels[ch][2] for ch in RGB]
        tail = {row[0]: row for row in case.tails}
        center = int(np.argmin(np.abs(axis)))
        rgb_kernel = case.clear_rgb_kernel if clear else case.rgb_kernel
        rgb = case.clear_rgb_psf if clear else case.rgb_psf
        if rgb_kernel is not None:
            rgb = rgb_psf_png(rgb_kernel, rgb_reference)
        return dict(
            id=key, name=name, color=color, selected=True, reference=clear, dash=dash,
            psf=psf_url(normalized, 1), absolutePsf=psf_url(crop, 1),
            psfExtent=[float(axis[keep[0]] - step/2), float(axis[keep[-1]] + step/2)],
            rgbPsf=png_url(rgb) if rgb is not None else None,
            rgbLegacy=rgb is not None and rgb_kernel is None,
            rgbExtent=[float(axis[0] - step/2), float(axis[-1] + step/2)],
            maskPreview=None if clear else mask_payload(case.mask_previews),
            pupilPreview=pupil_payload(case.pupil_preview, clear),
            relative=[num / den if den > 0 else None for num, den in zip(captured, denominators)],
            captured=captured, total=[1 if clear else case.channels[ch][0] for ch in RGB],
            freq=freq, mtf=curve, bandMin=band_min, mtf300=mtf_300,
            tailFraction=[tail[r][3 if clear else 1] for r in radii],
            tailAbsolute=[tail[r][4 if clear else 2] for r in radii],
            center=float(kernel[center, center] / energy) if energy else None,
            scene=png_url(case.clear_scene if clear else case.scene),
            mean=[None, 1 if clear else case.conditions.get('scene_mean_g'), None],
        )

    records = []
    if show_clear:
        for i, (name, case) in enumerate(reference_groups(cases)):
            records.append(make_record(case, name, f'clear-{i + 1}', REFERENCE_COLORS[i % len(REFERENCE_COLORS)],
                                       clear=True))
    records.extend(make_record(case, case.name, f'case-{i + 1}', CASE_COLORS[i % len(CASE_COLORS)],
                               dash=(7, 5) if i % 2 else ())
                   for i, case in enumerate(cases))
    caution = len(cases) > 1 and any(evaluate(cases, metric)[0] != '可直接比较' for metric in RULES)
    return render_page({'records': records, 'palettes': psf_palettes(), 'floorDb': PSF_FLOOR_DB,
                        'workspace': {'fullRadius': full_radius, 'zoomRadius': min(25., full_radius),
                                      'band': band, 'tailRadii': radii, 'caution': caution}})
