"""Read UDC V2 results and exchange compact, versioned case packages."""

import argparse
from dataclasses import dataclass
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

import h5py
import numpy as np
from PIL import Image

from comparison_data import green_mtf, read_csv
from mask_previews import MASK_LAYERS, MAX_SIDE, MaskPreviews, read_mask_previews


OUTPUT = Path(__file__).resolve().parents[2] / 'Output'
RESULT_FILES = (
    'udc_results.mat', 'udc_channel_metrics.csv', 'udc_mtf_directional.csv',
    'udc_psf_tail.csv', 'udc_fixed_exposure.png', 'udc_clear.png',
)
MAX_AXIS = 2049
MAX_PACKAGE_MIB = 224
MAX_PACKAGE = MAX_PACKAGE_MIB * 1024 * 1024
MAX_IMAGE = 8 * 1024 * 1024
RGB = 'RGB'


@dataclass
class Case:
    name: str
    conditions: dict
    axis_um: np.ndarray
    kernel: np.ndarray
    clear_kernel: np.ndarray
    freq: np.ndarray
    mtf: np.ndarray
    clear_mtf: np.ndarray
    channels: dict
    tails: list
    scene: bytes
    clear_scene: bytes
    rgb_psf: bytes | None = None
    clear_rgb_psf: bytes | None = None
    rgb_kernel: np.ndarray | None = None
    clear_rgb_kernel: np.ndarray | None = None
    mask_previews: MaskPreviews | None = None


def fingerprint(*arrays):
    digest = sha256()
    for array in arrays:
        array = np.asarray(array, dtype='<f8', order='C')
        digest.update(repr(array.shape).encode('ascii'))
        digest.update(array.tobytes())
    return digest.hexdigest()


def discover(root):
    paths = (*root.glob('*/udc_results.mat'), *root.glob('*/*/udc_results.mat'))
    return sorted(path.parent for path in paths
                  if all((path.parent / name).is_file() for name in RESULT_FILES))


def source_stamp(folder):
    return tuple((p.stat().st_mtime_ns, p.stat().st_size) for p in
                 (folder / name for name in RESULT_FILES))


def preview_fingerprint(preview):
    if len(preview) > MAX_IMAGE or not preview.startswith(b'\x89PNG\r\n\x1a\n'):
        raise ValueError('场景预览必须为有效大小的 PNG。')
    with Image.open(BytesIO(preview)) as image:
        if image.width * image.height > 16_000_000:
            raise ValueError('场景预览超过 1600 万像素。')
        pixels = np.asarray(image.convert('RGB'))
    return sha256(repr(pixels.shape).encode('ascii') + pixels.tobytes()).hexdigest()


def read_result(folder, name):
    with h5py.File(folder / 'udc_results.mat') as f:
        cfg = f['cfg']

        def scalar(key, factor=1):
            value = float(cfg[key][()].item()) * factor if key in cfg else None
            return value if value is not None and np.isfinite(value) else None

        def spectral_signature():
            keys = ('wavelengths', 'spectralWeights', 'channelResponse')
            if all(key in cfg for key in keys):
                wavelengths, weights, response = (cfg[key][()] for key in keys)
                combined = response * weights.ravel()
                normalized = combined / combined.sum(axis=1, keepdims=True)
                return fingerprint(wavelengths.ravel(), np.round(normalized, 14))

        saved_config = f.get('projectConfig')
        project = None
        if saved_config is not None and saved_config.attrs.get('MATLAB_class') == b'char':
            project = json.loads(saved_config[()].astype('<u2').tobytes().decode('utf-16-le'))
        project_cfg = project['cfg'] if project else {}
        receiver = project_cfg.get('receiverMethod', 'midpoint') if project else None
        all_kernel = f['optics/kernel'][()]
        kernel = all_kernel[1].T
        all_clear = f['optics/clearKernel'][()]
        input_scene = f.get('images/inputLinear')
        fixed_scene = f.get('images/fixedExposure')
        clear_scene = f.get('images/clear')
        axis = f['optics/coordinate'][()].ravel() * 1e6
        previews = read_mask_previews(f)
        wavelengths = cfg['wavelengths'][()].ravel() * 1e9 if 'wavelengths' in cfg else []
        pitch = scalar('pitch', 1e6)
        samples = scalar('samplesPerPeriod')
        periods = cfg['periods'][()].ravel().astype(int).tolist() if 'periods' in cfg else None
        conditions = {
            'spectrum': spectral_signature(),
            'wavelength_range': [round(float(wavelengths[0]), 6),
                                 round(float(wavelengths[-1]), 6), len(wavelengths)] if len(wavelengths) else None,
            'pupil_mm': ('无限入瞳' if 'pupilRadius' in cfg and np.isposinf(cfg['pupilRadius'][()].item())
                         else scalar('pupilRadius', 2e3)),
            'focus_mm': scalar('focalLength', 1e3),
            'pixel_um': scalar('sensorPitch', 1e6),
            'object_step_um': pitch / samples if pitch and samples else None,
            'periods': periods,
            'receiver': receiver,
            'receiver_samples': scalar('receiverFftPaddingFactor' if receiver == 'fftBins' else
                                       'receiverBinSamples') if receiver else None,
            'grid': fingerprint(np.round(axis, 9)),
            'window_um': [round(float(axis[0]), 9), round(float(axis[-1]), 9), len(axis)],
            'energy_basis': 'V2 / incident pupil energy' if project else None,
            'clear': fingerprint(all_clear),
            'scene': fingerprint(input_scene[()]) if input_scene is not None else None,
            'scene_shape': list(fixed_scene.shape) if fixed_scene is not None else None,
            'exposure': 'V2 fixedExposure / unit gain' if project and fixed_scene is not None else None,
            'mapping': 'V2 / linear clamp 0-1 to sRGB' if project and fixed_scene is not None else None,
        }
        if fixed_scene is not None and clear_scene is not None:
            denominator = float(clear_scene[1].mean())
            if denominator > 0:
                conditions['scene_mean_g'] = float(fixed_scene[1].mean()) / denominator
    freq, mtf, clear_mtf = green_mtf(read_csv(folder, 'udc_mtf_directional.csv'))
    defined = np.isfinite(freq) & np.isfinite(mtf) & np.isfinite(clear_mtf)
    channels = {}
    for row in read_csv(folder, 'udc_channel_metrics.csv'):
        channels[row['Channel']] = [float(row[key]) for key in
                                    ('Throughput', 'CapturedThroughput', 'ClearCapturedThroughput')]
    tails = [[float(row[key]) for key in ('RadiusUm', 'TailFractionWithinWindow', 'TailThroughput',
                                         'ClearTailFractionWithinWindow', 'ClearTailThroughput')]
             for row in read_csv(folder, 'udc_psf_tail.csv')
             if row['Channel'] == 'G' and row['RadiusFitsWindow'] == '1']
    return validate(Case(name, conditions, axis, kernel, all_clear[1].T,
                         freq[defined], mtf[defined], clear_mtf[defined], channels, tails,
                         (folder / 'udc_fixed_exposure.png').read_bytes(),
                         (folder / 'udc_clear.png').read_bytes(),
                         rgb_kernel=all_kernel, clear_rgb_kernel=all_clear, mask_previews=previews))


def validate(case):
    n = len(case.axis_um)
    if (not isinstance(case.name, str) or not case.name.strip() or len(case.name) > 120
            or not isinstance(case.conditions, dict)):
        raise ValueError('方案名称或条件清单无效。')
    for key in ('clear', 'grid', 'scene', 'spectrum', 'receiver', 'exposure', 'mapping', 'energy_basis'):
        if case.conditions.get(key) is not None and not isinstance(case.conditions[key], str):
            raise ValueError(f'条件字段 {key} 的类型无效。')
    if (case.axis_um.ndim != 1 or not 3 <= n <= MAX_AXIS
            or not np.all(np.isfinite(case.axis_um)) or not np.all(np.diff(case.axis_um) > 0)):
        raise ValueError('PSF 坐标必须为递增的一维有限数组。')
    for kernel in (case.kernel, case.clear_kernel):
        if kernel.shape != (n, n) or not np.all(np.isfinite(kernel)) or np.any(kernel < 0):
            raise ValueError('PSF 能量网格与坐标不一致或含有无效值。')
    for rgb, green in ((case.rgb_kernel, case.kernel), (case.clear_rgb_kernel, case.clear_kernel)):
        if rgb is not None:
            if rgb.shape != (3, n, n) or not np.all(np.isfinite(rgb)) or np.any(rgb < 0):
                raise ValueError('RGB PSF 能量与接收网格不一致或含有无效值。')
            np.testing.assert_array_equal(rgb[1].T, green, err_msg='RGB 能量中的 G 通道与 PSF 不一致')
    if (case.freq.ndim != 1 or len(case.freq) < 2 or not np.all(np.diff(case.freq) > 0)
            or not np.all(np.isfinite(case.freq)) or case.freq[0] < 0
            or any(v.shape != case.freq.shape or not np.all(np.isfinite(v))
                   for v in (case.mtf, case.clear_mtf))):
        raise ValueError('MTF 曲线或频率坐标无效。')
    for channel in RGB:
        metrics = np.asarray(case.channels[channel], dtype=float)
        if metrics.shape != (3,) or not np.all(np.isfinite(metrics)) or np.any(metrics < 0):
            raise ValueError(f'{channel} 通量数值无效。')
    np.testing.assert_allclose([case.kernel.sum(), case.clear_kernel.sum()], case.channels['G'][1:],
                               atol=1e-12, rtol=1e-8, err_msg='G 通道通量与接收能量网格不一致')
    for tail in case.tails:
        values = np.asarray(tail, dtype=float)
        if values.shape != (5,) or not np.all(np.isfinite(values)) or np.any(values < 0):
            raise ValueError('拖尾数据无效。')
    preview_fingerprint(case.scene)
    case.conditions['clear_preview'] = preview_fingerprint(case.clear_scene)
    for rgb in (case.rgb_psf, case.clear_rgb_psf):
        if rgb is not None:
            if len(rgb) > MAX_IMAGE or not rgb.startswith(b'\x89PNG\r\n\x1a\n'):
                raise ValueError('RGB PSF 预览必须为有效大小的 PNG。')
            with Image.open(BytesIO(rgb)) as image:
                if image.size != (n, n) or image.mode != 'RGB':
                    raise ValueError('RGB PSF 图片尺寸或通道与接收网格不一致。')
                image.verify()
    return case


def write_package(case, target):
    manifest = dict(format='udc-case', version=1, name=case.name, conditions=case.conditions,
                    axis_um=case.axis_um.tolist(), freq=case.freq.tolist(), mtf=case.mtf.tolist(),
                    clear_mtf=case.clear_mtf.tolist(), channels=case.channels, tails=case.tails)
    previews = case.mask_previews
    if previews is not None:
        manifest['mask_preview'] = dict(width_um=previews.width_um, height_um=previews.height_um,
                                        tiled=previews.tiled, subsampled=previews.subsampled)
    with ZipFile(target, 'w', ZIP_DEFLATED) as archive:
        archive.writestr('case.json', json.dumps(manifest, ensure_ascii=False, allow_nan=False))
        for key in ('kernel', 'clear_kernel'):
            archive.writestr(f'{key}.f64', np.asarray(getattr(case, key), dtype='<f8').tobytes())
        # Reuse the existing G plane in v1 packages and store only extra R/B planes.
        for key in ('rgb_kernel', 'clear_rgb_kernel'):
            rgb = getattr(case, key)
            if rgb is not None:
                archive.writestr(f'{key}_rb.f64', np.asarray(rgb[[0, 2]], dtype='<f8').tobytes())
        archive.writestr('scene.png', case.scene)
        archive.writestr('clear.png', case.clear_scene)
        if case.rgb_psf is not None:
            archive.writestr('rgb-psf.png', case.rgb_psf)
        if case.clear_rgb_psf is not None:
            archive.writestr('clear-rgb-psf.png', case.clear_rgb_psf)
        if previews is not None:
            for name, png in previews.layers.items():
                archive.writestr(f'masks/{name}.png', png)


def read_package(data):
    with ZipFile(BytesIO(data)) as archive:
        if sum(info.file_size for info in archive.infolist()) > MAX_PACKAGE:
            raise ValueError(f'解压后的 Case 数据包超过 {MAX_PACKAGE_MIB} MiB。')
        meta = json.loads(archive.read('case.json'))
        if meta.get('format') != 'udc-case' or meta.get('version') != 1:
            raise ValueError('需要由打包脚本生成的 UDC Case 数据包 v1。')
        axis = np.asarray(meta['axis_um'], dtype=float)
        n = len(axis)
        if not 3 <= n <= MAX_AXIS:
            raise ValueError('数据包的 PSF 网格尺寸超出支持范围。')

        def kernel(key):
            buf = archive.read(f'{key}.f64')
            if len(buf) != n * n * 8:
                raise ValueError('数据包的 PSF 数组尺寸错误。')
            return np.frombuffer(buf, dtype='<f8').reshape(n, n)

        names = set(archive.namelist())

        def optional(name):
            return archive.read(name) if name in names else None

        def rgb_kernel(name, green):
            data = optional(f'{name}_rb.f64')
            if data is None:
                return None
            if len(data) != 2 * n * n * 8:
                raise ValueError('数据包的 RGB PSF 数组尺寸错误。')
            red, blue = np.frombuffer(data, dtype='<f8').reshape(2, n, n)
            return np.stack((red, green.T, blue))

        green, clear_green = kernel('kernel'), kernel('clear_kernel')
        previews = None
        if (mask_meta := meta.get('mask_preview')) is not None:
            width, height = [float(mask_meta[key]) for key in ('width_um', 'height_um')]
            if not all(np.isfinite(v) and v > 0 for v in (width, height)):
                raise ValueError('数据包的掩膜显示尺寸无效。')
            layers = {}
            for name in MASK_LAYERS:
                png = optional(f'masks/{name}.png')
                if png is None:
                    continue
                if len(png) > MAX_IMAGE or not png.startswith(b'\x89PNG\r\n\x1a\n'):
                    raise ValueError('掩膜预览必须为有效大小的 PNG。')
                with Image.open(BytesIO(png)) as image:
                    if max(image.size) > MAX_SIDE:
                        raise ValueError('数据包的掩膜预览超出支持的采样尺寸。')
                    image.verify()
                layers[name] = png
            previews = MaskPreviews(width, height, bool(mask_meta['tiled']), bool(mask_meta['subsampled']), layers)
        return validate(Case(meta['name'], meta['conditions'], axis, green, clear_green,
                             np.asarray(meta['freq'], dtype=float),
                             np.asarray(meta['mtf'], dtype=float), np.asarray(meta['clear_mtf'], dtype=float),
                             meta['channels'], meta['tails'], archive.read('scene.png'), archive.read('clear.png'),
                             optional('rgb-psf.png'), optional('clear-rgb-psf.png'),
                             rgb_kernel('rgb_kernel', green), rgb_kernel('clear_rgb_kernel', clear_green),
                             mask_previews=previews))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('result_folder', type=Path, help='包含 udc_results.mat 和相关 CSV/PNG 的目录')
    parser.add_argument('--name', help='方案显示名称')
    parser.add_argument('--output', required=True, type=Path, help='输出的 .udc.zip 文件路径')
    args = parser.parse_args()
    folder = args.result_folder.resolve()
    case = read_result(folder, args.name or f'{folder.parent.name} · {folder.name}')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_package(case, args.output)
    print(args.output.resolve())


if __name__ == '__main__':
    main()
