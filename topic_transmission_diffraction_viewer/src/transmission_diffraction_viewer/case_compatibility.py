"""Compare measurement conditions without treating design changes as input errors."""

from itertools import combinations
import math


FIELDS = {
    'spectrum': '波长及 RGB 光谱权重',
    'pupil_mm': '入瞳直径',
    'focus_mm': '焦距',
    'image_distance_mm': '像距',
    'pixel_um': '像元尺寸',
    'object_step_um': '物面采样间距',
    'periods': '物面阵列周期数',
    'receiver': '接收积分算法',
    'receiver_samples': '接收积分采样设置',
    'grid': '接收网格 / 窗口',
    'energy_basis': '入射能量归一化',
    'clear': '无遮挡基准',
    'scene': '输入场景',
    'scene_mode': '场景模式',
    'scene_shape': '图像尺寸',
    'exposure': '曝光定义',
    'mapping': '预览映射',
}
OPTICS = ('pupil_mm', 'focus_mm', 'image_distance_mm', 'spectrum', 'pixel_um', 'object_step_um',
          'periods', 'receiver', 'receiver_samples', 'grid')
RULES = {
    '归一化 PSF': OPTICS,
    '绝对能量 PSF': (*OPTICS, 'energy_basis'),
    'MTF': OPTICS,
    'RGB 通量': (*OPTICS, 'energy_basis'),
    '相对无遮挡通量': (*OPTICS, 'energy_basis', 'clear'),
    'G 通道拖尾': (*OPTICS, 'energy_basis'),
    '同曝光图': (*OPTICS, 'energy_basis', 'scene_mode', 'scene', 'scene_shape', 'exposure', 'mapping'),
}


def same_value(a, b):
    if isinstance(a, (float, int)) and isinstance(b, (float, int)):
        return math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-9)
    return a == b


def evaluate(cases, metric):
    differences = []
    missing = []
    for key in RULES[metric]:
        vals = [case.conditions.get(key) for case in cases]
        known = [value for value in vals if value is not None]
        if len(known) != len(vals):
            missing.append(FIELDS[key])
        if any(not same_value(a, b) for a, b in combinations(known, 2)):
            differences.append(FIELDS[key])
    if metric == 'MTF' and common_band(cases) is None:
        return '无共同频段', differences, missing
    if metric == 'G 通道拖尾' and not common_tail_radii(cases):
        return '无共同有效半径', differences, missing
    state = ('可比性待确认' if missing else '需注明条件' if differences else '可直接比较')
    return state, differences, missing


def common_band(cases):
    low = max(case.freq[0] for case in cases)
    high = min(case.freq[-1] for case in cases)
    return (low, high) if low < high else None


def common_tail_radii(cases):
    return sorted(set.intersection(*({row[0] for row in case.tails} for case in cases)))
