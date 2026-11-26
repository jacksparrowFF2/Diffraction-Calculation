"""Local workspace for comparing current and future UDC simulations."""

from collections import Counter
from hashlib import sha256
from io import BytesIO
from pathlib import Path
from zipfile import BadZipFile

import streamlit as st

from case_compatibility import RULES, evaluate
from case_workspace import (MAX_PACKAGE_MIB, OUTPUT, discover, read_package, read_result,
                            source_stamp, write_package)
from workspace_figure import workspace_page


READ_ERRORS = (OSError, ValueError, KeyError, TypeError, IndexError, AssertionError, BadZipFile)


@st.cache_data(max_entries=12, show_spinner=False)
def load_local(path, stamp, name):
    return read_result(path, name=name)


@st.cache_data(max_entries=8, show_spinner=False)
def load_upload(data):
    return read_package(data)


def show_package_download(root, paths, names):
    with st.sidebar.expander('打包 MATLAB 仿真结果'):
        with st.container(height=200, border=False):
            key = st.radio('选择要打包的本地结果', list(paths), index=None, format_func=names.get)
        if key is None:
            return
        folder = paths[key]
        try:
            case = load_local(folder, source_stamp(folder), names[key])
        except READ_ERRORS as error:
            st.error(f'无法打包 {names[key]}。{error}')
            return

        def make_zip():
            with BytesIO() as buffer:
                write_package(case, buffer)
                return buffer.getvalue()

        filename = '-'.join(folder.relative_to(root).parts) + '.udc.zip'
        st.download_button('生成并下载 Case 包', data=make_zip, file_name=filename,
                           mime='application/zip', on_click='ignore',
                           icon=':material/download:', width='stretch')


def catalog():
    st.sidebar.header('方案来源')
    directory = st.sidebar.text_input('本地结果目录', str(OUTPUT))
    packs = st.sidebar.file_uploader('导入 UDC Case 数据包', type='zip',
                                    accept_multiple_files=True, max_upload_size=MAX_PACKAGE_MIB)
    root = Path(directory.strip().strip('"')).expanduser()
    root_exists = bool(directory.strip()) and root.is_dir()
    if directory.strip() and not root_exists:
        st.sidebar.warning('找不到结果目录。可以继续导入数据包。')
    names, paths, imported = {}, {}, {}
    for folder in discover(root) if root_exists else []:
        relative = folder.relative_to(root)
        uid = str(folder.resolve())
        names[uid] = ' · '.join(relative.parts)
        paths[uid] = folder
    if paths:
        show_package_download(root, paths, names)
    for pack in packs:
        data = pack.getvalue()
        uid = 'pack:' + sha256(data).hexdigest()
        try:
            case = load_upload(data)
        except READ_ERRORS as error:
            st.sidebar.error(f'{pack.name}：{error}')
            continue
        imported[uid] = case
        names[uid] = f'{case.name} · 导入 ({uid[5:11]})'
    if not names:
        st.info('请选择包含仿真输出的本地目录，或导入 Case 数据包。')
        st.stop()
    preferred = [key for key, name in names.items() if name in ('case7 · nk', 'case8 · nk')]
    defaults = preferred or list(names)[:2]
    selected = st.sidebar.multiselect('参与对比的方案', names, default=defaults,
                                     format_func=names.get, key='selected_cases')
    loaded = []
    with st.sidebar.expander('修改方案显示名称'):
        for uid in selected:
            display = st.text_input(names[uid], key=f'name:{uid}', value=names[uid], max_chars=120).strip()
            try:
                if uid in imported:
                    case = imported[uid]
                else:
                    folder = paths[uid]
                    case = load_local(folder, source_stamp(folder), names[uid])
            except READ_ERRORS as error:
                st.error(f'{names[uid]}：{error}')
                continue
            case.name = display or names[uid]
            loaded.append(case)
    return loaded


def show_compatibility(cases):
    st.subheader('可比性判定')
    rows = []
    for metric in RULES:
        status, differences, missing = evaluate(cases, metric)
        rows.append({'指标': metric, '判定': status, '差异项': '、'.join(differences) or '—',
                     '缺失项': '、'.join(missing) or '—'})
    counts = Counter(row['判定'] for row in rows)
    st.caption(' · '.join(f'{count} 项{status}' for status, count in counts.items()))
    with st.expander('查看逐项判定和各方案的仿真条件'):
        st.dataframe(rows, hide_index=True, alt='各指标的可比性、不同设置及缺失字段')
        st.dataframe(
            [{'方案': case.name, '入瞳直径 (mm)': case.conditions.get('pupil_mm'),
              '焦距 (mm)': case.conditions.get('focus_mm'), '像元 (μm)': case.conditions.get('pixel_um'),
              '波长范围 (nm) / 点数': str(case.conditions.get('wavelength_range') or '未知'),
              '接收窗口 x/y (μm) / 点数': str(case.conditions.get('window_um') or '未知'),
              '物面步长 (μm)': case.conditions.get('object_step_um'),
              '接收积分': case.conditions.get('receiver'),
              '无遮挡基准 ID': (case.conditions.get('clear') or '未知')[:12],
              '场景 ID': (case.conditions.get('scene') or '未知')[:12]} for case in cases],
            hide_index=True, alt='各方案用于比较的仿真参数')


def main():
    st.set_page_config(page_title='UDC 仿真对比', layout='wide')
    st.title('UDC 仿真对比')
    st.caption('选择本地仿真结果或导入 Case 数据包，先检查各指标的比较条件。')
    cases = catalog()
    if not cases:
        st.info('请选择至少一个有效方案。')
        return
    show_clear = st.sidebar.checkbox('在图中显示无遮挡参考', value=True)
    if len(cases) > 1:
        show_compatibility(cases)
    try:
        page = workspace_page(cases, show_clear)
    except ValueError as error:
        st.error(str(error))
        return
    st.iframe(page, height=760, alt='所选方案的 BM 和 Cath 掩膜、PSF、MTF、通量、拖尾与同曝光成像总图')


main()
