"""Local workspace for comparing current and future UDC simulations."""

from collections import Counter
from hashlib import sha256
from pathlib import Path
from subprocess import CalledProcessError
from zipfile import BadZipFile
from zlib import error as DecompressionError

import streamlit as st

from transmission_diffraction_viewer.case_compatibility import RULES, evaluate
from transmission_diffraction_viewer.case_workspace import (
    MAX_PACKAGE_MIB, default_results_root, discover, read_package, read_result, source_stamp,
)
from transmission_diffraction_viewer.folder_picker import pick_folder
from transmission_diffraction_viewer.local_packages import (
    current_local_package, ensure_local_package, upgrade_pupil_package,
)
from transmission_diffraction_viewer.workspace_figure import workspace_page


READ_ERRORS = (OSError, ValueError, KeyError, TypeError, IndexError, AssertionError, BadZipFile, DecompressionError)
RESULTS_DIR_KEY = 'results_directory'
SELECTED_KEY = 'selected_cases'


@st.cache_data(max_entries=12, show_spinner=False)
def load_local(path, stamp, name):
    return read_result(path, name=name)


@st.cache_data(max_entries=8, show_spinner=False)
def load_upload(data):
    return read_package(data)


@st.cache_data(max_entries=12, show_spinner=False)
def load_local_package(package, stamp):
    return read_package(package.read_bytes())


def read_cached_case(package):
    stat = package.stat()
    return load_local_package(package, (stat.st_mtime_ns, stat.st_size))


def read_local_case(folder, name):
    stamp = source_stamp(folder)
    if package := current_local_package(folder, stamp):
        try:
            return read_cached_case(package)
        except READ_ERRORS:
            pass
    return load_local(folder, stamp, name)


def read_imported_local_case(folder, name, package):
    case = read_cached_case(package)
    if case.pupil_preview is None:
        try:
            if upgrade_pupil_package(folder, name, package):
                return read_cached_case(package)
        except READ_ERRORS:
            pass
    return case


def browse_results():
    try:
        folder = pick_folder(st.session_state[RESULTS_DIR_KEY].strip().strip('"'))
    except (OSError, CalledProcessError):
        st.session_state['results_picker_error'] = '无法打开本机目录选择器，请在路径框中输入目录。'
        return
    if folder:
        st.session_state[RESULTS_DIR_KEY] = folder


def show_local_import(paths, names, local_imports):
    messages = {'created': '已打包并加入对比', 'updated': '已更新数据包和对比结果',
                'reused': '已复用已有数据包并加入对比'}
    with st.sidebar.expander('打包并导入本地结果'):
        pending = st.multiselect('选择要导入的本地结果', list(paths), format_func=names.get)
        if st.button('导入所选结果', disabled=not pending, type='primary', width='stretch'):
            selected = list(st.session_state.get(SELECTED_KEY, ()))
            with st.spinner('正在检查本地结果…'):
                for uid in pending:
                    name = names[uid]
                    try:
                        package, state = ensure_local_package(paths[uid], name)
                        try:
                            read_cached_case(package)
                        except READ_ERRORS:
                            if state != 'reused':
                                raise
                            package, state = ensure_local_package(paths[uid], name, rebuild=True)
                            read_cached_case(package)
                    except READ_ERRORS as error:
                        st.error(f'{name} · 导入失败 · {error}')
                        continue
                    if state == 'reused' and uid in selected and uid in local_imports:
                        st.info(f'{name} · 结果无更新，已在对比中')
                        continue
                    local_imports[uid] = (name, package)
                    if uid not in selected:
                        selected.append(uid)
                    st.success(f'{name} · {messages[state]}')
            st.session_state[SELECTED_KEY] = selected
        for uid in pending:
            if uid in local_imports:
                name, package = local_imports[uid]
                if package.is_file():
                    st.download_button(f'下载 {name} 的 Case 包', data=package.read_bytes,
                                       file_name=f"{name.replace(' · ', '-')}.udc.zip",
                                       mime='application/zip', on_click='ignore', key=f'local-download:{uid}')


def catalog():
    st.sidebar.header('方案来源')
    if RESULTS_DIR_KEY not in st.session_state:
        st.session_state[RESULTS_DIR_KEY] = str(default_results_root())
    path_field, browse_button = st.sidebar.columns([3, 1], vertical_alignment='bottom')
    directory = path_field.text_input('本地结果目录', key=RESULTS_DIR_KEY)
    browse_button.button('浏览…', on_click=browse_results, width='stretch')
    if error := st.session_state.pop('results_picker_error', None):
        st.sidebar.error(error)
    packs = st.sidebar.file_uploader('导入 UDC Case 数据包', type='zip',
                                    accept_multiple_files=True, max_upload_size=MAX_PACKAGE_MIB)
    root = Path(directory.strip().strip('"')).expanduser()
    root_exists = bool(directory.strip()) and root.is_dir()
    if directory.strip() and not root_exists:
        st.sidebar.warning('找不到结果目录。可以继续导入数据包。')
    names, paths, imported = {}, {}, {}
    for folder in discover(root) if root_exists else []:
        uid = str(folder.resolve())
        parts = folder.relative_to(root).parts or (folder.parent.name, folder.name)
        names[uid] = ' · '.join(parts)
        paths[uid] = folder
    local_imports = st.session_state.setdefault('local_result_imports', {})
    if paths:
        show_local_import(paths, names, local_imports)
    elif root_exists:
        st.sidebar.info('此目录没有完整的 UDC V2 结果。')
    for uid, (name, _) in local_imports.items():
        names.setdefault(uid, name)
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
    previous = st.session_state.get(SELECTED_KEY, defaults)
    st.session_state[SELECTED_KEY] = [uid for uid in previous if uid in names]
    selected = st.sidebar.multiselect('参与对比的方案', names, format_func=names.get, key=SELECTED_KEY)
    loaded = []
    with st.sidebar.expander('修改方案显示名称'):
        for uid in selected:
            display = st.text_input(names[uid], key=f'name:{uid}', value=names[uid], max_chars=120).strip()
            try:
                if uid in imported:
                    case = imported[uid]
                elif uid in local_imports:
                    name, package = local_imports[uid]
                    case = read_imported_local_case(Path(uid), name, package)
                else:
                    case = read_local_case(paths[uid], names[uid])
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
    st.set_page_config(page_title='透射衍射对比工作台', layout='wide')
    st.title('透射衍射对比工作台')
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
    st.iframe(page, height=760, alt='所选方案的 BM 和 Cath 掩膜、物面入瞳、PSF、MTF、通量、拖尾与同曝光成像总图')


main()
