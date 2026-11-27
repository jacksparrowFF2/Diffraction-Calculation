# 透射衍射对比工作台

工作台读取 UDC V2 的 MATLAB 仿真结果，可以检查方案可比性、并排查看掩膜和光学指标、导出 PNG。离线 HTML 提供固定的 Case1–8 对比。MATLAB 求解入口见 [UDC V2 使用说明](../../topic_udc_imagingV2/docs/README.md)。

## 启动工作台

Windows 下双击[启动器](../Start_Transmission_Diffraction_Viewer.cmd)。启动器会查找装有依赖的项目虚拟环境、系统 Python 和常见的 Anaconda/Miniconda 环境，选中的 Python 路径会显示在窗口中。工作台默认从 `8501` 开始找可用端口，并在浏览器打开。关闭启动器的 CMD 窗口即可停止服务并释放端口。

需要 Python 3.10 或以上。首次安装依赖及手动启动时，在本 topic 目录运行以下命令。

```powershell
python -m pip install -r requirements.txt
python -m streamlit run src/streamlit_app.py --server.address 127.0.0.1 --server.headless true --browser.gatherUsageStats false
```

在侧栏「本地结果目录」旁点「浏览…」，通过本机文件夹窗口选择目录，也可直接输入路径。工作台优先发现本 topic 的 `Input/` 中完整的结果，其次读取同级 `topic_udc_imagingV2/Output/`。独立复制本 topic 后，可以选择任意结果目录，也可以直接导入 `.udc.zip`。

## 结果目录和 Case 包

可选择结果根目录，也可直接选中一个包含 `udc_results.mat`、三个指标 CSV（通道指标、方向 MTF、PSF 拖尾）和同曝光、无遮挡两个 PNG 的方案目录。结果根目录支持 `<方案>/<透射模式>/`、`<方案>/object/<透射模式>/` 和 `<方案>/`。

侧栏「场景模式」可以筛选 sensor 或 object。选「全部」并载入同一 Case 的两套结果后，成像采用上排 sensor、下排 object，光学图表每个 Case 显示一套 object 结果。侧栏「成像显示」可切换同曝光与亮度匹配。Case9–12 的 object 模式使用 50×50 cm 目标、50 cm 物距和 4.8×3.6 mm 画幅。

展开「打包并导入本地结果」，勾选方案后点击「导入所选结果」。工作台自动打包并将方案加入当前对比，切换结果目录后已导入方案仍可参与对比。数据包保存在 `Output/case_packages/`，也可以从同一面板下载。工作台根据源文件的修改时间和大小检查更新，源文件未变时复用已有包，已加入的方案不会重复添加。再次打开工作台时，直接选择已有有效数据包的本地方案也会复用该包。源文件更新后，再点一次导入即可更新原方案。

数据包保留线性 RGB PSF、BM/Cath 掩膜及物面入瞳预览。已加入对比的本地旧包缺少物面预览时，工作台可用同一批原始结果自动补齐。单独上传的旧 `.udc.zip` 若缺少预览，可选择原结果目录重新导入。

也可用命令行打包，例如从本 topic 目录运行：

```powershell
python src/case_workspace.py ../topic_udc_imagingV2/Output/case8/nk --name Case8 --output Output/Case8.udc.zip
```

## 生成离线 HTML

```powershell
python src/build_interactive_comparison.py
```

脚本在 `Input/` 或同级 UDC V2 的 `Output/` 中查找完整的 4 mm 入瞳 Case1–8，核对仿真参数、场景和共同无遮挡参考，生成 `Output/udc_comparison.html`。可用 `--results` 指定其它结果根目录，用 `--output` 指定 HTML 保存路径。双击生成的文件即可离线使用。

运行 `python src/build_interactive_comparison.py --dual-scenes` 可生成 Case9–12 的 `Output/udc_case9-12_dual_scene.html`。页内默认显示两排成像，可用场景筛选切换成单排。

## 比较结果

- 方案列头固定在顶部，无遮挡参考在最左侧。可编辑方案名称、开关参考，并从工具栏下载采用当前设置的 PNG。切换到「读图说明」可用 Case7 和 Case8 学习如何结合各个指标。返回时保留已选方案和图表设置。
- BM/Cath 掩膜来自 MAT 中实际参与仿真的覆盖数据。叠加时白色为共同开口、蓝色为仅 BM、橙色为仅 Cath、紫色为两层重叠。图下标注所示的物面范围。
- 「物面 · 有效振幅入瞳」采用与 MATLAB 报告相同的几何入瞳与透射振幅，选用最接近 550 nm 的仿真波长。各列共用物面坐标和 0–1 振幅灰度。MTF 与拖尾曲线旁各有方案图例，名称和线型随当前对比同步更新。新图及图例会随总图一起导出。
- G 通道归一化 PSF 以各方案的窗口能量归一为 1。固定离线页的 0 dB 为全部内置方案及无遮挡的共同峰值，工作台的 0 dB 为单像元占全部窗口能量。绝对能量模式以入射能量为分母。两个模式共用 −70 至 0 dB 的显示范围，中心放大沿用同一色标。色图含 MATLAB 预定义色图、Magma 与 Viridis，可用 [MATLAB 导出脚本](../src/transmission_diffraction_viewer/export_matlab_colormaps.m)重新导出色表。
- RGB PSF 对线性 RGB 能量应用 COE 的对数亮度映射 `[log1p(60 × Y) / log1p(60)]^0.72`，以共同峰值确定曝光，再编码为 sRGB。亮度映射用于显示。只有旧版 RGB PNG 的 Case 包会标为「旧版映射」。
- 工作台按指标列出有差异和缺失的仿真条件。PSF 显示共同接收视场，MTF 显示共同有效频段，拖尾按共同有效半径比较。窗口捕获量可选以入射能量或各自无遮挡捕获量为分母。
