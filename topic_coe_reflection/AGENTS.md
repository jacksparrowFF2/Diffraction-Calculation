# COE 模块维护约定

- 用户运行入口为 `src/COEReflection_main.m`，集中呈现输入、参数、计算、展示及保存流程，并保留必要中文注释。
- `src/` 顶层仅保留主脚本，以及已供外部代码调用的公共接口 `coeReflectionDefaults.m`、`runCOEReflection.m`。未经明确需求，不新增公共接口。
- 新增内部计算、输入读取、绘图与保存辅助函数时，直接创建在 `src/private/`，不要先放在主脚本同级再要求用户手动移动。
- 内部函数及其文件使用 lowerCamelCase（小驼峰），例如 `prepareCoeStack.m`、`composeCoeReflectionRgb.m`；文件名必须与首函数名一致。沿用的两个公共接口名称和主脚本名称保留兼容。
- 保留 MATLAB 特殊目录名 `private`，不要对其单独 `addpath`，不要使用 `genpath` 展开私有目录。模块只需添加 `src/`。
- 移动函数时检查基于 `mfilename('fullpath')` 的输入路径，并同步修改所有调用、测试、说明与打包脚本。通过主脚本和独立复制模块的回归验证真实私有函数解析。
- 新的功能测试优先通过公共入口验证。现有数值算法测试的临时副本由 `tests/private/addCoeAlgorithmFixture.m` 管理，不为测试暴露生产私有函数。
- 实际版图按完整物面导入，不隐式缩放、平铺或改变阳极反射范围；额外 `heightMap` 默认保持 0 nm。
- 用户图片集按 `Input/项目名称/` 管理，主脚本使用 `projectName` 切换；输出对应 `Output/项目名称/`。CIE 数据共用 `Input/Colorimetry/`。不要在日常运行时覆盖用户输入；示例重新生成必须由显式开关控制。

- 正式仿真不提供像素轮廓调整，用户输入阳极、BM、RGB CF 图片即可；不可重新暴露圆角、倾角、菱形尺寸等参数。参考图中实线是 PDL（本模型阳极），虚线是 Mask，不是 BM。预生成示例仅供开发维护。
- 常用参数按项目保存为 `Input/项目名称/config.txt`，在光谱和几何计算前加载；配置覆盖默认值。解析文本不得执行 MATLAB 表达式；未知/重复/类型错误应报错。结果中保留配置记录。

- 使用说明、输入/配置规则、模型与指标解释、迁移交接统一维护在 `docs/README.md`；不再新增平行功能说明。Input 中的简短项目清单和 CIE 数据署名/许可保留在原位。

- 主脚本统一按stack.layerNames/stack.layerLevels构建膜堆，不再提供simple/multilayer选择。每项method为nk有限膜或reflectance不透明有效边界；RGB CF同层互斥。无基底膜用medium出射介质，COE外部默认absorber终端，不能把悬空膜未反射能量直接称为吸收。

- heightMap全零时只运行一次传播、单次色度和角域积分，并使用单状态图表及导出；非零高度保留两状态比较。历史baseline等数组字段可复用作兼容，不得重新计算相同结果。膜厚本身的相位不可因关闭额外heightMap而删除。

- 照明支持plane、有限距离/尺寸point和近轴TEM00 laser；有限源点按强度叠加，勿相干叠加。光谱与空间模型独立，支持D65/自定义SPD/单色谱线，不增加相机预览。入射功率按物面平均辐照度标定。
- 默认IJP采用全白图片，与其它图层同尺寸同坐标，n=1.5、k=0、厚10μm，在Anode/BM之间。透明不代表无界面反射或无相位。厚膜需报告光谱欠采样风险。默认角域±6°，输出801；只裁切显示不等于提高计算分辨率。
