# F192Base：实际版图输入

- 主脚本设置 `projectName="F192Base"`，保持 `useProjectConfig=true`。
- 六张 PNG 为 1000×1000、4×4 Pitch，原生250像素/Pitch；Pitch为59.06 μm，当前按64点/Pitch在内存中采样，不改写图片。
- 当前配置：4×4 Pitch、接收面25 cm、±5°、平面波、无独立入瞳，D65 380:10:780 nm。
- anode、CF_R/B/G白色表示材料；BM、IJP黑色表示材料且invert=true。层序为Anode→IJP→BM→CF_R→CF_B→CF_G。
- Materials中保存BM与RGB CF的三列波长/n/k数据；当前config中的nkFile为空，仍使用配置常数，需显式填写相对路径后启用。CAD为 `../Cad/F192/(20261103)F192Base.dwg`。
- 模型、材料启用方式及输入规则见[统一使用说明](../../docs/README.md)。
