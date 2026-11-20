# UDC 衍射仿真与评价文献

## 版图、PSF 与方向性 MTF

Anqi Yang、Aswin C. Sankaranarayanan，*Designing Display Pixel Layouts for Under-Panel Cameras*，IEEE TPAMI，DOI 10.1109/TPAMI.2021.3075978。第 3 节说明共面入瞳模型、周期版图衍射和低 MTF 区域。图 4 的四列依次是有效孔径图样、色调映射的 RGB PSF、对数显示的 G 通道 PSF，以及采用孔径自相关尺度的 PSF 傅里叶结果。图中的亮点和光刺是计算 PSF 的结构。V2 的四联图采用相同的阅读顺序，最后一列直接以传感器空间频率 cycles/mm 表示二维 MTF。 [论文全文](https://imagesci.ece.cmu.edu/files/paper/2021/UDC_PAMI21.pdf) · [补充材料](https://imagesci.ece.cmu.edu/files/paper/2021/UDC_PAMI21-Supp.pdf)

## 透过率、参考成像与共面模型

Yuqian Zhou、David Ren、Neil Emerton、Sehoon Lim、Timothy Large，*Image Restoration for Under-Display Camera*，CVPR 2021。第 3 节把信号衰减与衍射核分开处理，并在近入瞳假设下由透射图求 PSF。V2 用同一入瞳的无遮挡核作为参照分别评价能量损失和 PSF 形状。论文补充材料还讨论了屏幕与镜头间存在间距时的可见条纹。 [论文与实验材料](https://arxiv.org/abs/2003.04857) · [补充材料](https://openaccess.thecvf.com/content/CVPR2021/supplemental/Zhou_Image_Restoration_for_CVPR_2021_supplemental.pdf)

## HDR 点光源和长尾

Ruicheng Feng 等，*Removing Diffraction Image Artifacts in Under-Display Camera via Dynamic Skip Connection Networks*，CVPR 2021。作者使用手机原型测得真实 PSF，并用 HDR 场景生成衍射降质图像。V2 的 HDR 点光源预览对应暗场中单个亮点的理想光学响应，用同一曝光与显示映射比较峰值、彩色光刺与拖尾。绝对尾部透过率和窗口内能量占比是 V2 的对照指标。 [项目页、论文和数据](https://jnjaby.github.io/projects/UDC/)

## 膜层相位与版图优化

Anqi Yang、Eunhee Kang、Hyong-Euk Lee、Aswin C. Sankaranarayanan，*Designing Phase Masks for Under-Display Cameras*，ICCV 2023。论文从入瞳复场出发讨论相位与颜色，以及用 PSF 可逆性和通光量评价设计，适合作为 V2 中空间变化相位和薄膜方案比较的参考。 [论文与实验图](https://anqiyang.github.io/ICCV23_UDC/index.html)
