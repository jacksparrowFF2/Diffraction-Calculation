from pathlib import Path
import argparse
import csv

import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
from PIL import Image

from comparison_data import PSF_COLORMAP, PSF_FLOOR_DB, green_mtf, read_csv

module = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description="Compare UDC results against one clear pupil.")
parser.add_argument("--apertures", action="store_true", help="Compare Case1 and Case5–8.")
parser.add_argument("--results", type=Path, default=module / "Output",
                    help="Directory containing the MATLAB case results.")
parser.add_argument("--output", type=Path, help="Directory for the generated figure and CSV.")
args = parser.parse_args()
apertures = args.apertures
results = args.results
out_name = ("case1_case5_case6_case7_case8_4mm_500pitch" if apertures
            else "case1_case2_case3_case4_4mm_500pitch")
out = args.output or results / out_name
out.mkdir(parents=True, exist_ok=True)
if apertures:
    cases = ('case1','case5','case6','case7','case8')
    labels = ('Case1 · BM 59 μm','Case5 · BM 59 / 阴极 54 μm','Case6 · BM 54 μm',
              'Case7 · BM 59 μm + 连续阴极','Case8 · Case7 + 50% 衰减')
    colors = ('#255b8a','#b54531','#28805f','#7953a1','#b06b20')
else:
    cases = ('case1','case2','case3','case4')
    labels = ('Case1 · BM','Case2 · BM + 图案阴极','Case3 · BM','Case4 · BM + 连续阴极')
    colors = ('#255b8a','#b54531','#28805f','#7953a1')
channels = 'RGB'
mask_path = 'optics/layerGeometry/layerMasks'

def verify_full_cathode(data, case):
    cath, bm = (data[f'{mask_path}/{name}'] for name in ('Cath','BM'))
    pupil = data['optics/pupil']
    region = data['optics/regionMap']
    np.testing.assert_array_equal(data['optics/regionPatterns'][()].T,
                                  [[1, 0], [1, 1]])
    grid_count = pupil_count = opening_count = 0
    with h5py.File(results/'case1/nk/udc_results.mat') as case1:
        for start in range(0, cath.shape[0], 256):
            rows = slice(start, start+256)
            assert np.all(cath[rows]), f'{case} has uncovered cathode pixels'
            np.testing.assert_array_equal(bm[rows], case1[f'{mask_path}/BM'][rows])
            np.testing.assert_array_equal(region[rows], bm[rows])
            inside = pupil[rows] != 0
            grid_count += inside.size
            pupil_count += np.count_nonzero(inside)
            opening_count += np.count_nonzero(inside & (bm[rows] == 0))
    assert 0 < opening_count < pupil_count
    print(f'{case} cathode coverage: {grid_count:,} grid pixels, '
          f'{pupil_count:,} pupil pixels, {opening_count:,} pupil pixels under BM openings')

def verify_sampled_openings(data, case):
    layer = 'Cath' if case == 'case5' else 'BM'
    with Image.open(module/'Input'/case/f'{layer}.png') as image:
        opening = np.asarray(image.convert('L'))[1::2,1::2].T == 255
    mask = data[f'{mask_path}/{layer}']
    pupil = data['optics/pupil']
    pupil_open = pupil_size = 0
    cols = np.arange(mask.shape[1]) % opening.shape[1]
    with h5py.File(results/'case1/nk/udc_results.mat') as base:
        for start in range(0, mask.shape[0], 256):
            stop = min(start+256, mask.shape[0])
            rows = np.arange(start,stop) % opening.shape[0]
            expected = opening[np.ix_(rows,cols)]
            np.testing.assert_array_equal(mask[start:stop] == 0, expected)
            bm = base[f'{mask_path}/BM'][start:stop]
            assert not np.any(expected & (bm != 0)), '54 um opening extends beyond the Case1 BM'
            if case == 'case5':
                np.testing.assert_array_equal(data[f'{mask_path}/BM'][start:stop],bm)
            inside = pupil[start:stop] != 0
            pupil_open += np.count_nonzero(inside & expected)
            pupil_size += np.count_nonzero(inside)
    print(f'{case} {layer} open in pupil: {pupil_open:,}/{pupil_size:,} ({pupil_open/pupil_size:.6%})')

def verify_scaled_result(data, case):
    reference, scale = ('case4', 1) if case == 'case7' else ('case7', 0.5)
    with h5py.File(results/reference/'nk/udc_results.mat') as base:
        for key in ('optics/kernel','optics/spectralEnergy','optics/quality/hdrPoint/withMask',
                    'optics/channelThroughput','optics/channelCapturedThroughput','images/fixedExposure'):
            np.testing.assert_array_equal(data[key][()],base[key][()]*scale)
        for key in ('optics/psf','optics/mtf','optics/clearKernel',
                    'optics/quality/hdrPoint/clear','images/clear','images/brightnessMatched'):
            np.testing.assert_array_equal(data[key][()],base[key][()])
        lookup = 'optics/transmissionLookup'
        for part in ('real','imag'):
            np.testing.assert_allclose(data[lookup][part],base[lookup][part]*np.sqrt(scale),
                                       rtol=1e-15,atol=0)

records = []
clear_reference = None
for case in cases:
    folder = results/case/'nk'
    channel = read_csv(folder,'udc_channel_metrics.csv')
    tail = read_csv(folder,'udc_psf_tail.csv')
    mtf = read_csv(folder,'udc_mtf_directional_summary.csv')
    profiles = read_csv(folder,'udc_mtf_directional.csv')
    with h5py.File(folder/'udc_results.mat') as data:
        np.testing.assert_array_equal(data['cfg/periods'][()].ravel(), [17,17])
        for key,value in dict(pitch=236.2e-6, samplesPerPeriod=500, pupilRadius=0.002,
                              focalLength=0.003, sensorPitch=1e-6, outputSize=201,
                              **{'quality/hdrDisplayHalfWidthUm':50}).items():
            np.testing.assert_allclose(data[f'cfg/{key}'][()].item(),value,rtol=1e-12,atol=0)
        np.testing.assert_allclose(data['cfg/wavelengths'][()].ravel(),np.arange(380,781,10)*1e-9)
        shared_inputs = {key:data[key][()] for key in
                         ('cfg/spectralWeights','cfg/channelResponse','images/inputLinear')}
        x = data['optics/coordinate'][()].ravel()*1e6
        psf = data['optics/psf'][1].T
        clear_psf = data['optics/quality/clearPsf'][1].T
        region = data['optics/regionMap']
        assert region.shape == (8500,8500)
        assert psf.shape == (201,201)
        assert abs(psf.sum()-1) < 1e-10
        linear_mean = data['images/fixedExposure'][()].mean(axis=(1,2))
        clear_mean = data['images/clear'][()].mean(axis=(1,2))
        if case in ('case4','case7','case8'):
            verify_full_cathode(data,case)
        if case in ('case5','case6'):
            verify_sampled_openings(data,case)
        if case in ('case7','case8'):
            verify_scaled_result(data,case)
    g_tail = {float(row['RadiusUm']):float(row['TailFractionWithinWindow']) for row in tail if row['Channel']=='G'}
    g_clear_tail = {float(row['RadiusUm']):float(row['ClearTailFractionWithinWindow']) for row in tail if row['Channel']=='G'}
    freq, mtf_worst, clear_worst = green_mtf(profiles)
    relative = np.array([float(row['RelativeToClear']) for row in channel])
    absolute = np.array([float(row['CapturedThroughput']) for row in channel])
    total = np.array([float(row['Throughput']) for row in channel])
    worst_band = min(float(row['BandMinMTF']) for row in mtf if row['Channel']=='G')
    worst_300 = np.interp(300,freq,mtf_worst)
    with Image.open(folder/'udc_fixed_exposure.png') as image:
        scene = np.asarray(image.convert('RGB'))
    with Image.open(folder/'udc_clear.png') as image:
        clear_scene = np.asarray(image.convert('RGB'))
    if clear_reference is None:
        clear_reference = dict(x=x,psf=clear_psf,tail=g_clear_tail,scene=clear_scene,
                               inputs=shared_inputs, clear_mtf=clear_worst, freq=freq,
                               absolute=np.array([float(r['ClearCapturedThroughput']) for r in channel]),
                               band=min(float(r['ClearBandMinMTF']) for r in mtf if r['Channel']=='G'))
    else:
        np.testing.assert_array_equal(clear_reference['scene'],clear_scene)
        np.testing.assert_allclose(clear_reference['psf'],clear_psf,rtol=0,atol=1e-12)
        np.testing.assert_array_equal(clear_reference['freq'],freq)
        np.testing.assert_allclose(clear_reference['clear_mtf'],clear_worst,rtol=0,atol=1e-12)
        for key,values in shared_inputs.items():
            np.testing.assert_array_equal(clear_reference['inputs'][key],values)
    records.append(dict(x=x,psf=psf,relative=relative,absolute=absolute,tail=g_tail,
                        clear_tail=g_clear_tail,freq=freq,mtf=mtf_worst,clear=clear_worst,
                        band=worst_band,mtf300=worst_300,scene=scene,total=total,
                        scene_mean_relative=linear_mean/clear_mean))
    print(case,'RGB capture relative:',relative,'tail:',g_tail,
          'clear tail:',g_clear_tail,'min MTF 10–100:',worst_band,'MTF300:',worst_300)

with (out/'summary.csv').open('w',encoding='utf-8',newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Case','RelativeCaptureR','RelativeCaptureG','RelativeCaptureB',
                     'GMinMTF10To100','GWorstDirectionMTFAt300',
                     'GTailFractionOutside10Um','GTailFractionOutside20Um','GTailFractionOutside50Um',
                     'GClearTailFractionOutside20Um','GAbsoluteTailThroughputOutside20Um',
                     'TotalThroughputR','TotalThroughputG','TotalThroughputB',
                     'GPSFCenterPixelFraction',
                     'RelativeLinearSceneMeanR','RelativeLinearSceneMeanG','RelativeLinearSceneMeanB'])
    for case,row in zip(cases,records):
        writer.writerow([case,*row['relative'],row['band'],row['mtf300'],
                         *[row['tail'][r] for r in (10,20,50)],row['clear_tail'][20],
                         row['tail'][20]*row['absolute'][1], *row['total'],
                         row['psf'][100,100], *row['scene_mean_relative']])
    if apertures:
        clear = clear_reference
        writer.writerow(['clear',1,1,1,clear['band'],np.interp(300,clear['freq'],clear['clear_mtf']),
                         *[clear['tail'][r] for r in (10,20,50)],clear['tail'][20],
                         clear['tail'][20]*clear['absolute'][1],1,1,1,clear['psf'][100,100],1,1,1])

font_path = 'C:/Windows/Fonts/msyh.ttc'
font_manager.fontManager.addfont(font_path)
plt.rcParams.update({'font.family':font_manager.FontProperties(fname=font_path).get_name(),
                     'axes.unicode_minus':False,'axes.spines.top':False,'axes.spines.right':False,
                     'font.size':10,'text.color':'#172a42','axes.labelcolor':'#465b70'})
visual_labels = ('无遮挡 · 仅 4 mm 入瞳', *labels)
visual_records = [clear_reference, *records]
n_panels = len(visual_records)
fig = plt.figure(figsize=(4*n_panels,13.5),facecolor='white')
grid = fig.add_gridspec(4,n_panels+1,width_ratios=[1]*n_panels+[0.045],height_ratios=(1,1,0.85,1),
                       left=0.075,right=0.928,bottom=0.035,top=0.9,wspace=0.29,hspace=0.56)
title = 'Case1 与 Case5–8' if apertures else '四种 UDC 方案'
fig.text(0.075,0.973,f'4 mm 入瞳 · 500 点/Pitch · {title}对比',size=23,weight='bold',va='top')
fig.text(0.075,0.935,'G 通道 PSF、RGB 进光量、四方向 MTF 与同曝光成像 · 3 mm 焦距 · D65',
         size=10.8,color='#4b637b',va='top')
peak = max(row['psf'].max() for row in visual_records)
for i,(label,row) in enumerate(zip(visual_labels,visual_records)):
    for plot_row,limit in enumerate((100,25)):
        ax = fig.add_subplot(grid[plot_row,i])
        inside = np.abs(row['x'])<=limit+1e-6
        crop = row['psf'][np.ix_(inside,inside)]
        db = 10*np.log10(np.maximum(crop/peak,10**(PSF_FLOOR_DB/10)))
        edge = limit+0.5
        im = ax.imshow(db,extent=[-edge,edge]*2,origin='lower',cmap=PSF_COLORMAP,vmin=PSF_FLOOR_DB,vmax=0,
                       interpolation='nearest')
        ax.set(xlabel='传感器 x (μm)',ylabel='传感器 y (μm)')
        ax.set_title(label if plot_row==0 else f'中心放大 · 20 μm 外能量 {row["tail"][20]:.1%}',
                     loc='left',weight='bold',fontsize=10)
fig.colorbar(im,cax=fig.add_subplot(grid[:2,-1]),label='归一化 PSF，相对共同峰值 (dB)')

ax = fig.add_subplot(grid[2,:3])
for case,label,color,row in zip(cases,labels,colors,records):
    ax.plot(row['freq'],row['mtf'],color=color,lw=1.8,label=label,
            ls='--' if case == 'case8' else '-')
ax.plot(records[0]['freq'],records[0]['clear'],ls='--',color='#727e8d',label='无遮挡')
ax.set(xlabel='空间频率 (cycles/mm)',ylabel='四方向最低 MTF',xlim=(0,500),ylim=(0,1.03))
ax.grid(alpha=0.18)
ax.set_title('G 通道 MTF · 逐频率取四方向最小值',loc='left',weight='bold',pad=8)
ax.legend(loc="upper right",frameon=False,fontsize=8,ncol=2)

ax = fig.add_subplot(grid[2,3:n_panels])
xp = np.arange(len(cases))
rgb = ('#b64549','#368153','#365dab')
for i,(channel,color) in enumerate(zip(channels,rgb)):
    ax.bar(xp+(i-1)*0.22,[100*row['relative'][i] for row in records],
           width=0.20,color=color,label=channel)
ax.set(ylabel='相对无遮挡 (%)',xticks=xp,xticklabels=cases)
ax.set_title('窗口捕获能量',loc='left',weight='bold',pad=8)
ax.grid(axis='y',alpha=0.18)
ax.legend(frameon=False,fontsize=8,ncol=3)

for i,(label,row) in enumerate(zip(visual_labels,visual_records)):
    ax = fig.add_subplot(grid[3,i])
    ax.imshow(row['scene'])
    ax.set_axis_off()
    ax.set_title(f'同曝光成像 · {label}',loc='left',weight='bold',fontsize=10)
with_reference = out/'comparison.png'
fig.savefig(with_reference,dpi=155,facecolor='white')
print(with_reference.resolve())
