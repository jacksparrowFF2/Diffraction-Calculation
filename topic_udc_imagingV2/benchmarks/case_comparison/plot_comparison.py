from pathlib import Path
import csv

import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
from PIL import Image

module = Path(__file__).resolve().parents[2]
out = module / "Output" / "case1_case2_case3_4mm_500pitch"
out.mkdir(parents=True, exist_ok=True)
cases = ('case1','case2','case3')
labels = ('Case1 · BM','Case2 · BM + Cath','Case3 · BM')
colors = ('#255b8a','#b54531','#28805f')
channels = 'RGB'
directions = ('Horizontal','Vertical','Diagonal45','Diagonal135')

def read_csv(folder, name):
    with (folder/name).open(encoding='utf-8') as stream:
        return list(csv.DictReader(stream))

records = []
clear_reference = None
for case in cases:
    folder = out.parent/case/'nk'
    channel = read_csv(folder,'udc_channel_metrics.csv')
    tail = read_csv(folder,'udc_psf_tail.csv')
    mtf = read_csv(folder,'udc_mtf_directional_summary.csv')
    profiles = read_csv(folder,'udc_mtf_directional.csv')
    with h5py.File(folder/'udc_results.mat') as data:
        np.testing.assert_array_equal(data['cfg/periods'][()].ravel(), [17,17])
        np.testing.assert_allclose(data['cfg/pupilRadius'][()].item(),0.002)
        x = data['optics/coordinate'][()].ravel()*1e6
        psf = data['optics/psf'][1].T
        clear_psf = data['optics/quality/clearPsf'][1].T
        region = data['optics/regionMap']
        assert region.shape == (8500,8500)
        assert psf.shape == (201,201)
        assert abs(psf.sum()-1) < 1e-10
    g_tail = {float(row['RadiusUm']):float(row['TailFractionWithinWindow']) for row in tail if row['Channel']=='G'}
    g_clear_tail = {float(row['RadiusUm']):float(row['ClearTailFractionWithinWindow']) for row in tail if row['Channel']=='G'}
    paths = [[r for r in profiles if r['Channel']=='G' and r['Direction']==direction] for direction in directions]
    freq = np.array([float(r['CyclesPerMm']) for r in paths[0]])
    mtf_worst = np.min([[float(r['MTF']) for r in path] for path in paths],axis=0)
    clear_worst = np.min([[float(r['ClearMTF']) for r in path] for path in paths],axis=0)
    relative = np.array([float(row['RelativeToClear']) for row in channel])
    absolute = np.array([float(row['CapturedThroughput']) for row in channel])
    worst_band = min(float(row['BandMinMTF']) for row in mtf if row['Channel']=='G')
    worst_300 = np.interp(300,freq,mtf_worst)
    with Image.open(folder/'udc_fixed_exposure.png') as image:
        scene = np.asarray(image.convert('RGB'))
    with Image.open(folder/'udc_clear.png') as image:
        clear_scene = np.asarray(image.convert('RGB'))
    if clear_reference is None:
        clear_reference = dict(x=x,psf=clear_psf,tail=g_clear_tail,scene=clear_scene)
    else:
        np.testing.assert_array_equal(clear_reference['scene'],clear_scene)
        np.testing.assert_allclose(clear_reference['psf'],clear_psf,rtol=0,atol=1e-12)
    records.append(dict(x=x,psf=psf,relative=relative,absolute=absolute,tail=g_tail,
                        clear_tail=g_clear_tail,freq=freq,mtf=mtf_worst,clear=clear_worst,
                        band=worst_band,mtf300=worst_300,scene=scene))
    print(case,'RGB capture relative:',relative,'tail:',g_tail,
          'clear tail:',g_clear_tail,'min MTF 10–100:',worst_band,'MTF300:',worst_300)

with (out/'summary.csv').open('w',encoding='utf-8',newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Case','RelativeCaptureR','RelativeCaptureG','RelativeCaptureB',
                     'GMinMTF10To100','GWorstDirectionMTFAt300',
                     'GTailFractionOutside10Um','GTailFractionOutside20Um','GTailFractionOutside50Um',
                     'GClearTailFractionOutside20Um','GAbsoluteTailThroughputOutside20Um'])
    for case,row in zip(cases,records):
        writer.writerow([case,*row['relative'],row['band'],row['mtf300'],
                         *[row['tail'][r] for r in (10,20,50)],row['clear_tail'][20],
                         row['tail'][20]*row['absolute'][1]])

font_path = 'C:/Windows/Fonts/msyh.ttc'
font_manager.fontManager.addfont(font_path)
plt.rcParams.update({'font.family':font_manager.FontProperties(fname=font_path).get_name(),
                     'axes.unicode_minus':False,'axes.spines.top':False,'axes.spines.right':False,
                     'font.size':10,'text.color':'#172a42','axes.labelcolor':'#465b70'})
visual_labels = ('无遮挡 · 仅 4 mm 入瞳', *labels)
visual_records = [clear_reference, *records]
fig = plt.figure(figsize=(16.8,13.5),facecolor='white')
grid = fig.add_gridspec(4,5,width_ratios=(1,1,1,1,0.045),height_ratios=(1,1,0.85,1),
                       left=0.075,right=0.928,bottom=0.035,top=0.9,wspace=0.29,hspace=0.56)
fig.text(0.075,0.973,'4 mm 入瞳 · 500 点/Pitch 三种 UDC 图样对比',size=23,weight='bold',va='top')
fig.text(0.075,0.935,'G 通道 PSF、RGB 进光量、四方向 MTF 与同曝光成像 · 3 mm 焦距 · D65',
         size=10.8,color='#4b637b',va='top')
peak = max(row['psf'].max() for row in visual_records)
for i,(label,row) in enumerate(zip(visual_labels,visual_records)):
    for plot_row,limit in enumerate((100,25)):
        ax = fig.add_subplot(grid[plot_row,i])
        inside = np.abs(row['x'])<=limit+1e-6
        crop = row['psf'][np.ix_(inside,inside)]
        db = 10*np.log10(np.maximum(crop/peak,1e-7))
        edge = limit+0.5
        im = ax.imshow(db,extent=[-edge,edge]*2,origin='lower',cmap='magma',vmin=-70,vmax=0,
                       interpolation='nearest')
        ax.set(xlabel='传感器 x (μm)',ylabel='传感器 y (μm)')
        ax.set_title(label if plot_row==0 else f'中心放大 · 20 μm 外能量 {row["tail"][20]:.1%}',
                     loc='left',weight='bold',fontsize=10)
fig.colorbar(im,cax=fig.add_subplot(grid[:2,4]),label='归一化 PSF，相对共同峰值 (dB)')

ax = fig.add_subplot(grid[2,:3])
for label,color,row in zip(labels,colors,records):
    ax.plot(row['freq'],row['mtf'],color=color,lw=1.8,label=label)
ax.plot(records[0]['freq'],records[0]['clear'],ls='--',color='#727e8d',label='无遮挡')
ax.set(xlabel='空间频率 (cycles/mm)',ylabel='四方向最低 MTF',xlim=(0,500),ylim=(0,1.03))
ax.grid(alpha=0.18)
ax.set_title('G 通道 MTF · 逐频率取四方向最小值',loc='left',weight='bold',pad=8)
ax.legend(loc="upper right",frameon=False,fontsize=8,ncol=2)

ax = fig.add_subplot(grid[2,3])
xp = np.arange(3)
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
