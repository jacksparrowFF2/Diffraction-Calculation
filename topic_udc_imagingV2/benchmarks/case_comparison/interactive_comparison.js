const {records, palettes, floorDb, workspace, rgbDisplay} = JSON.parse(document.getElementById('data').textContent);
const byId = id => document.getElementById(id);
const figure = byId('figure');
const empty = byId('empty');
const status = byId('status');
const exportButton = byId('export');
const psfMode = byId('psf-mode');
const psfPalette = byId('psf-palette');
const tailMode = byId('tail-mode');
const fluxMode = byId('flux-mode');
const tailRadius = byId('tail-radius');
const radii = workspace?.tailRadii ?? [10, 20, 50];
const ink = '#172638';
const muted = '#516277';
const gridColor = '#e0e8ee';
const rgb = ['#bd4956', '#31805b', '#416bb5'];
const radians = Math.PI * 2;
const scale = 2;
let ready = false;

const percent = value => Number.isFinite(value) ? `${(100 * value).toFixed(2)}%` : '—';
const number = (value, digits) => Number.isFinite(value) ? value.toFixed(digits) : '—';
const nameOf = record => record.nameInput?.value.trim() || record.name;
const isReference = record => record.reference || record.id === 'clear';
const selected = () => workspace ? records : records.filter(record => record.checkbox.checked);

function addControls(record) {
  const row = document.createElement('div');
  row.className = 'control-row';
  const checkbox = document.createElement('input');
  checkbox.type = 'checkbox';
  checkbox.id = `select-${record.id}`;
  checkbox.checked = record.selected;
  checkbox.addEventListener('change', update);
  const fields = document.createElement('div');
  const label = document.createElement('label');
  label.htmlFor = checkbox.id;
  const idLabel = record.id === 'clear' ? '无遮挡参考' : record.id.replace('case', 'Case');
  const swatch = document.createElement('span');
  swatch.className = 'swatch';
  swatch.style.backgroundColor = record.color;
  label.append(swatch, idLabel);
  const nameInput = document.createElement('input');
  nameInput.type = 'text';
  nameInput.maxLength = 48;
  nameInput.value = record.name;
  nameInput.setAttribute('aria-label', `${idLabel}显示名称`);
  nameInput.addEventListener('input', update);
  fields.append(label, nameInput);
  row.append(checkbox, fields);
  byId(record.id === 'clear' ? 'reference-control' : 'case-controls').append(row);
  Object.assign(record, {checkbox, nameInput});
}

function setAll(checked) {
  for (const record of records) {
    if (record.id !== 'clear') record.checkbox.checked = checked;
  }
  update();
}

function reset() {
  for (const record of records) {
    record.checkbox.checked = record.selected;
    record.nameInput.value = record.name;
  }
  psfMode.value = 'normalized';
  psfPalette.value = 'magma';
  tailMode.value = 'tailFraction';
  update();
}

function loadImage(url) {
  return new Promise((resolve, reject) => {
    const image = new Image();
    image.onload = () => resolve(image);
    image.onerror = reject;
    image.src = url;
  });
}

async function loadPsfRaster(url) {
  const image = await loadImage(url);
  const canvas = document.createElement('canvas');
  canvas.width = image.naturalWidth;
  canvas.height = image.naturalHeight;
  const ctx = canvas.getContext('2d', {willReadFrequently: true});
  ctx.drawImage(image, 0, 0);
  const pixels = ctx.getImageData(0, 0, canvas.width, canvas.height);
  const levels = new Uint8Array(canvas.width * canvas.height);
  for (let i = 0; i < levels.length; i++) levels[i] = pixels.data[4 * i];
  return {canvas, ctx, pixels, levels, paletteName: null};
}

function colorPsf(raster) {
  const name = psfPalette.value;
  if (raster.paletteName === name) return raster.canvas;
  const palette = palettes[name];
  const {data} = raster.pixels;
  for (let i = 0; i < raster.levels.length; i++) {
    const [r, g, b] = palette[raster.levels[i]];
    const offset = 4 * i;
    data[offset] = r;
    data[offset + 1] = g;
    data[offset + 2] = b;
  }
  raster.ctx.putImageData(raster.pixels, 0, 0);
  raster.paletteName = name;
  return raster.canvas;
}

function label(ctx, value, x, y, {size = 12, weight = 400, color = ink, align = 'left', max} = {}) {
  ctx.font = `${weight} ${size}px system-ui, "Microsoft YaHei", sans-serif`;
  ctx.textBaseline = 'alphabetic';
  ctx.textAlign = align;
  ctx.fillStyle = color;
  if (max) ctx.fillText(value, x, y, max);
  else ctx.fillText(value, x, y);
}

function sectionTitle(ctx, title, note, y) {
  label(ctx, title, 28, y, {size: 17, weight: 680});
  label(ctx, note, 28, y + 18, {color: muted});
}

function drawLegend(ctx, cases, width) {
  const cols = Math.min(cases.length, Math.max(1, Math.floor((width - 56) / 225)));
  const cell = (width - 56) / cols;
  cases.forEach((record, i) => {
    const x = 28 + (i % cols) * cell;
    const y = 91 + Math.floor(i / cols) * 24;
    lineStyle(ctx, record);
    ctx.beginPath();
    ctx.moveTo(x, y - 4);
    ctx.lineTo(x + 21, y - 4);
    ctx.stroke();
    ctx.setLineDash([]);
    label(ctx, nameOf(record), x + 29, y, {size: 12, max: cell - 36});
  });
  return 113 + (Math.ceil(cases.length / cols) - 1) * 24;
}

function lineStyle(ctx, record) {
  ctx.strokeStyle = record.color;
  ctx.lineWidth = 2.2;
  ctx.setLineDash(isReference(record) ? [3, 5] : record.dash ?? (record.id === 'case8' ? [7, 5] : []));
}

function drawHeatScale(ctx, x, y, side) {
  const grad = ctx.createLinearGradient(0, y + side, 0, y);
  palettes[psfPalette.value].forEach(([r, g, b], i) => grad.addColorStop(i / 255, `rgb(${r},${g},${b})`));
  ctx.fillStyle = grad;
  ctx.fillRect(x, y, 12, side);
  for (const db of [0, -20, -40, -60, floorDb]) {
    const top = y + side * db / floorDb;
    label(ctx, String(db), x + 17, top + 4, {size: 10, color: muted});
  }
  label(ctx, 'dB', x + 4, y - 7, {size: 10, color: muted});
}

function drawPsf(ctx, cases, width, top, radius, side, mode) {
  const rgbPsf = mode === 'rgb';
  const left = 70;
  const cell = (width - left - 64) / cases.length;
  cases.forEach((record, i) => {
    const center = left + cell * (i + .5);
    const x = center - side / 2;
    const titleY = top + 40;
    const y = top + 50;
    label(ctx, nameOf(record), center, titleY, {align: 'center', weight: 660, max: cell - 10});
    ctx.imageSmoothingEnabled = false;
    const image = rgbPsf ? record.rgbPsfImage : colorPsf(record.psfRasters[mode]);
    const extent = rgbPsf ? record.rgbExtent : record.psfExtent;
    if (!image) {
      ctx.fillStyle = '#f2f6f8';
      ctx.fillRect(x, y, side, side);
      label(ctx, '缺少 RGB PSF 数据', center, y + side / 2, {align: 'center', color: muted});
    }
    const count = 2 * radius + 1;
    if (image && extent) {
      const [lo, hi] = extent;
      const widthPerUm = image.width / (hi - lo);
      const heightPerUm = image.height / (hi - lo);
      ctx.drawImage(image, (-radius - lo) * widthPerUm, (hi - radius) * heightPerUm,
                    2 * radius * widthPerUm, 2 * radius * heightPerUm, x, y, side, side);
    } else if (image) {
      ctx.drawImage(image, 100 - radius, 100 - radius, count, count, x, y, side, side);
    }
    for (const tick of [-radius, 0, radius]) {
      const frac = tick / (extent ? 2 * radius : count);
      label(ctx, String(tick), center + side * frac, y + side + 14,
            {size: 10, color: muted, align: 'center'});
      if (i === 0) {
        label(ctx, String(tick), x - 5, y + side * (.5 - frac) + 3,
              {size: 10, color: muted, align: 'right'});
      }
    }
    label(ctx, 'x (μm)', center, y + side + 29, {size: 10, color: muted, align: 'center'});
    if (i === 0) label(ctx, 'y (μm)', x - 4, y - 7, {size: 10, color: muted, align: 'right'});
    if (rgbPsf && record.rgbLegacy) {
      label(ctx, '旧版映射 · 曝光尺度独立', center, y + side + 43, {size: 10, color: muted, align: 'center'});
    }
  });
  if (!rgbPsf) drawHeatScale(ctx, width - 52, top + 50, side);
}

function chartAxes(ctx, x, y, width, height, ticks, maximum, format) {
  ctx.lineWidth = 1;
  ctx.strokeStyle = gridColor;
  ticks.forEach(value => {
    const top = y + height * (1 - value / maximum);
    ctx.beginPath();
    ctx.moveTo(x, top);
    ctx.lineTo(x + width, top);
    ctx.stroke();
    label(ctx, format(value), x - 7, top + 3, {size: 10, color: muted, align: 'right'});
  });
}

function drawCurves(ctx, cases, box, axis, values) {
  const {x, y, width, height} = box;
  const min = axis.xMin ?? 0;
  const toX = value => x + (value - min) / (axis.xMax - min) * width;
  chartAxes(ctx, x, y, width, height, axis.yTicks, axis.yMax, axis.format);
  axis.xTicks.forEach(value => {
    const left = toX(value);
    label(ctx, String(value), left, y + height + 16, {size: 10, color: muted, align: 'center'});
  });
  cases.forEach(record => {
    const [xs, ys] = values(record);
    const points = xs.map((value, i) => [toX(value),
      y + height * (1 - ys[i] / axis.yMax)]);
    lineStyle(ctx, record);
    ctx.beginPath();
    points.forEach(([px, py], i) => i ? ctx.lineTo(px, py) : ctx.moveTo(px, py));
    ctx.stroke();
    ctx.setLineDash([]);
    if (axis.markers) {
      ctx.fillStyle = record.color;
      points.forEach(([px, py]) => {
        ctx.beginPath();
        ctx.arc(px, py, 3, 0, radians);
        ctx.fill();
      });
    }
  });
}

function drawThroughput(ctx, cases, x, y, width, height) {
  const field = fluxMode.value;
  const max = Math.max(1, ...cases.flatMap(record => record[field].filter(Number.isFinite)));
  chartAxes(ctx, x, y, width, height, [0, .25, .5, .75, 1].map(t => t * max), max,
            n => `${Math.round(100 * n)}%`);
  const cell = width / cases.length;
  const bar = Math.min(22, (cell - 7) / 4);
  cases.forEach((record, index) => {
    const center = x + (index + .5) * cell;
    record[field].forEach((value, channel) => {
      if (!Number.isFinite(value)) return;
      ctx.fillStyle = rgb[channel];
      ctx.fillRect(center + (channel - 1.5) * bar, y + (1 - value / max) * height,
                   bar - 2, value / max * height);
    });
    ctx.fillStyle = record.color;
    ctx.beginPath();
    ctx.arc(center, y + height + 10, 3.5, 0, radians);
    ctx.fill();
    const name = nameOf(record);
    const style = {size: 11, align: 'center', max: cell - 6};
    ctx.font = '400 11px system-ui, "Microsoft YaHei", sans-serif';
    if (ctx.measureText(name).width > style.max) {
      const chars = Array.from(name);
      const middle = Math.ceil(chars.length / 2);
      label(ctx, chars.slice(0, middle).join('').trimEnd(), center, y + height + 25, style);
      label(ctx, chars.slice(middle).join('').trimStart(), center, y + height + 38, style);
    } else {
      label(ctx, name, center, y + height + 28, style);
    }
  });
  rgb.forEach((color, channel) => {
    const left = x + width - 83 + channel * 28;
    ctx.fillStyle = color;
    ctx.fillRect(left, y - 16, 8, 8);
    label(ctx, 'RGB'[channel], left + 11, y - 8, {size: 10, color: muted});
  });
}

function drawScenes(ctx, cases, width, top, side) {
  const left = 70;
  const cell = (width - left - 64) / cases.length;
  ctx.imageSmoothingEnabled = true;
  ctx.imageSmoothingQuality = 'high';
  cases.forEach((record, i) => {
    const center = left + (i + .5) * cell;
    label(ctx, nameOf(record), center, top + 41, {align: 'center', weight: 660, max: cell - 10});
    const image = record.sceneImage;
    const fit = side / Math.max(image.naturalWidth, image.naturalHeight);
    const w = image.naturalWidth * fit;
    const h = image.naturalHeight * fit;
    ctx.drawImage(image, center - w / 2, top + 52 + (side - h) / 2, w, h);
    label(ctx, `G 均值 ${percent(record.mean[1])}`, center, top + side + 69,
          {size: 11, color: muted, align: 'center'});
  });
}

function drawFigure(cases) {
  const width = Math.max(900, 138 * cases.length + 152, figure.parentElement.clientWidth);
  const side = Math.min(202, Math.floor((width - 134) / cases.length) - 28);
  const legendColumns = Math.min(cases.length, Math.floor((width - 56) / 225));
  const height = 961 + 4 * side + (Math.ceil(cases.length / legendColumns) - 1) * 24;
  figure.width = width * scale;
  figure.height = height * scale;
  figure.style.minWidth = `${Math.max(900, 138 * cases.length + 152)}px`;
  const ctx = figure.getContext('2d', {alpha: false});
  ctx.scale(scale, scale);
  ctx.fillStyle = '#ffffff';
  ctx.fillRect(0, 0, width, height);
  label(ctx, 'UDC 衍射成像对比', 28, 36, {size: 25, weight: 700});
  const subtitle = workspace
    ? (workspace.caution ? '所选方案存在条件差异或缺失项 · 详见工作台可比性判定' : '选中方案 · 固定色标与物理坐标')
    : '4 mm 入瞳 · 500 点/Pitch · D65 · USAF500';
  label(ctx, subtitle, 28, 58, {color: muted});
  let top = drawLegend(ctx, cases, width);
  const mode = psfMode.value;
  const modeLabel = psfMode.selectedOptions[0].text;
  const paletteLabel = psfPalette.selectedOptions[0].text;
  const reference = mode === 'absolute'
    ? '像元能量 / 入射能量，0 dB = 100% 入射能量'
    : (workspace ? '各方案窗口能量归一为 1，0 dB = 单像元占全部窗口能量'
                 : '各方案窗口能量归一为 1，0 dB = 全部方案的共同峰值');
  const full = workspace?.fullRadius ?? 100;
  const zoom = workspace?.zoomRadius ?? 25;
  sectionTitle(ctx, `G 通道 ${modeLabel} · ${workspace ? '共同视场' : '全视场'}`,
               `±${full} μm · ${reference} · −70 至 0 dB · ${paletteLabel}`, top);
  drawPsf(ctx, cases, width, top, full, side, mode);
  top += side + 102;
  sectionTitle(ctx, `G 通道 ${modeLabel} · 中心放大`, `±${zoom} μm · 沿用同一色标`, top);
  drawPsf(ctx, cases, width, top, zoom, side, mode);
  top += side + 108;
  sectionTitle(ctx, 'RGB 衍射 PSF · 对数色调映射',
               `±${full} μm · 共同曝光 · 亮度对数强度 ${rgbDisplay.toneStrength} · γ ${rgbDisplay.gamma} · sRGB`, top);
  drawPsf(ctx, cases, width, top, full, side, 'rgb');
  top += side + 108;
  const mid = width / 2;
  const box = {x: 67, y: top + 44, width: mid - 97, height: 148};
  const band = workspace ? workspace.band : [0, 500];
  const note = workspace && band ? ` · 共同频段 ${band[0].toFixed(2)}–${band[1].toFixed(2)}` : '';
  sectionTitle(ctx, 'G 通道 MTF', `逐频率取四方向最低值${note} · cycles/mm`, top);
  if (band) {
    const [xMin, xMax] = band;
    const xTicks = Array.from({length: 6}, (_, i) => Number((xMin + i * (xMax - xMin) / 5).toFixed(1)));
    drawCurves(ctx, cases, box, {xMin, xMax, xTicks, yMax: 1,
        yTicks: [0, .25, .5, .75, 1], format: n => n.toFixed(2)},
      record => [record.freq, record.mtf]);
  } else {
    label(ctx, '无共同有效频段', box.x, box.y + 30, {color: muted});
  }
  label(ctx, 'RGB 窗口捕获量', mid + 24, top, {size: 17, weight: 680});
  const fluxNote = fluxMode.value === 'captured' ? '占入射能量' :
    (workspace ? '相对各自无遮挡参考' : '201×201 接收格 · 相对无遮挡');
  label(ctx, fluxNote, mid + 24, top + 18, {color: muted});
  drawThroughput(ctx, cases, mid + 51, box.y, mid - 79, box.height);
  top += 243;
  const absolute = tailMode.value === 'tailAbsolute';
  sectionTitle(ctx, 'G 通道 PSF 拖尾', `半径之外、计算窗口以内的能量 · μm · 占${absolute ? '入射' : '窗口'}能量`, top);
  if (radii.length) {
    const peak = Math.max(...cases.flatMap(record => record[tailMode.value]));
    const maxTail = Math.max(absolute ? .12 : .42, peak);
    const tickStep = maxTail > (absolute ? .12 : .42) ? maxTail / 4 : (absolute ? .03 : .1);
    const xTicks = workspace ? [0, ...radii.filter(radius => radius > 0)] : [0, 10, 20, 30, 40, 50];
    drawCurves(ctx, cases, {x: 67, y: top + 41, width: width - 98, height: 112},
      {xMax: Math.max(1, radii.at(-1)), xTicks, yMax: maxTail,
        yTicks: [0, tickStep, 2*tickStep, 3*tickStep, 4*tickStep],
        markers: true, format: n => `${Math.round(100 * n)}%`},
      record => [radii, record[tailMode.value]]);
  } else {
    label(ctx, '无共同有效半径', 67, top + 86, {color: muted});
  }
  top += 194;
  sectionTitle(ctx, workspace ? '同曝光成像' : 'USAF500 · 同曝光成像',
               '使用仿真输出的同曝光预览 · G 通道线性均值以各自无遮挡为 100%', top);
  drawScenes(ctx, cases, width, top, side);
  figure.setAttribute('aria-label',
    `${modeLabel}（${paletteLabel}）、RGB 衍射 PSF、MTF、RGB 通量、拖尾及同曝光对比｜${cases.map(nameOf).join('、')}`);
}

function update() {
  if (!ready) return;
  const cases = selected();
  const references = cases.filter(isReference).length;
  const count = cases.length - references;
  status.textContent = `已选 ${count} 个方案${references ? ' · 已显示无遮挡参考' : ' · 已隐藏无遮挡参考'}`;
  figure.hidden = !cases.length;
  empty.hidden = cases.length > 0;
  empty.textContent = '请勾选至少一个方案或无遮挡参考。';
  exportButton.disabled = !cases.length;
  const tailIndex = workspace ? radii.indexOf(Number(tailRadius.value)) : 1;
  if (workspace) {
    byId('tail-heading').textContent = tailIndex >= 0 ? `${radii[tailIndex]} μm 外拖尾` : '外拖尾';
    const basis = fluxMode.value === 'captured' ? '入射能量' : '各自无遮挡捕获量';
    byId('metrics-note').textContent = `RGB 捕获量以${basis}为分母。MTF 和拖尾取 G 通道。`;
  }
  const rows = cases.map(record => {
    const row = document.createElement('tr');
    const values = [nameOf(record), record[fluxMode.value].map(percent).join(' / '),
      record.total.map(percent).join(' / '), percent(record.center), number(record.bandMin, 5),
      number(record.mtf300, 4), [record.tailFraction[tailIndex], record.tailAbsolute[tailIndex]].map(percent).join(' / '),
      percent(record.mean[1])];
    values.forEach((value, i) => {
      const cell = document.createElement(i ? 'td' : 'th');
      if (!i) cell.scope = 'row';
      cell.textContent = value;
      row.append(cell);
    });
    return row;
  });
  byId('metrics').replaceChildren(...rows);
  if (cases.length) drawFigure(cases);
}

exportButton.addEventListener('click', () => {
  const ids = selected().map(record => record.id).join('-');
  const suffix = psfPalette.value === 'magma' ? '' : `-${psfPalette.value}`;
  const filename = `UDC-${ids}-psf-${psfMode.value}${suffix}.png`;
  figure.toBlob(blob => {
    if (!blob) {
      status.textContent = '无法导出 PNG，请重试。';
      return;
    }
    const link = document.createElement('a');
    link.href = URL.createObjectURL(blob);
    link.download = filename;
    document.body.append(link);
    link.click();
    link.remove();
    setTimeout(() => URL.revokeObjectURL(link.href), 30000);
  }, 'image/png');
});

if (workspace) {
  document.body.classList.add('embedded');
  byId('flux-label').hidden = fluxMode.hidden = false;
  if (radii.length) {
    byId('radius-label').hidden = tailRadius.hidden = false;
    for (const radius of radii) tailRadius.add(new Option(`${radius} μm`, radius));
    if (radii.includes(20)) tailRadius.value = '20';
  }
} else {
  byId('all').addEventListener('click', () => setAll(true));
  byId('none').addEventListener('click', () => setAll(false));
  byId('reset').addEventListener('click', reset);
  records.forEach(addControls);
}
const paletteLabels = {jet: 'Jet · 蓝–红', turbo: 'Turbo', parula: 'Parula'};
for (const name of Object.keys(palettes)) {
  if (!Array.from(psfPalette.options).some(option => option.value === name)) {
    psfPalette.add(new Option(paletteLabels[name] ?? name, name));
  }
}
for (const control of [psfMode, psfPalette, tailMode, fluxMode, tailRadius]) {
  const storageKey = `udc-workspace:${control.id}`;
  if (workspace) {
    try {
      const saved = sessionStorage.getItem(storageKey);
      if (Array.from(control.options).some(option => option.value === saved)) control.value = saved;
    } catch {}
  }
  control.addEventListener('change', () => {
    if (workspace) {
      try { sessionStorage.setItem(storageKey, control.value); } catch {}
    }
    update();
  });
}
window.addEventListener('resize', update);
Promise.all(records.map(async record => {
  const [normalized, absolute, scene, rgbPsf] = await Promise.all([
    loadPsfRaster(record.psf), loadPsfRaster(record.absolutePsf), loadImage(record.scene),
    record.rgbPsf ? loadImage(record.rgbPsf) : null
  ]);
  record.psfRasters = {normalized, absolute};
  record.sceneImage = scene;
  record.rgbPsfImage = rgbPsf;
})).then(() => {
  ready = true;
  update();
}).catch(error => {
  status.textContent = '本地仿真图像加载失败。';
  figure.hidden = true;
  empty.hidden = false;
  empty.textContent = '请重新生成完整的 HTML 文件。';
  exportButton.disabled = true;
  console.error(error);
});
