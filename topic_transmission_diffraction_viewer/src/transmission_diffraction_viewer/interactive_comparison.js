const {records, palettes, floorDb, workspace, rgbDisplay} = JSON.parse(document.getElementById('data').textContent);
const byId = id => document.getElementById(id);
const comparisonView = byId('comparison-view');
const guideView = byId('guide-view');
const comparisonLink = byId('comparison-link');
const guideLink = byId('guide-link');
const figure = byId('figure');
const plots = byId('plots');
const caseHeader = byId('case-header');
const caseTrack = byId('case-track');
const plotScrollers = [caseHeader, ...plots.querySelectorAll('.plot-window')];
const empty = byId('empty');
const status = byId('status');
const exportButton = byId('export');
const maskMode = byId('mask-mode');
const psfMode = byId('psf-mode');
const psfPalette = byId('psf-palette');
const tailMode = byId('tail-mode');
const fluxMode = byId('flux-mode');
const tailRadius = byId('tail-radius');
const sceneFilter = byId('scene-filter');
const embedded = workspace && !workspace.standalone;
const radii = workspace?.tailRadii ?? [10, 20, 50];
const ink = '#172638';
const muted = '#516277';
const gridColor = '#e0e8ee';
const rgb = ['#bd4956', '#31805b', '#416bb5'];
const radians = Math.PI * 2;
const scale = 2;
const columnLeft = 70;
const columnRight = 64;
const columnGutters = columnLeft + columnRight;
const legendRowHeight = 22;
const maskLayers = ['BM', 'Cath'];
// Palette indices encode coverage: bit 0 = BM, bit 1 = Cath.
const maskColors = [[255, 255, 255], [38, 117, 172], [243, 164, 59], [116, 81, 154]];
let ready = false;
let comparisonScrollY = 0;

const percent = value => Number.isFinite(value) ? `${(100 * value).toFixed(2)}%` : '—';
const number = (value, digits) => Number.isFinite(value) ? value.toFixed(digits) : '—';
const nameOf = record => record.nameInput?.value.trim() || record.name;
const isReference = record => record.reference || record.id === 'clear';
const matchesScene = record => !workspace?.standalone || sceneFilter.value === 'all' || record.sceneMode === sceneFilter.value;
const groups = new Map();
for (const record of records) {
  const key = record.scenePair || record.id;
  if (!groups.has(key)) groups.set(key, []);
  groups.get(key).push(record);
}
const opticalRecord = variants => variants.find(record => record.sceneMode === 'object') ?? variants[0];
const activeGroups = () => Array.from(groups.values(), variants =>
  variants.filter(record => matchesScene(record) && (embedded || record.checkbox.checked))).filter(group => group.length);
const selected = () => activeGroups().map(opticalRecord);

function addControls(variants) {
  const record = opticalRecord(variants);
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
  const idLabel = workspace?.standalone ? record.name :
                  record.id === 'clear' ? '无遮挡参考' : record.id.replace('case', 'Case');
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
  byId(isReference(record) ? 'reference-control' : 'case-controls').append(row);
  for (const variant of variants) Object.assign(variant, {checkbox, nameInput, controlRow: row});
}

function setAll(checked) {
  for (const variants of groups.values()) {
    const record = opticalRecord(variants);
    if (!isReference(record) && variants.some(matchesScene)) record.checkbox.checked = checked;
  }
  update();
}

function reset() {
  for (const variants of groups.values()) {
    const record = opticalRecord(variants);
    record.checkbox.checked = record.selected;
    record.nameInput.value = record.name;
  }
  psfMode.value = 'normalized';
  psfPalette.value = 'magma';
  tailMode.value = 'tailFraction';
  maskMode.value = 'layers';
  if (workspace?.standalone) sceneFilter.value = 'all';
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

function legendColumns(cases, width) {
  return Math.min(cases.length, Math.max(1, Math.floor(width / 180)));
}

function legendHeight(cases, width) {
  return legendRowHeight * Math.ceil(cases.length / legendColumns(cases, width));
}

function drawLegend(ctx, cases, left, top, width) {
  const cols = legendColumns(cases, width);
  const cell = width / cols;
  cases.forEach((record, i) => {
    const x = left + (i % cols) * cell;
    const y = top + Math.floor(i / cols) * legendRowHeight;
    lineStyle(ctx, record);
    ctx.beginPath();
    ctx.moveTo(x, y - 4);
    ctx.lineTo(x + 21, y - 4);
    ctx.stroke();
    ctx.setLineDash([]);
    label(ctx, nameOf(record), x + 29, y, {size: 12, max: cell - 36});
  });
}

function caseDash(record) {
  return isReference(record) ? [3, 5] : record.dash ?? (record.id === 'case8' ? [7, 5] : []);
}

function lineStyle(ctx, record) {
  ctx.strokeStyle = record.color;
  ctx.lineWidth = 2.2;
  ctx.setLineDash(caseDash(record));
}

function syncScroll(left) {
  for (const scroller of plotScrollers) {
    if (scroller.scrollLeft !== left) scroller.scrollLeft = left;
  }
}

function drawCaseHeader(cases, width) {
  caseTrack.style.width = `${width}px`;
  caseTrack.style.paddingInline = `${columnLeft}px ${columnRight}px`;
  caseTrack.style.gridTemplateColumns = `repeat(${cases.length}, minmax(0, 1fr))`;
  caseHeader.tabIndex = width > caseHeader.clientWidth ? 0 : -1;
  caseTrack.replaceChildren(...cases.map(record => {
    const heading = document.createElement('div');
    heading.className = 'case-heading';
    heading.title = nameOf(record);
    const mark = document.createElement('span');
    mark.className = 'case-mark';
    mark.style.borderColor = record.color;
    if (caseDash(record).length) mark.style.borderBlockStartStyle = 'dashed';
    const name = document.createElement('span');
    name.textContent = heading.title;
    heading.append(mark, name);
    return heading;
  }));
}

function drawPlot(id, width, start, end) {
  const canvas = byId(id);
  const height = (end - start) * scale;
  canvas.width = width * scale;
  canvas.height = height;
  canvas.style.width = `${width}px`;
  // Use the same finished raster for the plots on screen and the PNG export.
  canvas.getContext('2d', {alpha: false}).drawImage(figure, 0, start * scale, canvas.width, height,
                                                   0, 0, canvas.width, height);
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

function composeMaskOverlay({BM, Cath}) {
  const source = BM ?? Cath;
  if (!source) return {image: null};
  const {naturalWidth: width, naturalHeight: height} = source;
  if (BM && Cath && (Cath.naturalWidth !== width || Cath.naturalHeight !== height)) {
    return {image: null, message: '掩膜网格尺寸不同'};
  }

  const canvas = document.createElement('canvas');
  canvas.width = width;
  canvas.height = height;
  const ctx = canvas.getContext('2d', {willReadFrequently: true});
  const readCoverage = image => {
    if (!image) return null;
    ctx.clearRect(0, 0, width, height);
    ctx.drawImage(image, 0, 0);
    return ctx.getImageData(0, 0, width, height).data;
  };
  const bm = readCoverage(BM);
  const cath = readCoverage(Cath);
  const pixels = ctx.createImageData(width, height);
  for (let i = 0; i < pixels.data.length; i += 4) {
    // Coverage previews store dark material and white openings on the same grid.
    const index = (bm && bm[i] < 128 ? 1 : 0) + (cath && cath[i] < 128 ? 2 : 0);
    const [r, g, b] = maskColors[index];
    pixels.data[i] = r;
    pixels.data[i + 1] = g;
    pixels.data[i + 2] = b;
    pixels.data[i + 3] = 255;
  }
  ctx.putImageData(pixels, 0, 0);
  return {image: canvas};
}

function drawMaskLegend(ctx, y) {
  const names = ['共同开口', '仅 BM', '仅 Cath', 'BM + Cath'];
  let x = 28;
  names.forEach((name, i) => {
    ctx.fillStyle = `rgb(${maskColors[i].join(',')})`;
    ctx.fillRect(x, y - 11, 12, 12);
    ctx.strokeStyle = gridColor;
    ctx.lineWidth = 1;
    ctx.strokeRect(x, y - 11, 12, 12);
    label(ctx, name, x + 18, y, {color: muted});
    x += ctx.measureText(name).width + 35;
  });
}

function drawMasks(ctx, cases, width, top, side, overlay) {
  const cell = (width - columnGutters) / cases.length;
  const rows = overlay ? ['叠加'] : maskLayers;
  const rowY = index => top + 50 + index * (side + 12);
  if (overlay) {
    label(ctx, 'BM / Cath 掩膜叠加', 28, top, {size: 17, weight: 680});
    drawMaskLegend(ctx, top + 18);
  } else {
    sectionTitle(ctx, 'BM / Cath 掩膜', '物面覆盖 · 深色表示材料覆盖，白色表示开口', top);
  }
  for (const [index, layer] of rows.entries()) {
    label(ctx, layer, 28, rowY(index) + side / 2 + 4, {weight: 660});
  }
  ctx.imageSmoothingEnabled = false;
  cases.forEach((record, i) => {
    const center = columnLeft + cell * (i + .5);
    const preview = record.maskPreview;
    label(ctx, nameOf(record), center, top + 40, {align: 'center', weight: 660, max: cell - 10});
    if (overlay) record.maskOverlay ??= composeMaskOverlay(record.maskImages);
    let drawn = false;
    for (const [index, layer] of rows.entries()) {
      const middle = rowY(index) + side / 2;
      const image = overlay ? record.maskOverlay.image : record.maskImages[layer];
      if (image) {
        drawn = true;
        const fit = side / Math.max(preview.widthUm, preview.heightUm);
        const w = preview.widthUm * fit;
        const h = preview.heightUm * fit;
        ctx.drawImage(image, center - w / 2, middle - h / 2, w, h);
        ctx.strokeStyle = gridColor;
        ctx.lineWidth = 1;
        ctx.strokeRect(center - w / 2, middle - h / 2, w, h);
      } else {
        ctx.fillStyle = '#f2f6f8';
        ctx.fillRect(center - side / 2, rowY(index), side, side);
        const absent = overlay ? '未设置 BM / Cath' : `未设置 ${layer}`;
        const message = (overlay && record.maskOverlay.message) ||
          (preview || isReference(record) ? absent : '缺少掩膜数据');
        label(ctx, message, center, middle + 4, {size: 11, color: muted, align: 'center', max: side - 8});
      }
    }
    if (drawn) {
      const captionY = rowY(rows.length - 1) + side + 17;
      const extent = `${Number(preview.widthUm.toPrecision(6))} × ${Number(preview.heightUm.toPrecision(6))} μm`;
      const scope = preview.tiled ? '平铺单元' : '完整物面';
      const caption = `${scope}${preview.subsampled ? '（抽样）' : ''} · ${extent}`;
      const style = {size: 10, color: muted, align: 'center', max: cell - 8};
      label(ctx, caption, center, captionY, style);
      const missingLayer = overlay && maskLayers.find(layer => !record.maskImages[layer]);
      if (missingLayer) label(ctx, `未设置 ${missingLayer}`, center, captionY + 15, style);
    }
  });
}

function drawPsf(ctx, cases, width, top, radius, side, mode) {
  const rgbPsf = mode === 'rgb';
  const cell = (width - columnGutters) / cases.length;
  cases.forEach((record, i) => {
    const center = columnLeft + cell * (i + .5);
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

function drawPupils(ctx, cases, width, top, side) {
  const cell = (width - columnGutters) / cases.length;
  const extent = Math.max(0, ...cases.flatMap(record => record.pupilPreview?.extent.map(Math.abs) ?? []));
  const rounding = 10 ** (Math.floor(Math.log10(extent || 1)) - 1);
  const radius = Math.ceil(extent / rounding) * rounding;
  const firstImage = cases.findIndex(record => record.pupilImage);
  sectionTitle(ctx, '物面 · 有效振幅入瞳', '选用最接近 550 nm 的仿真波长 · 共同物面坐标和 0–1 振幅色标', top);
  const y = top + 50;
  ctx.imageSmoothingEnabled = true;
  ctx.imageSmoothingQuality = 'high';
  cases.forEach((record, i) => {
    const center = columnLeft + (i + .5) * cell;
    const x = center - side / 2;
    label(ctx, nameOf(record), center, top + 40, {align: 'center', weight: 660, max: cell - 10});
    ctx.fillStyle = '#f2f6f8';
    ctx.fillRect(x, y, side, side);
    if (!record.pupilImage) {
      label(ctx, '缺少物面数据', center, y + side / 2 + 4,
            {size: 11, color: muted, align: 'center', max: side - 8});
      return;
    }

    const [left, right, bottom, upper] = record.pupilPreview.extent;
    const px = side / (2 * radius);
    ctx.drawImage(record.pupilImage, center + left * px, y + (radius - upper) * px,
                  (right - left) * px, (upper - bottom) * px);
    ctx.strokeStyle = gridColor;
    ctx.lineWidth = 1;
    ctx.strokeRect(x, y, side, side);
    for (const tick of [-radius, 0, radius]) {
      const text = String(Number(tick.toPrecision(5)));
      label(ctx, text, center + tick * px, y + side + 14, {size: 10, color: muted, align: 'center'});
      if (i === firstImage) {
        label(ctx, text, x - 5, y + side / 2 - tick * px + 3,
              {size: 10, color: muted, align: 'right'});
      }
    }
    if (i === firstImage) {
      label(ctx, 'y (μm)', x - 4, y - 7, {size: 10, color: muted, align: 'right'});
    }
    const detail = isReference(record) ? '几何入瞳' : `${Number(record.pupilPreview.wavelengthNm.toFixed(2))} nm`;
    label(ctx, `x (μm) · ${detail}`, center, y + side + 29, {size: 10, color: muted, align: 'center'});
  });
  if (radius > 0) {
    const x = width - 52;
    const gradient = ctx.createLinearGradient(0, y, 0, y + side);
    gradient.addColorStop(0, '#fff');
    gradient.addColorStop(1, '#000');
    ctx.fillStyle = gradient;
    ctx.fillRect(x, y, 12, side);
    ctx.strokeStyle = gridColor;
    ctx.strokeRect(x, y, 12, side);
    for (const value of [1, .5, 0]) {
      label(ctx, String(value), x + 16, y + (1 - value) * side + 4, {size: 10, color: muted});
    }
  }
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

function drawScenes(ctx, columns, width, top, side, stacked) {
  const cell = (width - columnGutters) / columns.length;
  const rowModes = stacked ? ['sensor', 'object'] : [null];
  ctx.imageSmoothingEnabled = true;
  ctx.imageSmoothingQuality = 'high';
  columns.forEach((variants, i) => {
    const center = columnLeft + (i + .5) * cell;
    label(ctx, nameOf(opticalRecord(variants)), center, top + 41, {align: 'center', weight: 660, max: cell - 10});
    rowModes.forEach((mode, row) => {
      const record = stacked ? variants.find(item => item.sceneMode === mode)
        ?? (row === 0 && variants.find(item => !item.sceneMode)) : variants[0];
      const rowTop = top + 52 + row * (side + 35);
      if (stacked && i === 0) label(ctx, mode, 24, rowTop + side / 2, {size: 12, weight: 660, color: muted});
      if (!record) return;
      const image = record.sceneImage;
      const fit = side / Math.max(image.naturalWidth, image.naturalHeight);
      const w = image.naturalWidth * fit;
      const h = image.naturalHeight * fit;
      ctx.drawImage(image, center - w / 2, rowTop + (side - h) / 2, w, h);
      const note = workspace?.brightnessMatched ? record.sceneNote : `G 均值 ${percent(record.mean[1])}`;
      if (note) label(ctx, note, center, rowTop + side + 17, {size: 11, color: muted, align: 'center'});
    });
  });
}

function drawFigure(cases, sceneColumns) {
  const width = Math.max(900, 138 * cases.length + 152, plots.clientWidth);
  const cellWidth = Math.floor((width - columnGutters) / cases.length);
  const side = Math.min(202, cellWidth - 28);
  const overlay = maskMode.value === 'overlay';
  const maskSide = overlay ? side : Math.min(152, side);
  const maskHeight = overlay ? maskSide + 110 : 2 * maskSide + 108;
  const pupilSide = Math.min(202, cellWidth - 40);
  const pupilHeight = pupilSide + 105;
  const mid = width / 2;
  const band = workspace ? workspace.band : [0, 500];
  const responseHeight = band ? 240 + legendHeight(cases, mid - 97) : 257;
  const tailHeight = radii.length ? 200 + legendHeight(cases, width - 98) : 194;
  const maskTop = 91 + legendHeight(cases, width - 56);
  const stackedScenes = sceneColumns.some(variants => variants.length > 1);
  const height = maskTop + maskHeight + pupilHeight + (3 * side + 318) + responseHeight + tailHeight
    + side + 93 + (stackedScenes ? side + 35 : 0);
  figure.width = width * scale;
  figure.height = height * scale;
  const ctx = figure.getContext('2d', {alpha: false});
  ctx.scale(scale, scale);
  ctx.fillStyle = '#ffffff';
  ctx.fillRect(0, 0, width, height);
  label(ctx, 'UDC 衍射成像对比', 28, 36, {size: 25, weight: 700});
  const subtitle = workspace
    ? (stackedScenes ? '光学图表采用 object · 成像上排 sensor、下排 object' :
       !workspace.standalone && workspace.caution ? '所选方案存在条件差异或缺失项 · 详见工作台可比性判定' : '选中方案 · 固定色标与物理坐标')
    : '4 mm 入瞳 · 500 点/Pitch · D65 · USAF500';
  label(ctx, subtitle, 28, 58, {color: muted});
  drawLegend(ctx, cases, 28, 91, width - 56);
  drawMasks(ctx, cases, width, maskTop, maskSide, overlay);
  const pupilTop = maskTop + maskHeight;
  drawPupils(ctx, cases, width, pupilTop, pupilSide);
  let top = pupilTop + pupilHeight;
  const psfTop = top;
  const mode = psfMode.value;
  const modeLabel = psfMode.selectedOptions[0].text;
  const paletteLabel = psfPalette.selectedOptions[0].text;
  const reference = mode === 'absolute'
    ? '像元能量 / 入射能量，0 dB = 100% 入射能量'
    : (workspace ? '各方案窗口能量归一为 1，0 dB = 单像元占全部窗口能量'
                 : '各方案窗口能量归一为 1，0 dB = 全部方案的共同峰值');
  const full = workspace?.fullRadius ?? 100;
  const zoom = workspace?.zoomRadius ?? 25;
  const psfTitle = `G 通道 ${modeLabel} · ${workspace ? '共同视场' : '全视场'}`;
  sectionTitle(ctx, psfTitle, `±${full} μm · ${reference} · −70 至 0 dB · ${paletteLabel}`, top);
  drawPsf(ctx, cases, width, top, full, side, mode);
  top += side + 102;
  sectionTitle(ctx, `G 通道 ${modeLabel} · 中心放大`, `±${zoom} μm · 沿用同一色标`, top);
  drawPsf(ctx, cases, width, top, zoom, side, mode);
  top += side + 108;
  sectionTitle(ctx, 'RGB 衍射 PSF · 对数色调映射',
               `±${full} μm · 共同曝光 · 亮度对数强度 ${rgbDisplay.toneStrength} · γ ${rgbDisplay.gamma} · sRGB`, top);
  drawPsf(ctx, cases, width, top, full, side, 'rgb');
  top += side + 108;
  const responseTop = top;
  const box = {x: 67, y: top + 44, width: mid - 97, height: 148};
  const note = workspace && band ? ` · 共同频段 ${band[0].toFixed(2)}–${band[1].toFixed(2)}` : '';
  sectionTitle(ctx, 'G 通道 MTF', `逐频率取四方向最低值${note} · cycles/mm`, top);
  if (band) {
    const [xMin, xMax] = band;
    const xTicks = Array.from({length: 6}, (_, i) => Number((xMin + i * (xMax - xMin) / 5).toFixed(1)));
    drawCurves(ctx, cases, box, {xMin, xMax, xTicks, yMax: 1,
        yTicks: [0, .25, .5, .75, 1], format: n => n.toFixed(2)},
      record => [record.freq, record.mtf]);
    drawLegend(ctx, cases, box.x, top + 232, box.width);
  } else {
    label(ctx, '无共同有效频段', box.x, box.y + 30, {color: muted});
  }
  label(ctx, 'RGB 窗口捕获量', mid + 24, top, {size: 17, weight: 680});
  const fluxNote = fluxMode.value === 'captured' ? '占入射能量' :
    (workspace ? '相对各自无遮挡参考' : '201×201 接收格 · 相对无遮挡');
  label(ctx, fluxNote, mid + 24, top + 18, {color: muted});
  drawThroughput(ctx, cases, mid + 51, box.y, mid - 79, box.height);
  top += responseHeight;
  const tailTop = top;
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
    drawLegend(ctx, cases, 67, top + 193, width - 98);
  } else {
    label(ctx, '无共同有效半径', 67, top + 86, {color: muted});
  }
  top += tailHeight;
  const sceneTitle = workspace?.brightnessMatched ? '亮度匹配成像' :
                     workspace ? '同曝光成像' : 'USAF500 · 同曝光成像';
  const sceneCaption = stackedScenes ? '同一 Case 共用一列 · sensor 在上 · object 在下' :
                       workspace?.brightnessMatched ? '使用仿真输出的亮度匹配预览' :
                       '使用仿真输出的同曝光预览 · G 通道线性均值以各自无遮挡为 100%';
  sectionTitle(ctx, sceneTitle, sceneCaption, top);
  drawScenes(ctx, sceneColumns, width, top, side, stackedScenes);

  byId('mask-title').textContent = overlay ? 'BM / Cath 掩膜叠加' : 'BM / Cath 掩膜';
  byId('psf-title').textContent = psfTitle;
  byId('scene-title').textContent = sceneTitle;
  byId('scene-plot').setAttribute('aria-label', stackedScenes ? `${sceneTitle}，sensor 在上、object 在下` : sceneTitle);
  const left = caseHeader.scrollLeft;
  drawCaseHeader(cases, width);
  // HTML titles replace the first canvas title in sections with a header bar.
  drawPlot('mask-plot', width, maskTop + 3, pupilTop - 18);
  drawPlot('pupil-plot', width, pupilTop + 3, psfTop - 18);
  drawPlot('psf-plot', width, psfTop + 3, responseTop - 18);
  drawPlot('response-plot', width, responseTop - 18, tailTop - 18);
  drawPlot('tail-plot', width, tailTop + 3, top - 18);
  drawPlot('scene-plot', width, top + 3, height);
  syncScroll(Math.min(left, width - caseHeader.clientWidth));
}

function update() {
  if (!ready || comparisonView.hidden) return;
  if (workspace?.standalone) {
    for (const variants of groups.values()) opticalRecord(variants).controlRow.hidden = !variants.some(matchesScene);
  }
  const sceneColumns = activeGroups();
  const cases = sceneColumns.map(opticalRecord);
  const references = cases.filter(isReference).length;
  const count = cases.length - references;
  status.textContent = `已选 ${count} 个方案${references ? ' · 已显示无遮挡参考' : ' · 已隐藏无遮挡参考'}`;
  plots.hidden = caseHeader.hidden = !cases.length;
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
  if (cases.length) drawFigure(cases, sceneColumns);
}

function showView(showGuide) {
  if (showGuide === comparisonView.hidden) return;
  if (showGuide) comparisonScrollY = window.scrollY;
  comparisonView.hidden = showGuide;
  guideView.hidden = !showGuide;
  const activeLink = showGuide ? guideLink : comparisonLink;
  const inactiveLink = showGuide ? comparisonLink : guideLink;
  activeLink.setAttribute('aria-current', 'page');
  inactiveLink.removeAttribute('aria-current');
  if (!showGuide) update();
  const focusTarget = showGuide ? byId('guide-title') : comparisonLink;
  focusTarget.focus({preventScroll: true});
  window.scrollTo(0, showGuide ? 0 : comparisonScrollY);
}

exportButton.addEventListener('click', () => {
  const ids = selected().map(record => record.id).join('-');
  const suffix = psfPalette.value === 'magma' ? '' : `-${psfPalette.value}`;
  const maskSuffix = maskMode.value === 'overlay' ? '-mask-overlay' : '';
  const filename = `UDC-${ids}-psf-${psfMode.value}${suffix}${maskSuffix}.png`;
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
  byId('flux-controls').hidden = false;
  if (radii.length) {
    byId('radius-label').hidden = tailRadius.hidden = false;
    for (const radius of radii) tailRadius.add(new Option(`${radius} μm`, radius));
    if (radii.includes(20)) tailRadius.value = '20';
  }
}
if (embedded) {
  document.body.classList.add('embedded');
  // Fragment links in srcdoc resolve against the host page. Switch views in place.
  for (const link of document.querySelectorAll('a[href="#guide"], a[href="#comparison"]')) {
    link.addEventListener('click', event => {
      event.preventDefault();
      showView(link === guideLink);
    });
  }
} else {
  byId('all').addEventListener('click', () => setAll(true));
  byId('none').addEventListener('click', () => setAll(false));
  byId('reset').addEventListener('click', reset);
  groups.forEach(addControls);
  if (workspace?.standalone) {
    byId('comparison-summary').textContent = 'Case9–12 · sensor / object · 亮度匹配成像';
    byId('scene-filter-control').hidden = false;
    sceneFilter.addEventListener('change', update);
  }
  const syncView = () => showView(window.location.hash === '#guide');
  window.addEventListener('hashchange', syncView);
  syncView();
}
const paletteLabels = {jet: 'Jet · 蓝–红', turbo: 'Turbo', parula: 'Parula'};
for (const name of Object.keys(palettes)) {
  if (!Array.from(psfPalette.options).some(option => option.value === name)) {
    psfPalette.add(new Option(paletteLabels[name] ?? name, name));
  }
}
for (const control of [maskMode, psfMode, psfPalette, tailMode, fluxMode, tailRadius]) {
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
for (const scroller of plotScrollers) {
  scroller.addEventListener('scroll', () => syncScroll(scroller.scrollLeft));
}
window.addEventListener('resize', update);
Promise.all(records.map(async record => {
  const [normalized, absolute, scene, rgbPsf, masks, pupil] = await Promise.all([
    loadPsfRaster(record.psf), loadPsfRaster(record.absolutePsf), loadImage(record.scene),
    record.rgbPsf ? loadImage(record.rgbPsf) : null,
    Promise.all(maskLayers.map(layer => {
      const url = record.maskPreview?.layers[layer];
      return url ? loadImage(url) : null;
    })),
    record.pupilPreview ? loadImage(record.pupilPreview.png) : null,
  ]);
  record.psfRasters = {normalized, absolute};
  record.sceneImage = scene;
  record.rgbPsfImage = rgbPsf;
  record.maskImages = Object.fromEntries(maskLayers.map((layer, i) => [layer, masks[i]]));
  record.pupilImage = pupil;
})).then(() => {
  ready = true;
  update();
}).catch(error => {
  status.textContent = '本地仿真图像加载失败。';
  plots.hidden = caseHeader.hidden = true;
  empty.hidden = false;
  empty.textContent = '请重新生成完整的 HTML 文件。';
  exportButton.disabled = true;
  console.error(error);
});
