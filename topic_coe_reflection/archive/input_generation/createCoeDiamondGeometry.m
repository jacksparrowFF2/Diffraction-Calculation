function geometry = createCoeDiamondGeometry(cfg, stack)
%createCoeDiamondGeometry 参考 PDL 实线生成阳极：R/B 交错，G 位于四邻中心。
% Pitch 为 R/B 网格间距；每个 Pitch×Pitch 区域含一个 R 或 B、一个 G。
% 参考图虚线是 Mask，不作为 BM；BM/CF 按约定独立从阳极外扩。
    % 仅用于预生成示例和测试的固定轮廓，不属于用户仿真配置。
    demo.anodeHalfDiagonal = [10, 7, 14]*1e-6;
    demo.cornerRadius = [2, 2, 2]*1e-6;
    demo.greenAspectRatio = 0.65;
    demo.greenTiltDeg = 45;
    demo.cfExpansion = 7e-6;
    validateattributes(demo.anodeHalfDiagonal, {'numeric'}, {'row', 'numel', 3, 'positive', 'finite'});
    validateattributes(demo.greenAspectRatio, {'numeric'}, {'scalar', 'positive', '<=', 1});
    validateattributes(demo.cornerRadius, {'numeric'}, {'row', 'numel', 3, 'nonnegative', 'finite'});
    validateattributes(demo.greenTiltDeg, {'numeric'}, {'scalar', 'real', 'finite'});
    validateattributes(demo.cfExpansion, {'numeric'}, {'scalar', 'real', 'finite', '>', stack.bmExpansion});
    assert(all(mod(cfg.periods, 2) == 0), 'COEReflection:DiamondPeriods', 'Diamond 晶胞周期数必须为偶数。');
    n = cfg.samplesPerPeriod;
    dx = cfg.pitch/n;
    % 在完整 2×2 周期晶胞中栅格化，并包含越界邻居，避免边缘子像素被错误截断。
    coordinate = (0:2*n-1)*dx;
    [xx, yy] = meshgrid(coordinate);
    masks = false(2*n, 2*n, 3);
    for row = -1:2
        for col = -1:2
            color = 1 + 2*mod(row+col, 2); % R/B 棋盘交错
            masks(:,:,color) = masks(:,:,color) | roundedPixel(xx, yy, ...
                (col+0.25)*cfg.pitch, (row+0.25)*cfg.pitch, ...
                demo.anodeHalfDiagonal(color), 1, 45, demo.cornerRadius(color));
            % 图像行向下，绿色左右倾斜方向随列交替；所有行均有绿色。
            tilt = demo.greenTiltDeg*(-1)^col;
            masks(:,:,2) = masks(:,:,2) | roundedPixel(xx, yy, ...
                (col+0.75)*cfg.pitch, (row+0.75)*cfg.pitch, ...
                demo.anodeHalfDiagonal(2), demo.greenAspectRatio, tilt, demo.cornerRadius(2));
        end
    end
    % 外扩采用周期延拓；先在三倍晶胞上处理，再裁出中央晶胞。
    bmCells = expandPeriodic(any(masks, 3), stack.bmExpansion, dx);
    cfCells = false(size(masks));
    for color = 1:3
        cfCells(:,:,color) = expandPeriodic(masks(:,:,color), demo.cfExpansion, dx);
    end
    repeats = cfg.periods([2,1])/2;
    geometry.anodeMask = repmat(any(masks,3), repeats);
    geometry.bmOpenMask = repmat(bmCells, repeats);
    geometry.cfRMask = repmat(cfCells(:,:,1), repeats);
    geometry.cfGMask = repmat(cfCells(:,:,2), repeats);
    geometry.cfBMask = repmat(cfCells(:,:,3), repeats);
end

function mask = roundedPixel(xx, yy, cx, cy, halfDiagonal, aspect, angle, radius)
    % R/B 为旋转 45° 的圆角正方形；G 为可倾斜的圆角长方形。
    halfWidth = halfDiagonal/sqrt(2);
    halfHeight = halfWidth*aspect;
    assert(radius <= min(halfWidth, halfHeight), 'COEReflection:CornerRadius', '圆角半径不能超过局部短半轴。');
    x = (xx-cx)*cosd(angle) + (yy-cy)*sind(angle);
    y = -(xx-cx)*sind(angle) + (yy-cy)*cosd(angle);
    qx = abs(x) - (halfWidth-radius);
    qy = abs(y) - (halfHeight-radius);
    distance = hypot(max(qx,0),max(qy,0)) + min(max(qx,qy),0) - radius;
    mask = distance <= 1e-15;
end

function expanded = expandPeriodic(mask, distance, dx)
    radius = ceil(distance/dx);
    [xx, yy] = meshgrid(-radius:radius);
    disk = hypot(xx,yy)*dx <= distance+1e-15;
    tiled = conv2(double(repmat(mask,3,3)), double(disk), 'same') > 0;
    shape = size(mask);
    expanded = tiled(shape(1)+(1:shape(1)), shape(2)+(1:shape(2)));
end
