function Aperture = CustomPattern(angle,x,y,a,b,offsetx,offsety)
    % 生成椭圆或圆的离散化矩阵
    % 椭圆旋转角度
    theta = deg2rad(angle); 
    % 中心点是否偏移
    Xs = x - offsetx;
    Ys = y - offsety;
    Xr = Xs*cos(theta) + Ys*sin(theta);
    Yr = -Xs*sin(theta) + Ys*cos(theta);
    
    Aperture = double(Xr.^2/a^2 + Yr.^2/b^2 <= 1);
end