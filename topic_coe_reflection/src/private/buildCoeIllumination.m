function [fields, weights, info] = buildCoeIllumination(cfg, xx, yy, pupil, lambda)
%buildCoeIllumination 各互不相干源点的物面复场；各自场内保留空间相干。
% 总入射功率统一标定到 Iinc*几何光瞳面积，不猜测未知的LED总功率。
    source = coeSourceDefaults();
    if isfield(cfg,'source'), source=cfg.source; end
    mustBeMember(source.type,["plane","point","laser"]);
    validateattributes(source.centerXY,{'numeric'},{'size',[1,2],'real','finite'});
    weights=1; fields={complex(ones(size(xx)))};
    info=struct('type',source.type,'modeCount',1,'normalizationScale',1, ...
        'maxIncidenceAngleDeg',0,'beamRadius',NaN,'divergenceHalfAngleDeg',NaN);
    if source.type=="point"
        validateattributes(source.distance,{'numeric'},{'scalar','real','finite','positive'});
        validateattributes(source.diameter,{'numeric'},{'scalar','real','finite','nonnegative'});
        validateattributes(source.samplesPerAxis,{'numeric'},{'scalar','integer','positive','odd'});
        validateattributes(source.beamHalfAngleDeg,{'numeric'},{'scalar','real','finite','>',0,'<',90});
        positions=[0,0];
        if source.diameter>0
            assert(source.samplesPerAxis>=3,'COEReflection:SourceSampling', ...
                '有限尺寸光源需至少3×3采样；请进一步增加采样检查收敛。');
            axis=((0:source.samplesPerAxis-1)+0.5)/source.samplesPerAxis-0.5;
            [sx,sy]=meshgrid(axis*source.diameter);
            inside=hypot(sx,sy)<=source.diameter/2;
            positions=[sx(inside),sy(inside)];
        end
        weights=ones(1,size(positions,1))/size(positions,1);
        fields=cell(1,numel(weights));
        for j=1:numel(weights)
            radius=hypot(xx-source.centerXY(1)-positions(j,1),yy-source.centerXY(2)-positions(j,2));
            range=hypot(source.distance,radius);
            angle=atan2d(radius,source.distance);
            % 路径差的稳定写法，避免 R-z 在远距离下的消减误差。
            path=radius.^2./(range+source.distance);
            amplitude=(source.distance./range).^1.5 ...
                .*exp(-0.5*log(2)*(angle/source.beamHalfAngleDeg).^2);
            fields{j}=amplitude.*exp(2i*pi*path/lambda);
            info.maxIncidenceAngleDeg=max(info.maxIncidenceAngleDeg,max(angle(pupil>0)));
        end
    elseif source.type=="laser"
        mustBeMember(source.laserSpecification,["waist","divergence"]);
        validateattributes(source.waistDistance,{'numeric'},{'scalar','real','finite'});
        if source.laserSpecification=="waist"
            validateattributes(source.waistRadius,{'numeric'},{'scalar','real','finite','positive'});
            waist=source.waistRadius;
        else
            validateattributes(source.divergenceHalfAngleDeg,{'numeric'},{'scalar','real','finite','>',0,'<=',5});
            waist=lambda/(pi*deg2rad(source.divergenceHalfAngleDeg));
        end
        rayleigh=pi*waist^2/lambda;
        distance=source.waistDistance;
        width=waist*sqrt(1+(distance/rayleigh)^2);
        radiusSquared=(xx-source.centerXY(1)).^2+(yy-source.centerXY(2)).^2;
        inverseCurvature=distance/(distance^2+rayleigh^2);
        fields={exp(-radiusSquared/width^2).*exp(1i*pi*radiusSquared*inverseCurvature/lambda)};
        info.beamRadius=width;
        info.divergenceHalfAngleDeg=rad2deg(lambda/(pi*waist));
        assert(info.divergenceHalfAngleDeg<=5,'COEReflection:GaussianParaxial', ...
            '当前基模近轴高斯模型要求发散半角不超过5°。');
        info.maxIncidenceAngleDeg=rad2deg(max(sqrt(radiusSquared(pupil>0)))*abs(inverseCurvature));
    end
    info.modeCount=numel(weights);
    power=0;
    for j=1:numel(weights), power=power+weights(j)*sum(abs(fields{j}.*pupil).^2,'all'); end
    assert(power>0 && isfinite(power),'COEReflection:SourcePower','光源在物面上的入射能量为零或非有限值。');
    scale=sqrt(sum(pupil(:).^2)/power);
    for j=1:numel(weights), fields{j}=scale*fields{j}; end
    info.normalizationScale=scale;
end
