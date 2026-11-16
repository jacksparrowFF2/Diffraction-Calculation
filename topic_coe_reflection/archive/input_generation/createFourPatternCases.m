function createFourPatternCases()
%createFourPatternCases 按用户指定四图样一次性建立case01~case04；不覆盖既有输入。
% 图中毫米结构仅保留比例，以59.06μm间距缩放，便于与真实版图的角域对比。
    moduleDir=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    names=["case01","case02","case03","case04"];
    labels=["单圆","等面积椭圆","四圆阵列","四椭圆阵列（不同取向）"];
    for name=names
        assert(~isfolder(fullfile(moduleDir,'Input',name)), ...
            'COEReflection:ExistingInput','输入已存在，拒绝覆盖：%s',name);
    end
    pitch=59.06e-6; samples=250; radius=pitch/4; ratio=1.15;
    axis=((0:2*samples-1)-samples)*pitch/samples;
    [xx,yy]=meshgrid(axis);
    % 截图未给出各椭圆精确角度，此处明确采用示意角度，不声称反推得到。
    angles=[60,30;-30,-60];
    patterns=false(2*samples,2*samples,4);
    patterns(:,:,1)=xx.^2+yy.^2<=radius^2;
    patterns(:,:,2)=ellipseMask(xx,yy,radius,ratio,0);
    for row=1:2
        for col=1:2
            cx=(col-1.5)*pitch;
            cy=(1.5-row)*pitch;
            patterns(:,:,3)=patterns(:,:,3) | (xx-cx).^2+(yy-cy).^2<=radius^2;
            patterns(:,:,4)=patterns(:,:,4) | ellipseMask(xx-cx,yy-cy,radius,ratio,angles(row,col));
        end
    end
    sourceConfig=fileread(fullfile(moduleDir,'Input','F192_AnodeBM_PlaneRGB','config.txt'));
    sourceConfig=strrep(sourceConfig,'cfg.periods = [16, 16]','cfg.periods = [8, 8]');
    for j=1:4
        folder=fullfile(moduleDir,'Input',names(j)); mkdir(folder);
        imwrite(uint8(patterns(:,:,j))*255,fullfile(folder,'anode.png'));
        imwrite(uint8(patterns(:,:,j))*255,fullfile(folder,'BM.png'));
        file=fopen(fullfile(folder,'config.txt'),'w','n','UTF-8'); assert(file>=0);
        fprintf(file,'# %s：%s；图样按截图比例示意，不是产品真实版图\n%s',names(j),labels(j),sourceConfig);
        fclose(file);
        file=fopen(fullfile(folder,'README.md'),'w','n','UTF-8'); assert(file>=0);
        fprintf(file,['# %s：%s\n\n用户指定四种理想图样之一，非真实F192版图。' ...
            '500×500图片覆盖2×2 Pitch，每Pitch250点，Pitch=59.06 μm。' ...
            '圆半径=Pitch/4；椭圆长短半轴比1.15、连续面积与圆相等。' ...
            'case04左上/右上/左下/右下角度为60°/30°/−30°/−60°（示意假设）。' ...
            'BM白色为开口，invert=true。输入一次性建立，正式运行不生成或修改图片。\n\n' ...
            '完整说明见../../docs/README.md四图样对比章节。\n'],names(j),labels(j));
        fclose(file);
    end
    fprintf('已建立case01~case04，原有真实版图未改动。\n');
end

function mask=ellipseMask(x,y,radius,ratio,angle)
    u=x*cosd(angle)+y*sind(angle);
    v=-x*sind(angle)+y*cosd(angle);
    a=radius*sqrt(ratio); b=radius/sqrt(ratio);
    mask=(u/a).^2+(v/b).^2<=1;
end
