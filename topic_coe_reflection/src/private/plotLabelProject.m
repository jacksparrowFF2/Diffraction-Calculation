function plotLabelProject(fig,projectName)
%plotLabelProject 给整图总标题添加项目名，保留原有参数与TeX格式。
    if ~isgraphics(fig) || strlength(string(projectName))==0, return; end
    projectName=string(projectName);
    label="项目："+replace(projectName,"_","\_");
    objects=findall(fig);
    layouts=objects(arrayfun(@(item)isa(item,'matlab.graphics.layout.TiledChartLayout'),objects));
    if isempty(layouts)
        axesList=findall(fig,'Type','axes');
        targets=arrayfun(@(item)item.Title,axesList,'UniformOutput',false);
    else
        targets=arrayfun(@(item)item.Title,layouts,'UniformOutput',false);
    end
    for j=1:numel(targets)
        target=targets{j}; old=string(target.String);
        if any(startsWith(old(:),label)), continue; end
        target.String=cellstr([label;old(:)]);
        target.Interpreter='tex';
    end
    if ~startsWith(string(fig.Name),projectName+" · "), fig.Name=projectName+" · "+string(fig.Name); end
end
