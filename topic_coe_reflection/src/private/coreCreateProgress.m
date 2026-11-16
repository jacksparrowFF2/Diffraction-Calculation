function [updateProgress, closeProgress] = coreCreateProgress(enabled)
%coreCreateProgress 创建中文进度窗口；关闭窗口只隐藏进度，不中断计算。
% 进度按计算任务划分，并非剩余时间预测；回调不写入 cfg 或物理结果。
    progressFigure = gobjects(1); % 无效的标量句柄：关闭开关/无桌面时更新直接返回
    started = tic;
    lastDraw = -Inf;
    if enabled && usejava('desktop')
        progressFigure = waitbar(0, '准备仿真…', 'Name', 'COE 仿真进度', ...
            'WindowStyle', 'normal', 'HandleVisibility', 'off', ...
            'CloseRequestFcn', @(window,~) delete(window));
        set(findall(progressFigure, 'Type', 'text'), 'Interpreter', 'none');
    end
    updateProgress = @update;
    closeProgress = @finish;

    function update(fraction, message)
        if ~isgraphics(progressFigure), return; end
        elapsed = toc(started);
        % 限制刷新频率，避免多源点计算时 UI 开销影响传播速度。
        if elapsed-lastDraw < 0.15 && fraction > 0 && fraction < 1, return; end
        fraction = min(1, max(0, fraction));
        waitbar(fraction, progressFigure, ...
            {sprintf('总进度 %.1f%%  |  已用 %.1f 秒', 100*fraction, elapsed), char(message)});
        lastDraw = elapsed;
        drawnow limitrate;
    end

    function finish()
        if isgraphics(progressFigure), delete(progressFigure); end
    end
end
