% 主函数 - 光谱分析实时处理系统
function spectral_analysis_real_time()
    % 全局变量用于存储绘图数据和图形句柄
    global processed_files_data current_spectrum ax_handles all_spectra_data
    
    % 初始化数据存储
    processed_files_data = struct();
    current_spectrum = [];
    ax_handles = [];
    all_spectra_data = struct(); 
    
    % =========================================================================
    % --- 用户配置区域：报警阈值设置 ---
    % =========================================================================
    
    % 设置目标值 (Target) 和 允许偏差 (Tolerance)
    % 格式：[目标值, 偏差值]
    
    % PP (Peak Position)
    cfg_PP_target = 669.4;   
    cfg_PP_tol    = 0.5;   
    
    % PPC (Peak Position Center)
    cfg_PPC_target = 690.7;  
    cfg_PPC_tol    = 0.5;  
    
    % PI (Peak Intensity)
    cfg_PI_target = 10262; 
    cfg_PI_tol    = 500;   
    
    % =========================================================================
    
    % 创建图形窗口 
    fig = figure('Name', '实时光谱分析系统 (Real-time Spectral Analysis)', 'NumberTitle', 'off', ...
                 'Position', [50, 20, 1200, 1200], 'Color', 'w');
    
    % 输入参数 - 文件夹路径
    input_directory = "*****"; 
    
    % 检查输入目录是否存在
    if ~exist(input_directory, 'dir')
        error('输入目录不存在: %s', input_directory);
    end
    
    % CSV文件名生成
    folders = strsplit(input_directory, '\');
    last_folder = folders{end};
    base_name = sprintf('spectral_analysis_results_%s.csv', last_folder);
    output_csv = get_new_csv_filename(input_directory, base_name);
    
    fprintf('输出文件将保存为: %s\n', output_csv);
    fprintf('报警条件设置:\n');
    fprintf('  PP:  %.1f ± %.1f nm\n', cfg_PP_target, cfg_PP_tol);
    fprintf('  PPC: %.1f ± %.1f nm\n', cfg_PPC_target, cfg_PPC_tol);
    fprintf('  PI:  %.0f ± %.0f\n', cfg_PI_target, cfg_PI_tol);
    
    processed_files = {}; 
    
    % 定义颜色
    peak_color = [0.85, 0.1, 0.1];      % 深红色
    avg_color = [0.1, 0.4, 0.8];        % 深蓝色
    intensity_color = [0.95, 0.6, 0.0]; % 深橙色
    
    % 初始化图形布局
    initialize_plots();
    
    try
        while true
            files = dir(fullfile(input_directory, '*.txt'));
            
            for i = 1:length(files)
                file = files(i).name;
                file_path = fullfile(input_directory, file);
                
                if ~ismember(file, processed_files)
                    fprintf('正在处理 %s...\n', file);
                    data = get_spectral_data(file_path);
                    
                    % 降噪处理
                    denoised_data = denoise_data(data);
                    current_spectrum = denoised_data;
                    
                    % 存储所有光谱数据用于2D图
                    all_spectra_data.(matlab.lang.makeValidName(file)) = denoised_data;
                    
                    % 使用多项式拟合来检测峰值
                    [fitted_peak_wavelength, fitted_max_intensity] = fit_peak_enhanced(denoised_data);
                    
                    % 分析光谱数据
                    [peak_wavelength, average_wavelength, max_intensity, left_wavelength, right_wavelength] = analyze_spectrum(denoised_data);
                    
                    % --- 报警逻辑判断 ---
                    if isnan(cfg_PP_target)
                        is_pp_ok = true;
                    else
                        is_pp_ok = abs(fitted_peak_wavelength - cfg_PP_target) <= cfg_PP_tol;
                    end
                    
                    if isnan(cfg_PPC_target)
                        is_ppc_ok = true;
                    else
                        is_ppc_ok = abs(average_wavelength - cfg_PPC_target) <= cfg_PPC_tol;
                    end
                    
                    if isnan(cfg_PI_target)
                        is_pi_ok = true;
                    else
                        is_pi_ok = abs(max_intensity - cfg_PI_target) <= cfg_PI_tol;
                    end
                    
                    is_alarm_triggered = is_pp_ok && is_ppc_ok && is_pi_ok;
                    
                    if is_alarm_triggered
                        fprintf('>>> 警报触发！数据在目标范围内 <<<\n');
                        play_alarm_sound(); 
                    end
                    
                    % 存储数据
                    processed_files_data.(matlab.lang.makeValidName(file)) = struct(...
                        'peak', fitted_peak_wavelength, ...
                        'avg', average_wavelength, ...
                        'intensity', max_intensity, ...
                        'is_alarm', is_alarm_triggered);
                    
                    % 写入CSV
                    row_data = {file, fitted_peak_wavelength, average_wavelength, max_intensity, left_wavelength, right_wavelength};
                    write_csv_data(output_csv, row_data);
                    
                    fprintf('文件: %s\n', file);
                    fprintf('PP:  %.1f (目标: %.1f ± %.1f) -> %s\n', fitted_peak_wavelength, cfg_PP_target, cfg_PP_tol, string_status(is_pp_ok));
                    fprintf('PPC: %.1f (目标: %.1f ± %.1f) -> %s\n', average_wavelength, cfg_PPC_target, cfg_PPC_tol, string_status(is_ppc_ok));
                    fprintf('PI:  %d (目标: %d ± %d) -> %s\n', max_intensity, cfg_PI_target, cfg_PI_tol, string_status(is_pi_ok));
                    fprintf('----------------------------------------------------\n');
                    
                    % 更新图表
                    update_plot_data(processed_files_data, current_spectrum, all_spectra_data, peak_color, avg_color, intensity_color);
                    drawnow;
                    
                    processed_files{end+1} = file;
                end
            end
            pause(1); 
        end
    catch ME
        if strcmp(ME.identifier, 'MATLAB:operatingSystem:Interrupt')
            fprintf('\n程序被用户中断\n');
        else
            rethrow(ME);
        end
    end
end

% 辅助函数：打印状态字符串
function s = string_status(bool_val)
    if bool_val, s = 'OK'; else, s = 'Wait'; end
end

% 播放报警音效
function play_alarm_sound()
    fs = 8192; T = 0.25; t = 0:1/fs:T;
    y = 0.5 * sin(2*pi*880*t) + 0.5 * sin(2*pi*1100*t);
    envelope = [linspace(0,1,100), ones(1,length(t)-200), linspace(1,0,100)];
    sound(y .* envelope, fs);
end

% 初始化图形布局
function initialize_plots()
    global ax_handles
    clf;
    
    % 设置字体和样式
    set(gcf, 'DefaultAxesFontSize', 12);        
    set(gcf, 'DefaultAxesFontName', 'Microsoft YaHei'); 
    set(gcf, 'DefaultAxesLineWidth', 1.5);      
    set(gcf, 'DefaultLineLineWidth', 2.0);      
    
    % 调整间距参数，适应纵向布局
    gap = [0.06, 0.08];    % [垂直间隙, 水平间隙]
    marg_h = [0.05, 0.05]; % [下边距, 上边距]
    marg_w = [0.10, 0.05]; % [左边距, 右边距]
    
    ax_handles = gobjects(8, 1);
    
    % 1. 当前光谱
    ax_handles(1) = subtightplot(4, 2, 1, gap, marg_h, marg_w);
    title('当前光谱', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('波长 (nm)'); ylabel('强度 (a.u.)');
    grid on; ax = ax_handles(1); ax.GridAlpha = 0.4; ax.MinorGridAlpha = 0.2;
    set(ax, 'TickDir', 'out', 'Box', 'off'); 
    hold on;
    
    % 2. 2D 热图
    ax_handles(2) = subtightplot(4, 2, 2, gap, marg_h, marg_w);
    title('全谱演变', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('波长 (nm)');
    set(ax_handles(2), 'TickDir', 'out', 'Box', 'off');
    hold on;
    
    % 3-8. 趋势图设置
    titles = {'PP 趋势 (全时段)', 'PP 趋势 (最近20点)', ...
              'PPC 趋势 (全时段)', 'PPC 趋势 (最近20点)', ...
              'PI 趋势 (全时段)', 'PI 趋势 (最近20点)'};
    ylabels = {'PP (nm)', '', 'PPC (nm)', '', 'PI (a.u.)', ''};
    
    for i = 3:8
        ax_handles(i) = subtightplot(4, 2, i, gap, marg_h, marg_w);
        grid on; ax = ax_handles(i); ax.GridAlpha = 0.4; ax.MinorGridAlpha = 0.2;
        
        set(ax, 'TickDir', 'out', 'Box', 'off'); 
        hold on;
        
        title(titles{i-2}, 'FontSize', 12, 'FontWeight', 'normal');
        
        if ~isempty(ylabels{i-2}), ylabel(ylabels{i-2}, 'FontWeight', 'bold'); end
        
        if i >= 7
            xlabel('反应时间 (min)', 'FontWeight', 'bold');
        end
    end
    
    drawnow;
end

% 更新数据和绘图
function update_plot_data(processed_files_data, current_spectrum, all_spectra_data, peak_color, avg_color, intensity_color)
    global ax_handles
    
    % 1. 更新当前光谱图
    if ~isempty(current_spectrum)
        cla(ax_handles(1));
        wl = current_spectrum(:, 1);
        int = current_spectrum(:, 2);
        
        % 灰色填充区域
        fill_color = [0.2, 0.2, 0.2]; 
        area(ax_handles(1), wl, int, 'FaceColor', fill_color, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
        hold(ax_handles(1), 'on');
        
        plot(ax_handles(1), wl, int, 'Color', 'k', 'LineWidth', 2);
        
        [~, avg_wl, max_int, ~, ~] = analyze_spectrum(current_spectrum);
        [fitted_peak_wl, ~] = fit_peak_enhanced(current_spectrum);
        
        xline(ax_handles(1), fitted_peak_wl, '--', 'Color', peak_color, 'LineWidth', 1.5, 'Alpha', 0.8);
        xline(ax_handles(1), avg_wl, '--', 'Color', avg_color, 'LineWidth', 1.5, 'Alpha', 0.8);
        yline(ax_handles(1), max_int, '--', 'Color', intensity_color, 'LineWidth', 1.5, 'Alpha', 0.8);

        
        text(ax_handles(1), 0.96, 0.94, sprintf('PP: %.1f', fitted_peak_wl), ...
            'Units', 'normalized', 'HorizontalAlignment', 'right', 'FontSize', 16, 'Color', peak_color, 'FontWeight', 'bold');
        text(ax_handles(1), 0.96, 0.82, sprintf('PPC: %.1f', avg_wl), ...
            'Units', 'normalized', 'HorizontalAlignment', 'right', 'FontSize', 16, 'Color', avg_color, 'FontWeight', 'bold');
        text(ax_handles(1), 0.96, 0.70, sprintf('PI: %d', max_int), ...
            'Units', 'normalized', 'HorizontalAlignment', 'right', 'FontSize', 16, 'Color', intensity_color, 'FontWeight', 'bold');
            
        ylim(ax_handles(1), [0, max(int) * 1.15]);
        xlim(ax_handles(1), [400, 1000]);
    end
    
    % 2. 更新2D热图
    if ~isempty(fieldnames(all_spectra_data))
        cla(ax_handles(2));
        files = fieldnames(all_spectra_data);
        num = length(files);
        
        file1 = all_spectra_data.(files{1});
        wl_base = file1(:, 1);
        z_matrix = zeros(num, length(wl_base));
        
        for i = 1:num
            dat = all_spectra_data.(files{i});
            if length(dat) == length(wl_base)
                z_matrix(i, :) = dat(:, 2);
            else
                z_matrix(i, :) = interp1(dat(:,1), dat(:,2), wl_base, 'linear', 0);
            end
        end
        
        imagesc(ax_handles(2), wl_base, 1:num, z_matrix);
        set(ax_handles(2), 'YDir', 'normal');
        colormap(ax_handles(2), 'hot');
        
        xlim(ax_handles(2), [400, 1000]);
        ylim(ax_handles(2), [0.5, num + 0.5]);
        set(ax_handles(2), 'YTickLabel', []);
        
        title(ax_handles(2), sprintf('反应时间: %d min', num), 'FontWeight', 'bold');
    end
    
    % 3. 更新趋势图
    if ~isempty(fieldnames(processed_files_data))
        files = fieldnames(processed_files_data);
        num = length(files);
        
        pp = zeros(num, 1); ppc = zeros(num, 1); pi_val = zeros(num, 1);
        alarms = false(num, 1);
        
        for i = 1:num
            d = processed_files_data.(files{i});
            pp(i) = d.peak; ppc(i) = d.avg; pi_val(i) = d.intensity;
            if isfield(d, 'is_alarm'), alarms(i) = d.is_alarm; end
        end
        alarm_idx = find(alarms);
        
        if num > 1, start_recent = max(1, num-19); else, start_recent = 1; end
        idx_full = 1:num;
        idx_recent = start_recent:num;
        
        plot_beautified(ax_handles(3), idx_full, pp, peak_color, alarm_idx);
        xlim(ax_handles(3), [0.5, num+0.5]); 
        set(ax_handles(3), 'XTickLabel', []);
        
        plot_beautified(ax_handles(4), idx_recent, pp(idx_recent), peak_color, alarm_idx);
        xlim(ax_handles(4), [idx_recent(1)-0.5, idx_recent(end)+0.5]);
        set(ax_handles(4), 'XTick', idx_recent(1):2:idx_recent(end));
        set(ax_handles(4), 'XTickLabel', []);
        
        plot_beautified(ax_handles(5), idx_full, ppc, avg_color, alarm_idx);
        xlim(ax_handles(5), [0.5, num+0.5]);
        set(ax_handles(5), 'XTickLabel', []);
        
        plot_beautified(ax_handles(6), idx_recent, ppc(idx_recent), avg_color, alarm_idx);
        xlim(ax_handles(6), [idx_recent(1)-0.5, idx_recent(end)+0.5]);
        set(ax_handles(6), 'XTick', idx_recent(1):2:idx_recent(end));
        set(ax_handles(6), 'XTickLabel', []);
        
        plot_beautified(ax_handles(7), idx_full, pi_val, intensity_color, alarm_idx);
        xlim(ax_handles(7), [0.5, num+0.5]);
        
        plot_beautified(ax_handles(8), idx_recent, pi_val(idx_recent), intensity_color, alarm_idx);
        xlim(ax_handles(8), [idx_recent(1)-0.5, idx_recent(end)+0.5]);
        set(ax_handles(8), 'XTick', idx_recent(1):2:idx_recent(end));
    end
end

% === 绘图核心函数 ===
function plot_beautified(ax, x, y, main_color, alarm_global_indices)
    cla(ax);
    
    % 1. 绘制连线
    plot(ax, x, y, '-', 'Color', main_color, 'LineWidth', 2.0);
    hold(ax, 'on');
    
    is_alarm = ismember(x, alarm_global_indices);
    is_normal = ~is_alarm;
    
    % 2. 绘制普通点：空心圆
    if any(is_normal)
        plot(ax, x(is_normal), y(is_normal), 'o', ...
             'MarkerSize', 6, ...
             'LineWidth', 1.5, ...
             'MarkerEdgeColor', main_color, ...
             'MarkerFaceColor', 'w'); 
    end
    
    % 3. 绘制报警点：实心圆
    if any(is_alarm)
        plot(ax, x(is_alarm), y(is_alarm), 'o', ...
             'MarkerSize', 8, ...          
             'LineWidth', 1.0, ...
             'MarkerEdgeColor', 'k', ...    
             'MarkerFaceColor', 'k');      
    end
end

% --- 数据处理辅助函数 ---
function data = get_spectral_data(file_path)
    fid = fopen(file_path, 'r', 'n', 'GBK');
    if fid == -1, error('无法打开'); end
    data = []; found = false;
    while ~feof(fid)
        line = fgetl(fid);
        if contains(line, '>>>>>Begin Spectral Data<<<<<'), found = true;
        elseif found
            tmp = str2double(strsplit(line, '\t'));
            if length(tmp) >= 2, data = [data; tmp(1), tmp(2)]; end
        end
    end
    fclose(fid);
end

function d = denoise_data(data)
    d = [data(:,1), sgolayfilt(data(:,2), 3, 11)];
end

function [pk_wl, pk_int] = fit_peak_enhanced(data)
    wl = data(:,1); int = data(:,2);
    sm = sgolayfilt(int, 3, 15);
    [max_i, ~] = max(int);
    [pks, locs] = findpeaks(sm, 'MinPeakHeight', max_i*0.8, 'MinPeakDistance', 15);
    
    if length(pks) >= 2
        wl_pks = wl(locs);
        [min_wl, idx] = min(wl_pks);
        mask = wl >= (min_wl-12) & wl <= (min_wl+12);
        if sum(mask) > 4
            [p, S, mu] = polyfit(wl(mask), int(mask), 3);
            x_ft = linspace(min(wl(mask)), max(wl(mask)), 100);
            y_ft = polyval(p, x_ft, S, mu);
            [pk_int, mi] = max(y_ft); pk_wl = x_ft(mi);
        else
            pk_wl = min_wl; pk_int = pks(idx);
        end
    else
        mask = int >= max_i*0.95;
        [p, S, mu] = polyfit(wl(mask), int(mask), 3);
        x_ft = linspace(min(wl(mask)), max(wl(mask)), 100);
        y_ft = polyval(p, x_ft, S, mu);
        [pk_int, mi] = max(y_ft); pk_wl = x_ft(mi);
    end
    pk_wl = round(pk_wl, 1);
end

function [pk, avg, mx, l, r] = analyze_spectrum(data)
    wl = data(:,1); int = data(:,2);
    [mx, idx] = max(int); pk = wl(idx);
    half = mx/2;
    l_idx = find(int(1:idx) < half, 1, 'last'); if isempty(l_idx), l_idx=1; end
    r_idx = find(int(idx:end) < half, 1, 'first') + idx - 1; if isempty(r_idx), r_idx=length(int); end
    l = wl(l_idx); r = wl(r_idx); avg = (l+r)/2;
    pk = round(pk,1); avg = round(avg,1); mx = round(mx);
end

function f = get_new_csv_filename(d, b)
    [~, n, e] = fileparts(b); c = 1;
    f = fullfile(d, sprintf('%s_%d%s', n, c, e));
    while exist(f, 'file'), c=c+1; f = fullfile(d, sprintf('%s_%d%s', n, c, e)); end
end

function write_csv_data(fn, rd)
    mode = 'a'; if ~exist(fn, 'file'), mode = 'w'; end
    fid = fopen(fn, mode, 'n', 'UTF-8');
    for i=1:length(rd)
        if i>1, fprintf(fid, ','); end
        if ischar(rd{i}), fprintf(fid, '"%s"', rd{i});
        elseif isinteger(rd{i}), fprintf(fid, '%d', rd{i});
        else, fprintf(fid, '%.1f', rd{i}); end
    end
    fprintf(fid, '\n'); fclose(fid);
end

function h = subtightplot(m, n, p, gap, mh, mw)
    if isscalar(gap), gap=[gap gap]; end
    if isscalar(mh), mh=[mh mh]; end
    if isscalar(mw), mw=[mw mw]; end
    w = (1-sum(mw)-(n-1)*gap(2))/n; h = (1-sum(mh)-(m-1)*gap(1))/m;
    if isscalar(p)
        r = m-floor((p-1)/n); c = mod(p-1,n)+1;
        pos = [mw(1)+(c-1)*(w+gap(2)), mh(1)+(r-1)*(h+gap(1)), w, h];
    else
        pos = [0.1 0.1 0.8 0.8]; 
    end
    h = axes('Position', pos);
end