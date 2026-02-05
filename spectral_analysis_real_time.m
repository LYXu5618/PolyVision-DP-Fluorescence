function spectral_analysis_real_time()
% 全集成的光谱分析与报警系统
% 功能：自动触发硬件采集 -> 监测新文件 -> 降噪与拟合 -> 目标比对报警

    % 全局变量用于存储绘图数据和图形句柄
    global processed_files_data current_spectrum ax_handles all_spectra_data
    
    % --- 初始化数据存储 ---
    processed_files_data = struct();
    current_spectrum = [];
    ax_handles = [];
    all_spectra_data = struct(); 
    
    % =========================================================================
    % --- 用户配置区域：报警阈值设置 (目标值 ± 允许偏差) ---
    % =========================================================================
    cfg_PP_target  = 669.4;  cfg_PP_tol  = 0.5;   % 峰值波长 (Peak Position)
    cfg_PPC_target = 690.7;  cfg_PPC_tol = 0.5;   % 峰值中心 (Peak Position Center)
    cfg_PI_target  = 10262;  cfg_PI_tol  = 500;   % 峰值强度 (Peak Intensity)
    % =========================================================================
    
    % --- 运行参数配置 ---
    input_directory = '*****'; % 请修改为你的实际路径
    acquisition_interval = 60; % 自动采集间隔 (秒)
    
    % 检查路径
    if ~exist(input_directory, 'dir')
        error('输入目录不存在，请检查路径设置: %s', input_directory);
    end
    
    % 创建图形窗口
    fig = figure('Name', 'Real-time Spectral Analysis System', 'NumberTitle', 'off', ...
                 'Position', [50, 50, 1100, 950], 'Color', 'w');
    
    % 生成 CSV 文件名
    folders = strsplit(input_directory, '\');
    last_folder = folders{end};
    base_name = sprintf('spectral_analysis_results_%s.csv', last_folder);
    output_csv = get_new_csv_filename(input_directory, base_name);
    
    processed_files = {}; 
    last_acquisition_time = now; 
    
    % 定义绘图颜色
    peak_color = [0.85, 0.1, 0.1];      % 红色
    avg_color = [0.1, 0.4, 0.8];       % 蓝色
    intensity_color = [0.95, 0.6, 0];   % 橙色
    
    % 初始化绘图布局
    initialize_plots();
    
    fprintf('实时光谱系统已启动\n');
    fprintf('目标设定: PP=%.1f, PPC=%.1f, PI=%.0f\n', cfg_PP_target, cfg_PPC_target, cfg_PI_target);
    
    % 执行初始采集
    fprintf('执行启动初始硬件采集...\n');
    hardware_acquisition_cycle();
    last_acquisition_time = now;
    
    try
        while true
            % --- 1. 硬件采集逻辑 ---
            current_time = now;
            time_since = (current_time - last_acquisition_time) * 86400; % 转换秒
            
            if time_since >= acquisition_interval
                fprintf('达到采集时间点，开始采集...\n');
                hardware_acquisition_cycle();
                last_acquisition_time = current_time;
                pause(1.5); % 等待文件写入磁盘
            end
            
            % --- 2. 文件轮询逻辑 ---
            files = dir(fullfile(input_directory, '*.txt'));
            for i = 1:length(files)
                file = files(i).name;
                if ~ismember(file, processed_files)
                    fprintf('正在分析新文件: %s\n', file);
                    try
                        % 读取与降噪
                        data = get_spectral_data(fullfile(input_directory, file));
                        denoised = denoise_data(data);
                        current_spectrum = denoised;
                        all_spectra_data.(matlab.lang.makeValidName(file)) = denoised;
                        
                        % 参数计算 (拟合 + 分析)
                        [f_pp, f_pi] = fit_peak_enhanced(denoised);
                        [~, f_ppc, ~, l_wl, r_wl] = analyze_spectrum(denoised);
                        
                        % 报警逻辑判断
                        is_pp_ok = abs(f_pp - cfg_PP_target) <= cfg_PP_tol;
                        is_ppc_ok = abs(f_ppc - cfg_PPC_target) <= cfg_PPC_tol;
                        is_pi_ok = abs(f_pi - cfg_PI_target) <= cfg_PI_tol;
                        
                        is_alarm = is_pp_ok && is_ppc_ok && is_pi_ok;
                        
                        if is_alarm
                            fprintf('数据已达标! 正在播放提示音...\n');
                            play_alarm_sound();
                        end
                        
                        % 存储结果
                        processed_files_data.(matlab.lang.makeValidName(file)) = struct(...
                            'peak', f_pp, 'avg', f_ppc, 'intensity', f_pi, 'is_alarm', is_alarm);
                        
                        % 写入 CSV
                        write_csv_data(output_csv, {file, f_pp, f_ppc, f_pi, l_wl, r_wl});
                        
                        % 控制台打印状态
                        fprintf('PP: %.1f [%s] | PPC: %.1f [%s] | PI: %.0f [%s]\n', ...
                            f_pp, string_status(is_pp_ok), f_ppc, string_status(is_ppc_ok), f_pi, string_status(is_pi_ok));
                        
                        % 更新绘图
                        update_plot_data(processed_files_data, current_spectrum, all_spectra_data, ...
                            peak_color, avg_color, intensity_color);
                        drawnow;
                        
                        processed_files{end+1} = file;
                    catch ME
                        fprintf('文件处理出错: %s\n', ME.message);
                    end
                end
            end
            
            % 倒计时更新显示
            wait_time = max(0, acquisition_interval - (now - last_acquisition_time)*86400);
            fprintf('距离下次采集还剩: %.0fs | 已分析: %d 个文件\r', wait_time, length(processed_files));
            pause(1);
        end
    catch ME
        if ~strcmp(ME.identifier, 'MATLAB:operatingSystem:Interrupt'), rethrow(ME); end
        fprintf('\n[STOP] 程序已被用户手动停止。\n');
    end
end

% --- 硬件控制子函数 ---
function hardware_acquisition_cycle()
    s1 = []; s = [];
    try
        % 控制快门 (COM23)
        s1 = serialport('COM23', 9600);
        configureTerminator(s1, "CR");
        write(s1, uint8(hex2dec({'A0'; '01'; '01'; 'A2'})), "uint8"); % 打开光源
        pause(0.5);
        
        % 控制光谱仪 (COM22)
        s = serialport('COM22', 115200);
        configureTerminator(s, "CR");
        write(s, uint8(hex2dec({'3A'; '01'; '01'})), "uint8"); % 发送采集触发
        pause(1.0); % 等待曝光完成
        
        write(s, uint8(hex2dec({'3A'; '01'; '00'})), "uint8"); % 停止触发信号
        write(s1, uint8(hex2dec({'A0'; '01'; '00'; 'A1'})), "uint8"); % 关闭光源
        fprintf('硬件采集指令执行成功\n');
    catch ME
        fprintf('硬件串口连接失败: %s\n', ME.message);
    end
    if ~isempty(s), delete(s); end
    if ~isempty(s1), delete(s1); end
end

% --- 绘图初始化 ---
function initialize_plots()
    global ax_handles
    clf;
    set(gcf, 'DefaultAxesFontSize', 10, 'DefaultAxesFontName', 'Microsoft YaHei');
    gap = [0.07, 0.07]; marg_h = [0.06, 0.05]; marg_w = [0.08, 0.05];
    ax_handles = gobjects(8, 1);
    
    titles = {'Current Spectrum', 'Spectral Evolution (2D Heatmap)', ...
              'PP Trend (Full)', 'PP Trend (Recent 20)', ...
              'PPC Trend (Full)', 'PPC Trend (Recent 20)', ...
              'PI Trend (Full)', 'PI Trend (Recent 20)'};
    
    for i = 1:8
        ax_handles(i) = subtightplot(4, 2, i, gap, marg_h, marg_w);
        title(titles{i}); grid on; hold on;
        if i == 7, xlabel('Time (Counts)'); ylabel('PI (a.u.)'); end
        if i == 5, ylabel('PPC (nm)'); end
        if i == 3, ylabel('PP (nm)'); end
    end
end

% --- 数据绘图更新 ---
function update_plot_data(processed_files_data, current_spectrum, all_spectra_data, p_col, a_col, i_col)
    global ax_handles
    
    % 1. 更新当前光谱
    if ~isempty(current_spectrum)
        cla(ax_handles(1));
        wl = current_spectrum(:,1); int = current_spectrum(:,2);
        area(ax_handles(1), wl, int, 'FaceColor', [0.8 0.8 0.8], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(ax_handles(1), wl, int, 'k', 'LineWidth', 1.5);
        
        [f_pp, f_pi] = fit_peak_enhanced(current_spectrum);
        [~, f_ppc] = analyze_spectrum(current_spectrum);
        
        xline(ax_handles(1), f_pp, '--', 'Color', p_col);
        xline(ax_handles(1), f_ppc, '--', 'Color', a_col);
        
        text(ax_handles(1), 0.95, 0.9, sprintf('PP: %.1f', f_pp), 'Units', 'normalized', 'Color', p_col, 'FontWeight', 'bold', 'FontSize', 12, 'HorizontalAlignment', 'right');
        text(ax_handles(1), 0.95, 0.8, sprintf('PPC: %.1f', f_ppc), 'Units', 'normalized', 'Color', a_col, 'FontWeight', 'bold', 'FontSize', 12, 'HorizontalAlignment', 'right');
        text(ax_handles(1), 0.95, 0.7, sprintf('PI: %.0f', f_pi), 'Units', 'normalized', 'Color', i_col, 'FontWeight', 'bold', 'FontSize', 12, 'HorizontalAlignment', 'right');
        xlim(ax_handles(1), [400, 900]); ylim(ax_handles(1), [0, max(int)*1.2]);
    end
    
    % 2. 更新2D热图
    if ~isempty(fieldnames(all_spectra_data))
        cla(ax_handles(2));
        fnames = fieldnames(all_spectra_data);
        n_files = length(fnames);
        wl_base = all_spectra_data.(fnames{1})(:,1);
        z_matrix = zeros(n_files, length(wl_base));
        for i = 1:n_files
            dat = all_spectra_data.(fnames{i});
            z_matrix(i,:) = interp1(dat(:,1), dat(:,2), wl_base, 'linear', 0);
        end
        imagesc(ax_handles(2), wl_base, 1:n_files, z_matrix);
        colormap(ax_handles(2), 'hot'); set(ax_handles(2), 'YDir', 'normal');
        xlim(ax_handles(2), [400, 900]); ylabel(ax_handles(2), 'File Index');
    end
    
    % 3. 趋势图处理
    if ~isempty(fieldnames(processed_files_data))
        f_names = fieldnames(processed_files_data);
        n = length(f_names);
        pp = zeros(n,1); ppc = zeros(n,1); pi_v = zeros(n,1); alms = false(n,1);
        for k=1:n
            d = processed_files_data.(f_names{k});
            pp(k)=d.peak; ppc(k)=d.avg; pi_v(k)=d.intensity; alms(k)=d.is_alarm;
        end
        
        alm_idx = find(alms);
        idx_full = 1:n;
        idx_rec = max(1, n-19):n;
        
        plot_beautified(ax_handles(3), idx_full, pp, p_col, alm_idx);
        plot_beautified(ax_handles(4), idx_rec, pp(idx_rec), p_col, alm_idx);
        plot_beautified(ax_handles(5), idx_full, ppc, a_col, alm_idx);
        plot_beautified(ax_handles(6), idx_rec, ppc(idx_rec), a_col, alm_idx);
        plot_beautified(ax_handles(7), idx_full, pi_v, i_col, alm_idx);
        plot_beautified(ax_handles(8), idx_rec, pi_v(idx_rec), i_col, alm_idx);
    end
end

% --- 趋势绘图美化 (报警点显示黑色实心) ---
function plot_beautified(ax, x, y, col, alm_global_idx)
    cla(ax);
    plot(ax, x, y, '-', 'Color', col, 'LineWidth', 1); hold(ax, 'on');
    is_alm = ismember(x, alm_global_idx);
    if any(~is_alm)
        plot(ax, x(~is_alm), y(~is_alm), 'o', 'MarkerEdgeColor', col, 'MarkerFaceColor', 'w', 'MarkerSize', 4);
    end
    if any(is_alm)
        plot(ax, x(is_alm), y(is_alm), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    end
end

% --- 数据读取辅助函数 ---
function data = get_spectral_data(file_path)
    fid = fopen(file_path, 'r', 'n', 'GBK');
    if fid == -1, error('无法打开文件'); end
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
    [max_i, midx] = max(int);
    % 选取峰值附近区域进行三阶多项式拟合
    mask = wl >= (wl(midx)-10) & wl <= (wl(midx)+10);
    if sum(mask) > 5
        p = polyfit(wl(mask), int(mask), 3);
        x_ft = linspace(min(wl(mask)), max(wl(mask)), 200);
        y_ft = polyval(p, x_ft);
        [pk_int, mi] = max(y_ft); pk_wl = x_ft(mi);
    else
        pk_wl = wl(midx); pk_int = max_i;
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
    % 子图间距控制函数
    if isscalar(gap), gap=[gap gap]; end
    if isscalar(mh), mh=[mh mh]; end
    if isscalar(mw), mw=[mw mw]; end
    w = (1-sum(mw)-(n-1)*gap(2))/n; h = (1-sum(mh)-(m-1)*gap(1))/m;
    row = m-floor((p-1)/n); col = mod(p-1,n)+1;
    pos = [mw(1)+(col-1)*(w+gap(2)), mh(1)+(row-1)*(h+gap(1)), w, h];
    h = axes('Position', pos);
end

function s = string_status(val)
    if val, s = 'OK'; else, s = 'Wait'; end
end

function play_alarm_sound()
    % 播放一个简单的 880Hz 提示音
    fs = 8192; t = 0:1/fs:0.15;
    y = 0.2 * sin(2*pi*880*t);
    sound(y, fs);
end