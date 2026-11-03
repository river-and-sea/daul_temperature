function plot_modes()
% plot_modes -- 调用 main 绘制各类井孔导波（n 与 is_st 不同）随 dT 的对比
%
%   导波分类（其余导波只改 n，is_st 传 0）：
%     n = 0, is_st = 1 : 斯通利波
%     n = 0, is_st = 0 : 伪瑞利波
%     n = 1, is_st = 0 : 弯曲波
%     n = 2, is_st = 0 : 螺旋波
%
%   每张图 2x2：相速度 / 井壁声压激发 / 井壁温度激发 / 衰减
%   频率轴 2e3 ~ 4e4 Hz
%
%   约定（AGENTS.md）：不归一化；不用浅色；dT 用深色 + 线型区分；图内文字中文；线性轴。

here = fileparts(mfilename('fullpath'));
addpath(here);
outdir = fullfile(here, 'fig');
if exist(outdir, 'dir') ~= 7, mkdir(outdir); end

dT_list = [0, 20, 40];
nDT     = numel(dT_list);
F_LO    = 2e3;
F_HI    = 4e4;
NF      = 256;

CLR = [0.000 0.150 0.550;
       0.750 0.100 0.100;
       0.100 0.450 0.100];
LS  = {'-', '--', '-.'};
LW  = 0.9;
FS  = 10.5;
set(0, 'DefaultAxesFontSize', FS);
set(0, 'DefaultTextFontSize', FS);
set(0, 'DefaultAxesLineWidth', 0.6);

dT_lab = cellfun(@(x) sprintf('%+gK', x), num2cell(dT_list), 'UniformOutput', false);

% ============ 图 1：斯通利波 + 伪瑞利波（n=0）============
Rst = cell(1, nDT);   Rpr = cell(1, nDT);
for k = 1:nDT
    fprintf('[1/4] n=0 is_st=1 (Stoneley) dT=%+g ...\n', dT_list(k));
    Rst{k} = main([], [], dT_list(k), 0, [], 1, F_LO, F_HI, NF);
    fprintf('[1/4] n=0 is_st=0 (pseudo-Rayleigh) dT=%+g ...\n', dT_list(k));
    Rpr{k} = main([], [], dT_list(k), 0, [], 0, F_LO, F_HI, NF);
end
info = sprintf('pc = %g MPa, pf = %g MPa, n = 0', Rst{1}.pc/1e6, Rst{1}.pf/1e6);

fig1 = figure('Visible', 'off', 'Position', [80 80 1200 900]);
Q1   = {'Vph', 'Vg', 'ap_w', 'at_w', 'atten_mode'};
ttl1 = {'(a) 相速度', '(b) 群速度', '(c) 井壁声压激发', '(d) 井壁温度激发', '(e) 衰减'};
for iq = 1:5
    subplot(2, 3, iq); hold on; box on; grid on; grid minor;
    for k = 1:nDT
        plot(Rst{k}.freqs, abs(Rst{k}.(Q1{iq})), 'LineStyle', LS{k}, 'Color', CLR(1,:), 'LineWidth', LW);
    end
    for k = 1:nDT
        plot(Rpr{k}.freqs, abs(Rpr{k}.(Q1{iq})), 'LineStyle', LS{k}, 'Color', CLR(2,:), 'LineWidth', LW);
    end
    xlabel('频率 f [Hz]'); ylabel(Q1{iq}); title(ttl1{iq});
    xlim([Rst{1}.freqs(1), Rst{1}.freqs(end)]);      % 只画有数据覆盖的频率范围
    hold off;
    if iq == 1
        leg = [cellfun(@(s) ['斯通利 ' s], dT_lab, 'UniformOutput', false), ...
               cellfun(@(s) ['伪瑞利 ' s], dT_lab, 'UniformOutput', false)];
        legend(leg, 'Location', 'southoutside', 'Orientation', 'horizontal', ...
               'NumColumns', 2, 'FontSize', FS-3);
    end
end
sgtitle(info);
saveas_print(fig1, fullfile(outdir, 'fig_导波频散对比'));

% ============ 图 2~4：其他导波（只改 n，is_st=0）============
cfg = { 1, '弯曲波', 'fig_弯曲波频散'; ...
        2, '螺旋波', 'fig_螺旋波频散' };

for ic = 1:size(cfg,1)
    nn    = cfg{ic,1};
    cname = cfg{ic,2};
    fname = cfg{ic,3};

    RR = cell(1, nDT);
    for k = 1:nDT
        fprintf('[%d/4] n=%d is_st=0 dT=%+g ...\n', ic+1, nn, dT_list(k));
        RR{k} = main([], [], dT_list(k), nn, [], 0, F_LO, F_HI, NF);
    end

    fig = figure('Visible', 'off', 'Position', [80 80 1200 900]);
    Q   = {'Vph', 'Vg', 'ap_w', 'at_w', 'atten_mode'};
    ttl = {'(a) 相速度', '(b) 群速度', '(c) 井壁声压激发', '(d) 井壁温度激发', '(e) 衰减'};
    for iq = 1:5
        subplot(2, 3, iq); hold on; box on; grid on; grid minor;
        for k = 1:nDT
            v = abs(RR{k}.(Q{iq}));
            if all(~isfinite(v)), continue; end
            plot(RR{k}.freqs, v, 'LineStyle', LS{k}, 'Color', CLR(k,:), 'LineWidth', LW);
        end
        xlabel('频率 f [Hz]'); ylabel(Q{iq}); title(ttl{iq});
        xlim([RR{1}.freqs(1), RR{1}.freqs(end)]);
        hold off;
        if iq == 1
            legend(cellfun(@(s) ['\DeltaT=' s], dT_lab, 'UniformOutput', false), ...
                   'Location', 'southoutside', 'Orientation', 'horizontal', ...
                   'NumColumns', 3, 'FontSize', FS-3);
        end
    end
    sgtitle(sprintf('pc = %g MPa, pf = %g MPa, n = %d', RR{1}.pc/1e6, RR{1}.pf/1e6, nn));
    saveas_print(fig, fullfile(outdir, fname));
end

fprintf('done: %s\n', outdir);
end

% ======================================================================
function saveas_print(fig, base)
set(fig, 'PaperPositionMode', 'auto');
print(fig, [base '.png'], '-dpng', '-r300');
print(fig, [base '.pdf'], '-dpdf', '-bestfit');
close(fig);
fprintf('  saved: %s.png / .pdf\n', base);
end