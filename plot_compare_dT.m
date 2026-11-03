function plot_compare_dT()
% plot_compare_dT -- 不同 dT 的导波对比 + 全波（调用 main / fullwave）
%
%   图 1  fig_导波频散对比.*    斯通利波（深色）+ 伪瑞利波（深红）× 3 组 dT
%   图 2  fig_伪瑞利波激发对比.*  伪瑞利波：井壁/井轴 × 声压/温度
%   图 3  fig_全波波形对比.*      z = 3 m 四通道波形
%
%   约定（AGENTS.md）：不归一化；深色 + 线型区分 dT；图内文字中文；线性轴；hold off。

here = fileparts(mfilename('fullpath'));
addpath(here);
outdir = fullfile(here, 'fig');
if exist(outdir, 'dir') ~= 7, mkdir(outdir); end

dT_list = [0, 20, 40];
nDT     = numel(dT_list);
F_LO = 2e3;  F_HI = 4e4;  NF = 256;
z_rec = 3;

CLR = [0.000 0.150 0.550;
       0.750 0.100 0.100;
       0.100 0.450 0.100];
LS  = {'-', '--', '-.'};
LW  = 0.9;  FS = 10.5;
set(0, 'DefaultAxesFontSize', FS);
set(0, 'DefaultTextFontSize', FS);
set(0, 'DefaultAxesLineWidth', 0.6);
dT_lab = cellfun(@(x) sprintf('%+gK', x), num2cell(dT_list), 'UniformOutput', false);

% ============ 1) 导波：斯通利 + 伪瑞利 ============
Rst = cell(1, nDT);  Rpr = cell(1, nDT);
for k = 1:nDT
    fprintf('[1/3] Stoneley dT=%+g ...\n', dT_list(k));
    Rst{k} = main([], [], dT_list(k), 0, [], 1, F_LO, F_HI, NF);
    fprintf('[1/3] pseudo-Rayleigh dT=%+g ...\n', dT_list(k));
    Rpr{k} = main([], [], dT_list(k), 0, [], 0, F_LO, F_HI, NF);
end
info = sprintf('pc = %g MPa, pf = %g MPa, n = 0', Rst{1}.pc/1e6, Rst{1}.pf/1e6);
f_ch = Rst{1}.freqs;

fig1 = figure('Visible', 'off', 'Position', [80 80 800 520]);
hold on; box on; grid on; grid minor;
for k = 1:nDT
    plot(f_ch, Rst{k}.Vph, 'LineStyle', LS{k}, 'Color', CLR(1,:), 'LineWidth', LW);
end
for k = 1:nDT
    plot(f_ch, Rpr{k}.Vph, 'LineStyle', LS{k}, 'Color', CLR(2,:), 'LineWidth', LW);
end
xlabel('频率 f [Hz]'); ylabel('相速度 [m/s]');
title({'(a) 斯通利波（深蓝）与伪瑞利波（深红）', info});
leg = [cellfun(@(s) ['斯通利 ' s], dT_lab, 'UniformOutput', false), ...
       cellfun(@(s) ['伪瑞利 ' s], dT_lab, 'UniformOutput', false)];
legend(leg, 'Location', 'northoutside', 'Orientation', 'horizontal', 'NumColumns', 3, 'FontSize', FS-3);
xlim([f_ch(1), f_ch(end)]);
hold off;
saveas_print(fig1, fullfile(outdir, 'fig_导波频散对比'));

% ============ 2) 伪瑞利波激发 ============
fig2 = figure('Visible', 'off', 'Position', [80 80 1200 820]);
tl2 = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
G = { Rst{1}.ap_w, Rst{2}.ap_w, Rst{3}.ap_w, '(a) 斯通利：井壁声压激发';
      Rst{1}.at_w, Rst{2}.at_w, Rst{3}.at_w, '(b) 斯通利：井壁温度激发';
      Rpr{1}.ap_w, Rpr{2}.ap_w, Rpr{3}.ap_w, '(c) 伪瑞利：井壁声压激发';
      Rpr{1}.at_w, Rpr{2}.at_w, Rpr{3}.at_w, '(d) 伪瑞利：井壁温度激发' };
for i = 1:4
    nexttile; hold on; box on; grid on; grid minor;
    for k = 1:nDT
        v = abs(G{i,k});
        if all(~isfinite(v)), continue; end
        plot(f_ch, v, 'LineStyle', LS{k}, 'Color', CLR(k,:), 'LineWidth', LW);
    end
    xlabel('频率 f [Hz]'); ylabel('|激发|'); title(G{i,4});
    xlim([f_ch(1), f_ch(end)]); hold off;
    if i == 1
        legend(cellfun(@(s) ['\DeltaT=' s], dT_lab, 'UniformOutput', false), ...
               'Location', 'northoutside', 'Orientation', 'horizontal', 'NumColumns', 3, 'FontSize', FS-3);
    end
end
title(tl2, info);
saveas_print(fig2, fullfile(outdir, 'fig_伪瑞利波激发对比'));

% ============ 3) 全波单道 ============
n = Rst{1}.n;
S = struct();
for k = 1:nDT
    fprintf('[3/3] fullwave dT=%+g ...\n', dT_list(k));
    o = fullwave([], [], dT_list(k), n, []);
    iz = find(o.z <= z_rec, 1, 'last');
    S.time = o.time;
    S.P_axis{k} = o.P_axis_r(:,iz) + o.P_axis_d(:,iz);
    S.P_wall{k} = o.P_wall_r(:,iz) + o.P_wall_d(:,iz);
    S.H_axis{k} = o.H_axis_r(:,iz) + o.H_axis_d(:,iz);
    S.H_wall{k} = o.H_wall_r(:,iz) + o.H_wall_d(:,iz);
    S.z_used = o.z(iz);  S.f_src = o.f_src;
    clear o
end

fig3 = figure('Visible', 'off', 'Position', [80 80 1000 820]);
tl3 = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
tt = S.time;                       % [s]（国际单位）
it = tt >= 1e-3 & tt <= 3e-3;     % 只看 1~3 ms 区间（放大显示）
names = {'P_axis', 'P_wall', 'H_axis', 'H_wall'};
labs  = {'P_{axis}', 'P_{wall}', 'H_{axis}', 'H_{wall}'};
for i = 1:4
    nexttile; hold on; box on; grid on; grid minor;
    for k = 1:nDT
        plot(tt(it), S.(names{i}){k}(it), 'LineStyle', LS{k}, 'Color', CLR(k,:), 'LineWidth', LW);
    end
    xlabel('时间 t [s]'); ylabel(labs{i});
    title(sprintf('(%c) %s，z = %.2f m', 'a'+i-1, labs{i}, S.z_used));
    hold off;
    if i == 1
        legend(cellfun(@(s) ['\DeltaT=' s], dT_lab, 'UniformOutput', false), ...
               'Location', 'northoutside', 'Orientation', 'horizontal', 'NumColumns', 3, 'FontSize', FS-3);
    end
end
title(tl3, sprintf('全波波形：%s，f_{src} = %g Hz', info, S.f_src));
saveas_print(fig3, fullfile(outdir, 'fig_全波波形对比'));

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