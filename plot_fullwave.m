function plot_fullwave(pc, pf, n_wave)
% plot_fullwave -- 全波波形随 dT 的对比（只画 fullwave 相关）
%
%   plot_fullwave()      % n = 0 / 1 / 2 各出一张
%   plot_fullwave(1)     % 只出 n = 1
%
%   输出：fig_全波波形_n<阶数>.png / .pdf
%
%   图 fig_全波波形对比.*：z = 3 m 处四通道（P_axis / P_wall / H_axis / H_wall）
%   三组 dT（0 / +20 / +40 K）叠画，t 窗 1e-3 ~ 3e-3 s（放大看主脉冲）
%
%   约定（AGENTS.md）：SI 单位（时间 s）；不归一化；深色 + 线型区分 dT；
%   图内文字中文；线性轴；四图等大（tiledlayout）；图例置顶。

if nargin < 1, pc = []; end
if nargin < 2, pf = []; end
if nargin < 3 || isempty(n_wave), n_wave = [0, 1, 2]; end

here = fileparts(mfilename('fullpath'));
addpath(here);
hdr0 = main(pc, pf, 0, 0, [], 1, 2e3, 2e3, 1);   % 取实际 pc/pf
pcv = hdr0.pc;  pfv = hdr0.pf;
outdir = fullfile(here, 'result', sprintf('pc%gMPa_pf%gMPa', pcv/1e6, pfv/1e6));
if exist(outdir, 'dir') ~= 7, mkdir(outdir); end

dT_list = [0, 5, 10];
nDT     = numel(dT_list);
z_rec   = 3;

CLR = [0.000 0.150 0.550;
       0.750 0.100 0.100;
       0.100 0.450 0.100];
LS  = {'-', '--', '-.'};
LW  = 0.9;  FS = 10.5;
set(0, 'DefaultAxesFontSize', FS);
set(0, 'DefaultTextFontSize', FS);
set(0, 'DefaultAxesLineWidth', 0.6);
dT_lab = cellfun(@(x) sprintf('%+gK', x), num2cell(dT_list), 'UniformOutput', false);

% ---------------- 算全波（按 n_wave 循环） ----------------
for n = n_wave
S = struct();
for k = 1:nDT
    fprintf('[n=%d %d/%d] fullwave dT=%+g ...\n', n, k, nDT, dT_list(k));
    o = fullwave(pc, pf, dT_list(k), n, []);
    iz = find(o.z <= z_rec, 1, 'last');
    S.time = o.time;
    S.P_axis{k} = o.P_axis_r(:,iz) + o.P_axis_d(:,iz);
    S.P_wall{k} = o.P_wall_r(:,iz) + o.P_wall_d(:,iz);
    S.H_axis{k} = o.H_axis_r(:,iz) + o.H_axis_d(:,iz);
    S.H_wall{k} = o.H_wall_r(:,iz) + o.H_wall_d(:,iz);
    S.z_used = o.z(iz);
    S.f_src  = o.f_src;
    S.pc = o.pc;  S.pf = o.pf;
    clear o
end

% ---------------- 绘图 ----------------
fig = figure('Visible', 'off', 'Position', [80 80 1000 820]);
tl  = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

tt = S.time;                         % [s]（国际单位）
it = tt >= 1e-3 & tt <= 4e-3;        % 只看 1~4 ms 区间

names = {'P_axis', 'P_wall', 'H_axis', 'H_wall'};
labs  = {'广义压力（井轴）', '井壁声压', '井轴温度', '井壁温度'};
ylabs = {'p_{axis}', 'p_{wall}', 'a_{axis}', 'a_{wall}'};
for i = 1:4
    nexttile; hold on; box on; grid on; grid minor;
    for k = 1:nDT
        plot(tt(it), S.(names{i}){k}(it), 'LineStyle', LS{k}, 'Color', CLR(k,:), 'LineWidth', LW);
    end
    xlabel('时间 t [s]'); ylabel(ylabs{i});
    title(sprintf('(%c) %s', 'a'+i-1, labs{i}));
    hold off;
    if i == 1
        legend(cellfun(@(s) ['\DeltaT=' s], dT_lab, 'UniformOutput', false), ...
               'Location', 'northoutside', 'Orientation', 'horizontal', ...
               'NumColumns', 3, 'FontSize', FS-3);
    end
end
title(tl, sprintf('全波波形：pc = %g MPa，pf = %g MPa，n = %d，f_{src} = %g Hz，z = %.2f m', ...
      S.pc/1e6, S.pf/1e6, n, S.f_src, S.z_used));

set(fig, 'PaperPositionMode', 'auto');
base = fullfile(outdir, sprintf('fig_全波波形_n%d', n));
print(fig, [base '.png'], '-dpng', '-r300');
print(fig, [base '.pdf'], '-dpdf', '-bestfit');
close(fig);
fprintf('saved: %s\n', [base '.png']);
end   % for n = n_wave
end