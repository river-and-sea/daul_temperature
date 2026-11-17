function plot_guidewave(pc, pf)
% plot_guidewave -- 各类井孔导波（不同 n 与 is_st）随 dT 的对比
%
%   plot_guidewave()            % 默认 pc=65 MPa, pf=5 MPa
%   plot_guidewave(45e6, 5e6)   % 指定围压/孔压（Pa）
%
%   输出目录：result/pc<pc>MPa_pf<pf>MPa/
%
%   每个 n 出两张图：
%     图 A  <模式>_频散：相速度 + 衰减（2 子图，图例放右下角）
%     图 B  <模式>_激发：广义压力（井轴）/ 井轴温度 / 井壁声压 / 井壁温度（2x2）
%
%   导波：n=0 is_st=1 斯通利；n=0 is_st=0 伪瑞利；n=1 弯曲；n=2 螺旋
%   约定：SI 单位；不归一化；深色 + 线型区分 dT；图内中文；线性轴；图例不挡数据。

if nargin < 1, pc = [0]; end
if nargin < 2, pf = [0]; end

here = fileparts(mfilename('fullpath'));
addpath(here);

% 输出目录（按 pc/pf 分子目录）
hdr0 = main(pc, pf, 0, 0, [], 1, 2e3, 2e3, 1);   % 只为取实际 pc/pf 数值
pcv = hdr0.pc;  pfv = hdr0.pf;
outdir = fullfile(here, 'result', sprintf('pc%gMPa_pf%gMPa', pcv/1e6, pfv/1e6));
if exist(outdir, 'dir') ~= 7, mkdir(outdir); end

dT_list = [0, 20, 40];
nDT     = numel(dT_list);
F_LO = 2e2;  F_HI = 4e4;  NF = 256;

CLR = [0.000 0.150 0.550;
       0.750 0.100 0.100;
       0.100 0.450 0.100];
LS  = {'-', '--', '-.'};
LW  = 0.9;  FS = 10.5;
set(0, 'DefaultAxesFontSize', FS);
set(0, 'DefaultTextFontSize', FS);
set(0, 'DefaultAxesLineWidth', 0.6);
dT_lab = cellfun(@(x) sprintf('%+gK', x), num2cell(dT_list), 'UniformOutput', false);

% 模式表：n, is_st, 中文名, 文件基名
modes = { 0, 1, '斯通利波', 'fig_斯通利波';
          0, 0, '伪瑞利波', 'fig_伪瑞利波';
          1, 0, '弯曲波',   'fig_弯曲波';
          2, 0, '螺旋波',   'fig_螺旋波' };

for im = 1:size(modes,1)
    nn    = modes{im,1};
    isst  = modes{im,2};
    cname = modes{im,3};
    base  = modes{im,4};

    R = cell(1, nDT);
    for k = 1:nDT
        fprintf('[%d/%d] n=%d is_st=%d dT=%+g ...\n', im, size(modes,1), nn, isst, dT_list(k));
        R{k} = main(pc, pf, dT_list(k), nn, [], isst, F_LO, F_HI, NF);
    end
    f_ch = R{1}.freqs;
    hdr  = sprintf('pc = %g MPa，pf = %g MPa，n = %d', R{1}.pc/1e6, R{1}.pf/1e6, nn);

    % ============ 图 A：相速度 + 衰减 ============
    figA = figure('Visible', 'off', 'Position', [80 80 1100 460]);
    tA   = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile; hold on; box on; grid on; grid minor;
    for k = 1:nDT
        plot(f_ch, R{k}.Vph, 'LineStyle', LS{k}, 'Color', CLR(k,:), 'LineWidth', LW);
    end
    for k = 1:nDT
        plot(f_ch, R{k}.Vg, 'LineStyle', LS{k}, 'Color', CLR(k,:), 'LineWidth', LW*2.4);
    end
    xlabel('频率 f [Hz]');
    title('(a) 相速度与群速度');
    xlim([f_ch(1), f_ch(end)]);
    if nn == 0 && isst == 1, ylim([800, 2000]); end     % 斯通利波纵轴固定
    hold off;

    legend(cellfun(@(s) ['\DeltaT=' s], dT_lab, 'UniformOutput', false), ...
           'Location', 'northoutside', 'Orientation', 'horizontal', ...
           'NumColumns', 3, 'FontSize', FS-3);

    nexttile; hold on; box on; grid on; grid minor;
    for k = 1:nDT
        plot(f_ch, abs(R{k}.atten_mode), 'LineStyle', LS{k}, 'Color', CLR(k,:), 'LineWidth', LW);
    end
    xlabel('频率 f [Hz]'); ylabel('$\mathrm{Im}(k_z)\,/\,\mathrm{Re}(k_z)$', 'Interpreter', 'latex'); title('(b) 衰减');
    xlim([f_ch(1), f_ch(end)]); hold off;

    title(tA, sprintf('%s：%s', cname, hdr));
    saveas_print(figA, fullfile(outdir, [base '_频散']));

    % ============ 图 B：四张激发子图 ============
    figB = figure('Visible', 'off', 'Position', [80 80 1100 900]);
    tB   = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    QQ = { 'ap_a', '(a) 轴·广义压力', 'p_{axis}';
           'at_a', '(b) 轴·温度',     'a_{axis}';
           'ap_w', '(c) 壁·声压',     'p_{wall}';
           'at_w', '(d) 壁·温度',     'a_{wall}' };

    for iq = 1:4
        nexttile; hold on; box on; grid on; grid minor;
        for k = 1:nDT
            v = abs(R{k}.(QQ{iq,1}));
            if all(~isfinite(v)), continue; end
            plot(f_ch, v, 'LineStyle', LS{k}, 'Color', CLR(k,:), 'LineWidth', LW);
        end
        xlabel('频率 f [Hz]'); ylabel(QQ{iq,3}); title(QQ{iq,2});
        xlim([f_ch(1), f_ch(end)]); hold off;
        if iq == 1
            legend(cellfun(@(s) ['\DeltaT=' s], dT_lab, 'UniformOutput', false), ...
                   'Location', 'northoutside', 'Orientation', 'horizontal', ...
                   'NumColumns', 3, 'FontSize', FS-3);
        end
    end
    title(tB, sprintf('%s：%s', cname, hdr));
    saveas_print(figB, fullfile(outdir, [base '_激发']));
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