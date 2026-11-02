function main()
% main —— 重构版主程序：参数 → 求根 → 沿支追踪 → 频率相关量 → 体波 + 井孔模式
%
%   param_sym() 先把符号算式解一遍（persistent 缓存），之后每个频点只做数值求值。

here = fileparts(mfilename('fullpath'));
addpath(here);

% 先解一遍：符号推导（结果缓存在 param_sym 里，只算一次）
param_sym();

% ---------------- 参数 ----------------
freqs   = logspace(3, 6, 16);  % 频率轴 [Hz]（快速验证用 16 点；正式跑改 128）
nFreq   = numel(freqs);
f0_seed = 1e4;                 % 播种频率：从这里向两端延拓
v1      = 1000;                % 粗扫相速度下限 [m/s]（上限由 fr.ks_val 给）

param = param_0();             % pc/pf/dT/n/r 用 param_0 里的默认值

[~, i0] = min(abs(freqs - f0_seed));
order   = [i0:-1:1, i0:nFreq];   % 双向延拓

Vp    = nan(4, nFreq);
Atten = nan(4, nFreq);
Vs    = nan(1, nFreq);
kz_all = nan(nFreq, 1);
ap_w = nan(nFreq,1);  at_w = nan(nFreq,1);
ap_a = nan(nFreq,1);  at_a = nan(nFreq,1);

lambda_pre = [];                        % 上一频点的 4 个根（首点空）
kz_pre     = NaN;                       % 上一频点的 kz（首点 NaN -> 触发粗扫）
kz_f0      = NaN;                       % f0 处的根（两条腿共用）

for ii = 1:numel(order)
    idx       = order(ii);
    omega     = 2*pi*freqs(idx);
    omega_pre = 2*pi*freqs(order(max(ii-1, 1)));        % 由 order 导出，首点取自身

    if idx == i0 && ~isnan(kz_f0), kz_pre = kz_f0;  omega_pre = omega;  end   % 第二条腿：复用 f0 的根

    fr         = ki_omega(omega, param, lambda_pre);
    lambda_pre = fr.lambda_current;

    % 体波（色散方程根直接给）
    Vp(:, idx)    = omega ./ real(fr.ki_vec);
    Atten(:, idx) = imag(fr.ki_vec) ./ real(fr.ki_vec);
    Vs(idx)       = omega / real(fr.ks_val);

    % 井孔模式（8x8 边界条件 det = 0）
    detf         = @(kz) M8_solve(fr, param, kz);
    [kz, kz_pre] = kz_solve(detf, omega, omega_pre, [v1, fr.ks_val], kz_pre);
    kz_all(idx)  = kz;

    % 极点处的激发留数（只在算到根时才有意义）
    if isfinite(kz)
        ex = excite(fr, param, kz, omega);
        ap_w(idx) = ex.ap_wall;  at_w(idx) = ex.at_wall;
        ap_a(idx) = ex.ap_axis;  at_a(idx) = ex.at_axis;
    end

    if idx == i0 && isnan(kz_f0), kz_f0 = kz;  end                            % 第一次到 f0：记下根
end

Vph = 2*pi*freqs(:) ./ real(kz_all);    % 井孔模式相速度
ok  = isfinite(Vph);

fprintf('P 波相速度   %8.1f ~ %8.1f m/s\n', min(Vp(:)), max(Vp(:)));
fprintf('S 波相速度   %8.1f ~ %8.1f m/s\n', min(Vs), max(Vs));
fprintf('P 波衰减     %8.3e ~ %8.3e\n', min(abs(Atten(:))), max(abs(Atten(:))));
if any(ok(:))
    fprintf('井孔模式相速度 %8.1f ~ %8.1f m/s（有效点 %d/%d）\n', ...
        min(Vph(ok)), max(Vph(ok)), sum(ok(:)), numel(ok));
end

% ---------------- 出图 ----------------
fig = figure('Visible', 'off', 'Position', [100 100 900 1200]);

subplot(5,1,1);
semilogx(freqs, Vp.', 'LineWidth', 1); hold on;
semilogx(freqs, Vs, 'k--', 'LineWidth', 1.2);
grid on; xlabel('f [Hz]'); ylabel('phase velocity [m/s]');
legend({'P_1','P_2','P_3','P_4','S'}, 'Location', 'best');
title('body-wave phase velocity');

subplot(5,1,2);
semilogx(freqs, Atten.', 'LineWidth', 1);
grid on; xlabel('f [Hz]'); ylabel('Im(k)/Re(k)');
title('attenuation');

subplot(5,1,3);
semilogx(freqs, Vph, 'LineWidth', 1.2);
grid on; xlabel('f [Hz]'); ylabel('kz phase velocity [m/s]');
title(sprintf('borehole mode (n = %d)', param.n));

subplot(5,1,4);
semilogx(freqs, abs(ap_w), 'LineWidth', 1.2); hold on;
semilogx(freqs, abs(ap_a), '--', 'LineWidth', 1.2);
grid on; xlabel('f [Hz]'); ylabel('|ap|');
legend({'ap\_wall','ap\_axis'}, 'Location', 'best');
title('pressure excitation at the poles');

subplot(5,1,5);
semilogx(freqs, abs(at_w), 'LineWidth', 1.2); hold on;
semilogx(freqs, abs(at_a), '--', 'LineWidth', 1.2);
grid on; xlabel('f [Hz]'); ylabel('|at|');
legend({'at\_wall','at\_axis'}, 'Location', 'best');
title('temperature excitation at the poles');

outdir = fullfile(here, 'fig');
if exist(outdir, 'dir') ~= 7, mkdir(outdir); end
saveas(fig, fullfile(outdir, 'fig_dispersion.png'));
close(fig);
fprintf('saved: %s\n', fullfile(outdir, 'fig_dispersion.png'));
end
