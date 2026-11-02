function out = main(pc, pf, dT, n, r)
% main —— 频散与激发：体波 + 井孔模式 + 极点激发留数（纯计算，不绘图）
%
%   out = main()                    % 五个算例参数全用 param_0 的默认值
%   out = main(pc, pf, dT, n, r)    % 围压[Pa] / 孔压[Pa] / 温度变化[K] / 阶数 / 接收半径[m]
%
%   输出：
%     out.freqs   1×nFreq  频率轴 [Hz]
%     out.Vp      4×nFreq  体波 P 相速度 [m/s]（4 支）
%     out.Vs      1×nFreq  S 波相速度 [m/s]
%     out.Atten   4×nFreq  体波衰减 Im(k)/Re(k)
%     out.Vph     1×nFreq  井孔模式相速度 [m/s]
%     out.kz_mode 1×nFreq  井孔模式波数 [1/m]
%     out.ap_w    1×nFreq  井壁声压激发
%     out.ap_a    1×nFreq  井轴声压激发
%     out.at_w    1×nFreq  井壁温度激发
%     out.at_a    1×nFreq  井轴温度激发
%     out.n / out.dT / out.pc  元数据
%
%   param_sym() 先把符号算式解一遍（persistent 缓存），之后每个频点只做数值求值。

here = fileparts(mfilename('fullpath'));
addpath(here);

if nargin < 1, pc = []; end
if nargin < 2, pf = []; end
if nargin < 3, dT = []; end
if nargin < 4, n  = []; end
if nargin < 5, r  = []; end

% 先解一遍：符号推导（结果缓存在 param_sym 里，只算一次）
param_sym();

% ---------------- 参数 ----------------
freqs   = logspace(log10(2e3), 6, 48);  % 频率轴 [Hz]（2000 Hz 起，48 点保光滑）
nFreq   = numel(freqs);
f0_seed = 1e4;                 % 播种频率：从这里向两端延拓
v1      = 1000;                % 粗扫相速度下限 [m/s]（上限由 fr.ks_val 给）

param = param_0(pc, pf, dT, n, r);

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

% ---------------- 打包 ----------------
out = struct( ...
    'freqs',   freqs(:).',   ...
    'Vp',      Vp,           ...
    'Vs',      Vs,           ...
    'Atten',   Atten,        ...
    'Vph',     Vph(:).',     ...
    'kz_mode', kz_all(:).',  ...
    'ap_w',    ap_w(:).',    ...
    'ap_a',    ap_a(:).',    ...
    'at_w',    at_w(:).',    ...
    'at_a',    at_a(:).',    ...
    'n',       param.n,      ...
    'dT',      param.dT,     ...
    'pc',      param.pc,     ...
    'pf',      param.pf);

fprintf('P 波相速度   %8.1f ~ %8.1f m/s\n', min(Vp(:)), max(Vp(:)));
fprintf('S 波相速度   %8.1f ~ %8.1f m/s\n', min(Vs), max(Vs));
fprintf('P 波衰减     %8.3e ~ %8.3e\n', min(abs(Atten(:))), max(abs(Atten(:))));
ok = isfinite(Vph);
if any(ok(:))
    fprintf('井孔模式相速度 %8.1f ~ %8.1f m/s（有效点 %d/%d）\n', ...
        min(Vph(ok)), max(Vph(ok)), sum(ok(:)), numel(ok));
end
end
