function out = main(pc, pf, dT, n, r, is_st, freqs_lo, freqs_hi, nFreq)
% main —— 频散与激发：体波 + 单个井孔导波模式（纯计算，不绘图）
%
%   out = main()                                   % 全用默认值
%   out = main(pc, pf, dT, n, r, is_st)            % 指定算例参数与分支
%   out = main(pc, pf, dT, n, r, is_st, f_lo, f_hi, nF)   % 再指定频率轴
%
%   输入
%     pc, pf, dT, n, r   围压[Pa] / 孔压[Pa] / 温度变化[K] / 方位阶数 / 接收半径[m]
%     is_st   导波分支标志（nMode 恒为 1）：
%              1 -> 斯通利型分支，相速度窗口 [0.8 v_f, v_f]
%              0 -> 伪瑞利型分支，相速度窗口 [v_f, V_S]
%              （v_f = 井内流体速度，V_S = 该频点 S 波速度）
%     freqs_lo, freqs_hi, nFreq   频率轴起点/终点 [Hz] 与点数（默认 2e3 / 2e4 / 128）
%
%   输出
%     out.freqs      1×nFreq     频率轴 [Hz]
%     out.Vp         4×nFreq     体波 P 相速度 [m/s]（4 支）
%     out.Vs         1×nFreq     S 波相速度 [m/s]
%     out.Atten      4×nFreq     体波衰减 Im(k)/Re(k)
%     out.Vph        1×nFreq     导波相速度 [m/s]（V_ph = V_S 的截止点保留，之后 NaN）
%     out.Vg         1×nFreq     导波群速度 [m/s]（dw/dk 数值梯度）
%     out.kz_mode    1×nFreq     导波波数 [1/m]
%     out.atten_mode 1×nFreq     导波衰减 Im(k)/Re(k)
%     out.ap_w / ap_a / at_w / at_a   1×nFreq  激发留数（壁/轴 × 声压/温度）
%     out.is_st / out.n / out.dT / out.pc / out.pf  元数据
%
%   其他导波（弯曲波 n=1、螺旋波 n=2 …）：只改 n，is_st 传 0。

if nargin < 1, pc = []; end
if nargin < 2, pf = []; end
if nargin < 3, dT = []; end
if nargin < 4, n  = []; end
if nargin < 5, r  = []; end
if nargin < 6 || isempty(is_st),    is_st    = 1;   end
if nargin < 7 || isempty(freqs_lo), freqs_lo = 2e3; end
if nargin < 8 || isempty(freqs_hi), freqs_hi = 2e4; end
if nargin < 9 || isempty(nFreq),    nFreq    = 128; end

here = fileparts(mfilename('fullpath'));
addpath(here);
param_sym();

% ---------------- 参数 ----------------
freqs   = logspace(log10(freqs_lo), log10(freqs_hi), nFreq);
f0_seed = (freqs_lo + freqs_hi) / 2;

param = param_0(pc, pf, dT, n, r);

[~, i0] = min(abs(freqs - f0_seed));
order   = [i0:-1:1, i0:nFreq];   % 双向延拓

Vp     = nan(4, nFreq);
Atten  = nan(4, nFreq);
Vs     = nan(1, nFreq);
kz_all = nan(1, nFreq);
ap_w = nan(1, nFreq);  at_w = nan(1, nFreq);
ap_a = nan(1, nFreq);  at_a = nan(1, nFreq);

lambda_pre = [];                        % 上一频点的 4 个体波根（ki_omega 沿支跟踪用）
omega_pre  = NaN;                       % 上一频点角频率
kz_pre     = NaN;                       % 上一频点导波根
kz_f0      = NaN;                       % f0 处导波根（两条腿共用）
active     = true;                      % 是否继续求解（V_ph = V_S 后停止，该点保留）

for ii = 1:numel(order)
    if ii == i0 + 1, active = true; end      % 第二腿（向高频）重新开始

    idx   = order(ii);
    omega = 2*pi*freqs(idx);

    fr         = ki_omega(omega, param, lambda_pre);
    lambda_pre = fr.lambda_current;

    % 体波（色散方程根直接给）
    Vp(:, idx)    = omega ./ real(fr.ki_vec);
    Atten(:, idx) = imag(fr.ki_vec) ./ real(fr.ki_vec);
    Vs(idx)       = omega / real(fr.ks_val);

    if idx == i0 && ~isnan(kz_f0), kz_pre = kz_f0;  omega_pre = omega; end
    if isnan(omega_pre), omega_pre = omega; end

    % ---- 单个井孔模式 ----
    if active
        detf = @(kz) M8_solve(fr, param, kz);
        vf   = omega / real(fr.kni_vec(1));     % 井内流体速度 [m/s]
        vS   = Vs(idx);                          % 本频点 S 波速度 [m/s]

        if is_st
            vwin = [0.8*vf, 0.99*vf];                 % 斯通利型分支
        else
            vwin = [1.05*vf, 0.95*vS];                % 伪瑞利型分支（上界略收，避开 V_S）
        end

        [kz, kz_pre] = kz_solve(detf, omega, omega_pre, vwin, kz_pre);   % 只在首点粗扫一次

        if isfinite(kz) && abs(omega/real(kz) - vS)/vS < 1e-3
            active = false;
            kz     = NaN;                        % 落到 V_S：此点也舍去，之后停止
        end
        kz_all(idx) = kz;

        if idx == i0 && isnan(kz_f0), kz_f0 = kz;  end
    end
    omega_pre = omega;

    % 极点处的激发留数
    if isfinite(kz_all(idx))
        ex = excite(fr, param, kz_all(idx), omega);
        ap_w(idx) = ex.ap_wall;  at_w(idx) = ex.at_wall;
        ap_a(idx) = ex.ap_axis;  at_a(idx) = ex.at_axis;
    end
end

Vph        = 2*pi*freqs(:).' ./ real(kz_all);
atten_mode = imag(kz_all) ./ real(kz_all);

% 群速度 vg = dw/dk（数值梯度；无效段自然传播 NaN）
wv  = 2*pi*freqs(:);
kv  = wv ./ Vph(:);
Vg  = (gradient(wv) ./ gradient(kv)).';

% ---------------- 打包 ----------------
out = struct( ...
    'freqs',      freqs(:).',   ...
    'Vp',         Vp,           ...
    'Vs',         Vs,           ...
    'Atten',      Atten,        ...
    'Vph',        Vph,          ...
    'Vg',         Vg,           ...
    'kz_mode',    kz_all,       ...
    'atten_mode', atten_mode,   ...
    'ap_w',       ap_w,         ...
    'ap_a',       ap_a,         ...
    'at_w',       at_w,         ...
    'at_a',       at_a,         ...
    'is_st',      is_st,        ...
    'n',          param.n,      ...
    'dT',         param.dT,     ...
    'pc',         param.pc,     ...
    'pf',         param.pf);

ok = isfinite(Vph);
if any(ok)
    fprintf('mode (n=%d, is_st=%d): %.1f ~ %.1f m/s, %d/%d pts\n', ...
        n, is_st, min(Vph(ok)), max(Vph(ok)), sum(ok), nFreq);
else
    fprintf('mode (n=%d, is_st=%d): no valid point\n', n, is_st);
end
end
