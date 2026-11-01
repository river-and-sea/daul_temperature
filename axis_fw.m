% axis_fw —— 井轴上的"广义声压"和温度全波列（张碧星等 1995，式 (9)）
%
% 脚本程序：直接运行本文件即可。
%   逐 (ω, kz) 求出式 (9) 第二项的被积函数，再做两次 FFT 得到
%   时域总场 = 直达项 + 反射项，最后画 P（声压）和 H（温度）的波形。
%
%   井内两支：M8_solve 给出 An（解级、按 I_n 归一）。
%   声压通道取 sum(An)，温度通道取 sum(An .* alphani)（照 excite 的口径）。

%% ============ 0. 参数 ============
param  = param_0();                          % 五参数用默认：pc=0, pf=0, theta0=0, R=0.1, n=0
n      = param.n;
n_rec  = 3;                                  % 接收道数

% ---- 网格常量（原 dt_config.fw_*）----
fw_tc        = 0.5e-3;             % 源脉冲特征时间 [s]
fw_fc_list   = [6000 5000 6000];   % 各 n 的源中心频率 [Hz]
fw_n_list    = [0 1 2];            % 允许的 n
fw_dt        = 10e-6;              % 时间步 [s]
fw_n2        = 2049;               % 频域点数
fw_n1        = 4097;               % k 域点数
fw_dz        = 0.5;                % 接收器间距 [m]
fw_ddz0      = 0.0208;             % 粗空间步 [m]
fw_zp        = 2;                  % 第一道源距 [m]
fw_f0i_ratio = 1/100;              % 复频偏移比例

% ============ 栅格 ============
% 两次 FFT 各需要一套栅格：
%   z ←→ k_z（空间-波数），t ←→ f（时间-频率）

% ---- 源 ----
f_src = fw_fc_list(fw_n_list == n);        % 源中心频率 [Hz]（n=0/1/2 → 6000/5000/6000）
tc    = fw_tc;                             % 源脉冲特征时间 [s]

% ---- 空间栅格 z ----
dt     = fw_dt;                            % 时间步 [s]
n1     = fw_n1;                            % z 域点数（4097）
n2     = fw_n2;                            % t 域点数（2049）
n_step = fix(fw_dz/fw_ddz0);               % 道间距被细分的步数：0.5/0.0208 → 24
d_z    = fw_dz/n_step;                     % 空间步 [m] ≈ 0.0208
m0     = round(fw_zp/d_z);                 % 第一道（z = 2 m）对应的点号
m_rec  = 1 + m0 + (0:(n_rec-1))*n_step;    % n_rec 个接收道的点号
z_rec  = m_rec*d_z;                        % 接收道位置 [m]

% ---- 波数栅格 k_z（只取正半轴，另一半由共轭镜像补）----
n_k    = fix(n1/2) + 1;                    % 正半轴点数 = 2049
dk     = 2*pi/(d_z*(n1-1));                % 波数步 [1/m]（与 d_z 的 DFT 对偶）
k_z    = ((1:n_k)-1)*dk + 1e-36;           % +1e-36 只是避开 k = 0

% ---- 时间栅格 t（只取前一半，后一半是 FFT 周期延拓）----
n_t    = fix(0.5*(fix(n2/2)+1));           % 输出时间点数 = 512
i_t    = 1:n_t;                            % 输出时间点下标

% ---- 频率栅格 f（只算源谱带宽内的频点，不是全部 2049 个）----
dfp    = 1/(dt*(n2-1));                    % 频率步 [Hz] ≈ 48.8
iw_f   = (round((f_src-2/tc)/dfp)+1) : round((f_src+2/tc)/dfp);   % f_src ± 2/tc 内的主 FFT 下标
n_w    = numel(iw_f);                      % 频点数 ≈ 164
freq   = iw_f*dfp;                         % 实频率 [Hz]
% 复频偏移：把积分路径从实轴抬到 w + i*sigma（sigma = 2*pi*f_src*ratio ≈ 377 rad/s）。
%   由 P(t) = e^{sigma*t} * Int[ P(w+i*sigma) e^{-i w t} ] 可知：
%   用复频算出的 FFT 结果是 e^{-sigma*t} P(t)，故时域要乘回 e^{+sigma*t}（见下方 decay）。
%   作用：① 压制 FFT 的周期回卷（窗尾被压掉几个量级）；② 复频上根/Bessel 更稳。
%   偏移本身不改变最终结果，只是数值手段。
freq_c = freq + 1i*f_src*fw_f0i_ratio;
om_c   = 2*pi*freq_c;                      % 复角频率 [rad/s]
F      = source_spec_fw(freq_c, f_src, tc);% 源谱 F(ω)

%% ============ 1. 逐 (ω, kz) 求式 (9) 的被积函数 + 直达项 ============
p_r  = zeros(n_w, n_k);                    % 声压：反射项核
h_r  = zeros(n_w, n_k);                    % 温度：反射项核
P_d  = zeros(n_w, n_rec);                  % 声压：直达项
H_d  = zeros(n_w, n_rec);                  % 温度：直达项
r0   = param.a0;
angf = 1;                                  % cos(n(θ-θ0))，θ=θ0 时 = 1
chengzi = angf * (1 + (n > 0)) * (r0/4)^n / (4*pi*factorial(n));   % eps_n = 2 - δ_{n0}

for iw = 1:n_w
    fr   = ki_omega(om_c(iw), param, []);  % 复频下的解
    k_in = fr.kni_vec(:);                  % 2x1 井内两支轴向波数（列向量，直达项要做外层积）
    aln  = fr.alphani(:);                  % 2x1 两支的温度/势比

    % ---- 直达项（式 (9) 第一项）----
    Dp_ = zeros(1, n_rec);                 % 声压：两支直接求和
    Dh_ = zeros(1, n_rec);                 % 温度：两支乘 alphani 后求和
    for m = 0:n
        term = exp(1i*(k_in*z_rec)) .* ((-2i*k_in).^m) ./ (z_rec.^(2*n-m+1));
        coef = nchoosek(n,m)*factorial(2*n-m);
        Dp_ = Dp_ + coef .* sum(term, 1);
        Dh_ = Dh_ + coef .* sum(term .* aln, 1);
    end
    P_d(iw,:) = chengzi .* Dp_;
    H_d(iw,:) = chengzi .* Dh_;

    % ---- 反射项的被积函数 ----
    for ik = 1:n_k
        kz  = k_z(ik);
        knr = sqrt(kz.^2 - k_in.^2);       % 2x1 井内两支径向波数

        [~, ~, ~, An] = M8_solve(fr, param, kz);   % 第 4 输出：解级、已按 I_n 归一（不算 Dp）
        An = An.';                                  % 1x2 -> 2x1

        w2n = knr.^(2*n);                  % 多极源的 n 阶径向因子
        p_r(iw,ik) = chengzi * sum(w2n .* An);
        h_r(iw,ik) = chengzi * sum(w2n .* An .* aln);
    end
end

%% ============ 2. 双重 FFT -> 时域 ============
% 反射项：kz -> z（integrate_kz 补负半轴镜像、取接收道），再乘源谱
P_r_z = integrate_kz(p_r, n1, n_k, m_rec) .* F(:);
H_r_z = integrate_kz(h_r, n1, n_k, m_rec) .* F(:);

% 直达项：本来就在 z 域，乘源谱及除采样时间间隔
P_d_z = P_d .* F(:) / dt;
H_d_z = H_d .* F(:) / dt;

% ω -> t（integrate_omega 补负频共轭镜像）
Tp = integrate_omega(P_r_z, n2, iw_f);
Td = integrate_omega(P_d_z, n2, iw_f);
Th = integrate_omega(H_r_z, n2, iw_f);
Tg = integrate_omega(H_d_z, n2, iw_f);

% 乘回 e^{+sigma*t}：复频逆变换自带 e^{-sigma*t}（见 freq_c 处），这一步精确还原真实信号
decay = exp(2*pi*f_src*fw_f0i_ratio*i_t*dt).';

P_r_t = real(Tp(i_t,:) .* decay);
P_d_t = real(Td(i_t,:) .* decay);
H_r_t = real(Th(i_t,:) .* decay);
H_d_t = real(Tg(i_t,:) .* decay);

P = (P_r_t + P_d_t).';                     % n_rec x n_t 声压
H = (H_r_t + H_d_t).';                     % n_rec x n_t 温度

fprintf('n=%d pc=%g MPa  z=%.4f..%.4f m\n', n, (param.pc)/1e6, z_rec(1), z_rec(end));
fprintf('  P: max|P_r|=%.4g  max|P_d|=%.4g  max|P|=%.4g\n', ...
    max(abs(P_r_t(:))), max(abs(P_d_t(:))), max(abs(P(:))));
fprintf('  H: max|H_r|=%.4g  max|H_d|=%.4g  max|H|=%.4g\n', ...
    max(abs(H_r_t(:))), max(abs(H_d_t(:))), max(abs(H(:))));

%% ============ 3. 画 P、H 的波形 ============
t_plot = i_t*dt*1000 - dt*1000;            % [ms]
fig = figure('Visible', 'off', 'Position', [100 100 900 800]);

subplot(2,1,1);
plot(t_plot, P.', 'LineWidth', 1);
grid on; xlabel('t [ms]'); ylabel('P');
leg = arrayfun(@(z) sprintf('z = %.3f m', z), z_rec, 'UniformOutput', false);
legend(leg, 'Location', 'best');
title(sprintf('axis pressure  (n = %d, pc = %g MPa)', n, param.pc/1e6));

subplot(2,1,2);
plot(t_plot, H.', 'LineWidth', 1);
grid on; xlabel('t [ms]'); ylabel('H');
legend(leg, 'Location', 'best');
title('axis temperature');

outdir = fullfile(fileparts(mfilename('fullpath')), 'fig');
if exist(outdir, 'dir') ~= 7, mkdir(outdir); end
saveas(fig, fullfile(outdir, 'fig_axis_fw.png'));
close(fig);
fprintf('saved: %s\n', fullfile(outdir, 'fig_axis_fw.png'));

% ======================================================================
function Zz = integrate_kz(Z, n1, n_k, m_rec)
% kz -> z：补负半轴镜像后做 FFT，再取 n_rec 个接收道
Zk = zeros(size(Z,1), n1, 'like', Z);
Zk(:, 1:n_k)          = Z;
Zk(:, n1+1-(2:n_k))   = Z(:, 2:n_k);          % 负半轴：复制式镜像
rec_z = fft(Zk(:, 1:n1-1).').';               % kz -> z
Zz = rec_z(:, m_rec);
end


% ======================================================================
function Zt = integrate_omega(Z, n2, iw_f)
% ω -> t：补负频共轭镜像后做 FFT
Xf = zeros(n2-1, size(Z,2), 'like', Z);
Xf(iw_f, :)         = Z;                      % 正频段
Xf(n2-iw_f+1, :)    = conj(Z);                % 负频共轭镜像
Zt = fft(Xf, [], 1);
end