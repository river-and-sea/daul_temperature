function out = fullwave(pc, pf, dT, n, r)
% fullwave —— 井轴 / 井壁全波：八个 (t,z) 场
%
%   out = fullwave()                    % 五个算例参数全用 param_0 的默认值
%   out = fullwave(pc, pf, dT, n, r)    % 围压[Pa] / 孔压[Pa] / 温度变化[K] / 阶数 / 接收半径[m]
%
%   输出（8 个场均 N_t×N_k = 512×2049，t 为行、z 为列）：
%     out.time        1×512    时间轴 [s]，0 ~ (N_t-1)·dt = 5.12 ms
%     out.z           1×2049   空间轴 [m]，0 ~ (N_k-1)·dz = 42.67 m
%     out.P_axis_r    N_t×N_k  井轴声压 反射
%     out.P_axis_d    N_t×N_k  井轴声压 直达
%     out.P_wall_r    N_t×N_k  井壁声压 反射
%     out.P_wall_d    N_t×N_k  井壁声压 直达
%     out.H_axis_r    N_t×N_k  井轴温度 反射
%     out.H_axis_d    N_t×N_k  井轴温度 直达
%     out.H_wall_r    N_t×N_k  井壁温度 反射
%     out.H_wall_d    N_t×N_k  井壁温度 直达
%     out.f_src       标量      源中心频率 [Hz]
%
%   物理口径：
%     井轴（A，轴压 / 张碧星式 (9)）：p = chengzi * sum( knr^(2n) * An )
%     井壁（W，全波 / 旧 zbx_fullwave）：p = angf * (src + x(5:6))
%
%   变换：
%     kz -> z：补负半轴（复制式）镜像 + fft，归一化 dk/(2π)
%     ω -> t ：乘源谱 F，补负频（共轭）镜像 + fft，归一化 dω/(2π)，
%              乘 e^{+σt} 还原复频偏移，取实部

if nargin < 1, pc = []; end
if nargin < 2, pf = []; end
if nargin < 3, dT = []; end
if nargin < 4, n  = []; end
if nargin < 5, r  = []; end

param = param_0(pc, pf, dT, n, r);
n = param.n;

%% ============ 0. 网格常量 ============
tc         = 0.5e-3;                  % 源脉冲特征时间 [s]
fw_fc_list = [6000 5000 6000];        % 各 n 的源中心频率 [Hz]
fw_n_list  = [0 1 2];                 % 允许的 n
dt         = 10e-6;                   % 时间步 [s]
N_t2       = 2049;                    % 频域点数
N_k2       = 4097;                    % k 域点数
f0_ratio   = 1/100;                   % 复频偏移比例
dz         = 1/48;                    % 空间步 [m]

f_src = fw_fc_list(fw_n_list == n);   % 源中心频率 [Hz]

N_k = fix(N_k2/2) + 1;                % 波数正半轴点数 = 2049
dk  = 2*pi/(dz*(N_k2-1));             % 波数步 [1/m]
kz  = linspace(0, (N_k-1)*dk, N_k) + 1e-36;

N_t  = fix(0.5*(fix(N_t2/2)+1));      % 输出时间点数 = 512
df   = 1/(dt*(N_t2-1));               % 频率步 [Hz]
i_fc = (round((f_src-2/tc)/df)+1) : round((f_src+2/tc)/df);
N_w  = numel(i_fc);
freq = i_fc*df;

freq_c = freq + 1i*f_src*f0_ratio;    % 复频偏移
om_c   = 2*pi*freq_c;
F      = source_spec_fw(freq_c, f_src, tc);

%% ============ 1. 逐 (ω,kz)：六个核 ============
pA_r = zeros(N_w, N_k);  hA_r = zeros(N_w, N_k);   % 井轴：声压/温度 反射核
pW_r = zeros(N_w, N_k);  pW_d = zeros(N_w, N_k);   % 井壁：声压 反射/直达核
hW_r = zeros(N_w, N_k);  hW_d = zeros(N_w, N_k);   % 井壁：温度 反射/直达核

k_in_all = zeros(N_w, 2);                          % 每频点井内两支轴向波数
aln_all  = zeros(N_w, 2);                          % 每频点两支温度/势比

r0      = param.r0;
angf    = cos(n*(param.th_rec - param.th_src));           % 角向因子 cos(nΔθ)
chengzi = param.s0*(1 + (n > 0)) * (r0/4)^n / (2*factorial(n));

for iw = 1:N_w
    omega = om_c(iw);
    fr   = ki_omega(omega, param, []);          % 复频下的解
    k_in = fr.kni_vec(:).';                        % 1x2
    aln  = fr.alphani(:);                          % 2x1
    k_in_all(iw,:) = k_in;
    aln_all(iw,:)  = aln.';

    for ik = 1:N_k
        [~, src_wall, An_wall, An_axis] = M8_solve(fr, param, kz(ik));

        pA_r(iw,ik) = angf * sum(An_axis);                 % 井轴反射
        hA_r(iw,ik) = angf * sum(An_axis .* aln.');
        pW_d(iw,ik) = angf * param.den_f_in * omega^2 * sum(src_wall);                % 井壁直达
        pW_r(iw,ik) = angf * param.den_f_in * omega^2 * sum(An_wall);                 % 井壁反射
        hW_d(iw,ik) = angf * sum(src_wall .* aln.');
        hW_r(iw,ik) = angf * sum(An_wall .* aln.');
    end
end

%% ============ 2. 八个场（N_t×N_k）============
z_all  = (0:N_k-1)*dz;                % z 轴 2049 点 [m]
iz_all = 1:N_k;                       % 全部 z 列
zw_all = z_all(:).';

% 井轴直达项（全 z 解析式，复用 §1 存下的 k_in/aln）
pAd = zeros(N_w, N_k);
hAd = zeros(N_w, N_k);
for iw = 1:N_w
    k_in = k_in_all(iw,:).';
    aln  = aln_all(iw,:).';
    pa = zeros(1, N_k);
    ha = zeros(1, N_k);
    for mm = 0:n
        term = exp(1i*(k_in*zw_all)) .* ((-2i*k_in).^mm) ./ (zw_all.^(2*n-mm+1));
        coef = nchoosek(n,mm)*factorial(2*n-mm);
        pa = pa + coef .* sum(term, 1);
        ha = ha + coef .* sum(term .* aln, 1);
    end
    pAd(iw,:) = chengzi .* pa;
    hAd(iw,:) = chengzi .* ha;
end

% 归一化 + 时间还原
N_full     = N_t2 - 1;
decay_full = exp(2*pi*f_src*f0_ratio*(1:N_full).'*dt);

% 两个变换器：kz 域核 -> (t,z)；z 域直达 -> (t,z)
to_tz = @(X) integrate_omega(integrate_kz(X, N_k2, N_k, iz_all, dz), N_t2, i_fc, dt, F, decay_full, N_t);
to_t  = @(X) integrate_omega(X, N_t2, i_fc, dt, F, decay_full, N_t);

out.P_axis_r = to_tz(pA_r);   out.P_axis_d = to_t(pAd);
out.P_wall_r = to_tz(pW_r);   out.P_wall_d = to_tz(pW_d);
out.H_axis_r = to_tz(hA_r);   out.H_axis_d = to_t(hAd);
out.H_wall_r = to_tz(hW_r);   out.H_wall_d = to_tz(hW_d);

out.time  = (0:N_t-1)*dt;             % 1×512  时间轴 [s]
out.z     = z_all;                    % 1×2049 空间轴 [m]
out.f_src = f_src;
out.n     = n;
out.pc    = param.pc;
out.pf    = param.pf;
end

% ======================================================================
function Zz = integrate_kz(Zk, N_k2, N_k, i_z, dz)
% kz -> z：补负半轴镜像（复制式）后做 FFT，再取 i_z 指定的列
%   Zk   : (N_w, N_k)  kz 域核，只含 kz >= 0 的半谱
%   N_k2 : z 域点数（fft 补零到这里，长度 N_k2-1）
%   i_z  : 要取的 z 列索引（向量则返回全部对应 z）
Zfull = zeros(size(Zk,1), N_k2, 'like', Zk);        % 补零到 N_k2 点
Zfull(:, 1:N_k)          = Zk;                       % 正半轴
Zfull(:, N_k2+1-(2:N_k)) = Zk(:, 2:N_k);             % 负半轴：复制式镜像
rec_z = fft(Zfull(:, 1:N_k2-1).').';                 % kz -> z（正指数效果）
Zz = rec_z(:, i_z) * (1/(dz*(N_k2-1)));              % 取列 + dk/(2π) 归一化
end

% ======================================================================
function Zt = integrate_omega(Zw, N_t2, i_fc, dt, F, decay, n_out)
% ω -> t：乘源谱 → 补负频共轭镜像 → fft → 乘衰减还原 → 取实部 → 截前 n_out 行
Zw  = Zw .* F(:);
Xf = zeros(N_t2-1, size(Zw,2), 'like', Zw);
Xf(i_fc, :)          = Zw;                     % 正频段
Xf(N_t2-i_fc+1, :)   = conj(Zw);               % 负频共轭镜像
Zt = real(fft(Xf, [], 1) * (1/(dt*(N_t2-1))) .* decay);   % dω/(2π) + e^{σt} 还原
if nargin >= 7 && ~isempty(n_out)
    Zt = Zt(1:n_out, :);                       % 只保留前 n_out 个时间点
end
end
