% axis_fw —— 井轴 / 井壁全波：八个 (t,z) 场
%
% 脚本程序：直接运行本文件即可（无输入、无绘图、不打包 struct）。
% 运行后工作区留下：
%     P_axis_r, P_axis_d, P_wall_r, P_wall_d   声压（t×z，2048×4096）
%     H_axis_r, H_axis_d, H_wall_r, H_wall_d   温度
%     t_all (1×2048)、z_all (1×4096)           轴
%
%   井轴（A，轴压 / 张碧星式 (9)，对应旧 zbx_axis_pressure）：
%       p = chengzi * sum( knr^(2n) * An )        An 已按 I_n 归一
%       直达项用 k_in 的解析式在全 z 上求值
%
%   井壁（W，全波，对应旧 zbx_fullwave）：
%       p = rho_f*om^2 * (src + x(5:6))           src 带 (knr^2*r0/4)^n*e^{-κ}
%
%   ── 命名规范 ──────────────────────────────────────────────────────────
%   核（ω,kz 域）：{p,h}{A,W}_{r,d}         小写 p/h，A=井轴 W=井壁，r=反射 d=直达
%   输出（t,z 域）：{P,H}_{axis,wall}[_r,_d]  大写、全名
%   kz->z 见 integrate_kz；ω->t（含乘源谱、衰减还原、取实部）见 integrate_omega

%% ============ 0. 参数 ============
param  = param_0();                          % 五参数用默认
n      = param.n;

% ---- 网格常量（原 dt_config.fw_*）----
tc        = 0.5e-3;             % 源脉冲特征时间 [s]
fw_fc_list   = [6000 5000 6000];   % 各 n 的源中心频率 [Hz]
fw_n_list    = [0 1 2];            % 允许的 n
dt        = 10e-6;              % 时间步 [s]
N_t2        = 2049;               % 频域点数
N_k2        = 4097;               % k 域点数
f0_ratio = 1/100;              % 复频偏移比例
d_z         = 1/48;              % 空间步 [m]（与 d_k 满足 d_k*d_z*(N_kz-1) = 2π）

% ============ 栅格 ============
% ---- 源 ----
f_src = fw_fc_list(fw_n_list == n);        % 源中心频率 [Hz]


% ---- 空间轴 / 波数轴（唯一约束：dk*d_z*(n1-1) = 2π）----
N_k = fix(N_k2/2) + 1;                       % 波数正半轴点数
d_k  = 2*pi/(d_z*(N_k2-1));                   % 波数步 [1/m]，由 d_z 的对偶关系定出
kz = linspace(0, (N_k-1)*d_k, N_k) + 1e-36;

% ---- 时间轴 ----
N_t = fix(0.5*(fix(N_t2/2)+1));             % 输出时间点数（N_t 是频域点数，别覆盖）
i_t  = 1:N_t;
time    = linspace(0, (N_t-1)*dt, N_t);

% ---- 频率轴 ----
df  = 1/(dt*(N_t2-1));
i_fc = (round((f_src-2/tc)/df)+1) : round((f_src+2/tc)/df);
N_w  = numel(i_fc);
freq = i_fc*df;

% 复频偏移：路径抬到 w + i*sigma；时域乘回 e^{+sigma*t}
freq_c = freq + 1i*f_src*f0_ratio;
om_c   = 2*pi*freq_c;
F      = source_spec_fw(freq_c, f_src, tc);

%% ============ 1. 逐 (ω, kz)：六个核 + 存 k_in/aln ============
pA_r = zeros(N_w, N_k);  hA_r = zeros(N_w, N_k);   % 井轴：声压/温度 反射核
pW_r = zeros(N_w, N_k);  pW_d = zeros(N_w, N_k);   % 井壁：声压 反射/直达核
hW_r = zeros(N_w, N_k);  hW_d = zeros(N_w, N_k);   % 井壁：温度 反射/直达核

k_in_all = zeros(N_w, 2);                          % 每频点井内两支轴向波数
aln_all  = zeros(N_w, 2);                          % 每频点两支温度/势比

r0   = param.r0;
angf = cos(n*(param.th_rec - param.th_src));          % 角向因子 cos(nΔθ)
chengzi = param.s0*(1 + (n > 0)) * (r0/4)^n / (2*factorial(n));

for iw = 1:N_w
    fr   = ki_omega(om_c(iw), param, []);          % 复频下的解
    k_in = fr.kni_vec(:).';                        % 1x2
    aln  = fr.alphani(:);                          % 2x1
    k_in_all(iw,:) = k_in;
    aln_all(iw,:)  = aln.';

    for ik = 1:N_k
        [~, src_wall, An_wall, An_axis] = M8_solve(fr, param, kz(ik));

        pA_r(iw,ik) = angf * sum(An_axis);                 % 井轴反射
        hA_r(iw,ik) = angf * sum(An_axis .* aln.');
        pW_d(iw,ik) = angf * sum(src_wall);                % 井壁直达
        pW_r(iw,ik) = angf * sum(An_wall);                 % 井壁反射
        hW_d(iw,ik) = angf * sum(src_wall .* aln.');
        hW_r(iw,ik) = angf * sum(An_wall .* aln.');
    end
end

%% ============ 2. 八个场：轴/壁 × 直/反 × P/H（2048×4096）============
z_all  = (0:N_k2-2)*d_z;                       % z 轴 [m]
iz_all = 1:(N_k2-1);                           % 全部 z 列
zw_all = z_all(:).';

% 井轴直达项（全 z 解析式，复用 §1 存下的 k_in/aln）
pAd = zeros(N_w, numel(z_all));
hAd = zeros(N_w, numel(z_all));
for iw = 1:N_w
    k_in = k_in_all(iw,:).';
    aln  = aln_all(iw,:).';
    pa = zeros(1, numel(z_all));
    ha = zeros(1, numel(z_all));
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
t_all      = (0:N_full-1)*dt;                  % 时间轴 [s]

% 八个场（t × z，2048×4096）
%   to_tz：kz 域核（反射项、井壁直达）→ kz->z -> 乘谱 -> ω->t
%   to_t ：z 域量（井轴直达，本来就是全 z 解析式）→ 乘谱 -> ω->t
to_tz = @(X) integrate_omega(integrate_kz(X, N_k2, N_k, iz_all, d_z), N_t2, i_fc, dt, F, decay_full);
to_t  = @(X) integrate_omega(X, N_t2, i_fc, dt, F, decay_full);

P_axis_r = to_tz(pA_r);   P_axis_d = to_t(pAd);
P_wall_r = to_tz(pW_r);   P_wall_d = to_tz(pW_d);
H_axis_r = to_tz(hA_r);   H_axis_d = to_t(hAd);
H_wall_r = to_tz(hW_r);   H_wall_d = to_tz(hW_d);

% ======================================================================
function Zz = integrate_kz(Zk, N_k2, N_k, i_z, d_z)
% kz -> z：补负半轴镜像（复制式）后做 FFT，再取 i_z 指定的列
%   Zk   : (N_w, N_k)  kz 域核，只含 kz >= 0 的半谱
%   N_k2 : z 域点数（fft 补零到这里，长度 N_k2-1）
%   i_z  : 要取的 z 列索引（向量则返回全部对应 z）
Zfull = zeros(size(Zk,1), N_k2, 'like', Zk);        % 补零到 N_k2 点
Zfull(:, 1:N_k)          = Zk;                       % 正半轴
Zfull(:, N_k2+1-(2:N_k)) = Zk(:, 2:N_k);             % 负半轴：复制式镜像
rec_z = fft(Zfull(:, 1:N_k2-1).').';                 % kz -> z（正指数效果）
Zz = rec_z(:, i_z) * (1/(d_z*(N_k2-1)));             % 取列 + dk/(2π) 归一化
end


% ======================================================================
function Zt = integrate_omega(Zw, N_t2, i_fc, dt, F, decay)
% ω -> t：乘源谱 → 补负频共轭镜像 → fft → 乘衰减还原 → 取实部
Zw  = Zw .* F(:);
Xf = zeros(N_t2-1, size(Zw,2), 'like', Zw);
Xf(i_fc, :)         = Zw;                      % 正频段
Xf(N_t2-i_fc+1, :)    = conj(Zw);                % 负频共轭镜像
Zt = real(fft(Xf, [], 1) * (1/(dt*(N_t2-1))) .* decay);   % dω/(2π) + e^{σt} 还原
end
