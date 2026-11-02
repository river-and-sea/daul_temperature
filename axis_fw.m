% axis_fw —— 井轴 / 井壁全波（合并 axis_fw 与 fullwave_fw）
%
% 脚本程序：直接运行本文件即可。单道（z），共用同一次 M8_solve。
%
%   井轴（A，轴压 / 张碧星式 (9)，对应旧 zbx_axis_pressure）：
%       p = chengzi * sum( knr^(2n) * An )        An 已按 I_n 归一
%       直达项用 k_in 的解析式直接在 z 上求值
%
%   井壁（B，全波，对应旧 zbx_fullwave）：
%       p = rho_f*om^2 * (src + x(5:6))           src 带 (knr^2*r0/4)^n*e^{-κ}
%
%   ── 命名规范 ──────────────────────────────────────────────────────────
%   核（ω,kz 域）：{p,h}{A,B}_{r,d}        小写 p/h，A=井轴 B=井壁，r=反射 d=直达
%   变换后      ：{P,H}{A,B}_{r,d}_z        kz->z
%   时域        ：{P,H}{A,B}_{r,d}_t        ω->t
%   输出 out    ：{P,H}_{axis,wall}[_r,_d]  总场与分量
%   两套口径都输出声压与温度；总场 = 反射 + 直达。

%% ============ 0. 参数 ============
param  = param_0();                          % 五参数用默认
n      = param.n;
z  = 3;                                  % 唯一的接收位置 [m]

% ---- 网格常量（原 dt_config.fw_*）----
tc        = 0.5e-3;             % 源脉冲特征时间 [s]
fw_fc_list   = [6000 5000 6000];   % 各 n 的源中心频率 [Hz]
fw_n_list    = [0 1 2];            % 允许的 n
dt        = 10e-6;              % 时间步 [s]
N_t        = 2049;               % 频域点数
N_kz        = 4097;               % k 域点数
f0_ratio = 1/100;              % 复频偏移比例
d_z         = 1/48;                              % 空间步 [m]（任取，只要 z 落在网格点）

% ============ 栅格 ============
% ---- 源 ----
f_src = fw_fc_list(fw_n_list == n);        % 源中心频率 [Hz]


% ---- 空间轴 / 波数轴（唯一约束：dk*d_z*(n1-1) = 2π）----
i_z = round(z/d_z);                      % z 落在第几个网格点
z     = i_z*d_z;                         % 对齐后的接收位置 [m]

N_k = fix(N_kz/2) + 1;                       % 波数正半轴点数
d_k  = 2*pi/(d_z*(N_kz-1));                   % 波数步 [1/m]，由 d_z 的对偶关系定出
kz = linspace(0, (N_k-1)*d_k, N_k) + 1e-36;

% ---- 时间轴 ----
N_tt = fix(0.5*(fix(N_t/2)+1));             % 输出时间点数（N_t 是频域点数，别覆盖）
i_t  = 1:N_tt;
time    = linspace(0, (N_tt-1)*dt, N_tt);

% ---- 频率轴 ----
df  = 1/(dt*(N_t-1));
iw_f = (round((f_src-2/tc)/df)+1) : round((f_src+2/tc)/df);
N_w  = numel(iw_f);
freq = iw_f*df;

% 复频偏移：路径抬到 w + i*sigma；时域乘回 e^{+sigma*t}
freq_c = freq + 1i*f_src*f0_ratio;
om_c   = 2*pi*freq_c;
F      = source_spec_fw(freq_c, f_src, tc);

%% ============ 1. 逐 (ω, kz)：两套口径一起算 ============
pA_r = zeros(N_w, N_k);   pA_d = zeros(N_w, 1);   % 井轴：声压 反射核 / 直达（z 域）
hA_r = zeros(N_w, N_k);   hA_d = zeros(N_w, 1);   % 井轴：温度
pW_r = zeros(N_w, N_k);   pW_d = zeros(N_w, N_k); % 井壁：声压 反射核 / 直达核
hB_r = zeros(N_w, N_k);   hB_d = zeros(N_w, N_k); % 井壁：温度

r0   = param.r0;
angf = 1;                                  % cos(n*theta_rec)，theta_rec = 0
chengzi = param.s0*(1 + (n > 0)) * (r0/4)^n / (2*factorial(n));   % 井轴口径前因子

for iw = 1:N_w
    fr    = ki_omega(om_c(iw), param, []); % 复频下的解
    om    = om_c(iw);
    k_in  = fr.kni_vec(:).';               % 1x2 井内两支轴向波数（行）
    aln   = fr.alphani(:);                 % 2x1 温度/势比

    % ---- 井轴直达项（式 (9) 第一项，z 域解析式）----
    p_acc = 0;
    h_acc = 0;
    for mm = 0:n
        term = exp(1i*(k_in*z)) .* ((-2i*k_in).^mm) ./ (z.^(2*n-mm+1));
        coef = nchoosek(n,mm)*factorial(2*n-mm);
        p_acc = p_acc + coef .* sum(term);
        h_acc = h_acc + coef .* sum(term .* aln.');
    end
    pA_d(iw) = chengzi .* p_acc;
    hA_d(iw) = chengzi .* h_acc;

    for ik = 1:N_k
        %   An_o = x(5:6).'（未归一），  An = An_o / I_n（已归一）
        %   src  = 源列本身（已含 (knr*r0/2)^n/n! 与 eps_n，第 2 支置 0）
        [~, src_wall, An_wall, An_axis] = M8_solve(fr, param, kz(ik));

        % ---- 井轴（A）：权重 × 归一系数 ----
        %   权重 wA = chengzi * knr^(2n)（不含 rho_f*om^2，声压/温度共用）
        pA_r(iw,ik) = angf *sum(An_axis);                      % 反射
        hA_r(iw,ik) = angf *sum(An_axis .* aln.');

        % ---- 井壁（B）：权重 × (源 + 反射) ----
        %   井壁场 = 直达 src + 反射 x(5:6)，两者同尺度，直接相加
        pW_d(iw,ik) = angf * sum(src_wall);
        pW_r(iw,ik) = angf * sum(An_wall);
        hB_d(iw,ik) = angf * sum(src_wall .* aln.');
        hB_r(iw,ik) = angf * sum(An_wall .* aln.');
    end
end

%% ============ 2. 双重 FFT -> 时域（单道）============
% kz -> z
PA_r_z = integrate_kz(pA_r, N_kz, N_k, i_z, d_z) .* F(:);
HA_r_z = integrate_kz(hA_r, N_kz, N_k, i_z, d_z) .* F(:);
PA_d_z = pA_d .* F(:);                % 井轴直达项本来就在 z 域
HA_d_z = hA_d .* F(:);

PW_r_z = integrate_kz(pW_r, N_kz, N_k, i_z, d_z) .* F(:);
PW_d_z = integrate_kz(pW_d, N_kz, N_k, i_z, d_z) .* F(:);
HB_r_z = integrate_kz(hB_r, N_kz, N_k, i_z, d_z) .* F(:);
HB_d_z = integrate_kz(hB_d, N_kz, N_k, i_z, d_z) .* F(:);

% ω -> t
PA_r_t = integrate_omega(PA_r_z, N_t, iw_f, dt);
PA_d_t = integrate_omega(PA_d_z, N_t, iw_f, dt);
HA_r_t = integrate_omega(HA_r_z, N_t, iw_f, dt);
HA_d_t = integrate_omega(HA_d_z, N_t, iw_f, dt);

PW_r_t = integrate_omega(PW_r_z, N_t, iw_f, dt);
PW_d_t = integrate_omega(PW_d_z, N_t, iw_f, dt);
HB_r_t = integrate_omega(HB_r_z, N_t, iw_f, dt);
HB_d_t = integrate_omega(HB_d_z, N_t, iw_f, dt);

decay = exp(2*pi*f_src*f0_ratio*i_t*dt).';   % 复频偏移的时域还原

% 取输出时段
PA_r_t = real(PA_r_t(i_t) .* decay);   PA_d_t = real(PA_d_t(i_t) .* decay);
HA_r_t = real(HA_r_t(i_t) .* decay);   HA_d_t = real(HA_d_t(i_t) .* decay);
PW_r_t = real(PW_r_t(i_t) .* decay);   PW_d_t = real(PW_d_t(i_t) .* decay);
HB_r_t = real(HB_r_t(i_t) .* decay);   HB_d_t = real(HB_d_t(i_t) .* decay);

% 总 = 反射 + 直达
P_axis = PA_r_t + PA_d_t;   H_axis = HA_r_t + HA_d_t;
P_wall = PW_r_t + PW_d_t;   H_wall = HB_r_t + HB_d_t;

out = struct( ...
    't',        time,       ...             % 时间轴 [s]
    'z',        z,     ...             % 接收位置 [m]
    'P_axis',   P_axis,    ...             % 井轴声压总场
    'P_axis_r', PA_r_t,    ...             % 井轴声压反射
    'P_axis_d', PA_d_t,    ...             % 井轴声压直达
    'H_axis',   H_axis,    ...             % 井轴温度总场
    'H_axis_r', HA_r_t,    ...
    'H_axis_d', HA_d_t,    ...
    'P_wall',   P_wall,    ...             % 井壁声压总场
    'P_wall_r', PW_r_t,    ...             % 井壁声压反射
    'P_wall_d', PW_d_t,    ...             % 井壁声压直达
    'H_wall',   H_wall,    ...             % 井壁温度总场
    'H_wall_r', HB_r_t,    ...
    'H_wall_d', HB_d_t,    ...
    'pA_r',     pA_r,      ...             % 井轴声压 (ω,kz) 核
    'pW_r',     pW_r,      ...             % 井壁声压 (ω,kz) 核
    'n',        n,         ...
    'pc',       param.pc,  ...
    'fc',       f_src,     ...
    'freq',     freq,      ...
    'kz',       kz);

fprintf('n=%d pc=%g MPa  z=%.4f m\n', n, param.pc/1e6, z);
fprintf('  井轴: max|P_r|=%.4g  max|P_d|=%.4g  max|P|=%.4g\n', ...
    max(abs(out.P_axis_r)), max(abs(out.P_axis_d)), max(abs(out.P_axis)));
fprintf('  井壁: max|P_r|=%.4g  max|P_d|=%.4g  max|P|=%.4g\n', ...
    max(abs(out.P_wall_r)), max(abs(out.P_wall_d)), max(abs(out.P_wall)));
fprintf('  温度: 井轴 max=%.4g  井壁 max=%.4g\n', max(abs(out.H_axis)), max(abs(out.H_wall)));

%% ============ 3. 画 12 张图：温压 × 轴壁 × (总/直/反) ============
t_plot = out.t;                            % [s]
fig = figure('Visible', 'off', 'Position', [100 100 1400 1100]);

ycol = {'total', 'direct', 'reflected'};
labs = {'P_{axis}', 'P_{wall}', 'H_{axis}', 'H_{wall}'};
dat  = { {out.P_axis, out.P_axis_d, out.P_axis_r}, ...
         {out.P_wall, out.P_wall_d, out.P_wall_r}, ...
         {out.H_axis, out.H_axis_d, out.H_axis_r}, ...
         {out.H_wall, out.H_wall_d, out.H_wall_r} };

for i = 1:4
    for j = 1:3
        subplot(4, 3, (i-1)*3 + j);
        plot(t_plot, dat{i}{j}, 'LineWidth', 1);
        grid on; xlabel('t [s]'); ylabel(labs{i});
        title(sprintf('%s  %s', labs{i}, ycol{j}));
    end
end
sgtitle(sprintf('n=%d, z=%.3f m, pc=%g MPa', n, z, param.pc/1e6));

outdir = fullfile(fileparts(mfilename('fullpath')), 'fig');
if exist(outdir, 'dir') ~= 7, mkdir(outdir); end
saveas(fig, fullfile(outdir, 'fig_axis_fw.png'));
close(fig);
fprintf('saved: %s\n', fullfile(outdir, 'fig_axis_fw.png'));

%% ============ 4. t-z 时空图（全波 = 直达 + 反射）============
z_all = (0:N_kz-2)*d_z;                          % fft 后的全 z 轴 [m]
iz    = find(z_all >= 0.5 & z_all <= 6);       % 观察窗
z_win = z_all(iz);
zw    = z_win(:).';                            % 1×nz

% ---- 井轴直达项：在 z 窗上解析求值（式 (9) 第一项）----
PAd_tz = zeros(N_w, numel(zw));
HAd_tz = zeros(N_w, numel(zw));
for iw = 1:N_w
    fr   = ki_omega(om_c(iw), param, []);
    k_in = fr.kni_vec(:).';                    % 1×2
    aln  = fr.alphani(:);                      % 2×1
    p_acc = zeros(1, numel(zw));
    h_acc = zeros(1, numel(zw));
    for mm = 0:n
        term = exp(1i*(k_in.'*zw)) .* ((-2i*k_in.').^mm) ./ (zw.^(2*n-mm+1));   % 2×nz
        coef = nchoosek(n,mm)*factorial(2*n-mm);
        p_acc = p_acc + coef .* sum(term, 1);
        h_acc = h_acc + coef .* sum(term .* aln, 1);
    end
    PAd_tz(iw,:) = chengzi .* p_acc;
    HAd_tz(iw,:) = chengzi .* h_acc;
end

% ---- 四通道反射项：kz -> z（取窗口列）-> 乘源谱 ----
PAr_tz = integrate_kz(pA_r, N_kz, N_k, iz, d_z) .* F(:);
PWr_tz = integrate_kz(pW_r, N_kz, N_k, iz, d_z) .* F(:);
HAr_tz = integrate_kz(hA_r, N_kz, N_k, iz, d_z) .* F(:);
HBr_tz = integrate_kz(hB_r, N_kz, N_k, iz, d_z) .* F(:);

% 井壁直达项：kz 域核 pW_d/hB_d 也走 kz -> z
PWd_tz = integrate_kz(pW_d, N_kz, N_k, iz, d_z) .* F(:);
HBd_tz = integrate_kz(hB_d, N_kz, N_k, iz, d_z) .* F(:);

% 井轴直达项：z 域直接乘源谱
PAd_tz = PAd_tz .* F(:);
HAd_tz = HAd_tz .* F(:);

% ---- ω -> t ----
tmp = integrate_omega(PAr_tz + PAd_tz, N_t, iw_f, dt);   PAw = real(tmp(i_t,:) .* decay);
tmp = integrate_omega(PWr_tz + PWd_tz, N_t, iw_f, dt);   PWw = real(tmp(i_t,:) .* decay);
tmp = integrate_omega(HAr_tz + HAd_tz, N_t, iw_f, dt);   HAw = real(tmp(i_t,:) .* decay);
tmp = integrate_omega(HBr_tz + HBd_tz, N_t, iw_f, dt);   HBw = real(tmp(i_t,:) .* decay);

fig2 = figure('Visible', 'off', 'Position', [100 100 1100 800]);
tzs = {abs(PAw), abs(PWw), abs(HAw), abs(HBw)};
ttl = {'|P_{axis}(t,z)|', '|P_{wall}(t,z)|', '|H_{axis}(t,z)|', '|H_{wall}(t,z)|'};
for k = 1:4
    subplot(2,2,k);
    imagesc(t_plot, z_win, tzs{k}.');
    axis xy; colorbar;
    xlabel('t [s]'); ylabel('z [m]');
    title(ttl{k});
end
saveas(fig2, fullfile(outdir, 'fig_axis_fw_tz.png'));
close(fig2);
fprintf('saved: %s\n', fullfile(outdir, 'fig_axis_fw_tz.png'));


% ======================================================================
function Zz = integrate_kz(Z, n1, n_k, i_z, d_z)
% kz -> z：补负半轴镜像后做 FFT，只取 m_rec 这一点
Zk = zeros(size(Z,1), n1, 'like', Z);
Zk(:, 1:n_k)          = Z;
Zk(:, n1+1-(2:n_k))   = Z(:, 2:n_k);          % 负半轴：复制式镜像
rec_z = fft(Zk(:, 1:n1-1).').';               % kz -> z
Zz = rec_z(:, i_z) * (1/(d_z*(n1-1)));      % 单道 + dk/(2π) 归一化
end


% ======================================================================
function Zt = integrate_omega(Z, n2, iw_f, dt)
% ω -> t：补负频共轭镜像后做 FFT
Xf = zeros(n2-1, size(Z,2), 'like', Z);
Xf(iw_f, :)         = Z;                      % 正频段
Xf(n2-iw_f+1, :)    = conj(Z);                % 负频共轭镜像
Zt = fft(Xf, [], 1) * (1/(dt*(n2-1)));        % dω/(2π) 归一化
end
