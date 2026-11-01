function e = excite(fr, param, kz, omega)
%   s0 内置为 1e-6（原 fw_src_amp）
% excite —— 极点处的激发留数（张碧星多极源；本次只加注释，代码未改动）
%
% ── 输入 ────────────────────────────────────────────────────────────────
%   fr     ki_omega 的输出，本函数只取两样：
%            fr.kni_vec   井内两支模式波数 k_{n,i}（2x1）
%            fr.alphani   井内两支温度/势比 alpha_{n,i}（2x1）
%   param  param_0 的静态量，本函数只取：
%            param.n      方位序数 n（Bessel 阶）
%            param.a0     多极源分布半径 a0 [m]
%            param.R      井半径 R [m]
%   kz     该支轴向波数（极点，由 kz_solve 求得）
%   s0     源强标定 fw_src_amp（1e-6），透传给 M8_solve 组源列
%   omega  角频率（当前代码里已不出现；ap_wall 的 rho_f*omega^2 因子待补）
%
% ── 数学脉络 ────────────────────────────────────────────────────────────
%   8x8 边界方程在极点：M8*x = b，克莱姆给出
%       Res(x_j) = N_j / D',   N_j = det(M8 第 j 列 <- b),  D' = d det(M8)/dk
%   M8_solve 一次组装给出 D1=det(M8)、N2=[N_5;N_6]、Dp=D'。
%   井内只有第 5、6 两支参与"声"，故本函数只用 N2/Dp。
%
%   井内是内部解，用 I_n（第一类 Bessel）展开：
%       轴上 r=0  ->  场按 (k_nr*r/2)^n 展开（n 阶导数因子）
%       井壁 r=R  ->  再乘 I_n(k_nr*R)
%   源（张碧星式 (5)）带多极因子 (k_nr*a0/2)^n/n! 与相位 e^{i kz z}。
%
% ── 输出（全部标量，两支已叠加）─────────────────────────────────────────
%   e.ap_wall  井壁声压激发      sum(A_wall)          （待乘 rho_f*omega^2）
%   e.at_wall  井壁温度激发      sum(A_wall.*alpha)
%   e.ap_axis  井轴声压激发      sum(A_axis)
%   e.at_axis  井轴温度激发      sum(A_axis.*alpha)
%   e.ok       留数是否可信（A 全有限）

% ── 1. 一次组装 8x8，取克莱姆分子与留数级反射系数 ───────────────────────
%   D1 = det(M8)：本函数不用（用 ~ 吃掉，ASGLU 抑制 lint）
%   N2 = [det(M8 第5列<-b); det(M8 第6列<-b)]，2x1
%   an = (N2/Dp).'/I_n(κ).*exp(-|Re κ|)，1x2（留数级反射系数）
[~, ~, ~, An, an_o, an] = M8_solve(fr, param, kz, 1e-6);   %#ok<ASGLU>  s0 内置 = fw_src_amp
%   An / an 归一化后的两支系数（1x2），见 §5 用法
% ── 2. 输出结构先占位（字段顺序固定；无效时保持空/ false）────────────────
e = struct('ap_wall', [], 'at_wall', [], 'ap_axis', [], 'at_axis', [], 'ok', false);

% ── 3. 局部量 ───────────────────────────────────────────────────────────
n   = param.n;                              % 方位序数（Bessel 阶）
z   = 3;                                    % 探子高度 [m]（原程序 fw_zp）
r   = param.R;                              % 井半径 [m]
knr = sqrt(kz.^2 - fr.kni_vec(:).^2).';     % 1x2 井内两支径向波数 k_nr=sqrt(kz^2-k_n^2)
aln = fr.alphani(:).';                      % 1x2 温度/势比，转行向量与 knr 同向

% ── 4. 井轴振幅 A（1x2，按支）───────────────────────────────────────────
%   ① an                  井内两支的反射系数留数（M8_solve 第 6 输出，1x2）
%   ② .* (knr*a0/2)^n/n!  张碧星式 (5) 的多极源因子（n=0 时为 1）
%   ③ .* exp( 1i*kz*z)     沿 z 传播相位 e^{i kz z}

% ── 5. 两种观测位置的场 ────────────────────────────────────────────────
A = an.* (knr*param.a0/2).^n / factorial(n) .* exp(1i*kz*z);
A_wall = A .* besseli(n, knr * r, 1).*exp(abs(real(knr * (r)))-abs(real(knr * 2*param.R)));   % 井壁：乘 I_n(knr*R)（轴到壁的传播）
A_axis = A .* (knr/2).^n .* exp(-abs(real(knr * 2*param.R)));   % 张碧星(5)

% ── 6. 四个输出：两支叠加成标量 ────────────────────────────────────────
e.ap_wall = sum(A_wall);                % 1x2，按支其实该*rho_f*omega^2
e.at_wall = sum(A_wall .* aln);         % 1x2，按支
e.ap_axis = sum(A_axis);
e.at_axis = sum(A_axis .* aln);

% ── 7. 有效性：极点处 Dp 为 0/NaN 时 an_o 会是 Inf/NaN，此处即 false ─────
e.ok = all(isfinite(an_o));

% ── 备注（未改动代码）──────────────────────────────────────────────────
%   ① ap_wall 尚缺 rho_f*omega^2（行内注释已标"其实该*"）。
%   ② 若探子不在井壁（r≠R），besseli 的宗量已随 r 变，无需另改；
%      但 z 目前写死 2 m，多道时需改成参数。
%   ③ 大 |knr*R| 时 besseli 未用缩放版（第三参 1），可能溢出。
end
