function [D1, src_wall, An_wall, An_axis, an_wall, an_axis] = M8_solve(fr, param, kz, z_obs)
% M8_solve —— 8x8 边界方程在 (ω, kz) 处的一次组装
%
%   fr    ki_omega 的输出；param 里取 R、n；kz 轴向波数
%   z_obs 观测距离 [m]（默认 3，只给 an_wall/an_axis 的相位用）
%   s0    源强标定从 param.s0 取
%
%   D1       det(M8)，M8 = M(:,1:8)；nargout == 1 时直接返回（牛顿只要这个）
%   src_wall 1x2 井壁源列：src.*K_n(knr*r)，r = 接收器径向位置
%   An_wall  1x2 解级井壁：An .* I_n(knr*r) .* e^{|Re knr r| - |Re 2knr R|}（无相位）
%   An_axis  1x2 解级井轴：An .* (knr/2)^n .* e^{-|Re 2knr R|}（无相位）
%   an_wall  1x2 留数级井壁：A .* I_n(knr*r) .* e^{|Re knr r| - |Re 2knr R|}
%   an_axis  1x2 留数级井轴：A .* (knr/2)^n .* e^{-|Re 2knr R|}
%            其中 A = 1i*an.*exp(1i*kz*z_obs)，an 为留数级归一系数
%
%   内部量（不输出）：N2、An_o = (N2/D1).'、An = An_o/I_n(κ)、an_o = (N2/Dp).'、an = an_o/I_n(κ)

M  = build_M(kz, fr, param);
M8 = M(:, 1:8);

D1 = det(M8);
if nargout < 2
    return
end

n  = param.n;
s0 = param.s0;

R = param.R;                                  % 井半径 [m]（只用于衰减因子 e^{-|Re 2knr R|}）
r = param.r;                                  % 接收器径向位置 [m]

knr = sqrt(kz.^2 - fr.kni_vec(:).^2).';      % 1x2 井内两支径向波数
kap = knr*R;                                  % κ = k_nr*R

src    = s0 * (1 + (n > 0)) .* (knr*param.r0/2).^n / factorial(n);
src(2) = 0;
src_wall = src .* besselk(n, knr*r);          % 井壁源列（第 2 支已置 0）
y        = src .* besselk(n, kap, 1);         % 源列（缩放版 K_n，供 b 用）
b        = -(M(:,9:10))*y(:);

N2 = zeros(2,1);
for j = 1:2
    Mj = M8;
    Mj(:, 4+j) = b;
    N2(j) = det(Mj);
end

An_o = (N2/D1).';
An   = An_o ./ besseli(n, kap, 1);

if nargout >= 3
    wall_fac = besseli(n, knr*r, 1) .* exp(abs(real(knr*r)) - abs(real(knr*2*R)));
    axis_fac = (knr/2).^n .* exp(-abs(real(knr*2*R)));

    An_wall = An .* wall_fac;
    An_axis = An .* axis_fac;
end
if nargout >= 5
    if nargin < 4, z_obs = 3; end

    % Dp = d det(M8)/dk —— 前向差分（一次额外组装）
    hz  = 1e-6*max(abs(kz), 1e-3);
    Mm  = build_M(kz - hz, fr, param);   Mm = Mm(:,1:8);
    Dp  = (D1 - det(Mm)) / hz;

    an_o = (N2/Dp).';
    an   = an_o ./ besseli(n, kap, 1);
    A    = 1i * an .* exp(1i*kz*z_obs);       % 留数 × 相位

    an_wall = A .* wall_fac;
    an_axis = A .* axis_fac;
end
end
