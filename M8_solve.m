function [D1, N2, An_o, An, an_o, an] = M8_solve(fr, param, kz, s0)
% M8_solve —— 8x8 边界方程在 (ω, kz) 处的一次组装
%
%   fr    ki_omega 的输出；param 里取 R、n；kz 轴向波数
%   s0    源强标定（默认 1e-6）
%
%   D1   det(M8)，M8 = M(:,1:8)
%        牛顿迭代只要这一个：nargout == 1 时直接返回，不算后面的
%   N2   2x1：把 M8 第 5、6 列分别换成源列 b 后的行列式
%        N2(j) = det( M8(:, 4+j) <- b )，于是 x(5:6) = N2/D1（克莱姆）
%   An_o 1x2  (N2/D1).'       解级、未按 I_n 归一
%   An   1x2  An_o ./ I_n(κ)  解级、已归一   ← 不需要 Dp，全波直接用这个
%   an_o 1x2  (N2/Dp).'       留数级、未归一 ← 从这里开始才算 Dp
%   an   1x2  an_o ./ I_n(κ)  留数级、已归一
%        （o 后缀 = original，未归一）

M  = build_M(kz, fr, param);
M8 = M(:, 1:8);

D1 = det(M8);
if nargout < 2
    return
end

if nargin < 4, s0 = 1e-6; end

n = param.n;
R = param.R;

knr = sqrt(kz.^2 - fr.kni_vec(:).^2).';      % 1x2 井内两支径向波数
kap = knr*R;                                  % κ = k_nr*R

src    = s0 .* besselk(n, kap, 1);
src(2) = 0;
b      = -(M(:,9:10))*src(:);

N2 = zeros(2,1);
for j = 1:2
    Mj = M8;
    Mj(:, 4+j) = b;
    N2(j) = det(Mj);
end

if nargout >= 3
    An_o = (N2/D1).';
end
if nargout >= 4
    An = An_o ./ besseli(n, kap, 1);          % 只需 D1，不算 Dp
end
if nargout >= 5
    % Dp = d det(M8)/dk —— 前向差分（一次额外组装）
    hz  = 1e-6*max(abs(kz), 1e-3);
    Mm  = build_M(kz - hz, fr, param);   Mm = Mm(:,1:8);
    Dp  = (D1 - det(Mm)) / hz;

    an_o = (N2/Dp).';
end
if nargout >= 6
    an = an_o ./ besseli(n, kap, 1);
end
end
