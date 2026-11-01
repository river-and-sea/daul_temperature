function F = source_spec_fw(fvec, fc, tc)
% source_spec_fw —— 源频谱 F(ω)，逐字照刘建宇 pf3.m + Sa.m
%
%   F(ω) = (tc/4)·[ Sa(cpt+p0t) + Sa(cpt-p0t)
%                   + 0.5·( Sa(cpt+a) + Sa(cpt-a) + Sa(cpt+b) + Sa(cpt-b) ) ]·exp(i·cpt)
%
%   cpt = π·tc·f（无量纲"频率×tc"）、p0t = π·tc·fc、a = p0t+π、b = p0t-π
%   Sa(x) = sin(x)/x（未归一化，Sa(0)=1）
%
% fvec 可为复数（复频口径）；本函数本身不必改，时域那一步补乘 e^{σt}。

cpt = pi*tc*fvec;
p0t = pi*tc*fc;
a   = p0t + pi;
b   = p0t - pi;
F = 0.25*tc*( Sa(cpt+p0t) + Sa(cpt-p0t) ...
    + 0.5*( Sa(cpt+a) + Sa(cpt-a) + Sa(cpt+b) + Sa(cpt-b) ) ) .* exp(1i*cpt);
end

% -------------------------------------------------------------------------
function y = Sa(x)
% Sa(x) = sin(x)/x，x=0 处取 1（未归一化；不等于 MATLAB 的 sinc(x)）
y = ones(size(x));
z = (x ~= 0);
y(z) = sin(x(z))./x(z);
end
