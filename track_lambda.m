function lambda_sorted = track_lambda(lambda_raw, prev_lambda)
% track_lambda —— 把 roots() 返回的乱序根配到固定的 4 条支上（沿支追踪）
%
%   lambda_raw   本频点 roots(C_current) 的原始输出（顺序任意，4x1）
%   prev_lambda  上一频点的 4 个根（即上一次本函数的输出）；首点传 []
%                首点无历史，退化为按 real(lambda) 排序
%
% 输出 lambda_sorted：与 prev_lambda 同一支序的 4 个根。
% 逐行照抄 run_dispersion.m 507-525 行，只把 model 换成 numel(prev_lambda)。

if isempty(prev_lambda)
    [~, idx] = sort(real(lambda_raw));
    lambda_sorted = lambda_raw(idx);
    return
end

n = numel(prev_lambda);
lambda_sorted = zeros(n,1);
used = false(n,1);
for jj = 1:n
    dist = abs(lambda_raw - prev_lambda(jj));   % 每个新根到"同支上一频点根"的距离
    dist(used) = inf;                           % 已配出去的不再参与
    [~, id] = min(dist);
    lambda_sorted(jj) = lambda_raw(id);
    used(id) = true;
end
end
