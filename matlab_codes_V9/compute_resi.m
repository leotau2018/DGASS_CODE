%% compute residuals
%     min  0.5*||x-y||^2,    s.t. l <= A*x <= u.
function resi = compute_resi(x, lambda, A, l, u)
Ax = A*x;
lambda_l = lambda+l-Ax;
lambda_u = lambda+u-Ax;
lambda_new = zeros(length(lambda),1);
lambda_new(lambda_l>0) = lambda_l(lambda_l>0);
lambda_new(lambda_u<0) = lambda_u(lambda_u<0);
resi = norm(lambda_new - lambda, inf);
end