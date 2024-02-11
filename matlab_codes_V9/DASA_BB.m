function [x, fvalue] =  DASA_BB(A, y, l, u, lambda0, epsilon)
% This function is for the Polyhedron Constrained Projction Problem:
% min || y - x ||_2^2
% s.t. l <= Ax <= b
%   x is a n-dimension vector
%   A is a m*n matrix
%   l and u are the lower and upper constrains vector with m-dimension.

%%
maxiter = 50;
m = size(A, 1);
n = size(A, 2);
y_norm2 = norm(y,2);
x0 = max(y + A' * lambda0, 0);
pvalue = norm(x0 - y, 2);
dvalue = D_funvalue(lambda0, A, y, l, u);
fprintf('\n iter   pvalue        dvalue            duality gap step-size \n');
fprintf(' %d     %5e   %5e     %5e   --\n',0, pvalue, dvalue, pvalue - dvalue)
%% preparation for main loop
lambdak_old = lambda0; k = 1;
xk_old = max(y + A'*lambdak_old, 0);
x_plus_indic = xk_old>0;
gamma1 = 1;
r = sqrt( 2*(pvalue - dvalue)/gamma1 );
A_norminf = norm(A, inf);
% only when the detection is safe(r is sufficient small -> the ragion is sufficient narrow)
if r < 0.5*min(u-l)/A_norminf
    T_p = A*xk_old + A_norminf * r;
    T_n = A*xk_old - A_norminf * r;
    AS_0 = (T_p<u) & (T_n>l);
    AS_p = (T_p>l) & (T_n<l);
    AS_n = (T_p>u) & (T_n<u);
    lambdak_old(AS_0) = 0; %screening the components to be zeros
end
dvalue_avg = dvalue; Qk= 1;
%
figure
% plot(k, sign(pvalue)*log10(pvalue), 'b*'); hold on;
% plot(k, sign(dvalue)*log10(abs(dvalue)), 'r.');
%% main loop
while k <= maxiter
%     gk_old = -A*xk_old;
    gk_old = zeros(m, 1);
    h_l_old = l - A*xk_old; h_u_old = u - A*xk_old;
    gk_old(lambdak_old>0)  = h_l_old(lambdak_old>0);
    gk_old(lambdak_old<=0) = h_u_old(lambdak_old<=0);
    if k >1
        alphakB = (A*xk_old - A*xk)' * (lambdak - lambdak_old) /norm(lambdak - lambdak_old, 2)^2;
        lambdak_old  = lambdak; xk_old = max(y + A'*lambdak_old, 0);
        h_l_old = l - A*xk_old; h_u_old = u - A*xk_old;
        % sparsa
        % linear search for sparsa
        alphak =alphakB; delta = 0.7; sigma = 0.4; xi = 0.9;
        dvalue_old = D_funvalue(lambdak_old, A, y, l, u);
        dvalue_avg = ( xi*Qk*dvalue_avg + dvalue_old )/(xi*Qk + 1);
        while 1
            lambdak_l = lambdak_old + h_l_old/alphak;
            lambdak_u = lambdak_old + h_u_old/alphak;
            lambdak = zeros(m,1);
            lambdak(lambdak_l>=0) = lambdak_l(lambdak_l>=0);
            lambdak(lambdak_u<=0) = lambdak_u(lambdak_u<=0);
            dvalue = D_funvalue(lambdak, A, y, l, u);
            if dvalue >=  dvalue_avg + 0.5 * sigma * alphak * norm(lambdak-lambdak_old, 2) %nonmonotone linear search
                Qk = xi*Qk + 1;
                break
            end
            alphak = 1/delta*alphak;
        end
    else
        dk = gk_old;
        alphak = 1 ;
        lambdak = lambdak_old + alphak * dk;
    end
    xk = max(y + A'*lambdak, 0);
    pvalue = norm(xk - y, 2);
    dvalue = D_funvalue(lambdak, A, y, l, u);
    %     r = sqrt( 2*(pvalue - dvalue)/gamma1 );
    plot(k, sign(pvalue)*log10(pvalue), 'b-*');hold on;
    plot(k, sign(dvalue)*log10(abs(dvalue)), 'r:.');
    legend({'pvalue','dvalue'},'Location','best')
    fprintf(' %d     %5e   %5e     %5e   %3e\n',k, pvalue, dvalue, pvalue - dvalue, alphak)
    %     if r < epsilon
    %         fprintf('=== convergence====  r:%5e \n', r); break; end
    % only when the detection is safe(r is sufficient small -> the ragion is sufficient narrow)
    if r < 0.5*min(u-l)/A_norminf && 0
        T_p = A*xk_old + A_norminf * r;
        T_n = A*xk_old - A_norminf * r;
        AS_0 = (T_p<u) & (T_n>l);
        AS_p = (T_p>l) & (T_n<l);
        AS_n = (T_p>u) & (T_n<u);
        lambdak_old(AS_0) = 0; %screening the components to be zeros
    end
    k = k+1;
end
hold off
x = xk;
fvalue.p = pvalue;
fvalue.d = dvalue;
%% Dual function value
function dvalue = D_funvalue(lambda, A, y, l, u)
x = max(y + A'*lambda, 0);
f   = 0.5*(y'*y - x'*x);
psi = l(lambda>0)'*lambda(lambda>0) + u(lambda<0)'*lambda(lambda<0);
dvalue = f + psi;

%% BB step point
function lambdak_plus = BBpoint(lambdak, s, w, gk, AS_p, AS_n)
alpha_min = 1e-30; alpha_max = 1e+30;
alphaBl = norm(s, 2)/(s'*w);
alphaBs = norm(s'*w)/norm(w, 2);
if alphaBs/alphaBl < 0.2
    alphaBk = min(alpha_max, max(alpha_min, alphaBs));
else
    alphaBk = min(alpha_max, max(alpha_min, alphaBl));
end
lambdak_plus = lambdak - alphaBk*gk;
lambdak_plus = max(lambdak_plus, 0).*AS_p + min(lambdak_plus, 0).*AS_n;

%% Truncation operator Trunc
function lambdak_new = Trunc(lambdak, AS_p, AS_n)
lambdak_new = zeros(size(lambdak));
lambdak_new(AS_p) = max(lambdak(AS_p ), 0);
lambdak_new(AS_n) = min(lambdak(AS_n ), 0);

