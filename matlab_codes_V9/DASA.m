function [x, fvalue] =  DASA(A, y, l, u, lambda0, epsilon, OPTION)
% This function is for the Polyhedron Constrained Projction Problem:
%   min || y - x ||_2^2
%   s.t. l <= Ax <= b
%       x is a n-dimension vector
%       A is a m*n matrix
%       l and u are the lower and upper constrains vector with m-dimension.
% Notation: This function is designed to solve the dual problem of the above
%   max 0.5(y'*y - x(lambda)'*x(lambda)) + z(lambda)'*lambda
%   s.t. x(lambda) = max(y + A'*lambda, 0);
%        z(lambda)_j = l_j if lambda_j >0 ; 
%        z(lambda)_j = u_j if lambda_j < 0;
%% trandport the imput from OPTION

alphakmin = 1e-30; alphakmax = 1e+30;

% the max iterations for solveing the first subproblem, default is 2000
if isfield(OPTION, 'maxiter')
    maxiter = OPTION.maxiter;
else
    maxiter = 2000;
end

% whether plot (default:don't plot)
% 0  don't plot, 1 plot point during iterations, 2 plot after iteration
if isfield(OPTION, 'plot')
    plot_true = OPTION.plot;
else 
    plot_true = 0;
end

% whether use BB step size (default: use long BB step size)
% 0 don't use BB step size when solving the first subproblem
% 1 use the long BBstep size
% 2 use the short BBstep size
% 3 use adaptive BBstep size
% if BBstep = 3, the parameter for choice BB_l or BB_s must be given
if isfield(OPTION, 'BBstep')
    BBstep = OPTION.BBstep;
else
    BBstep = 1;
end
if BBstep == 3
    if exist('OPTION.BBstep_kappa','var')
        BBstep_kappa = OPTION.BBstep_kappa;
    else
    error([ 'parameter (OPTION.BBstep_kappa) for choicing long or short BB step size is avaliable']);
    end
end
    

%%
m = size(A, 1);
n = size(A, 2);
y_norm2 = norm(y,2);
x0 = max(y + A' * lambda0, 0);
pvalue = 0.5*norm(x0 - y, 2)^2;
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
if plot_true == 1
    figure
    % plot(k, sign(pvalue)*log10(pvalue), 'b*'); hold on;
    % plot(k, sign(dvalue)*log10(abs(dvalue)), 'r.');
end
%% main loop
delta = 0.7; % for linear search scale
sigma = 0.4; % for linear search adjust factor
xi = 0.9;    % nonmonotune linear search; Qk = xi*Qk + 1
dvalue_avg = dvalue; Qk= 1; % nonmonotune linear search
while k <= maxiter
    gradk_old = -A*xk_old;
%     gk_old = zeros(m, 1);
%     h_l_old = l - A*xk_old; h_u_old = u - A*xk_old;
%     gk_old(lambdak_old>0)  = h_l_old(lambdak_old>0);
%     gk_old(lambdak_old<=0) = h_u_old(lambdak_old<=0);
    if k >1
        gradk_old = -A*xk_old;
        gradk = -A*xk;
        alphakB = (gradk_old - gradk)' * (lambdak - lambdak_old) /norm(lambdak - lambdak_old, 2)^2;
        % 为什么BB步长计算为负数？
%         alphak0 = min(alphakmax, max(alphakmin, alphakB));
        alphak0 = alphakB;
        lambdak_old  = lambdak;
        xk_old = max(y + A'*lambdak_old, 0);
        h_l_old = l - A*xk_old; 
        h_u_old = u - A*xk_old;
        % sparsa
        % linear search for sparsa
        alphak =alphak0; 
        
        dvalue_old = D_funvalue(lambdak_old, A, y, l, u);
        dvalue_avg = ( xi*Qk*dvalue_avg + dvalue_old )/(xi*Qk + 1);
        while 1  % nonmonotone linear search
            lambdak_l = lambdak_old + h_l_old/alphak;
            lambdak_u = lambdak_old + h_u_old/alphak;
            
            lambdak = zeros(m,1);
            lambdak(lambdak_l>=0) = lambdak_l(lambdak_l>=0);
            lambdak(lambdak_u<=0) = lambdak_u(lambdak_u<=0);
            dvalue = D_funvalue(lambdak, A, y, l, u);
            if dvalue >=  dvalue_avg + 0.5 * sigma * alphak * norm(lambdak-lambdak_old)^2 %nonmonotone linear search
                Qk = xi*Qk + 1;
                break
            end
            alphak = 1/delta*alphak;
        end
    else % the first setp: search direction is the gradient at current point
        dk = gradk_old;
        alphak = 1;
        while 1  % linear search
            lambdak = lambdak_old + alphak * dk;
            if dvalue >=  dvalue_avg + 0.5 * sigma * alphak * norm(lambdak-lambdak_old)^2 %nonmonotone linear search
                Qk = xi*Qk + 1;
                break
            end
            alphak = delta*alphak;
        end
    end
    xk = max(y + A'*lambdak, 0);
    pvalue  = 0.5*norm(xk - y,2)^2; % save object value of primal  and dual problems
    dvalue = D_funvalue(lambdak, A, y, l, u);
    Pvalue(k) = pvalue; Dvalue(k) = dvalue;
    r = sqrt( 2*(pvalue - dvalue)/gamma1 );
    if plot_true == 1
        plot(k, pvalue, 'b-.');hold on;
        plot(k, dvalue, 'r:.');
    end
    fprintf(' %d     %5e   %5e     %5e   %3e\n',k, pvalue, dvalue, pvalue - dvalue, alphak)
    if pvalue < dvalue
        pause
    end
    if pvalue - dvalue < epsilon && pvalue > dvalue
        convergence = 1;
        fprintf('=== convergence====  r=%.3e; gap=%.4e \n', r, pvalue-dvalue); 
        break;
    else 
        convergence = 0;
    end
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
    if mod(k,5) == 1
        pause(0.5)
    end
end
if plot_true == 2
    figure;
    plot(Pvalue, 'b-+');hold on;
    plot(Dvalue, 'r-d');
end
if  plot_true > 0 && convergence == 1
    %     plot([0,k],[dvalue,dvalue],'m-.')
    line([0,k],[dvalue,dvalue],'Color',[0.1,0.1,0.1],'LineStyle', '-.')
    title_text = ['\epsilon : ',num2str(epsilon)];
    title(title_text)
    legend({'pvalue','dvalue'},'Location','best');
    hold off;
end

%% output
x = xk;
fvalue.p = pvalue;
fvalue.d = dvalue;
%% Dual function value
function dvalue = D_funvalue(lambdak, A, y, l, u)
xk = max(y + A'*lambdak, 0);
f   = 0.5*(y'*y - xk'*xk);
psi = l(lambdak>0)'*lambdak(lambdak>0) + u(lambdak<0)'*lambdak(lambdak<0);

% f = 0.5*norm(y-x)^2;
% z = l;
% z(lambdak>0) = l(lambdak>0);
% z(lambdak<0) = u(lambdak<0);
% psi = lambdak'*(z - A*x);
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

