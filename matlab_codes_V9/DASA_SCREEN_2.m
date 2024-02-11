function [x, lambda, value] =  DASA_SCREEN_2(A, y, l, u, lambda0, varargin)
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
%
% IMPUT:
%       y: n-dim vector to be approximated in the primal problem
%       A: m*n matrix, which induce the polyhedron and the constraints problems
%       l: m-dim vector, represents the lower constrint bounds
%       u: m-dim vector, represents the upper constrint bounds
%       lambda0: starting-variable of the dual problem
% OPTION: struct object, transmits the optional varargins,
%       OPTION.epsilon:
%           the stop condition to terminate the iteration(default=1e-4)
%       OPTION.maiter:
%           max iteration number(default=2000)
%       OPTION.plot:
%           whether plot (default:don't plot)
%           0  don't plot, 1  plot point during iterations, 2  plot after iteration
%       OPTION.BBstep:
%           whether use BB step size (default: use long BB step size)
%           0  don't use BB step size when solving the first subproblem
%           1  use the long BBstep size
%           2  use the short BBstep size
%           3  use adaptive BBstep size
%           if BBstep = 3, the parameter for choice BB_l or BB_s is
%           better given
%       OPTION.BBstep_kappa:
%           the parameter for choice BB_l or BB_s, default = 0.2.
%       OPTION.active_method:
%           whether use active set method
%           0  default, don't use activ method
%           1  screening method, only screeing the doomed zero variables, AS0 set
%           2  screen the zero variables and at the same time gives the sign of other components

%% trandport the imput from OPTION

alphakmin = 1e-30; alphakmax = 1e+30;
if isempty(varargin)
    OPTION = [];
else
    OPTION = varargin{:};
end
% the stop condition to terminate the iteration, default is 1e-4
if isfield(OPTION, 'epsilon')
    epsilon = OPTION.epsilon;
else
    epsilon = 1e-4;
end
% the max iterations for solveing the first subproblem, default is 2000
if isfield(OPTION, 'maxiter')
    maxiter = OPTION.maxiter;
else
    maxiter = 2000;
end

% whether plot (default:don't plot)
% 0  don't plot, 1  plot point during iterations, 2  plot after iteration
if isfield(OPTION, 'plot')
    plot_true = OPTION.plot;
else
    plot_true = 0;
end

% whether use BB step size (default: use long BB step size)
% 0  don't use BB step size when solving the first subproblem
% 1  use the long BBstep size
% 2  use the short BBstep size
% 3  use adaptive BBstep size
% if BBstep = 3, the parameter for choice BB_l or BB_s is better given
if isfield(OPTION, 'BBstep')
    BBstep = OPTION.BBstep;
else
    BBstep = 1;
end
if BBstep == 3
    if exist('OPTION.BBstep_kappa','var')
        BBstep_kappa = OPTION.BBstep_kappa;
    else
        disp('parameter for choicing long or short BB step size is not given, default(0.2) is used')
        BBstep_kappa = 0.2;
    end
end

% whether use active set method
% 0, default, don't use activ method
% 1  only screen the doomed zero variables, AS0 set
% 2  screen the zero variables and given the sign of other components
if isfield(OPTION, 'active_method')
    active_method = OPTION.active_method;
    if isfield(OPTION, 'active_epoch')
        active_epoch = OPTION.active_epoch;
    else
        active_epoch = 5;
    end
else
    active_method = 0;
end




%% preparation for main loop
A_old = A; y_old = y; l_old = l; u_old = u;
m = size(A, 1);
n = size(A, 2);
y_norm2 = norm(y,2);
% A_norminf = norm(A, inf);
AAt = A*A';
A_rownorm = sqrt(diag(AAt));
x0 = max(y + A' * lambda0, 0);
pvalue = 0.5*norm(x0 - y, 2)^2;
dvalue = D_funvalue(lambda0, A, y, l, u);

lambdak_old = lambda0; k = 1;
xk_old = max(y + A'*lambdak_old, 0);
x_plus_indic = [1:n];
lambda_unzero_indic = [1:m];
gamma1 = 1;
pvalue_min = pvalue; dvalue_max = dvalue;
r = sqrt( 2*(pvalue_min - dvalue_max)/gamma1 );
% only when the detection is safe(r is sufficient small -> the ragion is sufficient narrow)
% if r < 0.5*min(u-l)./A_rownorm

if active_method == 1
    T_p = A*xk_old + A_rownorm * r;
    T_n = A*xk_old - A_rownorm * r;
    AS_0 = (T_p<u) & (T_n>l);
    AS_p = (T_p>l) & (T_n<l);
    AS_n = (T_p>u) & (T_n<u);
    %     lambdak_old(AS_0) = 0; %screening the components to be zeros
    %     lambdak(AS_0) = 0; %screening the components to be zeros
    lambda_unzero_indic(AS_0) = [];
    A(AS_0,:) = [] ;
    l(AS_0) = []; u(AS_0) = [];
    lambdak_old = lambdak_old(~AS_0);
    xk_old = max(y + A'*lambdak_old, 0);
else
    AS_0 = zeros(size(lambdak_old));
end
if plot_true == 1
    figure
end
fprintf('\n---------------------------\n')
fprintf(' row number of A: %d \n column number of A: %d \n [max,mean,min] of A''s row norm=[%0.2f, %0.2f, %0.2f]',m, n, full(max(A_rownorm)), full(mean(A_rownorm)), full(min(A_rownorm)));
fprintf('\n---------------------------')
fprintf('\n iter    pvalue        dvalue          step-size        r          duality-gap \n');
fprintf(' %2d     %.4e     %.4e       --          %.3e      %.3e',1, pvalue, dvalue, r, pvalue - dvalue)

test =1;
if test == 1
    while k<= maxiter
        lambdak = linsubp(A, y, l, u, x_plus_indic, lambdak_old);
        lambdak_old = lambdak;
        x = max( y + A'*lambdak, 0);
        x_plus_indic = (x>0);
        Axk=A*xk;
        h_l = max(0,l - Axk); h_u =min(0,u - Axk);
        %================== update xk to feasible area
        scaler = min( [l(h_l>0)./Axk(h_l>0); u(h_u<0)./Axk(h_u<0); 1] );
        xk = xk*scaler;
        %===================
        
%         pvalue  = 0.5*norm(xk - y,2)^2; % save object value of primal  and dual problems
%         dvalue = D_funvalue(lambdak, A, y, l, u);
%         Pvalue(k) = pvalue; Dvalue(k) = dvalue;
%         pvalue_min = min(pvalue_min, pvalue);
%         dvalue_max = max(dvalue_max, dvalue);
%         r = sqrt( 2*(pvalue - dvalue)/gamma1 );
    end
else
    %% main loop
    screening = 0;
    delta = 0.7; % for linear search scale
    sigma = 0.4; % for linear search adjust factor
    xi = 0.9;    % nonmonotune linear search; Qk = xi*Qk + 1
    dvalue_avg = dvalue; Qk= 1; % nonmonotune linear search
    while k <= maxiter
        gradk_old = -A*xk_old;
        if k >1
            alphak0 = 1;
            if BBstep > 0
                gradk_old = -A*xk_old;
                gradk = -A*xk;
                alphakB = (gradk_old - gradk)' * (lambdak - lambdak_old(~AS_0) ) /norm(lambdak - lambdak_old(~AS_0), 2)^2;
                alphak0 = min(alphakmax, max(alphakmin, alphakB));
            end
            
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
                
                lambdak = zeros(size(lambdak_old));
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
        %     xk = y + A'*lambdak;
        Axk=A*xk;
        h_l = max(0,l - Axk); h_u =min(0,u - Axk);
        %================== update xk to feasible area
        scaler = min( [l(h_l>0)./Axk(h_l>0); u(h_u<0)./Axk(h_u<0); 1] );
        xk = xk*scaler;
        %     xk = max(xk, 0);
        %===================
        sigma_1 = norm(lambdak, inf) + 1;
        sigma_bound_error =sigma_1*norm([h_l; h_u], 1);
        bound_error = sqrt( norm(h_l)^2 + norm(h_u)^2 );
        Sigma_bound_error(k) = sigma_bound_error;
        
        pvalue  = 0.5*norm(xk - y,2)^2; % save object value of primal  and dual problems
        dvalue = D_funvalue(lambdak, A, y, l, u);
        Pvalue(k) = pvalue; Dvalue(k) = dvalue;
        pvalue_min = min(pvalue_min, pvalue);
        dvalue_max = max(dvalue_max, dvalue);
        %     r = sqrt( 2*(pvalue + sigma_bound_error - dvalue)/gamma1 );
        r = sqrt( 2*(pvalue - dvalue)/gamma1 );
        % only when the detection is safe(r is sufficient small -> the ragion is sufficient narrow)
        %     if r < 0.01
        %         r;
        %     end
        if r < 1e-1
            screening = 1;
        end
        if screening == 1 && active_method == 1 && mod(k+1,active_epoch)==0
            A_rownorm = sqrt(sum(A.*A,2));
            T_p = A*xk + A_rownorm * r;
            T_n = A*xk - A_rownorm * r;
            AS_0 = (T_p<u) & (T_n>l);
            AS_p = (T_p>l) & (T_n<l);
            AS_n = (T_p>u) & (T_n<u);
            num_AS_0(k) = sum(AS_0);
            num_AS_p(k) = sum(AS_p);
            num_AS_n(k) = sum(AS_n);
            %lambdak(AS_0) = 0; %screening the components to be zeros
            lambda_unzero_indic(AS_0) = [];
            A(AS_0,:) = [] ;
            l(AS_0) = []; u(AS_0) = [];
            lambdak = lambdak(~AS_0);
            %xk = max(y + A'*lambdak, 0);
        else
            AS_0 = zeros(size(lambdak));
            AS_p = zeros(size(lambdak));
            AS_n = zeros(size(lambdak));
        end
        
        if plot_true == 1
            plot(k, pvalue, 'b-.');hold on;
            plot(k, dvalue, 'r:.');
            if mod(k,5) == 1
                pause(0.5)
            end
        end
        if  k > 1
            relative_error = abs((dvalue - Dvalue(k-1))/Dvalue(k-1));
            gap_error = pvalue - dvalue;
            Axk=A*xk;
            h_l = max(0,l - Axk); h_u =min(0,u - Axk);
            bound_error = sqrt( norm(h_l)^2 + norm(h_u)^2 );
            fprintf('\n k=%2d  p=%2.4e d=%2.4e alphak=%2.2e  bd_e=%0.3e  gapupper=%2.3e  #AS_0=%d  #AS_p=%d  #AS_n=%d',...
                k, pvalue, dvalue, alphak, bound_error, pvalue  - dvalue, sum(AS_0), sum(AS_p), sum(AS_n))
            if abs(gap_error) + bound_error < epsilon
                if pvalue >= dvalue || 1
                    convergence = 1;
                    fprintf('\n === convergence====  r=%.3e; gapceil=%.4e ', r, pvalue - dvalue);
                    %                 fprintf('\n === convergence====  relative_error=%.3e; \n', relative_error);
                    break;
                else
                    convergence = 0;
                    fprintf(' ==> OVERFLOW')
                    %                 continue
                end
            end
        else
            convergence = 0;
        end
        k = k+1;
    end
end
% end of the main loop
if plot_true == 2
    figure;
    marker_squence = max(floor(length(Pvalue)/30),1);
    plot(Pvalue, 'b-+', 'MarkerIndices',1:marker_squence:length(Pvalue) );hold on;
    plot(Dvalue, 'r-d', 'MarkerIndices',1:marker_squence:length(Pvalue) );
    %     plot(Sigma_bound_error, 'g:', 'MarkerIndices',1:marker_squence:length(Pvalue) );
end
if  plot_true > 0 %&& convergence == 1
    %     plot([0,k],[dvalue,dvalue],'m-.')
    line([0,k],[dvalue,dvalue],'Color',[0.1,0.1,0.1],'LineStyle', '-.')
    title_text = ['\epsilon : ',num2str(epsilon)];
    title(title_text)
    legend({'pvalue','dvalue'},'Location','best');
    hold off;
end
%% output
if active_method ==1
    lambda = zeros(m, 1);
    lambda(lambda_unzero_indic) = lambdak;
    x = max(y + A_old' * lambda, 0);
else
    lambda = lambdak;
    x = xk;
end
value.p = pvalue;
value.d = dvalue;
fprintf('\n>>>>>>>>>>>>>>>end of all the iterations<<<<<<<<<<<<<<')


%% Dual function value
function dvalue = D_funvalue(lambdak, A, y, l, u)
% used to compute the dual object
xk = max(y + A'*lambdak, 0);
f   = 0.5*(y'*y - xk'*xk);
psi = l(lambdak>0)'*lambdak(lambdak>0) + u(lambdak<0)'*lambdak(lambdak<0);
dvalue = f + psi;

%% BB step point
% unused yet
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
% unused yet
function lambdak_new = Trunc(lambdak, AS_p, AS_n)
lambdak_new = zeros(size(lambdak));
lambdak_new(AS_p) = max(lambdak(AS_p ), 0);
lambdak_new(AS_n) = min(lambdak(AS_n ), 0);

%% Sub problem--solving the linear systerm
function lambdak = linsubp(A, y, l, u, x_plus_indic, lambdak)
Asub = A(:,x_plus_indic);
AAtsub = Asub*Asub';
ysub = y(x_plus_indic);
z = zeros(size(lambdak));
z(lambdak>0) = l(lamndak>0);
z(lambdak<0) = u(lamndak<0);
lambdak = AAtsub\(z-Asub*ysub);


