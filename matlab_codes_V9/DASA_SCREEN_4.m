function [x, lambda, value, num_AS] =  DASA_SCREEN_4(A, y, l, u, lambda0, varargin)
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
%
% Copyright: USST, Yunlong Wang
% allow BB step size in SPARSA
% Version 4: add gradient information: allow remove zeros in AP_+ and AP_-.

%% get the input from OPTION

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
        BBstep_kappa = 0.01;
    end
end

% whether use active set method
% 0, default, don't use activ method
% 1  only screen the doomed zero variables, AS0 set
% 2  screen the zero variables and given the sign of other components
% if active_method > 0, active_epoch (default 5) and gradient_info (default 1) is better
% given
if isfield(OPTION, 'active_method')
    active_method = OPTION.active_method;
    if isfield(OPTION, 'active_epoch')
        active_epoch = OPTION.active_epoch;
    else
        active_epoch = 5;
    end
    if isfield(OPTION, 'gradient_info')
        gradient_info = OPTION.gradient_info;
    else
        gradient_info = 1;
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
A_rownorm = sqrt(sum(A.*A,2));

lambdak_old = lambda0; k = 1;
xk_old = max(y + A'*lambdak_old, 0);
%=vvv================= update xk to feasible area
Axk_old = A*xk_old;
h_l = max(0,l - Axk_old); h_u =min(0,u - Axk_old);
scaler = min( [l(h_l>0)./Axk_old(h_l>0); u(h_u<0)./Axk_old(h_u<0); 1] );
xk_scaler = xk_old*scaler;
pvalue_scaler = 0.5*norm(xk_scaler - y,2)^2;
dvalue_old = D_funvalue(lambdak_old, A, y, l, u);
%===================
x_plus_indic = [1:n];
lambda_nzero_indic = 1:m;
gamma1 = 1;
pvalue_min = pvalue_scaler; dvalue_max = dvalue_old;
r = sqrt( 2*(pvalue_min - dvalue_max)/gamma1 );
% only when the detection is safe(r is sufficient small -> the ragion is sufficient narrow)
% if r < 0.5*min(u-l)./A_rownorm
if active_method == 1 && r < 0.1*median(u-l)
    T_p = A*xk_old + A_rownorm * r;
    T_n = A*xk_old - A_rownorm * r;
    AS_0 = (T_p<u) & (T_n>l);
    AS_p = (T_p>l) & (T_n<l);
    AS_n = (T_p>u) & (T_n<u);
    
    AS_diff = AS_p - AS_n;
    AS_0 = find(AS_0 == 1);
    AS_p = find(AS_diff==1); % AS_diff(j)==1 means that l(j) is in the sphere while l(j) is not.
    AS_n = find(AS_diff==-1);% AS_diff(j)==-1 means that l(j) is not in the sphere while l(j) is.
% %     lambdak(AS_0) = 0; %screening the components to be zeros
%     lambda_nzero_indic(AS_0) = [];
%     A(AS_0,:) = [] ;
%     l(AS_0) = []; u(AS_0) = [];
%     lambdak_old = lambdak_old(~AS_0);
%     xk_old = max(y + A'*lambdak_old, 0);
else
    AS_0 = [];
    AS_p = []; AS_n = [];
end
if plot_true == 1
    figure
end
fprintf('\n---------------------------\n')
fprintf(' row number of A: %d \n column number of A: %d \n [max,mean,min] of A''s row norm=[%0.2f, %0.2f, %0.2f]',m, n, full(max(A_rownorm)), full(mean(A_rownorm)), full(min(A_rownorm)));
fprintf('\n---------------------------')
fprintf('\n iter    pvalue        dvalue          step-size        r          duality-gap \n');
fprintf(' %2d     %.4e     %.4e       --          %.3e      %.3e',1, pvalue_scaler, dvalue_old, r, pvalue_scaler - dvalue_old)

%% main loop
delta = 0.7; % for linear search scale
sigma = 0.04; % for linear search adjust factor
xi = 0.9;    % nonmonotune linear search; Qk = xi*Qk + 1
dvalue_avg = dvalue_old; Qk= 1; % nonmonotune linear search
t0 = clock;
screening = 0;
while k <= maxiter
    %======lambdak(AS_0) = 0; %screening the components to be zeros
    if screening == 1
%         lambdak(AS_p) = max(lambdak(AS_p), 0);
%         lambdak(AS_n) = min(lambdak(AS_n), 0);
        lambda_nzero_indic(AS_0) = []; 
        A(AS_0,:) = [] ;
        l(AS_0) = []; u(AS_0) = [];
        if k > 1
            lambdak(AS_0) = [];
            xk = max(y + A'*lambdak, 0);
            [AS_p, AS_n] = update_AS_pn(AS_0, AS_p, AS_n);
        end
        lambdak_old(AS_0) = []; % adjust the elements of lambdak_old / remove the zeros entries
    end
    %================
    
    gradk_old = -A*xk_old;
    if k >1
        alphak0 = 1;
        if BBstep > 0 %% BB-size step mode start
            gradk_old = -A*xk_old;
            gradk = -A*xk;
            if BBstep ==1 % use long BB step size
                alphakB = (gradk_old - gradk)' * (lambdak - lambdak_old ) /norm(lambdak - lambdak_old, 2)^2;
            elseif BBstep ==2 % use short BB step size
                alphakB = norm(gradk_old - gradk)^2 / ((gradk_old - gradk)' * (lambdak - lambdak_old ));
            elseif BBstep ==3 % use adaptive BB step size
                alphakB_l = (gradk_old - gradk)' * (lambdak - lambdak_old ) /norm(lambdak - lambdak_old, 2)^2;
                alphakB_s = norm(gradk_old - gradk)^2 / ((gradk_old - gradk)' * (lambdak - lambdak_old ));
                if alphakB_l/alphakB_s < BBstep_kappa
                    alphakB = alphakB_s;
                else
                    alphakB = alphakB_l;
                end
            end
            alphak0 = min(alphakmax, max(alphakmin, alphakB));
        end      %% BB-size step mode end
        
        lambdak_old  = lambdak;
        xk_old = max(y + A'*lambdak_old, 0);
        h_l_old = l - A*xk_old;
        h_u_old = u - A*xk_old;
        
        %%%%%%% SPARSA
        %%%% linear search for sparsa
        alphak =alphak0;
        dvalue_old = D_funvalue(lambdak_old, A, y, l, u);
        dvalue_avg = ( xi*Qk*dvalue_avg + dvalue_old )/(xi*Qk + 1);
        while 1  % nonmonotone linear search for sparsa
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
            if alphak > alphakmax
                break
            end
        end
    else % only the first setp: search direction is the gradient at current point
        dk = gradk_old; % graident direction
        alphak = 1;
        while 1  % linear search for sparsa
            lambdak = lambdak_old + alphak * dk;
            dvalue = D_funvalue(lambdak, A, y, l, u);
            if dvalue >=  dvalue_avg + 0.5 * sigma * alphak * norm(lambdak-lambdak_old)^2 %nonmonotone linear search
                Qk = xi*Qk + 1;
                break
            end
            alphak = delta*alphak;
        end
    end
    
    %=vvv============== second line search for the real step size
    if screening == 1, start_working = 1; end
    if screening == 1 || exist('start_working', 'var') && 0
        s = 1;
        ddk = lambdak  - lambdak_old;
        while 1
            lambdak = Trunc(lambdak_old + s * ddk, AS_p, AS_n) ;
            dvalue = D_funvalue(lambdak, A, y, l, u);
            if dvalue >=  dvalue_avg + 0.5 * sigma * s * norm(lambdak-lambdak_old)^2 %nonmonotone linear search
                Qk = xi*Qk + 1;
                break
            end
            s = delta*s;
            if s<1
                s
            end
            if s < alphakmin
                break
            end
        end
    end
    %============
    
    
    
    
    
    xk = max(y + A'*lambdak, 0);
    Axk=A*xk;
    h_l = max(0,l - Axk); h_u =min(0,u - Axk);
    %=vvv================= update xk to feasible area
    scaler = min( [l(h_l>0)./Axk(h_l>0); u(h_u<0)./Axk(h_u<0); 1] );
    
    xk_scaler = xk*scaler;
    pvalue_scaler = 0.5*norm(xk_scaler - y,2)^2;
    %===================
    
    
    sigma_1 = norm(lambdak, inf) + 1;
    sigma_p_feasi_error =sigma_1*norm([h_l; h_u], 1);
    p_feasi_error = sqrt( norm(h_l)^2 + norm(h_u)^2 );
    Sigma_p_feasi_error(k) = sigma_p_feasi_error;
    
    pvalue  = pvalue_scaler; % save object value of primal  and dual problems
    dvalue = D_funvalue(lambdak, A, y, l, u);
    Pvalue(k) = pvalue; Dvalue(k) = dvalue;
    pvalue_min = min(pvalue_min, pvalue_scaler);
    dvalue_max = max(dvalue_max, dvalue);
    r = sqrt( 2*(pvalue_scaler - dvalue_max)/gamma1 );
    % only when the detection is safe(r is sufficient small -> the ragion is sufficient narrow)
    if  active_method == 1
        if r < 0.01*median(u-l) && mod(k+1,active_epoch)==0
            screening = 1;
            fprintf('=> #:%d,screening:r=%.2e',length(lambdak),r)
            A_rownorm = sqrt(sum(A.*A,2));
            T_p = A*xk_scaler + A_rownorm * r;
            T_n = A*xk_scaler - A_rownorm * r;
            AS_0 = (T_p<u) & (T_n>l);
            AS_p = (T_p>l) & (T_n<l);
            AS_n = (T_p>u) & (T_n<u);
            AS_diff = AS_p - AS_n;
            AS_0 = find(AS_0 == 1);
            AS_p = find(AS_diff==1); % AS_diff(j)==1 means that l(j) is in the sphere while u(j) is not.
            AS_n = find(AS_diff==-1);% AS_diff(j)==-1 means that l(j) is not in the sphere while u(j) is.
            if  gradient_info == 1
                zer_p = intersect(find(lambdak ==0 & h_l == 0),AS_p); % find the the zeros set belonging AS_p
                zer_n = intersect(find(lambdak ==0 & h_u == 0),AS_n); % find the the zeros set belonging AS_n
                zer_pn = union(zer_p, zer_n); % the zero set in AS_p and AS_n, of which the variables are surely zeros
                if ~isempty(zer_p) || ~isempty(zer_n) || ~isempty(zer_pn)
                    fprintf('\n length(zer_p):  %d    length(zer_n):  %d    length(zer_pn):  %d', [length(zer_p), length(zer_n),length(zer_pn)])
                end
                AS_0 = union(AS_0, zer_pn);  % update the zero set (add zer_pn to AS_0)
                AS_p = setdiff(AS_p, zer_p); % update the AS_p (remove zer_p from AS_p)
                AS_n = setdiff(AS_n, zer_n); % update the AS_n (remove zer_n from AS_n)
            end
            num_AS_0(k) = length(AS_0);
            num_AS_p(k) = length(AS_p);
            num_AS_n(k) = length(AS_n);
        else
            screening = 0;
            AS_0 = [];
            AS_p = [];
            AS_n = [];
        end
    else
        num_AS_0 = 0;num_AS_p = 0;num_AS_n = 0;
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
        gap_error = pvalue_scaler - dvalue;
        Axk=A*xk_scaler;
        h_l = max(0,l - Axk); h_u =min(0,u - Axk);
        p_feasi_error = sqrt( norm(h_l)^2 + norm(h_u)^2 );
        time1 = etime(clock,t0); [hh, mm, ss] = mytime(time1);
        fprintf('\n k=%2d p=%2.4e d=%2.4e alphak=%2.2e  bd_e=%0.3e  gap_error=%2.3e  #AS_0=%d  #AS_p=%d  #AS_n=%d [%d:%d:%d]',...
            k, pvalue, dvalue, alphak, p_feasi_error, gap_error, length(AS_0), length(AS_p), length(AS_n), hh, mm, ss)
        
%         diffgap = pvalue + sigma_p_feasi_error - pvalue_scaler
%         if abs(gap_error) + sigma_p_feasi_error < epsilon
         if gap_error < epsilon
            if pvalue_scaler >= dvalue 
                convergence = 1;
                fprintf('\n === convergence====  r=%.3e; gap=%.4e; cpu time=%.3f', r,  gap_error, etime(clock,t0));
%                 fprintf('\n === convergence====  relative_error=%.3e; \n', relative_error);
                break;
            else
                convergence = 0;
                fprintf(' ==> OVERFLOW')
%                 continue
            end
         else
             convergence = 0;
         end
    end
    k = k+1;
end
% end of the main loop
if plot_true == 2
    figure;
    marker_squence = max(floor(length(Pvalue)/30),1);
    plot(Pvalue, 'b-+', 'MarkerIndices',1:marker_squence:length(Pvalue) );hold on;
    plot(Dvalue, 'r-d', 'MarkerIndices',1:marker_squence:length(Pvalue) );
%     plot(Sigma_p_feasi_error, 'g:', 'MarkerIndices',1:marker_squence:length(Pvalue) );
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
if active_method ==1
    lambda = zeros(m, 1);
    lambda(lambda_nzero_indic) = lambdak;
    x = xk_scaler;
else
    lambda = lambdak;
    x = xk_scaler;
end
value.p = pvalue;
value.d = dvalue;
num_AS.A0 = num_AS_0;
num_AS.Ap = num_AS_p;
num_AS.An = num_AS_n;
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
function lambdak_Trunc = Trunc(lambdak, AS_p, AS_n)
% lambdak_Trunc = zeros(length(lambda_nzero_indic), 1);
lambdak_Trunc = lambdak;
lambdak_Trunc(AS_p) = max(lambdak_Trunc(AS_p ), 0);
lambdak_Trunc(AS_n) = min(lambdak_Trunc(AS_n ), 0);
% lambdak_Trunc = lambdak_Trunc(lambda_nzero_indic);

%% time translate function
%%% To change the format of time
function [h,m,s] = mytime(t)
t = round(t);
h = floor(t/3600);
m = floor(rem(t,3600)/60);
s = rem(rem(t,60),60);

%% update sets AS_p and AS_n
function [AS_p_new, AS_n_new] = update_AS_pn(AS_0, AS_p, AS_n)
adj_value = zeros(length(AS_p), 1);
for i = 1:length(AS_p)
    adj_value(i) = sum(AS_0<AS_p(i));
end
AS_p_new = AS_p - adj_value;
adj_value = zeros(length(AS_n), 1);
for i = 1:length(AS_n)
    adj_value(i) = sum(AS_0<AS_n(i));
end
AS_n_new = AS_n - adj_value;
%%
%%% End of time.m

