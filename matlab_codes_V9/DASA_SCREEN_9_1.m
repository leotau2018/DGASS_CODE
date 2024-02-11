function [x, lambda, value, num_AS, k, time1, AS_0] =  DASA_SCREEN_8(A, y, l, u, lambda0, varargin)
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
%       OPTION.lbfgs_k:
%           whether use lbfgs method, zero or positive integer
%           0  don't use lbfgs method
%           k  positive integer  use lbfgs and set constant k as this positive integer(default=5);
%       OPTION.Newton:
%           whether use newton method
%           0 not use Newton method
%           1 use (modified) method. if ture, set OPTION.lbfgs_k=0



% Copyright: USST, Yunlong Wang
% allow BB step size in SPARSA
% Version 4: add gradient information: allow remove zeros in AP_+ and AP_-.
% Version 7: base on Version 4
%            change active sets to logical from index
% Version 8: add Newton method for solving(iterative solve the proximal method)
%            add LBFGS method(not complete)
%            add modified newton method(not complete)
%            ERROR: gradient_info fund to be not safe, will be fixed in Version 9

%% get the input from OPTION

alphakmin = 1e-8; alphakmax = 1e+8;
if isempty(varargin)
    OPTION = [];
else 
    OPTION = varargin{:};
end
% the stop condition to terminate the iteration, default is 1e-4
if isfield(OPTION, 'epsilon')
    epsilon = OPTION.epsilon;
else 
    epsilon = 1e-8;
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

% whether use lbfgs method
% 0  don't use lbfgs method
% positive integer  use lbfgs and set constant k as this positive integer(default=5);
if isfield(OPTION, 'lbfgs_k')
    lbfgs_k = OPTION.lbfgs_k;
else
    lbfgs_k = 5;
end

% whether use newton method
% 0 not use Newton method
% 1 use (modified) method. if ture, set OPTION.lbfgs_k=0
if isfield(OPTION, 'Newton')
    Newton = OPTION.Newton;
    lbfgs_k = 0;
else
    Newton = 0;
end
% Newton = 1;
sparsa_line_search = 1;
screen_open = 0;
%% preparation for main loop
A_raw = A; y_raw = y; l_raw = l; u_raw = u;
m = size(A, 1);
n = size(A, 2);
y_norm2 = norm(y,2);
A_rownorm_raw = sqrt(sum(A.*A,2));

lambdak_old = lambda0; k = 1;
lambdak_raw = lambdak_old;
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
% lambda_nzero_indic = 1:m;
gamma1 = 1;
pvalue_min = pvalue_scaler; dvalue_max = dvalue_old;
r = sqrt( 2*(pvalue_min - dvalue_max)/gamma1 );
% only when the detection is safe(r is sufficient small -> the ragion is sufficient narrow)
% if r < 0.5*min(u-l)./A_rownorm
if active_method == 1 && r < 0.1*median(u-l)
    T_p = A*xk_old + A_rownorm_raw * r;
    T_n = A*xk_old - A_rownorm_raw * r;
    AS_0 = (T_p<u) & (T_n>l);
    AS_p = (T_p>l) & (T_n<l);
    AS_n = (T_p>u) & (T_n<u);
else
    AS_0 = false(m,1);
    AS_p = false(m,1);
    AS_n = false(m,1);
end
if plot_true == 1
    figure
end
fprintf('\n---------------------------\n')
fprintf(' row number of A: %d \n column number of A: %d \n [max,mean,min] of A''s row norm=[%0.2f, %0.2f, %0.2f]',m, n, full(max(A_rownorm_raw)), full(mean(A_rownorm_raw)), full(min(A_rownorm_raw)));
fprintf('\n---------------------------')
fprintf('\n iter    pvalue        dvalue          step-size        r          duality-gap \n');
fprintf(' %2d     %.4e     %.4e       --          %.3e      %.3e',1, pvalue_scaler, dvalue_old, r, pvalue_scaler - dvalue_old)

%% main loop
delta = 0.5; % for linear search scale
sigma = 1e-2; % for linear search adjust factor
xi = 0.9;    % nonmonotune linear search; Qk = xi*Qk + 1
dvalue_avg = dvalue_old; Qk= 1; % nonmonotune linear search
t0 = clock;
screening = 0;
resi=inf;
convergence = 0;
while k <= maxiter
    %======lambdak_raw(AS_0) = 0; %screening the components to be zeros
    if screening == 1 % implement screening
        lambdak_raw(AS_0) = 0; % adjust the elements of lambda / remove the zeros entries
        A = A_raw(~AS_0,:);
        l = l_raw(~AS_0); u = u_raw(~AS_0);
        lambdak = lambdak_raw(~AS_0);
        lambdak_old  = lambdak_old_raw(~AS_0);
        xk = max(y + A'*lambdak, 0);
    end
    %================
    
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


        lambdak_old_raw = lambdak_raw;
        lambdak_old  = lambdak_raw(~AS_0); % lambdak_old_raw(~AS_0); it is equivalent
        xk_old = max(y + A'*lambdak_old, 0);
        
        J_x_indi = (xk_old >0);
        A_active = A(:,J_x_indi);
        Axk_old  = A_active * xk_old(J_x_indi);

        l_Axk_old = l - Axk_old;            % l_Axk_old is also the (sub)gradient corresponding to AS_p
        u_Axk_old = u - Axk_old;            % u_Axk_old is also the (sub)gradient corresponding to AS_n

        %%%%%%% SPARSA
        %%%% linear search for sparsa
        alphak = alphak0; % initial step size is BB-step size (if BB exist), or 1;
        dvalue_old = D_funvalue(lambdak_old, A, y, l, u);
        dvalue_avg = ( xi*Qk*dvalue_avg + dvalue_old )/(xi*Qk + 1);
        while 1  % nonmonotone linear search for sparsa / proximal gradient algorithm (PGA)

            lambdak_l = lambdak_old + l_Axk_old/alphak;
            lambdak_u = lambdak_old + u_Axk_old/alphak;
            %%%%% implement proximal operation
            lambdak = zeros(size(lambdak_old));
            lambdak(lambdak_l>=0) = lambdak_l(lambdak_l>=0);
            lambdak(lambdak_u<=0) = lambdak_u(lambdak_u<=0);
            resi=norm(lambdak-lambdak_old,inf)*alphak;
%             if resi < 1e-8
%                 convergence = 1;
%                 break
%             end
            sparsa_line_search = 1;
            if (Newton > 0) % && (resi<1e-1) %&& (screen_open == 1)                
                OPTION.newton_maxit=30;
                OPTION.newton_tol = min(1e-4, resi/2);
                [lambdak_n, newton_flag, newton_iter, newton_resi] = Newton_type(lambdak_old, Axk_old, l, u, A_active, alphak, OPTION);
                %                 d_k = lambdak_n-lambdak_old;
                if ( newton_flag == 1 )
                    sparsa_line_search = 0;
                    lambdak = lambdak_n;
                end
            else
                lambdak_l = lambdak_old + l_Axk_old/alphak;
                lambdak_u = lambdak_old + u_Axk_old/alphak;
                %%%%% implement proximal operation
                lambdak = zeros(size(lambdak_old));
                lambdak(lambdak_l>=0) = lambdak_l(lambdak_l>=0);
                lambdak(lambdak_u<=0) = lambdak_u(lambdak_u<=0);
            end

            
            
            if (sparsa_line_search == 1)
                dvalue = D_funvalue(lambdak, A, y, l, u);
                value_increament = 0.5 * sigma * alphak * norm(lambdak-lambdak_old)^2;
                if dvalue >=  dvalue_avg +  value_increament %nonmonotone linear search
                    Qk = xi*Qk + 1;
                    break
                end
                alphak = 1/delta*alphak;
                if alphak > alphakmax
                    break
                end
            elseif (sparsa_line_search == 0)
                ddk = lambdak_n  - lambdak_old;
                s = 1; % alphak;
                accept_step = 0;
                while 1
                    lambdak_n = lambdak_old + s \ ddk;
%                     lambdak_n = Trunc(lambdak_n, AS_p(~AS_0), AS_n(~AS_0)) ;
                    dvalue = D_funvalue(lambdak_n, A, y, l, u);
                    value_increament =  0.5 * sigma * s \ norm(lambdak_n - lambdak_old)^2;
                    if dvalue >=  dvalue_avg + value_increament %nonmonotone linear search
                        Qk = xi*Qk + 1;
                        accept_step = 1;
                        break
                    end
                    s = 1/delta*s;
                    if (s < alphakmin) || (s > alphakmax)
                        accept_step = 0;
                        break
                    end
                end
                if (accept_step == 1)
                    lambdak = lambdak_n;
                    fprintf("line_s, %.3e",s)
                    break
                end
                if (accept_step == 0)
                    lambdak = lambdak_n;
                    fprintf("unaccept step size, %.3e",s)
                    break
                end
            end
        end
    else % only the first setp: search direction is the gradient at current point
        gradk_old = -A*xk_old; % graident of f
        dk = gradk_old; % not exactly graident direction
        alphak = 1;
        while 1  % linear search for sparsa
            lambdak = lambdak_old + alphak \ dk;
            dvalue = D_funvalue(lambdak, A, y, l, u);
            value_increament = 0.5 * sigma * alphak * norm(lambdak-lambdak_old)^2;
            if dvalue >=  dvalue_avg + value_increament %nonmonotone linear search
                Qk = xi*Qk + 1;
                break
            end
            alphak = 1/delta*alphak;
        end
    end
    
    

    %=vvv============== second line search for the real step size, seems NOT Necessary[ignore]
%     if (screening == 1 || any(AS_p | AS_n)) && (sparsa_line_search == 0)
    if (sparsa_line_search == 0) && 0
        s = 1;
        ddk = lambdak  - lambdak_old;
        while 1
            lambdak = lambdak_old + s \ ddk;
            lambdak = Trunc(lambdak, AS_p(~AS_0), AS_n(~AS_0)) ;
            dvalue = D_funvalue(lambdak, A, y, l, u);
            value_increament =  0.5 * sigma * s * norm(lambdak - lambdak_old)^2;
            if dvalue >=  dvalue_avg + value_increament %nonmonotone linear search
                Qk = xi*Qk + 1;
                break
            end
            s = 1/delta*s;
            if s<1
                s;
            end
            if (s < alphakmin) || (s > alphakmax)
                break
            end
        end

        fprintf("line_s, %.3e",s)
    end
    %============
    
    
    lambdak_raw(~AS_0) = lambdak; % update the accepted lambdak


    xk = max(y + A'*lambdak, 0);
    Axk=A*xk;
    l_Axk = l-Axk; u_Axk = u-Axk;
    h_l = max(0,l_Axk); h_u =min(0,u_Axk);

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
        if (r < 0.1*median(u-l) || 0) && (mod(k+1,active_epoch)==0)
            fprintf('=> #:%d,screening:r=%.2e',sum(~AS_0),r)
            A_rownorm = A_rownorm_raw(~AS_0);
            T_p = A*xk_scaler + A_rownorm * r;
            T_n = A*xk_scaler - A_rownorm * r;
            indi_in = (T_p<u-1e-10) & (T_n>l+1e-10);
            indi_l  = (T_p>l) & (T_n<l);
            indi_u  = (T_p>u) & (T_n<u);
            indi_ll = indi_l & ~indi_u;         % indi_l & (~indi_u) means that l(j) is in the sphere while u(j) is not.
            indi_uu = indi_u & ~indi_l;         % indi_u & (~indi_l) means that u(j) is in the sphere while l(j) is not.
            
            if  gradient_info > 0  % gradient info is contained in h_l and h_u
                if gradient_info == 1 % use exact gradient information
                ZS_AS_p = (lambdak ==0 & l_Axk <= 0) & indi_ll; % find the the zeros set belonging AS_p
                ZS_AS_n = (lambdak ==0 & u_Axk >= 0) & indi_uu; % find the the zeros set belonging AS_n
                elseif gradient_info == 2 % use enhance gradient information (not exactly)
                    ZS_AS_p = (lambdak < 1e-6) & (lambdak + (l_Axk)/alphakB < -1e-9) & indi_ll; % find the the zeros set belonging AS_p
                    ZS_AS_n = (lambdak > 1e-6) & (lambdak + (u_Axk)/alphakB > 1e-9) & indi_uu; % find the the zeros set belonging AS_n
                end
                zer_pn = (ZS_AS_p | ZS_AS_n); % the zero set in AS_p and AS_n, of which the variables are surely zeros
                if any(ZS_AS_p) || any(ZS_AS_n) || any(zer_pn)
                    fprintf('\n length(zer_p):  %d    length(zer_n):  %d    length(zer_pn):  %d', [sum(ZS_AS_p), sum(ZS_AS_n), sum(zer_pn)])
                end
                indi_in = (indi_in | zer_pn);    % update the zero set (add zer_pn to AS_0)
                indi_ll  = (indi_ll & ~ZS_AS_p); % update the AS_p (remove zer_p from AS_p)
                indi_uu  = (indi_uu & ~ZS_AS_n); % update the AS_n (remove zer_n from AS_n)
            end
            % update AS_0, AS_p, AS_n
            AS_free = find(AS_0==0);
            if any (indi_in)                    % update AS_0
                screening = 1;
                screen_open = 1;
                AS_0(AS_free(indi_in)) = 1;
            else
                screening = 0;
            end
            if any(indi_ll)                     % update AS_p
                AS_p(AS_free(indi_ll)) = 1;
            end
            if any(indi_uu)                     % update AS_p
                AS_n(AS_free(indi_uu)) = 1;
            end
            num_AS_0(k) = sum(AS_0);
            num_AS_p(k) = sum(AS_p);
            num_AS_n(k) = sum(AS_n);
        else
            screening = 0;
        end
    else
        num_AS_0 = 0;
        num_AS_p = 0;
        num_AS_n = 0;
    end

    if plot_true == 1
        plot(k, pvalue, 'b:.');hold on;
        plot(k, dvalue, 'r:.');
        if mod(k,5) == 1
            pause(0.5)
        end
    end

    if  k > 1
        gap_error = pvalue_scaler - dvalue;
        relative_gap_error = abs((gap_error)/Dvalue(k-1));
        Axk=A*xk_scaler;
        h_l = max(0,l - Axk); h_u =min(0,u - Axk);
        p_feasi_error = sqrt( norm(h_l)^2 + norm(h_u)^2 );
        time1 = etime(clock,t0); [hh, mm, ss] = mytime(time1);
        fprintf('\n k=%2d p=%2.4e d=%2.4e alphak=%2.2e  bd_e=%0.3e  gap_error=%2.3e  #AS_0=%d  #AS_p=%d  #AS_n=%d [%d:%d:%d]',...
            k, pvalue, dvalue, alphak, p_feasi_error, gap_error, sum(AS_0), sum(AS_p), sum(AS_n), hh, mm, ss)
        
        lambdak_l = lambdak_old + l_Axk_old;
        lambdak_u = lambdak_old + u_Axk_old;
        %%%%% implement proximal operation
        lambdak_temp = zeros(size(lambdak_old));
        lambdak_temp(lambdak_l>=0) = lambdak_l(lambdak_l>=0);
        lambdak_temp(lambdak_u<=0) = lambdak_u(lambdak_u<=0);
        resi=norm(lambdak_temp-lambdak_old,inf);
         if (resi < OPTION.epsilon) % || (gap_error < epsilon)
            if pvalue_scaler >= dvalue 
                convergence = 1;
                break;
            else
                convergence = 0;
                fprintf(' ==> OVERFLOW')
            end
         else
             convergence = 0;
         end
    end
    k = k+1;
end

if convergence == 1
    fprintf('\n === convergence====  r=%.3e; gap=%.4e; cpu time=%.3f;', r,  gap_error, etime(clock,t0));
    fprintf('\n === convergence====  relative_gap_error=%.3e; ', relative_gap_error);
    fprintf('\n === convergence====  residual=%.3e; \n', resi);
end
% end of the main loop
% -------------------------------------------------------------------------
%% figure plot
if plot_true == 2
    figure;
    marker_squence = max(floor(length(Pvalue)/30),1);
    plot(Pvalue, 'b-+', 'MarkerIndices',1:marker_squence:length(Pvalue) );hold on;
    plot(Dvalue, 'r-d', 'MarkerIndices',1:marker_squence:length(Pvalue) );
%     plot(Sigma_p_feasi_error, 'g:', 'MarkerIndices',1:marker_squence:length(Pvalue) );
end
if  (plot_true > 0) && (convergence == 1)
    %     plot([0,k],[dvalue,dvalue],'m-.')
    line([0,k],[dvalue,dvalue],'Color',[0.1,0.1,0.1],'LineStyle', '-.')
    title_text = ['\epsilon : ',num2str(epsilon)];
    title(title_text)
    legend({'pvalue','dvalue'},'Location','best');
    hold off;
end
% plot obj error
if plot_true >=2
    figure
    semilogy(Pvalue-Pvalue(end), 'b-'); hold on;
    semilogy(abs(Dvalue-Dvalue(end)), 'r-.');
    title_text = ['Dual Gap Error and Primal Gap Error with \epsilon:', num2str(epsilon)];
    title(title_text)    
    legend({'p-gap-error','d-gap-error'},'Location','best');
end

%% output
lambda = zeros(m, 1);
lambda(~AS_0) = lambdak_raw(~AS_0);
x = xk_scaler;
value.p = pvalue;
value.d = dvalue;
value.gap = pvalue - dvalue;
value.resi = resi;
value.x = x;
value.lambda = lambda;
num_AS.A0 = num_AS_0;
num_AS.Ap = num_AS_p;
num_AS.An = num_AS_n;
fprintf('\n>>>>>>>>>>>>>>>end of all the iterations<<<<<<<<<<<<<<\n')
end

%% ================================= sub function ==================================
%% Dual function value
function dvalue = D_funvalue(lambdak, A, y, l, u)
% used to compute the dual object
xk = max(y + A'*lambdak, 0);
f   = 0.5*(y'*y - xk'*xk);
psi = l(lambdak>0)'*lambdak(lambdak>0) + u(lambdak<0)'*lambdak(lambdak<0);
dvalue = f + psi;
end

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
end

%% Truncation operator Trunc
function lambdak_Trunc = Trunc(lambdak, indi_ll, indi_uu)
% lambdak_Trunc = zeros(length(lambda_nzero_indic), 1);
lambdak_Trunc = lambdak;
lambdak_Trunc(indi_ll) = max(lambdak_Trunc(indi_ll ), 0);
lambdak_Trunc(indi_uu) = min(lambdak_Trunc(indi_uu ), 0);
% lambdak_Trunc = lambdak_Trunc(lambda_nzero_indic);
end

%% time translate function
%%% To change the format of time
function [h,m,s] = mytime(t)
t = round(t);
h = floor(t/3600);
m = floor(rem(t,3600)/60);
s = rem(rem(t,60),60);
end

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
end

