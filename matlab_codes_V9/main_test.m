%% main_test
% 本程序用来测试、记录并储存，算例包括Netlib LP数据集，和随机模拟数据（需要设定维数）
%% test function
DATA_netlib_data = 0; % 1 uses netlib lp data, 0 don't uses netlib lp data
DATA_random_data = 1; % 0 don't use random data, 1 uses generated random data, 2 uses my randomly generated data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
if DATA_netlib_data == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    clearvars -except  DATA_netlib_data  DATA_random_data
    
    %-------------------------------------------------
    OPTION.epsilon = 1e-8;
    % termination condition
    OPTION.maxiter = 50000;
    % the max iterations for solveing the first subproblem
    OPTION.maxtime = 300;
    % the max time for solveing the first subproblem
    OPTION.plot = 0;
    % 0 don't plot, 1 plot point during iterations, 2 plot after iteration
    OPTION.BBstep = 3;
    % 0 don't use BB step size when solving the first subproblem
    % 1 use the long BBstep size
    % 2 use the short BBstep size
    % 3 use adaptive BBstep size
    OPTION.active_method = 0;
    % 0 don't screenign, 1 screening
    OPTION.active_epoch = 5;
    %--------------------------------------------------
    
    
    
    usingDASA   = 1;
    usingSolver = 1;
    method = 2; % 1, only cplex, 2,only gurobi, 0, both cplex and gurobi;
    
    
    
    %---------------------------------------------------
    LP = load('LP_solution');
    run_prob = LP.run_prob;
    % run_prob=run_prob(1:10);
    tested_peoblem_list = find(run_prob==1);
    %Time_record = zeros(length(tested_peoblem_list),4);
    for i = 1:length(tested_peoblem_list)
        problem_number = tested_peoblem_list(i);
        [A,l,u,y,p_name] = rand_Netlib_data(problem_number, LP);
        [m, n] = size(A);
        lambda0 = ones(m, 1);
        %% using DASA
        if usingDASA == 1
            OPTION.active_method = 0;
            tic;
            [x30, lambda30, value30] =  DASA_SCREEN_9_3(A, y, l, u, lambda0, OPTION);
            time30 = toc;
            
            OPTION.active_method = 1;
            tic;
            [x31, lambda31, value31] =  DASA_SCREEN_9_3(A, y, l, u, lambda0, OPTION);
            time31 = toc;
            Time3 = [time30, time31];
        else
            Time3=[];
        end
        
        %% using solver
        if usingSolver ==1
            [X,OBJ,Solver_time]=yalmip_test(A,l,u,y,n,method);
        else
            Solver_time = [];
        end
        %%
        fprintf('\n========= %dth/%d problem, Name: %s==========',i, length(tested_peoblem_list), p_name)
        fprintf('\n========= SIZE:  [m=%d,n=%d]',size(A))
        Time_record(i,:) = [Time3, Solver_time];
        fprintf('\n========= Time:  %.3fs', Time_record(i,3:end))
        pause(1.5)
    end
    %sizeA = LP.sizeA(logical(run_prob),:);
    %prob_names = LP.prob_names(logical(run_prob));
    Test_Information = table([size(A),Time_record]) ;
    path = ['Netlib',datestr(date),'.xls'];
    writetable(Test_Information, path, 'Sheet',datestr(now,30));
   % 'VariableNames',string(prob_names),...
end


%%
if DATA_random_data == 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clearvars -except  DATA_netlib_data  DATA_random_data

    
    %-------------------------------------------------
    OPTION.epsilon = 1e-6;
    % termination condition
    OPTION.maxiter = 50000;
    % the max iterations for solveing the first subproblem
    OPTION.plot = 0;
    % 0 don't plot, 1 plot point during iterations, 2 plot after iteration
    OPTION.BBstep = 3;
    % 0 don't use BB step size when solving the first subproblem
    % 1 use the long BBstep size
    % 2 use the short BBstep size
    % 3 use adaptive BBstep size
    OPTION.active_method = 0;
    % 0 don't screenign, 1 screening
    OPTION.active_epoch = 5;
    %--------------------------------------------------
    
    
    
    usingDASA   = 1;
    usingSolver = 1;
    method = 0; % 1, only cplex, 2,only gurobi, 0, both cplex and gurobi;
    
    
    
    %---------------------------------------------------
    
    M=[50]*100; N=[5,10,20]*100;
    NM= [reshape(repmat(N,length(M),1),length(M)*length(N),1),repmat(M',length(N),1)];
    for i = 1:length(NM)
        n = NM(i,1); m = NM(i,2);
        [A,l,u,y, lambda0]=generate_random_data(n,m);
        %% using DASA
        if usingDASA == 1
            OPTION.active_method = 0;
            tic;
            [x30, lambda30, value30] =  DASA_SCREEN_9_3(A, y, l, u, lambda0, OPTION);
            time30 = toc;
            
            OPTION.active_method = 1;
            tic;
            [x31, lambda31, value31] =  DASA_SCREEN_9_3(A, y, l, u, lambda0, OPTION);
            time31 = toc;
            Time3 = [time30, time31];
        else
            Time3=[];
        end
        
        %% using solver
        if usingSolver ==1
            [X,OBJ,Solver_time]=yalmip_test(A,l,u,y,n,method);
        else
            Solver_time = [];
        end
        %%
        fprintf('\n========= %dth/%d problem==========',i, length(M)*length(N))
        fprintf('\n========= SIZE:  [m=%d,n=%d]',size(A))
        Time_record(i,:) = [Time3, Solver_time];
        fprintf('\n========= Time:  ');disp( Time_record(i,:))
        pause(2)
    end
    Test_Information = table([NM,Time_record]) ;
    path = ['Rand',datestr(date),'.xls'];
    writetable(Test_Information, path, 'Sheet',datestr(now,30));
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% End Main Function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%













%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------------------------------------------
%% sub functions--------bellow---------------------
%% ------------------------------------------------

% 44, 55, 130
% rand_Netlib_data(130, run_prob, LP) %本行可删
%% generate random data: A,l, u, y(infeasiblepoints), z(feasiblepoints)
%% generate Netlib data
function [A,l,u,y, p_name] = rand_Netlib_data(i, LP)

% run_prob = ones(size(LP.run_prob));
prob_names = LP.prob_names;
prob_solutions = LP.prob_solutions;
p_name = prob_names{i};
Problem = load(['LP_MATLAB\', p_name]);
rng(2019);
if isfield(Problem,'A')==0, Problem=Problem.Problem;end
A = Problem.A;
b = Problem.b;
[m, n] = size(A);
y = rand(n, 1); % rand vector to be approximated(infeasible points in th LP)
z = prob_solutions{i}; % feasible points of polyhedron constraints
Ay = A*y;
blu = [-0.05*norm(Ay, inf), 0.05*norm(Ay, inf)];
bl = blu(1);
bu = blu(2);

control = 0.1; % the smaller c, the wider the gap and smaller the noise
noise =  control*bl.*rand(m, 1); % noise is random variables from normal distribution
l = ones(m, 1)*bl + min(noise, -bl);
noise =  control*bu.*rand(m, 1); % noise is random variables from normal distribution
u = ones(m, 1)*bu + max(noise, -bu);

% % control = 100;boundgap1 = control*abs(b).*rand(m, 1); boundgap2 = control*abs(b).*rand(m, 1);
% control = 100; boundgap1 = control*mean(abs(b))*rand(m, 1); boundgap2 = control*mean(abs(b))*rand(m, 1);
% l = b - boundgap1;
% u = b + boundgap2;
% % M = 1e+20; l(b==0) = -M; u(b==0) = M;
end

%% generate random data
function [A,l,u,y, lambda0]=generate_random_data(n, m)
s = rng(2019,'twister');
Alu = [-10, 10]; % 2-dim variables, Alu(1)<Alu(2)  correaponding to the lower and upper bound of A
A = randi(Alu, [m, n]);
y = randn(n,1);Ay = A*y;
blu = [-0.05*norm(Ay, inf), 0.05*norm(Ay, inf)];
bl = blu(1);
bu = blu(2);

control = 0.1; % the smaller c, the wider the gap and smaller the noise
noise =  control*bl.*rand(m, 1); % noise is random variables from normal distribution
l = ones(m, 1)*bl + min(noise, -bl);
noise =  control*bu.*rand(m, 1); % noise is random variables from normal distribution
u = ones(m, 1)*bu + max(noise, -bu);

% lambda0 = ones(m, 1);
lambda0 = zeros(m, 1);
end




%% matlab_buildin_dolver, quadprog, fmincon
% [x,fval,exitflag,output]=matlab_buildin_dolver(A,l,u,y,n,method_buildin);

%% %%%%%%%%%%%%%%sub function end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%% End all %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
