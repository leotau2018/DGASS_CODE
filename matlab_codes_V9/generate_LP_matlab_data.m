% use netlib lp data set, down load from
% <http://users.clas.ufl.edu/hager/coap/format.html>
% use the polyhedron constraints matrix A
% then we generate the infeasibel points and feasible points ourself
% then generate the lower and upper bound of A*x
clear -DATA;
% working_path = 'D:\E\work\Projects\3.6 Dual Gap Active Set Strategy\numerical implement\matlab codes';
working_path =  pwd;
cd(working_path);
data_path= [working_path,'\LP_matlab\'];
% LP_matlab gives 146 instances in total
filenames = dir([data_path,'*.mat']);
for i = 1:length(filenames)
    prob_names{i} = filenames(i).name;
end
% generate the feasible point (solution) of the LP
options = optimoptions('linprog', 'Display', 'off');
run_prob = ones(length(prob_names), 1); % 0 1 , whether problem i is going to run
time_solve_LP = zeros(length(prob_names), 1); %  used time for solving problem i
prob_solutions = cell(length(prob_names), 1); % solution (vector) of problem i
fvals = time_solve_LP; exitflag = time_solve_LP; % Linear Program object values
sizeA = zeros(length(prob_names), 2); % each row is the dimision (size) of the matrix A of the problem
for i = 1:length(prob_names)
    load(['LP_MATLAB\', prob_names{i}]);
    sizeA(i,:) = size(A);
%     if sum(sizeA(i,:)) > 3e+4 || (sizeA(i,1)>1e+4 && sizeA(i,2) > 1e+4)
    if sum(sizeA(i,:)) > 2e+40
        run_prob(i) = 0; 
        fprintf('>>>prob number=%d  prob name:%7s   Dim:[%d,%d]***********\n', i, prob_names{i}, sizeA(i,: ));
        continue;
    else
        fprintf('>>>prob number=%d  prob name:%7s   Dim:[%d,%d]', i, prob_names{i}, sizeA(i,: ));
    end
    tic
    [prob_solutions{i}, fvals(i), exitflag(i)] = linprog(c,[],[],A,b,lo,hi, options);
    time1 = toc; 
    time_solve_LP(i) = time1;
    fprintf(' time:%.4f \n', time_solve_LP(i));
end
totaltime = sum(time_solve_LP);
num_successful_solved = sum(exitflag == 1);
fprintf('\n>>>>>>>>>>>>>>Informations:\n Total time:%4f \n solved/Total: %d/%d  \n', totaltime, num_successful_solved, length(prob_names))
% LPsolutions.feasible_points = feasible_points;
% LPsolutions.fvals = fvals;
% save('LPsolutions.mat', '-struct','LPsolutions');
% disp('Contents of LPsolutions.mat:');
% whos('-file','LPsolutions.mat');
% save('prob_solutions.mat', 'prob_solutions');
% save('prob_names.mat', 'prob_names');
% save('sazeA.mat', 'sizeA');
LP.prob_names = prob_names;
LP.run_prob = run_prob;
LP.prob_solutions = prob_solutions; 
LP.sizeA = sizeA;
save('LP_solution.mat','-struct','LP');
