% use miplib data set, down load from
% <http://miplib.zib.de/download.html>
% use the polyhedron constraints matrix A
% then we generate the infeasibel points and feasible points ourself
% then generate the lower and upper bound of A*x


% MPS_filePath = 'D:\E\work\Projects\3.6 Dual Gap Active Set Strategy\numerical implement\DATA_realData\MIPLIB2017_Benchmark\';
% mat_storePath= 'D:\E\work\Projects\3.6 Dual Gap Active Set Strategy\numerical implement\DATA_realData\MIPLIB2017_Benchmark_mat\';
MPS_filePath = 'D:\E\work\Projects\3.6 Dual Gap Active Set Strategy\numerical implement\DATA_realData\hard-v18\';
mat_storePath= 'D:\E\work\Projects\3.6 Dual Gap Active Set Strategy\numerical implement\DATA_realData\hard-v18_mat\';

data_path= MPS_filePath;
% LP_matlab gives 146 instances in total
filenames = dir([data_path,'*.mps']);
i=1;
%%
while i <= length(filenames)
    tic;
    try
        prob_name = filenames(i).name;
        prob_names{i} = prob_name(1:end-4);
        prob_MPS_filePath = [MPS_filePath,  prob_name];
        prob = mpsread(prob_MPS_filePath);
        prob.name = prob_name;
        sizeAineq(i,:) = [size(prob.Aineq), nnz(prob.Aineq)/numel(prob.Aineq)];
        sizeAeq(i,:) = [size(prob.Aeq), nnz(prob.Aeq)/numel(prob.Aeq)];
        time1 = toc; 
        time_used(i) = time1;
        fprintf('>>>prob number=%-2d   sizeAineq:[%-6d, %-6d] %-5.1e,   sizeAeq:[%-6d, %-6d] %-5.1e   prob name:%-25s \n', i,  sizeAineq(i,: ), sizeAeq(i,:), prob_names{i});
        save([mat_storePath, prob_names{i}, '.mat'], '-struct', 'prob')
    catch
        prob_names{i} = [];
    end
    i=i+1;
end
%%
totaltime = sum(time_used)
MIPLIB.prob_names = prob_names;
MIPLIB.time_used = time_used;
MIPLIB.sizeAineq = sizeAineq;
MIPLIB.sizeAeq = sizeAeq;
MIPLIB.solve_totaltime = totaltime;
save([mat_storePath, 'MIPLIB'], 'MIPLIB')

