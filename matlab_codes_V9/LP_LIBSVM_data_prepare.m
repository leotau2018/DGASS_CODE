% use libsvm lp data set, down load from
% <https://www.csie.ntu.edu.tw/~cjlin/libsvm/>
% use the polyhedron constraints matrix A
% then we generate the infeasibel points and feasible points ourself
% then generate the lower and upper bound of A*x
clear -DATA;
working_path = 'D:\E\work\Projects\3.6 Dual Gap Active Set Strategy\numerical implement\DATA_realData';
% cd(working_path);
data_path= [working_path,'\LIBSVM\'];
data_save_path = [working_path,'\LIBSVM_mat\'];
% LP_matlab gives 146 instances in total
filenames = dir([data_path]);
i = 0;
k = 0;
while i < length(filenames)
    i = i+1;
    if isequal(filenames(i).name, '.') || isequal(filenames(i).name, '..')
        continue
    end
    k = k+1;
    prob_names{k} = filenames(i).name;
    [label_vector, instance_matrix] = libsvmread([data_path, prob_names{k}]);
    prob_name = prob_names{k};
    if contains(prob_name, '.txt')
        prob_name = prob_name(1:end-4);
    end
    save([data_save_path,prob_name, '.mat'], 'instance_matrix')
    fprintf("=====i:%2d, Size:[%-5d, %-5d],  nnz:%-5d,  name:%10s\n", k, size(instance_matrix,1), size(instance_matrix,2), nnz(instance_matrix), prob_name)
end


