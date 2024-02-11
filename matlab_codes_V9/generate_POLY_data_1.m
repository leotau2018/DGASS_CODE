% use MIPLIB lp data set, down load from
%
% use the polyhedron constraints matrix A
% then we generate the infeasibel points and feasible points ourself
% then generate the lower and upper bound of A*x
function [y, A, l, u, m, n] = generate_POLY_data_1(prob_name, data_set_name, trans)
rng(2019);
mat_storePath= 'D:\E\work\Projects\3.6 Dual Gap Active Set Strategy\numerical implement\DATA_realData\poly\';
% filenames = dir([mat_storePath, '*.mat']);
% prob_name = filenames(1).name;
if strcmpi(data_set_name, 'poly')
    Problem_struct = load([mat_storePath, prob_name, '.mat']);
end
use_Aineq = 1;
if use_Aineq == 1
    A = Problem_struct.Aineq;
    if numel(A) == 0
        A = Problem_struct.Aeq;
    end
end
if use_Aineq == 0
    A = Problem_struct.Aeq;
    if numel(A) == 0
        A = Problem_struct.Aineq;
    end
end
if size(A,2) > 0.5 * size(A,1)
    trans = 't';
end
if strcmpi(trans, 't')
    A = A';
end

[m, n] = size(A);
y = 2*rand(n, 1)-1; % rand vector to be approximated(infeasible points in th LP)
Ay = A*y;
blu = [-0.05*norm(Ay, inf), 0.05*norm(Ay, inf)];
bl = blu(1);
bu = blu(2);


tight = 'loose';

if strcmpi(tight, 'tight')
    control = 0.1; % the smaller c, the wider the gap and smaller the noise
    noise =  control*bl.*randn(m, 1); % noise is random variables from normal distribution
    l = ones(m, 1)*bl + min(noise, -bl);
    noise =  control*bu.*randn(m, 1); % noise is random variables from normal distribution
    u = ones(m, 1)*bu + max(noise, -bu);

elseif strcmpi(tight, 'loose')
    control = 1; % the smaller c, the wider the gap and smaller the noise
    noise =  control*bl.*rand(m, 1); % noise is random variables from normal distribution
    l = ones(m, 1)*bl - noise;
    noise =  control*bu.*rand(m, 1); % noise is random variables from normal distribution
    u = ones(m, 1)*bu + noise;
else
    control = 0.1; % the smaller c, the wider the gap and smaller the noise
    noise =  control*bl.*rand(m, 1); % noise is random variables from normal distribution
    l = ones(m, 1)*bl - noise;
    noise =  control*bu.*rand(m, 1); % noise is random variables from normal distribution
    u = ones(m, 1)*bu + noise;
end
% % control = 100;boundgap1 = control*abs(b).*rand(m, 1); boundgap2 = control*abs(b).*rand(m, 1);
% control = 100; boundgap1 = control*mean(abs(b))*rand(m, 1); boundgap2 = control*mean(abs(b))*rand(m, 1);
% l = b - boundgap1;
% u = b + boundgap2;
% % M = 1e+20; l(b==0) = -M; u(b==0) = M;
end