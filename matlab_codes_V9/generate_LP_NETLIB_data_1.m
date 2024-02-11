% use netlib lp data set, down load from
% <http://users.clas.ufl.edu/hager/coap/format.html>
% use the polyhedron constraints matrix A
% then we generate the infeasibel points and feasible points ourself
% then generate the lower and upper bound of A*x
function [y, A, l, u, m, n] = generate_LP_NETLIB_data_1(prob_name, data_set_name, trans)
rng(2019);

if strcmpi(data_set_name, 'netlib')
    Problem_struct = load(['LP_MATLAB\', prob_name]);
end
if isfield(Problem_struct,'A')==0
    Problem_struct = Problem_struct.Problem;
end

A = Problem_struct.A;
% if size(A,2) > 0.5 * size(A,1)
%     trans = 't';
% end
if strcmpi(trans, 't')
    A = A';
end


[m, n] = size(A);
y = 2*rand(n, 1)-1; % rand vector to be approximated(infeasible points in th LP)
Ay = A*y;
blu = [-0.1*norm(Ay, inf), 0.1*norm(Ay, inf)];
bl = blu(1);
bu = blu(2);


tight = 'loose';

if strcmpi(tight, 'tight')
    control = 0.1; % the smaller c, the wider the gap and smaller the noise
    gauss_noise =  control*bl.*randn(m, 1); % noise is random variables from normal distribution
    l = ones(m, 1)*bl + min(gauss_noise, -bl);
    gauss_noise =  control*bu.*randn(m, 1); % noise is random variables from normal distribution
    u = ones(m, 1)*bu + max(gauss_noise, -bu);

elseif strcmpi(tight, 'loose')
    control = 1; % the smaller c, the wider the gap and smaller the noise
    noise =  control*bl.*rand(m, 1); % noise is random variables from normal distribution
    l = ones(m, 1)*bl + noise;
    noise =  control*bu.*rand(m, 1); % noise is random variables from normal distribution
    u = ones(m, 1)*bu + noise;
else
    l = bl*rand(m, 1);
    u = bu*rand(m, 1);
end

% control = 0.1; % the smaller c, the wider the gap and smaller the noise
% noise =  control*bl.*rand(m, 1); % noise is random variables from normal distribution
% l = ones(m, 1)*bl + min(noise, -bl);
% noise =  control*bu.*rand(m, 1); % noise is random variables from normal distribution
% u = ones(m, 1)*bu + max(noise, -bu);

end