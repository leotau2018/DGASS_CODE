function [y, A, l, u, m, n] = generate_random_data_1(m, n, tight)

Alu = [-10, 10]; % 2-dim variables, Alu(1)<Alu(2)  correaponding to the lower and upper bound of A
A = randi(Alu, [m, n]);
y = rand(n,1); y = y - prctile(y,50);
Ay = A*y;
blu = [-0.1*norm(Ay, inf), 0.1*norm(Ay, inf)];
bl = blu(1); % bl < 0;
bu = blu(2); % bu > 0;
%%
if ~exist('tight', 'var')
    tight = 'loose';
end
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
%     l = bl*rand(m, 1);
%     u = bu*rand(m, 1);
    l = min(Ay)-rand(m, 1);
    u = max(Ay)+rand(m, 1);
end

end