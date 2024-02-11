function [y, A, l, u, m, n, lambda0]  = geneate_CS_data()
image_name = 'Cameraman'; % Cameraman, Lena
image_mane_path = sprintf('D:/E/work/Projects/3.6 Dual Gap Active Set Strategy/numerical implement/DATA_realData/Images/%s.tif', image_name);
C = imread(image_mane_path, 'PixelRegion',{[50,200],[50,200]});
D = cast(C, 'single');
imagesc(D); %colormap(map);
num = numel(D);
x = reshape(cast(D, 'single'), num, 1);
N = num;
M = round(0.01*N);
Phi = randn(M, N);
Phi = Phi/sqrt(M); % Phi is the transformation matrix (Observation Matrix)
y = Phi*x;
y = y+0.01*randn(length(y),1);
psi = zeros(N);  % psi is the sparse representation matrix
for k=1:N
    t = zeros(N,1);
    t(k) = 1; 
    t = idct(t);
    psi(:,k) = t;
end
A = (Phi*psi)';
sigma = 1000;
l =-sigma*ones(N,1);
u = sigma*ones(N,1);
m = N;
n = M;
lambda0 = zeros(m,1);

end