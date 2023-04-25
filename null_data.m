%Input
T = readtable('fake_data.csv');
TR = 0.1;
M = table2array(T);
M  = M.';       %transpose (!!!)
Z = zscore(M);

%Training
addpath('.\MINDy-master\MINDy-master\MINDy_Base_v1.0');
out = MINDy_Simple(M, TR, "y"); 

%Get parameters
parameters = out.Param;
W = parameters{5};
D = parameters{6};
alpha = parameters{2};
b = parameters{3}(1);

%Generate null-data
Xt = propagate(W, D, alpha, b, TR, 5);
writematrix(Xt, "MINDy_null_data.csv")

%Define propagation
function Xt = propagate(W, D, alpha, b, TR, max_idx)
    n = size(W, 1);
    x1 = rand(n, 1);

    Xt = zeros(n, max_idx);
    Xt(:, 1) = x1;

    for i = 2:max_idx
        psi = sqrt(alpha.^2+(x1.*b+5).^2)-sqrt(alpha.^2+(x1.*b-5).^2);
    
        x2 = x1 + (W*psi - D.*x1).*TR;

        Xt(:, i) = x2;
        x1 = x2;
    end
end
