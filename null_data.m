%This file is meant to be used with functions imported by MINDy repository.
%Just copy-paste the MINDy_master repository in this folder before running. 

%Input
dataset = {};
data_path = ".\data";
file_list = dir(data_path);

n_skip = 0;
for i = 1:numel(file_list)
    file_name = file_list(i).name;
    
    %Skip directories (including "." and "..")
    if file_list(i).isdir
        n_skip = n_skip+1;
        continue;
    end

    %Try to load subjects data in dataset
    try
        subject = load(data_path + "\" + file_name);
        tseries = subject.tseries;
        tseries = zscore(tseries')'; 
        dataset{i-n_skip} = tseries;

         if i == 5
            disp(dataset(i-n_skip))
        end

    catch exception
        disp("Problem with " + file_name + " file")
    end
end



%Training with MINDy_Simple
addpath('.\MINDy-master\MINDy-master\MINDy_Base_v1.0');
TR = 0.7; %TR (s)
[out, Wfin, Dfin] = MINDy_Simple(dataset, TR, "y"); 

%Get parameters
parameters = out.Param;
W = parameters{5};
D = parameters{6};
alpha = parameters{2};
b = parameters{3}(1);
save(".\data_simulations\W.mat", "W")
save(".\data_simulations\D.mat", "D")
save(".\data_simulations\alpha.mat", "alpha")

%Generate null-data
max_idx = 1200;
Xt = propagate(W, D, alpha, b, TR, max_idx);
%writematrix(Xt, ".\data_simulation\MINDy_null_data.csv")



%Try MINDy_RAW_CV
addpath('.\MINDy-master\MINDy-master\MINDy_Filtering_and_Prediction_v1.0');

Dat1 = dataset{1};
Dat2 = dataset{2};
[out_2, Pred, Resid, R2] = MINDy_RAW_CV(Dat1, Dat2);
%[out_2, Pred, Resid, R2] = MINDy_RAW_CV(Dat1, []);

%Save MINDy_RAW_CV predictions
writematrix(Pred, ".\data_simulations\MINDy_null_data_2.csv")

%Get and save MINDy_RAW_CV parameters
parameters_2 = out_2.Param;
W_2 = parameters_2{5};
D_2 = parameters_2{6};
alpha_2 = parameters_2{2};
b_2 = parameters_2{3}(1);
save(".\data_simulations\W_2.mat", "W_2")
save(".\data_simulations\D_2.mat", "D_2")
save(".\data_simulations\alpha_2.mat", "alpha_2")

%Try using my propagation with MINDy_RAW_CV parameters
Xt = propagate(W_2, D_2, alpha_2, b_2, TR, max_idx);
%writematrix(Xt, ".\data_simulations\MINDy_null_data_3.csv")



%Define propagation
function Xt = propagate(W, D, alpha, b, TR, max_idx)
    %Generate random first point
    n = size(W, 1);
    x1 = rand(n, 1);

    %Initialize propagation matrix
    Xt = zeros(n, max_idx);
    Xt(:, 1) = x1;

    %Fill propagation matrix
    for i = 2:max_idx
        %Get psi(x_t) and x_{t+1}
        psi = sqrt(alpha.^2+(x1.*b+0.5).^2) - sqrt(alpha.^2+(x1.*b-0.5).^2);
        x2 = x1 + (W*psi - D.*x1).*TR;

        Xt(:, i) = x2;
        x1 = x2;
    end
end
