%This file is meant to be used with functions imported by MINDy repository.
%Just copy-paste the MINDy_master repository in this folder before running. 

%%% INPUT %%%

%Data path
dataset = {};
data_path = ".\data";
file_list = dir(data_path);

%Read and load data
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


%%% TRAINING %%%

%Define cells
Ws = {};
alphas = {};
Ds = {};

%MINDy path and hyperparameters
addpath('.\MINDy-master\MINDy-master\MINDy_Base_v1.0');
TR = 0.7; %TR (s)

%Train MINDy model to data
for i = 1:numel(dataset)
    %Training
    data = dataset{i};
    [out, Wfin, Dfin] = MINDy_Simple(data, TR, "y");
    
    %Get parameters
    parameters = out.Param;
    W = parameters{5};
    D = parameters{6};
    alpha = parameters{2};
    b = parameters{3}(1);
    
    %Store parameters
    Ws{i} = W;
    alphas{i} = alpha;
    Ds{i} = D;
end

%Save parameters
save(".\MINDy_parameters\W.mat", "Ws")
save(".\MINDy_parameters\alpha.mat", "alphas")
save(".\MINDy_parameters\D.mat", "Ds")


%%% SIMULATE DINAMICS %%%

%Define propagation
%This function is not used here but can be useful during debugging
function Xt = propagate(W, D, alpha, b, TR, max_idx, sigma)
    %Generate random first point
    n = size(W, 1);
    x1 = rand(1, 119) * 1.5 - 1;
    x1 = x1';
    
    %Initialize dynamics matrix
    Xt = zeros(n, max_idx);
    Xt(:, 1) = x1;

    %Fill dynamics matrix
    for i = 2:max_idx
        
        %Get psi(x_t) and x_{t+1}
        eps = normrnd(0, sigma, n, 1);
        psi = sqrt(alpha.^2+(x1.*b+0.5).^2) - sqrt(alpha.^2+(x1.*b-0.5).^2);
        x2 = x1 + (W*psi - D.*x1).*TR + eps;
        
        %Save x_{t+1}
        Xt(:, i) = x2;
        x1 = x2;
    end
end