function [ ward_error] = error_diffusedSignal_averageNetwork( ...
    diffusion_rate, signal_type, aggregate_network )
%Using the signal type (degree, closeness, eigenvector, etc) and the
%diffusion rate to compute the diffused signal signal * (I + \alpha L)
%^{-1}, where L is the average of the aggregated network of
%aggregate_network type (Eld, Chin, All). Use the diffused signals for all
%74 patients to do ward clustering and compute classification accuracy.

signal = load(['../ProcessedData/Signal/', signal_type, '.mat']);

% set the network
%
if exist(['../ProcessedData/Signal/', aggregate_network, '_all.mat'], 'file')
    load(['../ProcessedData/Signal/', aggregate_network, '_all.mat']);
else
    A_aggregate = zeros(116);
    switch aggregate_network
        case 'Eld'
            for i = 1:40
                loaded = load(['../RawData/Eld', num2str(i), '.mat']);
                A_aggregate = A_aggregate + loaded.matfile;
            end
        case 'Clin'
            for i = 1:34
                loaded = load(['../RawData/Clin', num2str(i), '.mat']);
                A_aggregate = A_aggregate + loaded.matfile;
            end
        case 'All'
            for i = 1:40
                loaded = load(['../RawData/Eld', num2str(i), '.mat']);
                A_aggregate = A_aggregate + loaded.matfile;
            end
            for i = 1:34
                loaded = load(['../RawData/Clin', num2str(i), '.mat']);
                A_aggregate = A_aggregate + loaded.matfile;
            end
    end
    save(['../ProcessedData/Signal/', aggregate_network, '_all.mat'], 'A_aggregate');
end
    
% diffused over network
L = diag(sum(A_aggregate)) - A_aggregate;
compute_diff = @(sig) norm(sig * inv(eye(size(L)) + diffusion_rate * L));

% ward_accuracy
%
data = signal.data;
N = size(data, 1);
dist_l2 = zeros(N);
for i = 1:N
    for j = i+1:N
        dist_l2(i, j) = compute_diff(data(i, :) - data(j, :));
        dist_l2(j, i) = dist_l2(i, j);
    end
end
Z = linkage(squareform(dist_l2), 'ward');
T = cluster(Z, 2);
temp1 = sum(T(1:40) ~= 2) + sum(T(41:74) ~= 1);
temp2 = sum(T(1:40) ~= 1) + sum(T(41:74) ~= 2);
ward_error = min(temp1, temp2);

end