function [ ward_error] = ...
    accuracy_diffused_signal( diffusion_rate, signal_type )
%Using the signal type (degree, closeness, eigenvector, etc) and the
%diffusion rate to compute the diffused signal signal * (I + \alpha L)
%^{-1}, and use the diffused signals for all 74 patients to do ward
%clustering as well as K-means and compute classification accuracy.

doi = [signal_type, '_', num2str(diffusion_rate)];

if exist(['../ProcessedData/Signal/', doi, '.mat'], 'file')
    load(['../ProcessedData/Signal/', doi, '.mat']);
else 
    % need to compute the diffused signal
    %
    data = zeros(74, 116);
    for i = 1:40
        loaded = load(['../RawData/Eld', num2str(i), '.mat']);
        data(i, :) = compute_diffused_signal(loaded.matfile,...
            diffusion_rate, signal_type);
    end
    for i = 1:34
        loaded = load(['../RawData/Clin', num2str(i), '.mat']);
        data(i + 40, :) = compute_diffused_signal(loaded.matfile,...
            diffusion_rate, signal_type);
    end
    % Store results for future references
    %
    save(['../ProcessedData/Signal/', doi, '.mat'], 'data');
end

% ward_accuracy
%
N = size(data, 1);
dist_l2 = zeros(N);
for i = 1:N
    for j = i+1:N
        dist_l2(i,j) = norm(data(i, :) - data(j, :));
        dist_l2(j,i) = dist_l2(i,j);
    end
end
Z = linkage(squareform(dist_l2), 'ward');
T = cluster(Z, 2);
temp1 = sum(T(1:40) ~= 2) + sum(T(41:74) ~= 1);
temp2 = sum(T(1:40) ~= 1) + sum(T(41:74) ~= 2);
ward_error = min(temp1, temp2);

end

function [diffused_signal] = compute_diffused_signal(input, ...
    diffusion_rate, signal_type)

switch signal_type
    case 'degree_centrality'
        signal = sum(input);
    case 'closeness_centrality'
        signal = closeness_centrality(input);
    case 'eigenvector_centrality'
        signal = eigen_centrality(input);
    case 'degree_centrality_nor'
        N = normalized_adjacencyMatrix( input );
        signal = sum(N);
    case 'closeness_centrality_nor'
        N = normalized_adjacencyMatrix( input );
        signal = closeness_centrality(N);
    case 'eigenvector_centrality_nor'
        N = normalized_adjacencyMatrix( input );
        signal = eigen_centrality(N);
end

diffused_signal = signal * inv(eye(size(input)) + diffusion_rate * input);
end