% explore with John's Clin & Eld data
%
% First we want to define the signal on the network using centrality
% measures where we try using degree centrality, closeness centrality,
% eigenvector centrality.

%% PART 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%  P  A  R  T  1 :   S  I  G  N  A  L  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Signal, degree centrality, using unnormalized adjacency matrix

data = zeros(74, 116);

for i = 1:40
    load(['../RawData/Eld', num2str(i), '.mat']);
    data(i, :) = sum(matfile);
end

for i = 1:34
    load(['../RawData/Clin', num2str(i), '.mat']);
    data(i + 40, :) = sum(matfile);
end

save('../ProcessedData/Signal/degree_centrality', 'data');

%% Signal, closeness centrality, using unnormalized adjacency matrix

data = zeros(74, 116);

for i = 1:40
    load(['../RawData/Eld', num2str(i), '.mat']);
    data(i, :) = closeness_centrality(matfile);
end

for i = 1:34
    load(['../RawData/Clin', num2str(i), '.mat']);
    data(i + 40, :) = closeness_centrality(matfile);
end

save('../ProcessedData/Signal/closeness_centrality', 'data');
    
%% Signal, eigenvector centrality, using unnormalized adjacency matrix

data = zeros(74, 116);

for i = 1:40
    load(['../RawData/Eld', num2str(i), '.mat']);
    data(i, :) = eigen_centrality(matfile);
end

for i = 1:34
    load(['../RawData/Clin', num2str(i), '.mat']);
    data(i + 40, :) = eigen_centrality(matfile);
end

save('../ProcessedData/Signal/eigenvector_centrality', 'data');
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  P A R T 1.A : N O R M A L I Z E D   S I G N A L  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Signal, degree centrality, using unnormalized adjacency matrix

data = zeros(74, 116);

for i = 1:40
    load(['../RawData/Eld', num2str(i), '.mat']);
    N = normalized_adjacencyMatrix( matfile );
    data(i, :) = sum(N)';
end

for i = 1:34
    load(['../RawData/Clin', num2str(i), '.mat']);
    N = normalized_adjacencyMatrix( matfile );
    data(i + 40, :) = sum(N)';
end

save('../ProcessedData/Signal/degree_centrality_nor', 'data');

%% Signal, closeness centrality, using unnormalized adjacency matrix

signal_type = 'degree_centrality';
for diffusion_rate = 1e-6:1e-6:1e-5

[ ward_error ] = ...
    accuracy_diffused_signal( diffusion_rate, signal_type );
    
fprintf('alpha = %.4f, ward = %d, %s\n',...
    diffusion_rate, ward_error, signal_type);
end
%% Signal, eigenvector centrality, using unnormalized adjacency matrix

data = zeros(74, 116);

for i = 1:40
    load(['../RawData/Eld', num2str(i), '.mat']);
    N = normalized_adjacencyMatrix( matfile );
    data(i, :) = eigen_centrality(N);
end

for i = 1:34
    load(['../RawData/Clin', num2str(i), '.mat']);
    N = normalized_adjacencyMatrix( matfile );
    data(i + 40, :) = eigen_centrality(N);
end

save('../ProcessedData/Signal/eigenvector_centrality_nor', 'data');

%% PART 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  P A R T 2 : D I F F U S E D   S I G N A L  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Signal, degree centrality, using unnormalized adjacency matrix

diffusion_rate = 1e-4;

data = zeros(74, 116);
for i = 1:40
    load(['../RawData/Eld', num2str(i), '.mat']);
    signal = sum(matfile);
    data(i, :) = signal * inv(eye(size(matfile)) + diffusion_rate * matfile);
end
for i = 1:34
    load(['../RawData/Clin', num2str(i), '.mat']);
    signal = sum(matfile);
    data(i + 40, :) = signal * inv(eye(size(matfile)) + diffusion_rate * matfile);
end

save(['../ProcessedData/Signal/degree_centrality', ...
    num2str(diffusion_rate), '.mat'], 'data');

%% Signal, closeness centrality, using unnormalized adjacency matrix

diffusion_rate = 1e-4;

data = zeros(74, 116);
for i = 1:40
    load(['../RawData/Eld', num2str(i), '.mat']);
    signal = closeness_centrality(matfile);
    data(i, :) = signal * inv(eye(size(matfile)) + diffusion_rate * matfile);
end
for i = 1:34
    load(['../RawData/Clin', num2str(i), '.mat']); 
    signal = closeness_centrality(matfile);
    data(i + 40, :) = signal * inv(eye(size(matfile)) + diffusion_rate * matfile);
end

save(['../ProcessedData/Signal/closeness_centrality', ...
    num2str(diffusion_rate), '.mat'], 'data');
    
%% Signal, eigenvector centrality, using unnormalized adjacency matrix

diffusion_rate = 1e-4;

data = zeros(74, 116);
for i = 1:40
    load(['../RawData/Eld', num2str(i), '.mat']);
    signal = eigen_centrality(matfile);
    data(i, :) = signal * inv(eye(size(matfile)) + diffusion_rate * matfile);
end
for i = 1:34
    load(['../RawData/Clin', num2str(i), '.mat']);
    signal = eigen_centrality(matfile);
    data(i + 40, :) = signal * inv(eye(size(matfile)) + diffusion_rate * matfile);
end

save(['../ProcessedData/Signal/eigenvector_centrality', ...
    num2str(diffusion_rate), '.mat'], 'data');


%% PART 3 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   P A R T 3  :   L 2   C O M P  %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doi = 'degree_centrality_0.0001';

load(['../ProcessedData/Signal/', doi, '.mat']);

N = size(data, 1);
dist_l2 = zeros(N);
for i = 1:N
    for j = i+1:N
        dist_l2(i,j) = norm(data(i, :) - data(j, :));
        dist_l2(j,i) = dist_l2(i,j);
    end
end

%% PART 3 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%   P A R T 3 A :   L 2   C O M P   F O R  %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Let the play start!
%
signal_type = 'degree_centrality';
for diffusion_rate = 1e-5:5e-6:5e-5

[ ward_error ] = ...
    accuracy_diffused_signal( diffusion_rate, signal_type );
    
fprintf('alpha = %.4f, ward = %d, %s\n',...
    diffusion_rate, ward_error, signal_type);
end

%% PART 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%  P  A  R  T  4 :   D I F F U S I O N   %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Diffused signal

signal_type = 'degree_centrality';
aggregate_network = 'All';
for diffusion_rate = 1e-6:1e-6:1e-5

[ward_error] = error_diffusedSignal_averageNetwork( ...
    diffusion_rate, signal_type, aggregate_network );

fprintf('alpha = %.4f, ward = %d, %s\n',...
    diffusion_rate, ward_error, signal_type);

end

%% PART5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  P  A  R  T  4 :   P  R  I  N  C  O  M  P   %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% principal component analysis

signal_type = 'degree_centrality';
signal = load(['../ProcessedData/Signal/', signal_type]);
data = signal.data;

[coeff, score] = princomp(data);

%% clinical and healthy subjects classification

% Define parameters here.
%
diffusion_rate_array = [0, logspace(-6, 1, 20)];
signal_type = 'DC+FA+GM';

addpath ../Scripts/libsvm-3.19/matlab/;
result = zeros(length(diffusion_rate_array), 1);
predictLabel = zeros(length(diffusion_rate_array), 74);

for diffusion_rate_idx = 1:length(diffusion_rate_array)
    
    diffusion_rate = diffusion_rate_array(diffusion_rate_idx);

    [ svm_accuracy, svm_label ] = accuracy_SVM_classify( diffusion_rate, signal_type);
    result(diffusion_rate_idx) = svm_accuracy;
    predictLabel(diffusion_rate_idx, :) = svm_label';
end

%% subgroup classification

signal_type = 'degree_centrality';
subgroup_type = 'Path';
diffusion_rate_array = [0, logspace(-6, 1, 20)];


addpath ../Scripts/libsvm-3.19/matlab/
result = zeros(length(diffusion_rate_array), 1);
predictLabel = zeros(length(diffusion_rate_array), 34);

for diffusion_rate_idx = 1:length(diffusion_rate_array)
    
    diffusion_rate = diffusion_rate_array(diffusion_rate_idx);

    [ svm_accuracy, svm_label ] = accuracy_SVM_subgroup( diffusion_rate, signal_type, subgroup_type );
    result(diffusion_rate_idx) = svm_accuracy;
    predictLabel(diffusion_rate_idx, :) = svm_label';
end

%% heatmap and pca to see FA and GM

load('../RawData_matlab/FA_all.mat', 'FA_all');
FA = reshape(cell2mat(FA_all), 116, 74)'; % For FA, each row is a subject.
[~, score] = princomp(FA);
figure;
imagesc(score(:, 1:3));

load('../RawData_matlab/GM_all.mat', 'GM_all');
GM = reshape(cell2mat(GM_all), 116, 74)'; % For GM, each row is a subject.
[~, score] = princomp(GM);
figure ;
imagesc(score(:, 1:3)); % we do see difference.

