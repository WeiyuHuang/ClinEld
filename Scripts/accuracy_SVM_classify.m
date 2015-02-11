function [ bestAccuracy, bestLabel] = accuracy_SVM_classify( diffusion_rate, signal_type )
%Using the signal type (degree, closeness, eigenvector, etc) and
%thediffusion rate to compute the diffused signal signal * (I + \alpha
%L)^{-1}, and use the diffused signals for all 74 patients to do ward
%clustering as well as K-means and compute classification accuracy.
%bestLabel is the assigned label of the parameter yileding best accuracy.

% Load brain network (structrual connectivity) since this will be used for
% any signal_type.
%
temp = load('../RawData_matlab/StrConn_all.mat');
network_all = temp.StrConn_all;

% If signal_type is FA or GM, load them for all subjects now. Then
% normlized them such that the maximum value in all FA or GM data is 1.
%
switch signal_type
    case 'FA'
        temp = load('../RawData_matlab/FA_all.mat');
        signal_all = temp.FA_all;
    case 'GM'
        temp = load('../RawData_matlab/GM_all.mat');
        signal_all = temp.GM_all;
    case {'DC+FA+GM', 'FA+GM'} % DC as shorthand for degree centrality
        temp = load('../RawData_matlab/FA_all.mat');
        FA_all = temp.FA_all;
        FA_all_mat = cell2mat(FA_all);
        FA_all = mat2cell(FA_all_mat / max(FA_all_mat), 116 * ones(74, 1));
        temp = load('../RawData_matlab/GM_all.mat');
        GM_all = temp.GM_all;
        GM_all_mat = cell2mat(GM_all);
        GM_all = mat2cell(GM_all_mat / max(GM_all_mat), 116 * ones(74, 1));
end

% Compute diffused signal DATA such that each row corresponds to a subject.
%
if strcmp(signal_type, 'DC+FA+GM')
    data = zeros(74, 3 * 116); % three features
    % Degree centrality needs to be normalized
    DC = cell(74, 1);
    DC_max = 0; 
    for i = 1:74
        network = network_all{i};
        DC{i} = sum( network );
        DC_max = max(DC_max, max(DC{i}));
    end % for i
    for i = 1:74
        DC{i} = DC{i} / DC_max;
    end
    % Fill in the data with features
    %
    for i = 1:74
        network = network_all{i};
        data(i, 1:116) = DC{i} * inv(eye(size(network)) + diffusion_rate * network); % DC
        data(i, 117:(2*116)) = FA_all{i}' * inv(eye(size(network)) + diffusion_rate * network); % FA
        data(i, (117+116):(3*116)) = GM_all{i}' * inv(eye(size(network)) + diffusion_rate * network); % GM
    end % for i
elseif strcmp(signal_type, 'FA+GM')
    data = zeros(74, 2 * 116); % three features
    for i = 1:74
        network = network_all{i};
        data(i, 1:116) = FA_all{i}' * inv(eye(size(network)) + diffusion_rate * network); % FA
        data(i, 117:(2*116)) = GM_all{i}' * inv(eye(size(network)) + diffusion_rate * network); % GM
    end % for i
else
    data = zeros(74, 116);
    signal = cell(74, 1);
    signal_max = 0;
    for i = 1:74
        network = network_all{i};
        switch signal_type
            case 'degree_centrality'
                signal{i} = sum( network );
            case 'closeness_centrality'
                signal{i} = closeness_centrality( network );
            case 'eigenvector_centrality'
                signal{i} = eigen_centrality( network );
            case 'degree_centrality_nor'
                N = normalized_adjacencyMatrix( network );
                signal{i} = sum(N);
            case 'closeness_centrality_nor'
                N = normalized_adjacencyMatrix( network );
                signal{i} = closeness_centrality(N);
            case 'eigenvector_centrality_nor'
                N = normalized_adjacencyMatrix( network );
                signal{i} = eigen_centrality(N);
            case {'FA', 'GM'}
                signal{i} = signal_all{i}';
        end % switch
        signal_max = max(signal_max, max(signal{i}));
    end % for i
    % Normalize the signal such that the maximum number is 1
    %
    for i = 1:74
        signal{i} = signal{i} / signal_max;
    end
    % Fill the data
    %
    for i = 1:74
        data(i, :) = signal{i} * inv(eye(size(network)) + diffusion_rate * network);
    end % for i
end % if DC+FA+GM

% Scaling such that the each feature is between 0 and 1.
%
for i = 1:size(data, 2)
    temp = data(:, i);
    data(:, i) = (temp - min(temp)) ./ (max(temp) - min(temp));
end

% Define labels Y.
%
label = [ones(40, 1); 2 * ones(34, 1)];

% Calculate svm accuracy.
%
C_power_array = -5:2:15;
GAMMA_power_array = -15:2:3;
cv_accuracy = zeros(length(C_power_array), length(GAMMA_power_array));
allLabel = cell(length(C_power_array), length(GAMMA_power_array));
for C_power = 1:11
    for GAMMA_power = 1:10
        
        temp = zeros(74, 1);
        C = 2.^ C_power_array(C_power);
        GAMMA = 2.^ GAMMA_power_array(GAMMA_power);
        CrtPredictLabel = zeros(74, 1);
        
        for i = 1:74
            if i == 1
                CRTtrain = data(2:74, :);
                CRTtrain_label = label(2:74);
                CRTtest = data(1, :);
                CRTtest_label = label(1);
            elseif i == 74
                CRTtrain = data(1:73, :);
                CRTtrain_label = label(1:73);
                CRTtest = data(74, :);
                CRTtest_label = label(74);
            else
                CRTtrain = data([1:i-1, i+1:74], :);
                CRTtrain_label = label([1:i-1, i+1:74]);
                CRTtest = data(i, :);
                CRTtest_label = label(i);
            end % if
            
            model = svmtrain(CRTtrain_label, CRTtrain, ['-c ',...
                num2str(C), ' -g ', num2str(GAMMA)]);
            
            [CrtPredictLabel(i), temp_accuracy, ~] = ...
                svmpredict(CRTtest_label, CRTtest, model);
            temp(i) = temp_accuracy(1);
            
        end % i
        
        cv_accuracy(C_power, GAMMA_power) = mean(temp);
        allLabel{C_power, GAMMA_power} = CrtPredictLabel;
    end % gamma_power
end % C_power

svm_accuracy = max(max(cv_accuracy));
[i, j] = find(cv_accuracy == svm_accuracy);
bestLabel = allLabel{i, j};
bestAccuracy = svm_accuracy * 74 / 100;

end