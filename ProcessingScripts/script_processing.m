% processing scripts

%% FA mean data processing - Clin

% Construct region matching algorithm.
%
temp = xlsread('../RawData/OASIS_MALF_labels_CBStracking.xlsx');
regionMatching = temp(:, 1:2);

% Create the list of files to be fetched about.
%
[num, temp] = xlsread('../RawData/behavior.xlsx', 1);
ID = temp(:, 2); % contains one more thing than Clin_ID
Clin_ID = num(:, 1);

for Clin_ID_idx = 1:length(Clin_ID)
% Read fixed width text file.
%
fid = fopen(['../RawData/FractionalAnisotropy/Clin/', ID{Clin_ID_idx + 1}, '_meanFAOASIS30.txt']);
g = textscan(fid, '%s', 'Delimiter', '\n');
fclose all;

S = g{1, 1};
FA_unMatched = zeros(length(S) ,2);
last = 1;
spacing = [8, 12];
for i = 2:length(S);
    Schar = (S{i, 1});
    for j = 1:2
        FA_unMatched(i, j) = str2double(Schar(last:last+spacing(j) - 1));
        last = last + spacing(j);
    end
    last = 1;
end

% Reorder the data into desired order, consistent with 116 brain regions.
%
FA = zeros(116, 1);
for i = 1:116
    idx = regionMatching(i, 2);
    FA(i) = FA_unMatched(FA_unMatched(:, 1) == idx, 2);
end

% Save the file.
%
save(['../RawData_matlab/FA/Clin', num2str(Clin_ID(Clin_ID_idx)), '.mat'], 'FA');

end

%% FA mean data processing - Eld

% Construct region matching algorithm.
%
temp = xlsread('../RawData/OASIS_MALF_labels_CBStracking.xlsx');
regionMatching = temp(:, 1:2);

% Create the list of files to be fetched about.
%
[num, temp] = xlsread('../RawData/behavior.xlsx', 2);
ID = temp(:, 2); % contains one more thing than Eld_ID
Eld_ID = num(:, 1);

for Eld_ID_idx = 1:length(Eld_ID)
% Read fixed width text file.
%
fid = fopen(['../RawData/FractionalAnisotropy/Eld/', ID{Eld_ID_idx + 1}, '_meanFAOASIS30.txt']);
g = textscan(fid, '%s', 'Delimiter', '\n');
fclose all;

S = g{1, 1};
FA_unMatched = zeros(length(S) ,2);
last = 1;
spacing = [8, 12];
for i = 2:length(S);
    Schar = (S{i, 1});
    for j = 1:2
        FA_unMatched(i, j) = str2double(Schar(last:last+spacing(j) - 1));
        last = last + spacing(j);
    end
    last = 1;
end

% Reorder the data into desired order, consistent with 116 brain regions.
%
FA = zeros(116, 1);
for i = 1:116
    idx = regionMatching(i, 2);
    FA(i) = FA_unMatched(FA_unMatched(:, 1) == idx, 2);
end

% Save the file.
%
save(['../RawData_matlab/FA/Eld', num2str(Eld_ID(Eld_ID_idx)), '.mat'], 'FA');

end


%% GM mean data processing - Clin

% Construct region matching algorithm.
%
temp = xlsread('../RawData/OASIS_MALF_labels_CBStracking.xlsx');
regionMatching = temp(:, 1:2);

% Create the list of files to be fetched about.
%
[num, temp] = xlsread('../RawData/behavior.xlsx', 1);
ID = temp(:, 2); % contains one more thing than Clin_ID
Clin_ID = num(:, 1);

for Clin_ID_idx = 1:length(Clin_ID)
% Read fixed width text file.
%
fid = fopen(['../RawData/meanGMThicknessFiles/Clin/', ID{Clin_ID_idx + 1}, '_meanThicknessOASIS30.txt']);
g = textscan(fid, '%s', 'Delimiter', '\n');
fclose all;

S = g{1, 1};
GM_unMatched = zeros(length(S) ,2);
last = 1;
spacing = [8, 12];
for i = 2:length(S);
    Schar = (S{i, 1});
    for j = 1:2
        GM_unMatched(i, j) = str2double(Schar(last:last+spacing(j) - 1));
        last = last + spacing(j);
    end
    last = 1;
end

% Reorder the data into desired order, consistent with 116 brain regions.
%
GM = zeros(116, 1);
for i = 1:116
    idx = regionMatching(i, 2);
    if any(GM_unMatched(:, 1) == idx)
        GM(i) = GM_unMatched(GM_unMatched(:, 1) == idx, 2);
    end
end

% Save the file.
%
save(['../RawData_matlab/GM/Clin', num2str(Clin_ID(Clin_ID_idx)), '.mat'], 'GM');

end

%% GM mean data processing - Eld

% Construct region matching algorithm.
%
temp = xlsread('../RawData/OASIS_MALF_labels_CBStracking.xlsx');
regionMatching = temp(:, 1:2);

% Create the list of files to be fetched about.
%
[num, temp] = xlsread('../RawData/behavior.xlsx', 2);
ID = temp(:, 2); % contains one more thing than Eld_ID
Eld_ID = num(:, 1);

for Eld_ID_idx = 1:length(Eld_ID)
% Read fixed width text file.
%
fid = fopen(['../RawData/meanGMThicknessFiles/Eld/', ID{Eld_ID_idx + 1}, '_meanThicknessOASIS30.txt']);
g = textscan(fid, '%s', 'Delimiter', '\n');
fclose all;

S = g{1, 1};
GM_unMatched = zeros(length(S) ,2);
last = 1;
spacing = [8, 12];
for i = 2:length(S);
    Schar = (S{i, 1});
    for j = 1:2
        GM_unMatched(i, j) = str2double(Schar(last:last+spacing(j) - 1));
        last = last + spacing(j);
    end
    last = 1;
end

% Reorder the data into desired order, consistent with 116 brain regions.
%
GM = zeros(116, 1);
for i = 1:116
    idx = regionMatching(i, 2);
    if any(GM_unMatched(:, 1) == idx)
        GM(i) = GM_unMatched(GM_unMatched(:, 1) == idx, 2);
    end
end

% Save the file.
%
save(['../RawData_matlab/GM/Eld', num2str(Eld_ID(Eld_ID_idx)), '.mat'], 'GM');

end

%% Group signal in a cell structure, first Eld, then Clin

% Group and save FA.
%
FA_all = cell(74, 1);
for i = 1:40
    data = load(['../RawData_matlab/FA/Eld', num2str(i), '.mat'], 'FA');
    FA_all{i} = data.FA;
end
for i = 1:34
    data = load(['../RawData_matlab/FA/Clin', num2str(i), '.mat'], 'FA');
    FA_all{i + 40} = data.FA;
end
save('../RawData_matlab/FA_all.mat', 'FA_all');

% Group and save GM.
%
GM_all = cell(74, 1);
for i = 1:40
    data = load(['../RawData_matlab/GM/Eld', num2str(i), '.mat'], 'GM');
    GM_all{i} = data.GM;
end
for i = 1:34
    data = load(['../RawData_matlab/GM/Clin', num2str(i), '.mat'], 'GM');
    GM_all{i + 40} = data.GM;
end
save('../RawData_matlab/GM_all.mat', 'GM_all');

% Group and save structural signal.
%
StrConn_all = cell(74, 1);
for i = 1:40
    data = load(['../RawData/StructualConnectivity/Eld', num2str(i), '.mat'], 'matfile');
    StrConn_all{i} = data.matfile;
end
for i = 1:34
    data = load(['../RawData/StructualConnectivity/Clin', num2str(i), '.mat'], 'matfile');
    StrConn_all{i + 40} = data.matfile;
end
save('../RawData_matlab/StrConn_all.mat', 'StrConn_all');




