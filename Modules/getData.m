% This function loads an array of structures corresponding to a cell array
% of N labels. Each structure contains S fields (one per stat), each of
% which is a 2D matrix of doubles representing spatial bins.

function [dataArray, dyeArray] = getData(labels, settings)
filenames = strrep(labels,'..',char(filesep));

%% Load data from hard drive
for i = length(labels):-1:1
    disp(['Getting Data: ' labels{i}])
    if ~exist([settings.thruData labels{i} '.mat'],'file')
        mkdir(settings.thruData);
        data = readTiff([settings.inData filenames{i} '.nd']);
        data = max(data, [], 4);
        dye = data(:,:,1);
        data = data(:,:,2:2:end);
        save([settings.thruData labels{i} '.mat'], 'data', 'dye');
    else
        load([settings.thruData labels{i} '.mat'], 'data', 'dye');
    end
    dataArray{i} = data;
    dyeArray{i} = dye;
end
end

