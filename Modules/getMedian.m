% This function loads an array of structures corresponding to a cell array
% of N labels. Each structure contains S fields (one per stat), each of
% which is a 2D matrix of doubles representing spatial bins.

function medProjections = getMedian(labels, settings)
%% Load data from hard drive
for i = length(labels):-1:1
    disp(['Getting Median: ' labels{i}])
    if ~exist([settings.thruMedian labels{i} '.mat'],'file')
        mkdir(settings.thruMedian);
        data = getData(labels(i), settings);
        medProj = median(data{1}, 3);
        save([settings.thruMedian labels{i} '.mat'], 'medProj');
    else
        load([settings.thruMedian labels{i} '.mat'], 'medProj');
    end
    medProjections(:,:,i) = medProj;
end
end

