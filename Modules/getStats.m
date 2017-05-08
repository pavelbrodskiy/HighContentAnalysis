% This function loads an array of structures corresponding to a cell array
% of N labels. Each structure contains S fields (one per stat), each of
% which is a 2D matrix of doubles representing spatial bins.

function statsArray = getStats(labels, settings, metadata)
%% Process data to make spatial maps if they are not yet on hard drive
thruGetStats(labels, settings, metadata);

%% Load spatial maps from hard drive
for i = 1:length(labels)
    disp(['Getting Statistics: ' labels{i}])
    load([settings.thruStats labels{i} '.mat']);
    statsArray(i) = stats;
end
end

