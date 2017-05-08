function [labels, metadata] = getLabels(settings)
% Obtain a list of the experiments that are in the analysis folder.
[numberTable, labelTabel] = xlsread(settings.labelTabel);
labels = labelTabel(2:end,1)';
for i = 2:size(labelTabel,2)
    metadata.(labelTabel{1, i}) = numberTable(1:end, i-1)';
end