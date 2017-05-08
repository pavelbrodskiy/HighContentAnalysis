% This function makes spatial maps for each raw tiff stack for which there
% is no rotated video.

function thruGetStats(labels, settings, metadata)
mkdir(settings.thruStats);
for i = 1:length(labels)
    if (~exist([settings.thruStats labels{i} '.mat'],'file'))||settings.force
        if ~exist('backgroundMedianProjection', 'var')
            medianProjections = getMedian(labels, settings);
            backgroundMedianProjection = median(medianProjections, 3);
            settings.background = backgroundMedianProjection;
            for j = length(labels):-1:1
                load([settings.thruData labels{i} '.mat'], 'dye');
                dyeArray(:,:,j) = dye;
            end
            settings.background2 = min(dyeArray, [], 3);
        end
        disp(['Extracting statistics: ' labels{i}])
        settings.timestep = metadata.Timestep(i);
        [data, dyeArray] = getData(labels(i), settings);
        stats = extractWellStatistics(data{1}, dyeArray{1}, settings);
        save([settings.thruStats labels{i} '.mat'], 'stats');
    end
end
end