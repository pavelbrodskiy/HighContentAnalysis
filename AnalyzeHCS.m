% %% Initialization
% clearvars
% close all
% settings = prepareWorkspace();
% 
% %% Declare constants
% 
% %% Load tables describing high and low frequency data
% tblLowFrequency = readtable(settings.tblLowFreq);
% [tblHighFrequency, tblPlateMap, tblPlateLegend] = getHighFrequencyTable(settings);
% processConfocalData(tblLowFrequency.Label, settings)
% 
% mkdir('MP4');
% 
% %% Make Low Frequency movies
% s = 32;
% metadata = load([settings.thruData 'LowFrequency_f1.mat']);
% 
% for j = 1:s
%     % Obtain condition information
%     str = metadata.metadata.global.get(['Stage' num2str(j)]);
%     coords = str(2:end-1);
%     
%     condition = cellfun(@(x) {num2str(x)},tblPlateMap(sub2ind(size(tblPlateMap),find(contains(settings.letters, coords(1)))+1,str2num(coords(2:end))+1)));
%     conditionIndex = cell2mat(cellfun(@(x,y) find(~cellfun('isempty',strfind(y,x))), ...
%         tblHighFrequency.Conditions, repmat({tblPlateLegend.Var1}, [1, length(tblHighFrequency.Conditions)])', 'UniformOutput',false));
%     
%     idx = find(contains(tblPlateLegend.Var1, condition));
%     
%     filename = ['MP4' filesep 'Cl8_ZB_' num2str(tblPlateLegend.Blebbistatin_uM(idx)) ...
%         '_uM_Bleb_' num2str(tblPlateLegend.DMSO_percent(idx)) '_pctDMSO_'...
%         coords '_LowFreq'];
%     
%     disp(['Processing: ' filename])
%     
%     if exist([filename, '.mp4'],'file')
%         continue
%     end
%     
%     % Load raw data
%     totalZProj = [];
%     for i = 1:length(tblLowFrequency.Label)
%         tmp = load([settings.thruData 'LowFrequency_s' num2str(j) '_f' num2str(i) '.mat'], 'zProj');
%         totalZProj = cat(3, totalZProj, tmp.zProj);
%     end
%     
%     % Convert into RGB shape
%     totalZProj = totalZProj(:,:,:,[2,1]);
%     totalZProj(:,:,:,3) = 0;
%     
%     writeMP4(totalZProj, filename, 15);
% end
% 
% %% Make high frequency movies
% for i = 1:length(tblHighFrequency.Label)
%     filename = ['MP4' filesep 'Cl8_ZB_' num2str(tblHighFrequency.Blebbistatin_uM(i)) ...
%         '_uM_Bleb_' num2str(tblHighFrequency.DMSO_percent(i)) '_pctDMSO_' ...
%         tblHighFrequency.PlateAddress{i} '_day_' num2str(round(tblHighFrequency.Days(i)*100)/100) '_HighFreq'];
%     
%     disp(['Processing: ' filename])
%     
%     if exist([filename, '.mp4'],'file')
%         continue
%     end
%     
%     if ~exist([settings.thruData tblHighFrequency.Label{i} '.mat'], 'file')
%         raw = bfopen([settings.inData 'HighFrequency' filesep tblHighFrequency.Label{i} '.tif']);
%         raw = cat(1,raw{:,1});
%         zStack = cat(3,raw{:,1});
%         mkdir(fileparts([settings.thruData tblHighFrequency.Label{i}]))
%         save([settings.thruData tblHighFrequency.Label{i} '.mat'], 'zStack');
%     else
%         load([settings.thruData tblHighFrequency.Label{i} '.mat'], 'zStack');
%     end
%     
%     writeMP4(zStack, filename, 60);
% end

%% Extract statistics
for i = 1:32
    coordArray{i} = metadata.metadata.global.get(['Stage' num2str(i)]);
    
    coords = coordArray{i}(2:end-1);
    condition = cellfun(@(x) {num2str(x)},tblPlateMap(sub2ind(size(tblPlateMap),find(contains(settings.letters, coords(1)))+1,str2num(coords(2:end))+1)));
    conditionIndex = cell2mat(cellfun(@(x,y) find(~cellfun('isempty',strfind(y,x))), ...
        tblHighFrequency.Conditions, repmat({tblPlateLegend.Var1}, [1, length(tblHighFrequency.Conditions)])', 'UniformOutput',false));
    
    idx(i) = find(contains(tblPlateLegend.Var1, condition));
    
end

coord = '"D3"';

s = find(contains(coordArray, coord));
conditions = tblPlateLegend(idx(s),:);
idxHighFreq = find(contains(tblHighFrequency.PlateAddress, coord(2:end-1)));


return






stats = getStats(labels, settings, metadata);

%% Make plots
imagingTime = metadata.Days * 24;
% zoFort = metadata.Zofort;

yArray = {[stats.meanAmp],[stats.meanfreq],[stats.cellNumber],[stats.meanWHM]};
ylabels = {'Amplitude','Frequency (mHz)','Cell Number','WHM (s)'};

figure(1)
clf

for i = 1:4
    subplot(2,2,i)
    ys = yArray{i};
    ys(isnan(ys)) = 0;
    %     scatter(imagingTime(logical(zoFort)), ys(logical(zoFort)),64,'.');
    plot(imagingTime(logical(zoFort)), ys(logical(zoFort)),'-o');
    hold on
    %     scatter(imagingTime(~zoFort), ys(~zoFort),64,'.');
    plot(imagingTime(~zoFort), ys(~zoFort),'-o');
    legend({'ZO fort','ZO'})
    xlabel('Time (hours)')
    ylabel(ylabels{i})
    
    axis([0,max(imagingTime),0,max(ys)*1.1]);
end

print('summaryGraph.png','-dpng','-r200')
saveas(gca,'summaryGraph.fig','fig')
return

%% Multiple regression
% predictors = table([stats.cellNumber]', [metadata.StartTimeDifference]', X', Y', FEX', [stats.meanAmp]', ...
%     'VariableNames', ...
%     {'CellNumber', 'Time', 'Column', 'Row', 'FEX', 'Amplitude'});
% modelAmp = fitlm(predictors)
%
% predictors = table([stats.cellNumber]', [metadata.StartTimeDifference]', X', Y', FEX', [stats.meanBasal]', ...
%     'VariableNames', ...
%     {'CellNumber', 'Time', 'Column', 'Row', 'FEX', 'Basal'});
% modelB = fitlm(predictors)
%
% predictors = table([stats.cellNumber]', [metadata.StartTimeDifference]', X', Y', FEX', [stats.meanWHM]', ...
%     'VariableNames', ...
%     {'CellNumber', 'Time', 'Column', 'Row', 'FEX', 'WHM'});
% modelWHM = fitlm(predictors)
%
% predictors = table([stats.cellNumber]', [metadata.StartTimeDifference]', X', Y', FEX', [stats.meanfreq]', ...
%     'VariableNames', ...
%     {'CellNumber', 'Time', 'Column', 'Row', 'FEX', 'Freq'});
% modelFreq = fitlm(predictors,'quadratic')
%
% predictors = table([metadata.StartTimeDifference]', X', Y', FEX', [stats.cellNumber]', ...
%     'VariableNames', ...
%     {'Time', 'Column', 'Row', 'FEX', 'CellNumber'});
% modelN = fitlm(predictors)


% bounds = [min(predictors1.WHM),max(predictors1.WHM)];
% WHMs = linspace(bounds(1),bounds(2),10000);
% subplot(2,2,1)
% plot(WHMs,WHMs)
% title('Accounting for: WHM, Amp,B,Freq,N,Time,X,Y')
% hold on
% scatter(predictors1.WHM,model1.feval(predictors1(:,1:end-1)));
% xlabel('Actual WHM')
% ylabel('Predicted WHM')
%
% subplot(2,2,2)
% plot(WHMs,WHMs)
% title('Accounting for: WHM, Amp,B,Freq')
% hold on
% scatter(predictors2.WHM,model2.feval(predictors2(:,1:end-1)));
% xlabel('Actual WHM')
% ylabel('Predicted WHM')
%
% subplot(2,2,3)
% plot(WHMs,WHMs)
% title('Accounting for: WHM, Amp, Freq')
% hold on
% scatter(predictors3.WHM,model3.feval(predictors3(:,1:end-1)));
% xlabel('Actual WHM')
% ylabel('Predicted WHM')
%
% subplot(2,2,4)
% plot(WHMs,WHMs)
% title('Accounting for: WHM')
% hold on
% scatter(predictors4.WHM,model4.feval(predictors4(:,1:end-1)));
% xlabel('Actual WHM')
% ylabel('Predicted WHM')

% % % % predictors1 = table([stats.meanAmp]', [stats.meanBasal]', [stats.meanWHM]', [stats.meanfreq]', ...
% % % %     [stats.cellNumber]', [metadata.StartTimeDifference]', X', Y', FEX', 'VariableNames', ...
% % % %     {'Amp','Basal','WHM','Freq','CellNumber', 'Time', 'Column', 'Row', 'FEX'});
% % % %
% % % % model1 = fitlm(predictors1)
% % % %
% % % % predictors2 = table([stats.meanAmp]', [stats.meanBasal]', [stats.meanWHM]', [stats.meanfreq]', ...
% % % %     FEX', 'VariableNames', ...
% % % %     {'Amp','Basal','WHM','Freq','FEX'});
% % % %
% % % % model2 = fitlm(predictors2)
% % % %
% % % % predictors3 = table([stats.meanAmp]', [stats.meanWHM]', [stats.meanfreq]', ...
% % % %     FEX', 'VariableNames', ...
% % % %     {'Amp','WHM','Freq','FEX'});
% % % %
% % % % model3 = fitlm(predictors3)
% % % %
% % % % predictors4 = table([stats.meanWHM]',...
% % % %     FEX', 'VariableNames', ...
% % % %     {'WHM','FEX'});
% % % %
% % % % model4 = fitlm(predictors4)
% % % %
% % % % bounds = [min(predictors1.WHM),max(predictors1.WHM)];
% % % % WHMs = linspace(bounds(1),bounds(2),10000);
% % % % subplot(2,2,1)
% % % % plot(WHMs,WHMs)
% % % % title('Accounting for: WHM, Amp,B,Freq,N,Time,X,Y')
% % % % hold on
% % % % scatter(predictors1.WHM,model1.feval(predictors1(:,1:end-1)));
% % % % xlabel('Actual WHM')
% % % % ylabel('Predicted WHM')
% % % %
% % % % subplot(2,2,2)
% % % % plot(WHMs,WHMs)
% % % % title('Accounting for: WHM, Amp,B,Freq')
% % % % hold on
% % % % scatter(predictors2.WHM,model2.feval(predictors2(:,1:end-1)));
% % % % xlabel('Actual WHM')
% % % % ylabel('Predicted WHM')
% % % %
% % % % subplot(2,2,3)
% % % % plot(WHMs,WHMs)
% % % % title('Accounting for: WHM, Amp, Freq')
% % % % hold on
% % % % scatter(predictors3.WHM,model3.feval(predictors3(:,1:end-1)));
% % % % xlabel('Actual WHM')
% % % % ylabel('Predicted WHM')
% % % %
% % % % subplot(2,2,4)
% % % % plot(WHMs,WHMs)
% % % % title('Accounting for: WHM')
% % % % hold on
% % % % scatter(predictors4.WHM,model4.feval(predictors4(:,1:end-1)));
% % % % xlabel('Actual WHM')
% % % % ylabel('Predicted WHM')
% % % %
% % % % return

%% Plot time vs summary stat
plot(metadata.StartTimeDifference, [stats.meanfreq])

%% Make plots
figure(1)
subplot(2,2,1)
ys = [stats.meanAmp];
ys(isnan(ys)) = 0;
scatter(FEX, ys);
hold on
xs = unique(FEX);
for i = 1:length(xs)
    meanVal(i) = mean(ys(xs(i) == FEX));
    stdVal(i) = mean(ys(xs(i) == FEX)) / sqrt(sum(xs(i)==FEX));
end
errorbar(xs, meanVal, stdVal);
xlabel('FEX Concentration (%)')
ylabel('Amplitude (\DeltaI/I_o)')

subplot(2,2,2)
ys = [stats.meanfreq];
ys(isnan(ys)) = 0;
scatter(FEX, ys);
hold on
xs = unique(FEX);
for i = 1:length(xs)
    meanVal(i) = mean(ys(xs(i) == FEX));
    stdVal(i) = mean(ys(xs(i) == FEX)) / sqrt(sum(xs(i)==FEX));
end
errorbar(xs, meanVal, stdVal);
xlabel('FEX Concentration (%)')
ylabel('Frequency (mHz)')

subplot(2,2,3)
ys = [stats.meanBasal];
ys(isnan(ys)) = 0;
scatter(FEX, ys);
hold on
xs = unique(FEX);
for i = 1:length(xs)
    meanVal(i) = mean(ys(xs(i) == FEX));
    stdVal(i) = mean(ys(xs(i) == FEX)) / sqrt(sum(xs(i)==FEX));
end
errorbar(xs, meanVal, stdVal);
xlabel('FEX Concentration (%)')
ylabel('Basal (a.u.)')

subplot(2,2,4)
ys = [stats.meanWHM];
ys(isnan(ys)) = 0;
scatter(FEX, ys);
hold on
xs = unique(FEX);
for i = 1:length(xs)
    meanVal(i) = mean(ys(xs(i) == FEX));
    stdVal(i) = mean(ys(xs(i) == FEX)) / sqrt(sum(xs(i)==FEX));
end
errorbar(xs, meanVal, stdVal);
xlabel('FEX Concentration (%)')
ylabel('WHM (s)')

%% Regression analysis
predictors = table([stats.cellNumber]', [metadata.StartTimeDifference]', X', Y', FEX', ...
    'VariableNames', ...
    {'CellNumber', 'Time', 'Column', 'Row', 'FEX'});

response = table([stats.meanAmp]', [stats.meanBasal]', [stats.meanWHM]', [stats.meanfreq]', ...
    'VariableNames', {'Amplitude','Basal','WHM','Frequency'});

X = table2array(predictors);
Y = table2array(response);

X(isnan(X)) = 0;
Y(isnan(Y)) = 0;

[beta,Sigma,E,CovB,logL] = mvregress(X,Y)






















