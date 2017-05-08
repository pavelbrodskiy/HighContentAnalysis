%% Initialization
clearvars
close all
settings = prepareWorkspace();

%% Declare constants
letters = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'};

%% Convert Low Frequency data to z-projected videos
tblLowFrequency = readtable('C:\ResearchTemp\Experiment7-Blebbistatin\Inputs\LowFreqLabels.xlsx');
processConfocalData(tblLowFrequency.Label, settings)

%% Processess inputs
tblHighFrequency = readtable('C:\ResearchTemp\Experiment7-Blebbistatin\Inputs\HighFreqLabels.xlsx');
[~,~,tblPlateMap] = xlsread('C:\ResearchTemp\Experiment7-Blebbistatin\Inputs\PlateLayout.xlsx');
tblPlateLegend = readtable('C:\ResearchTemp\Experiment7-Blebbistatin\Inputs\PlateLayoutLegend.xlsx');

fields = tblPlateLegend.Properties.VariableNames(2:end);

tblHighFrequency.Conditions = cellfun(@(x) {num2str(x)},tblPlateMap(sub2ind(size(tblPlateMap),tblHighFrequency.WellY+1,tblHighFrequency.WellX+1)));
conditionIndex = cell2mat(cellfun(@(x,y) find(~cellfun('isempty',strfind(y,x))), ...
    tblHighFrequency.Conditions, repmat({tblPlateLegend.Var1}, [1, length(tblHighFrequency.Conditions)])', 'UniformOutput',false));

for i = 1:length(fields)
    tblHighFrequency.(fields{i}) = tblPlateLegend.(fields{i})(conditionIndex);
end

tblHighFrequency.PlateAddress = arrayfun(@(x,y) [letters{x} num2str(y)], tblHighFrequency.WellY, tblHighFrequency.WellX, 'UniformOutput', 0);

mkdir('MP4');

%% Make Low Frequency movies
s = 32;
totalZProj = [];
for j = 1:s
    % Obtain condition information
    condition = cellfun(@(x) {num2str(x)},tblPlateMap(sub2ind(size(tblPlateMap),find(contains(letters, coords(1)))+1,str2num(coords(2:end))+1)));
    conditionIndex = cell2mat(cellfun(@(x,y) find(~cellfun('isempty',strfind(y,x))), ...
        tblHighFrequency.Conditions, repmat({tblPlateLegend.Var1}, [1, length(tblHighFrequency.Conditions)])', 'UniformOutput',false));
    
    for i = 1:length(fields)
        lowFreq.(fields{i}) = tblPlateLegend.(fields{i})(find(contains(tblPlateLegend.Var1, condition)));
    end
    
    filename = ['MP4' filesep 'Cl8_ZOf_' num2str(lowFreq.Blebbistatin_uM) ...
        '_uM_Bleb_' num2str(lowFreq.DMSO_percent) '_pctDMSO_'...
        coords '_LowFreq'];
    
    if exist([filename, '.mp4'],'file')
        continue
    end
    
    % Load raw data
    for i = 1:length(tblLowFrequency.Label)
        tmp = load([settings.thruData 'LowFrequency_s' num2str(j) '_f' num2str(i) '.mat'], 'zProj');
        totalZProj = cat(3, totalZProj, tmp.zProj);
    end
    
    % Convert into RGB shape
    totalZProj = totalZProj(:,:,:,[2,1]);
    totalZProj(:,:,:,3) = 0;
    
    writeMP4(permute(zProj, [1,2,4,3]), filename, 4);
end

return

for i = 1:size(raw,1)
    str = raw{i,2}.get(['Global Stage' num2str(i)]);
    coords = str(2:end-1);
    
    tmp = raw{i, 1};
    tmp2 = raw2{i, 1};
    zStack = cat(3, tmp{:,1}, tmp2{:,1});
    ts = length(zStack(:)) / 512 / 512 / zs / cs;
    hyperStack = reshape(zStack, 512, 512, zs, cs, ts);
    
    condition = cellfun(@(x) {num2str(x)},tblPlateMap(sub2ind(size(tblPlateMap),find(contains(letters, coords(1)))+1,str2num(coords(2:end))+1)));
    conditionIndex = cell2mat(cellfun(@(x,y) find(~cellfun('isempty',strfind(y,x))), ...
        tblHighFrequency.Conditions, repmat({tblPlateLegend.Var1}, [1, length(tblHighFrequency.Conditions)])', 'UniformOutput',false));
    
    for j = 1:length(fields)
        lowFreq.(fields{j}) = tblPlateLegend.(fields{j})(find(contains(tblPlateLegend.Var1, condition)));
    end
    
    zProj = squeeze(max(hyperStack, [], 3));
    zProj = zProj(:,:,[2,1],:);
    zProj(:,:,3,:) = 0;
    
    filename = ['MP4' filesep 'Cl8_ZOf_' num2str(lowFreq.Blebbistatin_uM) ...
        '_uM_Bleb_' num2str(lowFreq.DMSO_percent) '_pctDMSO_'...
        coords '_LowFreq'];
    
    if true %~exist([filename, '.mp4'],'file')
        writeMP4(permute(zProj, [1,2,4,3]), filename, 4);
    end
end

%% Make Low Frequency movies
% if ~exist([settings.thruData 'LowFreq.mat'], 'file')
%     raw = bfopen(beforeTreatmentDir);
%     raw2 = bfopen(afterTreatmentDir);
%     save([settings.thruData 'LowFreq.mat'], 'raw', 'raw2');
% else
%     load([settings.thruData 'LowFreq.mat'], 'raw', 'raw2');
% end
%
% zs = 11;
% cs = 2;
%
% for i = 1:size(raw,1)
%     str = raw{i,2}.get(['Global Stage' num2str(i)]);
%     coords = str(2:end-1);
%
%     tmp = raw{i, 1};
%     tmp2 = raw2{i, 1};
%     zStack = cat(3, tmp{:,1}, tmp2{:,1});
%     ts = length(zStack(:)) / 512 / 512 / zs / cs;
%     hyperStack = reshape(zStack, 512, 512, zs, cs, ts);
%
%     condition = cellfun(@(x) {num2str(x)},tblPlateMap(sub2ind(size(tblPlateMap),find(contains(letters, coords(1)))+1,str2num(coords(2:end))+1)));
%     conditionIndex = cell2mat(cellfun(@(x,y) find(~cellfun('isempty',strfind(y,x))), ...
%         tblHighFrequency.Conditions, repmat({tblPlateLegend.Var1}, [1, length(tblHighFrequency.Conditions)])', 'UniformOutput',false));
%
%     for j = 1:length(fields)
%         lowFreq.(fields{j}) = tblPlateLegend.(fields{j})(find(contains(tblPlateLegend.Var1, condition)));
%     end
%
%     zProj = squeeze(max(hyperStack, [], 3));
%     zProj = zProj(:,:,[2,1],:);
%     zProj(:,:,3,:) = 0;
%
%     filename = ['MP4' filesep 'Cl8_ZOf_' num2str(lowFreq.Blebbistatin_uM) ...
%         '_uM_Bleb_' num2str(lowFreq.DMSO_percent) '_pctDMSO_'...
%         coords '_LowFreq'];
%
%     if true %~exist([filename, '.mp4'],'file')
%         writeMP4(permute(zProj, [1,2,4,3]), filename, 4);
%     end
% end

%% Make high frequency movies
for i = 1:length(tblHighFrequency.Label)
    if ~exist([settings.thruData tblHighFrequency.Label{i} '.mat'], 'file')
        raw = bfopen([settings.inData 'HighFrequency' filesep tblHighFrequency.Label{i} '.tif']);
        raw = cat(1,raw{:,1});
        zStack = cat(3,raw{:,1});
        mkdir(fileparts([settings.thruData tblHighFrequency.Label{i}]))
        save([settings.thruData tblHighFrequency.Label{i} '.mat'], 'zStack');
    else
        load([settings.thruData tblHighFrequency.Label{i} '.mat'], 'zStack');
    end
    
    filename = ['MP4' filesep 'Cl8_ZOf_' num2str(tblHighFrequency.Blebbistatin_uM(i)) ...
        '_uM_Bleb_' num2str(tblHighFrequency.DMSO_percent(i)) '_pctDMSO_' ...
        tblHighFrequency.PlateAddress{i} '_day_' num2str(round(tblHighFrequency.Days(2)*100)/100) '_HighFreq'];
    
    if ~exist([filename, '.mp4'],'file')
        writeMP4(zStack, filename, 60);
    end
end

return

%% Extract statistics
settings.force = false;
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






















