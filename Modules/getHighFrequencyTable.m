function [tblHighFrequency, tblPlateMap, tblPlateLegend] = getHighFrequencyTable(settings)
tblHighFrequency = readtable(settings.tblHighFreq);
[~,~,tblPlateMap] = xlsread(settings.tblPlateMap);
tblPlateLegend = readtable(settings.tblPlateLegend);

fields = tblPlateLegend.Properties.VariableNames(2:end);

tblHighFrequency.Conditions = cellfun(@(x) {num2str(x)},tblPlateMap(sub2ind(size(tblPlateMap),tblHighFrequency.WellY+1,tblHighFrequency.WellX+1)));
conditionIndex = cell2mat(cellfun(@(x,y) find(~cellfun('isempty',strfind(y,x))), ...
    tblHighFrequency.Conditions, repmat({tblPlateLegend.Var1}, [1, length(tblHighFrequency.Conditions)])', 'UniformOutput',false));

for i = 1:length(fields)
    tblHighFrequency.(fields{i}) = tblPlateLegend.(fields{i})(conditionIndex);
end

tblHighFrequency.PlateAddress = arrayfun(@(x,y) [settings.letters{x} num2str(y)], tblHighFrequency.WellY, tblHighFrequency.WellX, 'UniformOutput', 0);

