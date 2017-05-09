% This function creates the settings structure used in the rest of the
% analysis tools. This makes it easier to keep track of all of the input
% and output files, as well as throughput data that is generated once and
% used many times.

function settings = getSettings()
%% Get root directories
[activeDir, ~, ~] = fileparts(mfilename('fullpath'));
settings.activeDir = [activeDir filesep];
settings.rootData = [activeDir filesep 'dataProcessing' filesep];
settings.output = [activeDir filesep 'Output' filesep];
settings.outMP4 = [settings.output 'MP4' filesep];

%% Get input directories
settings.inData = [settings.activeDir 'Data' filesep];
settings.inTables = [settings.activeDir 'Inputs' filesep];

settings.tblHighFreq = [settings.inTables 'HighFreqLabels.xlsx'];
settings.tblLowFreq = [settings.inTables 'LowFreqLabels.xlsx'];
settings.tblPlateMap = [settings.inTables 'PlateLayout.xlsx'];
settings.tblPlateLegend = [settings.inTables 'PlateLayoutLegend.xlsx'];

%% Get throughput directories
settings.thruData = [settings.rootData 'thruData' filesep];
settings.thruStats = [settings.rootData 'tempStatistics' filesep];
settings.thruMedian = [settings.rootData 'tempMedian' filesep];

%% Get output directories
settings.outRough = settings.output;
settings.outFinal = [settings.output 'outFigurePanels' filesep];

%% Get dependency directories
settings.depExt = [settings.activeDir 'Dependencies' filesep];
settings.depInt = [settings.activeDir 'Modules' filesep];

%% Add dependancy folders to path
addpath(genpath(settings.depExt))
addpath(genpath(settings.depInt))

%% Generate unique identifier for analysis
currentTime = now();
timeString = datestr(currentTime);
settings.uniqueIdentifier = strrep(timeString,':','-');

%% Set analysis settings
settings.force = false;             % Set true to rerun completed analysis
settings.scale.x20 = 0.7009;        % Scale of 20x objective (px/um)
settings.background = zeros(512);
settings.weiner = 12;
settings.minProminence = 0.1;
settings.timestep = 10;

settings.letters = {'A','B','C','D','E','F','G','H'};

