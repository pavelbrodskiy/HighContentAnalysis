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
% settings.outMP4 = [settings.output 'MP4' filesep];

%% Get input directories
settings.inData = [settings.activeDir 'Confocal_6PMto10AM' filesep];

settings.inData = [settings.activeDir 'Data' filesep];
% settings.inLowMag = [settings.inExperimentalData 'WholeDiscImages' filesep];
% settings.inExperimentalData = [settings.activeDir 'RawTiffStacks' filesep];
settings.inTables = [settings.activeDir 'Inputs' filesep];

settings.labelTabel = [settings.inTables 'lableTable.xlsx'];
% settings.rejectTabel = [settings.inTables filesep 'labelReject.xlsx'];
% settings.superLabelTabel = [settings.inTables filesep 'superCategoryReference.xlsx'];

%% Get throughput directories
settings.thruData = [settings.rootData 'thruData' filesep];
% settings.thruTime = [settings.rootData 'tempTime' filesep];
settings.thruStats = [settings.rootData 'tempStatistics' filesep];
% settings.thruPSD = [settings.rootData 'tempPSD' filesep];
settings.thruMedian = [settings.rootData 'tempMedian' filesep];
% settings.thruRot = [settings.rootData 'tempRotated' filesep];
% settings.thruRawMask = [settings.rootData 'tempRawMasks' filesep];
% settings.outArchive = [settings.rootData 'archive' filesep];
% settings.matProfiles = [settings.rootData 'tempProfiles.mat'];
% settings.matMasks = [settings.rootData 'tempManualMasks.mat'];
% settings.tmp = [settings.rootData 'tmp' filesep];

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
% settings.analysis.binSize = 4;      % Size of spatial bins in pixels
% settings.sizeStandard = 1;          % Standard scale of output images (pix in/pix out)
settings.force = false;             % Set true to rerun completed analysis
settings.scale.x20 = 0.7009;        % Scale of 20x objective (px/um)
settings.background = zeros(512);
settings.weiner = 12;
settings.minProminence = 0.1;
settings.timestep = 10;
% settings.minTime = 180;            % Minimum time for analysis (frames)

%% Default color maps for figures
% settings.colorMap.amp       = makeCMap('hot*');
% settings.colorMap.freq      = makeCMap('bone*');
% settings.colorMap.WHM       = makeCMap('copper*');
% settings.colorMap.dtyCyc    = parula(255);
% settings.colorMap.basal     = makeCMap('dusk*');
end