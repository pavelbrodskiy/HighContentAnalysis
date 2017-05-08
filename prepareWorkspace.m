% Run this at the beginning of each calcium wave analysis script.
% This identifies the dependancy folders, and sets folders for
% inputs and output, and creates them if they do not exist.
%
% prepareWorkspace()

function settings = prepareWorkspace(varargin)
% Clear off the workspace
tic
close all
fclose('all');

settings = getSettings();
cd(settings.activeDir);

end