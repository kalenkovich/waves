%% Start Brainstorm
run('brainstorm/brainstorm3/brainstorm');

%% Add external toolboxes
addpath(genpath('gptoolbox-master'));

%% Load the last session file
filenames = dir();
filenames = {filenames.name};
filenames = regexp(filenames ,'\w*.mat$', 'match');
filenames = sort([filenames{:}]);
load(filenames{end});
clear filenames;
%

%% open main.m for editing
edit main;

