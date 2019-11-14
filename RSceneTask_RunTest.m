function RSceneTask_RunTest(SubjNo, RunInfo, DataDir)

% 
% 
% 
% 

%% Startup

global Environment
global RealRun

RunType = 'Training';

RunNo = RunInfo(1);
NRuns = RunInfo(2);
ThisRunNo = RunInfo(3);

% Load settings
RST_Settings_General;
RST_Settings_Test;

%% Set diary and other files

LogDir = fullfile(RunDir,'Logs');
if ~exist(LogDir,'dir')
    mkdir(LogDir);
end

RunStart = datestr(now,'dd-mm-yyyy_HH-MM-SS');
diary(fullfile(LogDir,sprintf('Subj%02d_%s_%d_%s.txt', SubjNo, RunType, ThisRunNo, RunStart)));

DataFile = fullfile(RunDir, sprintf('Subj%02d_%s_%g.mat', SubjNo, RunType, ThisRunNo));

BUpDir = fullfile(RunDir,'Backup');
if ~exist(BUpDir,'dir')
    mkdir(BUpDir);
end
BUpFile = fullfile(BUpDir, sprintf('Subj%02d_%s_%g_%s.mat', SubjNo, RunType, ThisRunNo, RunStart));

%% Get and report flip interval, and compensate waiting times in settings

ifi = Screen('GetFlipInterval', w);
fprintf('\nFlip interval: %2.4g ms (should be %2.4g ms).\n\n', ifi*1000, (1000/RefRate));

if ifi > (1/RefRate)*1.5       % if flip interval is way too slow do not adapt timings
    warning( 'This is a serious timing problem! \n');
end

%% Message
% while loading images and waiting for scanner trigger

InstrTxt = sprintf(['Run %g/%g\n\nDoes the object turn\nclockwise or counterclockwise?'], RunNo, NRuns);

DrawFormattedText(w, InstrTxt, 'center', 'center', White);

Screen('Flip', w);

