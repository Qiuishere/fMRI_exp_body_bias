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
%RST_Settings_Test;