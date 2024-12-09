%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% fMRI study for probing the neural loci of biases in body representation %
%                                                                         %
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%                                                                         %
%    Qiu Han, Marco Gandolfo & Marius Peelen                              %
%    (DCCN host: Floris de Lange)                                         %
%    Partly based on Surya Gayet's code                                   %
%                                                                         %
%-------------------------------------------------------------------------%
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize

addpath(genpath('Functions'));
StartUp;

global Environment
Environment = 1; % 1 - office, 2 - dummy/beamer, 3 - skyra/prismafit (new screen)

global RealRun
RealRun = Environment == 2 || Environment == 3;

% Determine (main) data directory
if RealRun, DataDir = 'Data_Lab'; else, DataDir = 'Data_Demo'; end


SubNo = 4;


%% define experiment Runs

% experiment information
Runs = {'Train', ...
        'Test', ...
        'Test', ...
        'Localizer_Obj', ...
        'Localizer_BM', ...
        'Test', ...
        'Localizer_Obj', ...
        'Train', ...
        'Test', ...
        'Test'};

% get the current runno from data folder, as the default input
[RunNo, NOccur] = EstiRun(SubNo, Runs, DataDir); % NOccur returns the index 
% of each run inside their own category (e.g. the second training run)

% enquire about current run
prompt      =   {'Enter next run number:',...
    'Enter Participant number:',...
    'Which hand do you use the most? (left:1, right:2)'
    };
title       =   sprintf('Participant %g', SubNo);
defInput    =   {num2str( RunNo), num2str(SubNo), num2str(1)};
Input       =   inputdlg( prompt, title, [1 60], defInput);
RunNo       =   str2double(Input(1));

NextRun = Runs{RunNo};

RunInfo = [RunNo, numel(Runs), NOccur(RunNo)]; % which run, how many runs in total, which occurrence of run type

%% Run experiment

switch NextRun

    case 'Test'
        Body_RunRS(SubNo, RunInfo, DataDir);
    case 'Train'
        Body_RunTrain(SubNo, RunInfo, DataDir);
    case 'Localizer_Obj'
        Body_RunLocalizer_Obj(SubNo, RunInfo, DataDir);
    case 'Localizer_BM'
        Body_RunLocalizer_BM(SubNo, RunInfo, DataDir);
end

clearvars; sca; ShowCursor;
