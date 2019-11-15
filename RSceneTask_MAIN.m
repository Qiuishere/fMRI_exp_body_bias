%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%       Updating of visual templates from high-level scene context        %
%                                                                         %
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%                                                                         %
%    Giacomo Aldegheri, Surya Gayet & Marius Peelen                       %
%    (DCCN host: Floris de Lange)                                         %
%    Partly based on Surya Gayet's code                                   %
%                                                                         %
%-------------------------------------------------------------------------%
%                                                                         %
%  In this new version, the target always appears after occlusion.        %
%  We will compare decoding accuracies between expected and unexpected    %
%  trials.                                                                %
%                                                                         %
%                                                                         %
%                                                                         %
%                                                                         %
%                                                                         %
%                                                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize

addpath(genpath('Functions'));
StartUp;

global Environment
Environment = 1; % 1 - office, 2 - dummy, 3 - prismafit

global RealRun
RealRun = Environment == 2 || Environment == 3;

%% Determine (main) data directory

if RealRun, DataDir = 'Data_Lab'; else, DataDir = 'Data_Demo'; end

%% List of experiment runs
% Test, Training, Localizer, Demo

Runs = {'Training'};

%% Participant info

SubjNo = input('Please enter participant number:\n');

Confirm_subj = input(sprintf('\n\nParticipant no. %g. Confirm? [y / n]\n', SubjNo),'s');

if strcmp(Confirm_subj,'n')
    return;
elseif ~strcmp(Confirm_subj,'y')
    fprintf('Please enter y or n.\n');
    Confirm_subj = input(sprintf('\n\nParticipant no. %g. Confirm? [y / n]\n', SubjNo),'s');
end

fprintf('\nRUNS:\n----------------\n%s----------------\n', sprintf(' %s\n', Runs{:}));

% Estimate next run

[RunNo, NOccur] = EstiRun(SubjNo, Runs, DataDir);
NextRun = Runs{RunNo};

Confirm_run = input(sprintf('\n\nNext run: %g - %s %g. Confirm? [y / n]\n', RunNo, NextRun, NOccur(RunNo)),'s');

if strcmp(Confirm_run,'y')
elseif strcmp(Confirm_run,'n')
    RunNo = input('\nPlease enter next run number:\n');
    NextRun = Runs{RunNo};
end

RunInfo = [RunNo, numel(Runs), NOccur(RunNo)]; % which run, how many runs in total, which occurrence of run type

%% Run experiment

switch NextRun
    case 'Training'
        RSceneTask_RunTraining(SubjNo, RunInfo, DataDir);
    case 'Test'
        RSceneTask_RunTest(SubjNo, RunInfo, DataDir);
end
