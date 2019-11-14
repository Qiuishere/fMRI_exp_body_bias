
% Settings for TEST runs

%% Feedback
% 

GiveFB = 'Run'; % Run / Trial / None

%% Size of stuff

if RealRun
    
    FixOffset = 3.24; % 150
    
    ImH = 11.6448; % 540
    ImW = 20.5321; % 960
    
    XMargin = 1.08; % 50 % how many extra pixels around bounding box
    YMargin = 0.43; % 20
    
else
    
    FixOffset = 150;
    
    ImH = 540;
    ImW = 960;
    
    ImSize = [ImW ImH];
    
    XMargin = 50; % how many extra pixels around bounding box
    YMargin = 20;
    
end

%% Timing

% Trial:
T.Fix = 1.0;
T.Im = 0.5;
T.First = 2.0;
T.Final = 2.0;
T.Probe = 0.05;
T.ISI = 0.1;
%T.FB = 0.25; % timing of trial-by-trial feedback (if given)