
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
T.Fix = 1.5; % maximum
T.Seq(1) = 2.0;
T.Seq(2:4) = 0.5;
T.Seq(5) = 1.5;
T.Probe = 0.05;
T.ISI = 0.1;
T.PreRespDelay = 0.05; % very short pause before response screen
T.Resp = 1.5; % max. time to give a response

T.Delay = 15.0; % at the start and end of run
T.Count = 0.8;

T.RunFB = 2.0;
T.TrialFB = 0.25; % timing of trial-by-trial feedback (if given)

% in frames:

Frames = structfun(@(x) x * RefRate, T, 'UniformOutput', false);

TrialSequence = [T.Fix, T.Seq, T.Probe, T.ISI, T.Probe, T.PreRespDelay, T.Resp];
TrialSet = sum(TrialSequence); % how much a trial should last
             
%% Conversions (DVA to pixels)

if RealRun
    
    [FixOffset, ~, ~] = ang2pix3(FixOffset, ViewD, [ScreenResolution.width, ScreenResolution.height], ScreenSize);
    FixOffset = round(FixOffset(2));
    
    [ImSize, ~, ~] = ang2pix3([ImW ImH], ViewD, [ScreenResolution.width, ScreenResolution.height], ScreenSize);
    ImSize = round(ImSize);
    
    [XMargin, ~, ~] = ang2pix3(XMargin, ViewD, [ScreenResolution.width, ScreenResolution.height], ScreenSize);
    XMargin = round(XMargin(1));
    
    [YMargin, ~, ~] = ang2pix3(YMargin, ViewD, [ScreenResolution.width, ScreenResolution.height], ScreenSize);
    YMargin = round(YMargin(2));
    
end

%% Stimulus position (rect)

ImRect = CenterRectOnPoint([0 0 ImSize], xCenter, yCenter); 

%% Fixation

FixRct = CenterRectOnPoint([0 0 FixSize FixSize], xCenter, yCenter + FixOffset);

%% Stimulus directory

if IsOSX
    StimDir_main = '/Users/u010155/Documents/BlenderFiles/FinalRenders';
    StimDir_rot = '/Users/u010155/Documents/BlenderFiles/RotScene/Rotate_half';
else
    StimDir_main = fullfile('Stimuli', 'Test');
    StimDir_rot = fullfile(StimDir_main, 'Rotate');
end
