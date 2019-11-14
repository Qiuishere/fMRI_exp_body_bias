
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
T.Resp = 1.5; % max. time to give a response

T.Delay = 15.0; % at the start and end of run
switch GiveFB
    case 'Run'
        T.FB = 1.0;
    case 'Trial'
        %T.FB = 0.25; % timing of trial-by-trial feedback (if given)
end

% in frames:

Frames = structfun(@(x) x * RefRate, T, 'UniformOutput', false);

TrialSequence = [T.Fix, T.First, T.Im, T.Im, T.Im, ...
                 T.Final, T.Probe, T.ISI, T.Probe, T.Resp];
             
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


