
% Settings for TRAINING runs

%% Size of stuff

if RealRun
    
    FixOffset = 3.24;
    
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

T.On = 0.35;
T.Off = 0.40;
T.Fix = 6.75; % fixation after each 54 second miniblock
T.Delay = 15.0; % at the beginning and end of run
T.Count = 0.8;

TrialSet = round(1000*(T.On + T.Off)); % how much a trial should last
CompDur = round((T.Delay * 2) + (T.Count * 3) + ...
    ((T.On + T.Off) * RunTrials) + (4 * T.Fix)); % how much whole experiment should last

% in frames:

Frames = structfun(@(x) x * RefRate, T, 'UniformOutput', false);

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
    StimDir = '/Users/u010155/Documents/BlenderFiles/FinalRenders';
else
    StimDir = fullfile('Stimuli', 'Test');
end
