% GENERAL SETTINGS
%
%
%
%
%
%
%

if RealRun
    UseEyeTracker = 0;
else
    UseEyeTracker = 0;
    PsychDebugWindowConfiguration;
end

%% Visual context

% Fixation size

if RealRun
    FixSize = 0.2;
else
    FixSize = 15;
end

%% Keyboard

% control room keys

if IsWin
    EscKey      =   27;          % Escape
    SpaceBar      =   32;
    textfont   =   'Calibri';
elseif IsOSX
    EscKey       =   41;         % Escape
    SpaceBar       =   44;
    textfont    =   'Arial';
end

% participant keys

if RealRun
    RespKeys        =    [97 98 99];    %   INDEX - MIDDLE - RING (R)
    enabledKeys    =    [RespKeys EscKey SpaceBar];
else
    if IsWin
        RespKeys    =    [37 39];    %   (LR)
    elseif IsOSX
        RespKeys    =    [80 79];    %   (LR)
    end
    enabledKeys = [RespKeys EscKey SpaceBar];
end

% restrict keys that are read out by KbCheck (speeds up KbCheck)
RestrictKeysForKbCheck(enabledKeys);

%% Open screen & settings


Screen('Preference', 'SkipSyncTests', 2);
sca;

ScreenID = max(Screen('Screens'));

White = WhiteIndex(ScreenID);
Grey = White/2;
Black = BlackIndex(ScreenID);

ScreenResolution = Screen('Resolution', ScreenID);

[w, windowRect] = PsychImaging('OpenWindow', ScreenID, Grey, ...
    [], [], []);

ifi = Screen('GetFlipInterval', w);

% [screenXpixels, screenYpixels] = Screen('WindowSize',w);

[xCenter, yCenter] = RectCenter(windowRect);

topPriorityLevel = MaxPriority(w);
Priority(topPriorityLevel);

Screen('BlendFunction', w, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
if RealRun
    HideCursor(ScreenID);
end

% Text size

Screen('TextSize', w, 22);

% Screen settings

RefRate = 60;

if RealRun

    MyScreenRes = [1024 768];
    ScreenSize = [369 277]; % in mm
    ViewD = 955; % viewing distance in mm
    if (1 / ifi - RefRate) >= 1
        error('Framerate is not %g!', RefRate);
    end
    if ~isequal([ScreenResolution.width, ScreenResolution.height], MyScreenRes)
       error('Resolution is not %g x %g!', MyScreenRes(1), MyScreenRes(2))
    end

end

%% Make directories for participant (if needed)

RunDir = fullfile(DataDir,sprintf('Subj%02d', SubjNo), RunType);

if ~exist(RunDir, 'dir')
     mkdir(RunDir);
end

%% Scanner settings

% clear ports
if ~isempty(instrfind)
    fclose(instrfind);
    delete(instrfind);
end

% set usage of com-ports
if RealRun
    BitsiBB = bitsi_scanner('com2');
    BitsiScanner = bitsi_scanner('com3');
else
    BitsiBB = bitsi_scanner('');
end
timeout = 0.001;

%% Conversions (DVA to pixels)

if RealRun

    [FixSize, ~, ~] = ang2pix3(FixSize, ViewD, [ScreenResolution.width, ScreenResolution.height], ScreenSize);
    FixSize = round(FixSize(1));

end
