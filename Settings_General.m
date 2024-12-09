% GENERAL SETTINGS
%

%% Keyboard

KbName('UnifyKeyNames')
% participant keys
if RealRun
    Keys.name = {'Index', 'Middle'};
    Keys.number = [97, 98];
else
    Keys.name = {'F', 'J'};
    Keys.number = [KbName('F') KbName('J')];
end

% randomize key correspondence
if mod(SubNo,2)
    Keys = structfun(@fliplr,Keys,'UniformOutput',false);
end


RespKeys = Keys.number;
% restrict keys that are read out by KbCheck (speeds up KbCheck)
enabledKeys = [RespKeys KbName('ESCAPE') KbName('SPACE')];
RestrictKeysForKbCheck(enabledKeys);

%% Open screen & settings==================================================
sca;
AssertOpenGL;

ScreenID = max(Screen('Screens'));
monitor.Resolution = Screen('Resolution', ScreenID);

White = WhiteIndex(ScreenID);
Grey = White/2;
Black = BlackIndex(ScreenID);

text.Size = 24;
text.Color = White;
text.Font = 'Arial';

 
w = SetScreen('Window', ScreenID,'OpenGL', 1, 'BGColor',Grey, 'Blending', 'Transparent',...
     'FontSize',text.Size,'FontColor', text.Color, 'FontType',text.Font);

topPriorityLevel = MaxPriority(w.Number);
Priority(topPriorityLevel);

if RealRun
    HideCursor(ScreenID);
    
end
% monitor parameter

switch Environment
    
    case 1 % Office        
        correctRefRate = 60;
        monitor.Size = [527 296]; % in mm
        monitor.ViewDist = 1206;
        
    case 2 % Dummy scanner/projector        
        correctRefRate = 60;
        correctRes = [1280 720]; % This is the Resolution we are supposed to use
        monitor.Size = [369 277]; % in mm
        monitor.ViewDist = 955; % viewing distance in mm
        
    case 3 % Prismafit        
        % All details found in 'Optical_path_Skyra_BOLDscreen.pdf'
        correctRefRate = 60;
        correctRes = [1920 1080];
        monitor.Size = [698.4 392.9]; % in mm
        monitor.ViewDist = 1460; % viewing distance in mm (1340 from screen + 100 from mirror + 20 screen glass)
        
end

if RealRun
    if (w.RefreshRate - correctRefRate) >= 1
        error('Framerate is not %g!', correctRefRate);
    end
    if ~isequal([w.Width, w.Height], correctRes)
        error('Resolution is not %g x %g!', correctRes(1), correctRes(2))
    end
end



%% Make directories for participant (if needed)============================

RunDir = fullfile(DataDir, sprintf('Sub%02d', SubNo), RunType);

if ~exist(RunDir, 'dir')
    mkdir(RunDir);
end

%% Scanner settings========================================================

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


%% Visual content==========================================================

% Fixation size

if RealRun
    FixSize = 0.15;
else
    FixSize = 0.15;
end

monitor.Center = CenterRect([0 0 1 1],w.Rect);% this is only true when window is full screen

fix.Size = visAng2pix(FixSize, 0, monitor);
fix.Col =  [30 80 200];
fix.Col2 =  [200 0 0];

if RealRun
    time.before = 12*1000;
else
    time.before = 6*1000;
end

time.countdown = 3;