%
%
% adapted from 'Orient_Task_BlockFB_Demo' - added saving data
% 
% 
% 
% Adaptive staircase. The scenes have 2 different initial
% viewpoints (AllTrials.InitView = 1/2) and different final viewpoints (AllTrials.FinalView = 30/90)
% but in practice all that matters is whether they are consistent or not
% (AllTrials.Consistent = 1/0)
%
%
%

sca
clear all
close all
commandwindow
addpath(genpath('Functions'));

[seed, whichGen] = ClockRandSeed;
AssertOpenGL;

%% Settings

GiveFeedback = 1;

%% Environment

Environment = 'Office';

switch Environment
    case 'Office'
        PsychDebugWindowConfiguration;
        HideCursor;
        Screen('Preference','SkipSyncTests', 1);
        ScreenSize = [475.9 267.7];
        ViewD = 500;
        StimDir_main = '/Users/u010155/Documents/BlenderFiles/FinalRenders';
        StimDir_rot = '/Users/u010155/Documents/BlenderFiles/RotScene/Rotate_half';
    case 'Lab'
        % ...
        %PsychDebugWindowConfiguration;
        HideCursor;
        ScreenSize = [530 298];
        ViewD = 570;
        StimDir_main = fullfile('..', 'Stimuli');
        StimDir_rot = fullfile(StimDir_main, 'Rotate');
end

%% Experimental variables

TotNTrials = 10;

NScenes = 20;

%% Staircase settings

% Here the maximum difference is set at 20 (very reasonable)
% note that if you want to set it bigger you need to make sure that
% it doesn't choose non-existing images (>20 or <-20).

Diff_Start = 10;
Max_Diff = 20;
Min_Diff = 0.5;
Step_Size = 0.5;
NDown = 2;
OrientLims = [-10, 10];

%% AllTrials table

Views = {'A', 'B'};

Variables = {'Scene', 'InitView', 'FinalView', 'Diff', ...
    'Orients', 'Hit', 'RT'};

% Automatically if trials are evenly split between initial views, and then
% between final views, half will be wide and half narrow.

AllTrials = nan(TotNTrials, numel(Variables));

AllTrials(:, 2) = randi(2, TotNTrials, 1);
AllTrials(:, 3) = ones(TotNTrials, 1) * 30 + (rand(TotNTrials, 1)>0.5) * 60;

AllTrials = array2table(AllTrials, 'VariableNames', Variables);

%% Choose scenes and shuffle trials
% DO THE SCENES NEED TO BE ORDERED IN A SPECIFIC WAY? (to check for memory
% effects)
% Scenes here are assigned completely randomly:
AllTrials.Scene = randi(NScenes, TotNTrials, 1);

AllTrials = AllTrials(randperm(TotNTrials), :);

% Add sequence:

Seq = nan(TotNTrials, 5);

Seq(:, 1) = 0;
for i = 1:TotNTrials
    if AllTrials.FinalView(i)==30, Limit = 25; else Limit = 60; end
    Seq(i, 2:4) = sort(randsample(15:5:Limit, 3));
end
Seq(:, 5) = AllTrials.FinalView;

AllTrials.Sequence = num2cell(Seq, 2);

%% Assign initial difference

if rand > 0.5
    AllTrials.Diff(1) = Diff_Start;
else
    AllTrials.Diff(1) = -Diff_Start;
end

AllTrials.Orients = num2cell(AllTrials.Orients);

% Decide orientations in the first trial based on initial difference:
AllTrials.Orients{1}(1) = round(normrnd(0, 1));
AllTrials.Orients{1}(2) = max(min(AllTrials.Orients{1}(1) + AllTrials.Diff(1), ...
OrientLims(2)), OrientLims(1));

%% Data path

%% Keyboard

KbName('UnifyKeyNames');

escapeKey = KbName('ESCAPE');
KeyL = KbName('LeftArrow'); KeyR = KbName('RightArrow'); % consistent/inconsistent response keys

UsedKeys = [escapeKey, KeyL, KeyR];

%% Screen stuff

ScreenID = max(Screen('Screens'));

White = WhiteIndex(ScreenID);
Black = BlackIndex(ScreenID);
Grey = round((White+Black)/2);

[w, windowRect] = PsychImaging('OpenWindow', ScreenID, Grey, ...
    [], [], []);

ScreenResolution = Screen('Resolution', ScreenID);

ifi = Screen('GetFlipInterval', w);

[xCenter, yCenter] = RectCenter(windowRect);

topPriorityLevel = MaxPriority(w);
Priority(topPriorityLevel);

Screen('BlendFunction', w, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

%% Size of stuff

FixSize = 0.2;

FixOffset = 3.24; % 150

ImH = 11.6448; % 540
ImW = 20.5321; % 960

XMargin = 1.08; % 50 % how many extra pixels around bounding box
YMargin = 0.43; % 20

if ismac
    FixSize = FixSize/2;
    FixOffset = FixOffset/2;
    ImH = ImH/2;
    ImW = ImW/2;
    XMargin = XMargin/2;
    YMargin = YMargin/2;
end

%% Timing
%  (in seconds)

T.Fix = 0.5;
T.Seq(1) = 2.0;
T.Seq(2:4) = 0.5;
T.Seq(5) = 1.5;
T.Probe = 0.05;
T.ISI = 0.1;
T.FB = 0.25;

AllTrials.Jitter = round((T.Seq(5) - (rand(TotNTrials, 1) * 0.5))/ifi);

%% Conversions (DVA to pixels, seconds to frames)

[FixSize, ~, ~] = ang2pix3(FixSize, ViewD, [ScreenResolution.width, ScreenResolution.height], ScreenSize);
FixSize = round(FixSize(1));

[FixOffset, ~, ~] = ang2pix3(FixOffset, ViewD, [ScreenResolution.width, ScreenResolution.height], ScreenSize);
FixOffset = round(FixOffset(2));

[ImSize, ~, ~] = ang2pix3([ImW ImH], ViewD, [ScreenResolution.width, ScreenResolution.height], ScreenSize);
ImSize = round(ImSize);

[XMargin, ~, ~] = ang2pix3(XMargin, ViewD, [ScreenResolution.width, ScreenResolution.height], ScreenSize);
XMargin = round(XMargin(1));

[YMargin, ~, ~] = ang2pix3(YMargin, ViewD, [ScreenResolution.width, ScreenResolution.height], ScreenSize);
YMargin = round(YMargin(2));

T = structfun(@(x) round(x/ifi), T, 'UniformOutput', false);

%% How much has the image been scaled?

ImSize_orig = [960 540];
ImScaling = ImSize/ImSize_orig;

%% Rects

ImRect = CenterRectOnPoint([0 0 ImSize], xCenter, yCenter);

FixRct = CenterRectOnPoint([0 0 FixSize FixSize], xCenter, yCenter + FixOffset);

%% Initial text

WelcomeTxt = 'Welcome!\n Press a key to continue.';

Screen('TextSize', w, 30);

DrawFormattedText(w, WelcomeTxt, 'center', 'center', White);

Screen('Flip',w);
KbStrokeWait(-1);
Screen('Flip',w);
WaitSecs(0.2);

%% ======== TRIAL LOOP =========

for trial = 1:TotNTrials

    %% Load images and make textures

    ImgFile = sprintf('Scene%g%c', AllTrials.Scene(trial), Views{AllTrials.InitView(trial)});

    ThisSeq = AllTrials.Sequence(trial);
    ThisSeq = ThisSeq{1};

    Txtrs = nan(1, length(ThisSeq));

    for i = 1:length(ThisSeq)
       ThisImg = double(imread(fullfile(StimDir_main, ...
           [ImgFile '_' sprintf('%04d.png',  ThisSeq(i))])));
       ThisImg = ThisImg * 0.5 + 64;
       Txtrs(i) = Screen('MakeTexture', w, ThisImg);
       clear ThisImg
    end

    % Load bounding box

    BBoxes = dlmread(fullfile(StimDir_main, [ImgFile '_BBox.txt']));
    LargestBox = [min(BBoxes(:,1:2)), max(BBoxes(:,3:4))];
    LargestBox = LargestBox * ImScaling; % scale like image
    LargestBox(3:4) = LargestBox(1:2) + LargestBox(3:4); % convert width and height to x-end and y-end
    LargestBox(1) = LargestBox(1) - XMargin;
    LargestBox(2) = LargestBox(2) - YMargin;
    LargestBox(3) = LargestBox(3) + XMargin;
    LargestBox(4) = LargestBox(4) + YMargin;
    LargestBox([1,3]) = LargestBox([1,3]) + ImRect(1);
    LargestBox([2,4]) = LargestBox([2,4]) + ImRect(2);

    %% Load two probes

    ProbeView = Views{AllTrials.InitView(trial)};

    TheseOrients = AllTrials.Orients(trial);
    TheseOrients = TheseOrients{1};

    ProbeTxtrs = nan(1, 2);

    for i = 1:2

        %         if TheseOrients(i) == 0
        %             ProbeFile = fullfile(StimDir_main, sprintf('Scene%g%c_%04d.png', ...
        %                 AllTrials.Scene(trial), ProbeView, AllTrials.FinalView(trial)));
        %         else
        %             ProbeFile = fullfile(StimDir_rot, sprintf('Scene%g%c_%04d_%.1f.png', ...
        %                 AllTrials.Scene(trial), ProbeView, ...
        %                 AllTrials.FinalView(trial), TheseOrients(i)));
        %         end

        if TheseOrients(i)==0
            TheseOrients(i) = TheseOrients(i) + (rand>0.5)-0.5;
        end

        ProbeFile = fullfile(StimDir_rot, sprintf('Scene%g%c_%04d_%.1f.png', ...
            AllTrials.Scene(trial), ProbeView, ...
            AllTrials.FinalView(trial), TheseOrients(i)));

        ThisImg = double(imread(ProbeFile));
        ThisImg = ThisImg * 0.5 + 64;
        ProbeTxtrs(i) = Screen('MakeTexture', w, ThisImg);
        clear ThisImg

    end

    %% FINALLY SHOW STUFF

    % T.Fix = 0.5;
    % T.Seq(1) = 1.0;
    % T.Seq(2:4) = 0.5;
    % T.Seq(5) = 1.0;
    % T.Probe = 0.5;
    % T.ISI = 0.1;

    % Fixation

    for frame = 1:T.Fix

        Screen('FillOval', w, White, FixRct);
        Screen('FrameOval', w, Black, FixRct);
        Screen('Flip', w);

    end

    % Main sequence:

    for stim = 1:length(ThisSeq)

        if stim < 5
            StimFrames = T.Seq(stim);
        else
            StimFrames = AllTrials.Jitter(trial);
        end

        for frame = 1:StimFrames

            Screen('DrawTexture', w, Txtrs(stim), [], ImRect);
            if stim > 2
                Screen('FillRect', w,  Grey, LargestBox);
            end
            Screen('FillOval', w, White, FixRct);
            Screen('FrameOval', w, Black, FixRct);
            Screen('Flip', w);

        end

    end

    % Show the two probes:

    for pr = 1:2

        for frame = 1:T.Probe

            Screen('DrawTexture', w, ProbeTxtrs(pr), [], ImRect);
            Screen('FillOval', w, White, FixRct);
            Screen('FrameOval', w, Black, FixRct);
            Screen('Flip', w);

        end

        if pr==1

            for frame = 1:T.ISI

                Screen('FillOval', w, White, FixRct);
                Screen('FrameOval', w, Black, FixRct);
                Screen('Flip', w);

            end

        end

    end

    Screen('Flip', w);
    WaitSecs(0.05);

   %% Response

    while KbCheck(-1); end

    ResponseGiven = 0;
    Count = 1;
    t0 = GetSecs;

    while ~ResponseGiven

        DrawFormattedText(w, 'CW or CCW?', 'center', 'center', White);

        Screen('Flip', w);

        Count = Count + 1;

        RestrictKeysForKbCheck(UsedKeys);

        [keyIsDown, timeSecs, keyCode] = KbCheck(-1);

        if keyIsDown && ~ResponseGiven

            RespKey = KbName(keyCode);
            AllTrials.RT(trial) = timeSecs - t0;

            if keyCode(escapeKey)

                Priority(0);

                ShowCursor;
                Screen('CloseAll');

                return

            elseif keyCode(KeyL)||keyCode(KeyR)

                if keyCode(KeyL)
                    Response = 1;
                elseif keyCode(KeyR)
                    Response = 0;
                end

                AllTrials.Hit(trial) = Response==(TheseOrients(1)<TheseOrients(2));

                ResponseGiven = 1;
                RestrictKeysForKbCheck([]);

            end
        end
    end

    %% Feedback

    if GiveFeedback

       if AllTrials.Hit(trial)
           FBColor = [0 255 0];
       else
           FBColor = [255 0 0];
       end

       for frame = 1:T.FB

           Screen('FillOval', w, FBColor, FixRct);
           Screen('Flip', w);

       end

    end

    Screen('Close', [Txtrs ProbeTxtrs]);

    %% Choose difference intensity (and orientations) for next trial

    ThisInt = abs(AllTrials.Diff(trial));

    if AllTrials.Hit(trial) == 0 && (ThisInt + Step_Size <= Max_Diff)
        Int_next = ThisInt + Step_Size;
    elseif trial > NDown && sum(AllTrials.Hit(trial - (NDown-1):trial))==NDown && ...
            (ThisInt - Step_Size >= Min_Diff)
        Int_next = ThisInt - Step_Size;
    else
        Int_next = ThisInt;
    end

    if trial ~= TotNTrials

        if rand > 0.5
            AllTrials.Diff(trial + 1) = Int_next;
        else
            AllTrials.Diff(trial + 1) = -Int_next;
        end

        % Choose orientations:
        AllTrials.Orients{trial + 1}(1) = round(normrnd(0, 1));
        AllTrials.Orients{trial + 1}(2) = ...
        max(min(AllTrials.Orients{trial + 1}(1) + ...
        AllTrials.Diff(trial + 1), OrientLims(2)), OrientLims(1));

    end

end

%% END

MeanOrient = mean(abs(AllTrials.Diff(6:10)));

Endtxt = sprintf(['You were able to see\na difference of:\n', ...
    '%.2f degrees!\n\nReady for the experiment?'], MeanOrient);

DrawFormattedText(w, Endtxt, 'center', 'center', White);

Screen('Flip', w);

KbStrokeWait(-1);

Priority(0);
ShowCursor;
sca;
