function RSceneTask_RunDemo(SubjNo, RunInfo, DataDir)

%
%
%
%

%% Startup

global Environment
global RealRun

RunType = 'Demo';

RunNo = RunInfo(1);
NRuns = RunInfo(2);
ThisRunNo = RunInfo(3);

% Load settings
RST_Settings_General;
RST_Settings_Test;

%% Set diary and other files

LogDir = fullfile(RunDir,'Logs');
if ~exist(LogDir, 'dir')
    mkdir(LogDir);
end

RunStart = datestr(now, 'dd-mm-yyyy_HH-MM-SS');
diary(fullfile(LogDir,sprintf('Subj%02d_%s_%s.txt', SubjNo, RunType, RunStart)));

DataFile = fullfile(RunDir, sprintf('Subj%02d_%s.mat', SubjNo, RunType));

BUpDir = fullfile(RunDir, 'Backup');
if ~exist(BUpDir,'dir')
    mkdir(BUpDir);
end
BUpFile = fullfile(BUpDir, sprintf('Subj%02d_%s_%s.mat', SubjNo, RunType, RunStart));

%% Get and report flip interval, and compensate waiting times in settings

ifi = Screen('GetFlipInterval', w);
fprintf('\nFlip interval: %2.4g ms (should be %2.4g ms).\n\n', ifi*1000, (1000/RefRate));

if ifi > (1/RefRate)*1.5       % if flip interval is way too slow do not adapt timings
    warning( 'This is a serious timing problem! \n');
end

%% Message
% while loading images and waiting for scanner trigger

InstrTxt = 'Practice run\n\nDoes the object turn\nclockwise or counterclockwise?';

DrawFormattedText(w, InstrTxt, 'center', 'center', White);

Screen('Flip', w);

%% Experimental variables

RunTrials = 20;

Views = {'A', 'B'};

%% Staircase settings

% Here the maximum difference is set at 20 (very reasonable)
% note that if you want to set it bigger you need to make sure that
% it doesn't choose non-existing images (>20 or <-20).

Diff_Start = 5;
Max_Diff = 20;
Min_Diff = 0.5;
Step_Size = 0.5;
NDown = 2;
OrientLims = [-10, 10];

%% Create stupid trials table:

Variables = {'Scene', 'InitView', 'FinalView', 'Diff', ...
    'Orients', 'Hit', 'RT'};

% Automatically if trials are evenly split between initial views, and then
% between final views, half will be wide and half narrow.

AllTrials = nan(RunTrials, numel(Variables));

AllTrials(:, 2) = randi(2, RunTrials, 1);
AllTrials(:, 3) = ones(RunTrials, 1) * 30 + (rand(RunTrials, 1)>0.5) * 60;

AllTrials = array2table(AllTrials, 'VariableNames', Variables);

AllTrials.Scene = randi(NScenes, RunTrials, 1);

AllTrials = AllTrials(randperm(RunTrials), :);

% Add sequence:

Seq = nan(RunTrials, 5);

Seq(:, 1) = 0;
for i = 1:RunTrials
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

% Random jitter before target onset:
AllJitter = rand(RunTrials, 1) * 0.5;
AllTrials.Jitter = Frames.Seq(5) - round(AllJitter/ifi);

% Determine fixation time taking jitter into account
AllTrials.Fix = Frames.Fix - round((0.5 - AllJitter)/ifi);


%% Fixation

Screen('FillOval', w, White, FixRct);
Screen('FrameOval', w, Black, FixRct);
Screen('Flip', w);

fprintf('Waiting for the task to begin in %g seconds...\n', round(T.Delay));
WaitSecs(T.Delay);

% Check loading time to subtract from ITI
tic;

%% LOAD ARROWS FOR ILLUSTRATION

[Im_CW, ~, Alpha_CW] = imread('Arrow_CW_black.png');
[Im_CCW, ~, Alpha_CCW] = imread('Arrow_CCW_black.png');

ArrowTxt(1) = Screen('MakeTexture', w, cat(3, Im_CW, Alpha_CW));
ArrowTxt(2) = Screen('MakeTexture', w, cat(3, Im_CCW, Alpha_CCW));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRIAL LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for trial = 1:RunTrials

    fprintf('--- TRIAL %g/%g ---\n', trial, RunTrials);

    %% Prepare for stimulus presentation

    if RealRun
        BitsiBB.clearResponses();
    end

    % Start ITI
    Screen('FillOval', w, White, FixRct);
    Screen('FrameOval', w, Black, FixRct);
    TStamp.event(trial, 1) = Screen('Flip', w);

    %% Load images, create textures

    ThisView = Views{AllTrials.InitView(trial)};
    ThisSeq = AllTrials.Sequence{trial};

    ImgFile = sprintf('Scene%g%c', AllTrials.Scene(trial), ThisView);

    Txtrs = nan(1, length(ThisSeq));

    for i = 1:length(ThisSeq)
        ThisImg = double(imread(fullfile(StimDir_main, [ImgFile '_' sprintf('%04d.png', ThisSeq(i))])));
        ThisImg = ThisImg * 0.5 + 64;
        Txtrs(i) = Screen('MakeTexture', w, ThisImg);
        clear ThisImg
    end

    % Load bounding box

    BBoxes = dlmread(fullfile(StimDir_main,[ImgFile '_BBox.txt']));
    %BBoxes = BBoxes(Sequence/5+1,:);
    LargestBox = [min(BBoxes(:,1:2)), max(BBoxes(:,3:4))];
    LargestBox(3:4) = LargestBox(1:2) + LargestBox(3:4); % convert width and height to x-end and y-end
    LargestBox(1) = LargestBox(1) - XMargin;
    LargestBox(2) = LargestBox(2) - YMargin;
    LargestBox(3) = LargestBox(3) + XMargin;
    LargestBox(4) = LargestBox(4) + YMargin;
    LargestBox([1,3]) = LargestBox([1,3]) + ImRect(1);
    LargestBox([2,4]) = LargestBox([2,4]) + ImRect(2);

    %% Load two probes

    if AllTrials.Consistent(trial)
        ProbeView = Views{AllTrials.InitView(trial)};
    else
        ProbeView = Views{3 - AllTrials.InitView(trial)};
    end

    TheseOrients = AllTrials.Orients(trial);
    TheseOrients = TheseOrients{1};

    ProbeTxtrs = nan(1, 2);

    for i = 1:2

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

    %% ACTUALLY SHOW STUFF

    LoadDelay = round(toc/ifi);

    if LoadDelay < AllTrials.Fix(trial)

        for frame = 1:(AllTrials.Fix(trial) - LoadDelay)

            Screen('FillOval', w, White, FixRct);
            Screen('FrameOval', w, Black, FixRct);
            Screen('Flip', w);

        end
    end

    for stim = 1:length(ThisSeq)

        if stim < 5
            StimFrames = Frames.Seq(stim);
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
            tnow = Screen('Flip', w);

            if frame == 1
                TStamp.event(trial, stim + 1) = tnow;
            end

        end

    end

    % Show the two probes:

    for pr = 1:2

        for frame = 1:Frames.Probe

            Screen('DrawTexture', w, ProbeTxtrs(pr), [], ImRect);
            Screen('FillOval', w, White, FixRct);
            Screen('FrameOval', w, Black, FixRct);
            tnow = Screen('Flip', w);

            if frame == 1
                if pr == 1
                    TStamp.event(trial, 7) = tnow;
                elseif pr == 2
                    TStamp.event(trial, 9) = tnow;
                end
            end

        end

        if pr == 1

            for frame = 1:Frames.ISI

                Screen('FillOval', w, White, FixRct);
                Screen('FrameOval', w, Black, FixRct);
                tnow = Screen('Flip', w);

                if frame == 1
                    TStamp.event(trial, 8) = tnow;
                end

            end

        end

    end

    TStamp.event(trial, 10) = Screen('Flip', w);
    WaitSecs(T.PreRespDelay);

    %% Get response

    ButPres = 0;

    for frame = 1:Frames.Resp

        DrawFormattedText(w, 'CW or CCW?', 'center', 'center', White);
        sgtext(w, 'Clockwise', 
        tnow = Screen('Flip', w);

        if frame == 1
            TStamp.event(trial, 11) = tnow;
        end

        % Response from button box:

        if RealRun

            [wkey, timeStamp] = BitsiBB.getResponse(timeout, true);

            if ismember(wkey, RespKeys) && ~ButPres
                ButPres = 1;
                Response = find(wkey==RespKeys) - 1; % 0 for CW, 1 for CCW
                AllTrials.Hit(trial) = Response == (TheseOrients(1)<TheseOrients(2));
                AllTrials.RT(trial) = timeStamp - TStamp.event(trial, 11);
                TStamp.response(trial) = timeStamp;
            end

        end

        % Retrieve responses from keyboard (to abort run):

        if IsOSX
           [keydown, resptime, keyCode] = KbCheck(-1);
        else
           [keydown, resptime, keyCode] = KbCheck;
        end

        if keydown  && keyCode(EscKey)
            save(DataFile, 'AllTrials', 'TStamp');
            fprintf('\n\nExperiment terminated at %s, diary closed...\n', datestr(now));
            if RealRun
                close(BitsiScanner);
                fprintf('All serial ports closed.\n');
            end
            close(BitsiBB);
            delete(instrfind);
            ShowCursor;
            diary off
            Priority(0);
            sca; return
        elseif keydown && ~RealRun && any(keyCode(RespKeys)) && ~ButPres
            ButPres = 1;
            Response = find(keyCode(RespKeys)) - 1;
            AllTrials.Hit(trial) = Response == (TheseOrients(1)<TheseOrients(2));
            AllTrials.RT(trial) = resptime - TStamp.event(trial, 11);
        end % end of kbcheck

    end

    %% Feedback (if trial-by-trial)

    if strcmp(GiveFB, 'Trial')

        if AllTrials.Hit(trial) == 1
            FBColor = [0 255 0];
        elseif AllTrials.Hit(trial) == 0
            FBColor = [255 0 0];
        else % if response not given
            FBColor = Black;
        end

        for frame = 1:Frames.FB

            Screen('FillOval', w, FBColor, FixRct);
            Screen('Flip', w);

        end

    else

    end

    %% Save files & close textures

    tic;
    save(DataFile, 'AllTrials', 'TStamp');
    save(BUpFile); % save everything
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

    if trial ~= RunTrials

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

    %% END TRIAL LOOP
end

%% SAVE FILES

save(DataFile, 'AllTrials', 'TStamp');
save(BUpFile); % save everything

%% END OF RUN FEEDBACK

MeanOrient = mean(abs(AllTrials.Diff(trial-RunTrials+1:trial)));

EndTxt = sprintf(['End of run %01d/%01d.\n\nYou were able to', ...
    ' see\na difference of:\n%.2f degrees.\n\nGood job!'], ...
    RunNo, NRuns, MeanOrient);

DrawFormattedText(w, EndTxt, 'center', 'center', White);
Screen('Flip', w);

%% CLEAN UP

close(BitsiBB);
if RealRun
    close(BitsiScanner);
    wmsg = 'scanner and response box';
else
    wmsg = 'virtual response box';
end
delete( instrfind);
fprintf( '\nSerial ports (%s) CLOSED at %s. \n\n', wmsg, datestr( now));

% clean log
fprintf('Diary closed (%s)\n\n', datestr( now));
diary off

% clean memory
try %#ok<TRYNC>
    jheapcl;
end

warning('Press space to continue...')
waitforspace; waitfornokey;

Priority(0);
sca
ShowCursor;
