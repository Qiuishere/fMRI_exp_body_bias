function RSceneTask_RunTest(SubjNo, RunInfo, DataDir)

% 
% 
% 
% 

%% Startup

global Environment
global RealRun

RunType = 'Training';

RunNo = RunInfo(1);
NRuns = RunInfo(2);
ThisRunNo = RunInfo(3);

% Load settings
RST_Settings_General;
RST_Settings_Test;

%% Set diary and other files

LogDir = fullfile(RunDir,'Logs');
if ~exist(LogDir,'dir')
    mkdir(LogDir);
end

RunStart = datestr(now,'dd-mm-yyyy_HH-MM-SS');
diary(fullfile(LogDir,sprintf('Subj%02d_%s_%d_%s.txt', SubjNo, RunType, ThisRunNo, RunStart)));

DataFile = fullfile(RunDir, sprintf('Subj%02d_%s_%g.mat', SubjNo, RunType, ThisRunNo));

BUpDir = fullfile(RunDir,'Backup');
if ~exist(BUpDir,'dir')
    mkdir(BUpDir);
end
BUpFile = fullfile(BUpDir, sprintf('Subj%02d_%s_%g_%s.mat', SubjNo, RunType, ThisRunNo, RunStart));

%% Get and report flip interval, and compensate waiting times in settings

ifi = Screen('GetFlipInterval', w);
fprintf('\nFlip interval: %2.4g ms (should be %2.4g ms).\n\n', ifi*1000, (1000/RefRate));

if ifi > (1/RefRate)*1.5       % if flip interval is way too slow do not adapt timings
    warning( 'This is a serious timing problem! \n');
end

%% Message
% while loading images and waiting for scanner trigger

InstrTxt = sprintf(['Run %g/%g\n\nDoes the object turn\nclockwise or counterclockwise?'], RunNo, NRuns);

DrawFormattedText(w, InstrTxt, 'center', 'center', White);

Screen('Flip', w);

%% Experimental variables

RunTrials = 32;

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

%% Load or create trials table:

if exist(DataFile, 'file')
    load(DataFile);
else
    AllTrials = RSceneTask_CreateTestRunTable;
    save(DataFile, 'AllTrials');
end

TotNTrials = size(AllTrials, 1);

%% Assign initial difference (if not already done)

if ThisRunNo == 1
    
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
    AllJitter = rand(TotNTrials, 1) * 0.5;
    AllTrials.Jitter = Frames.Final - round(AllJitter/ifi);
    
    % Determine fixation time taking jitter into account
    AllTrials.Fix = Frames.Fix - round((0.5 - AllJitter)/ifi);
    
end

%% Create structure to save timestamps

TStamp = struct;

TStamp.trigger = NaN;
TStamp.event = nan(RunTrials, length(TrialSequence) + 1); % the '+ 1' is the response time

%% Scanner trigger

if RealRun

    warning('Waiting for scanner trigger...');

    BitsiScanner.clearResponses();
    firstScan = 0;

    while firstScan == 0
        while BitsiScanner.numberOfResponses() == 0
            WaitSecs(0.001);
        end
        [resp] = BitsiScanner.getResponse(0.001, true);
        if resp == 97
            firstScan = 1;
        end
    end

else
    warning('Sham run: press spacebar to continue.');
    waitforspace; waitfornokey;
end

% start (sham) scanning
TStamp.trigger = GetSecs;

fprintf('\n SCANNER TRIGGER (test run %g): %s \n', ThisRunNo, datestr(now,'dd-mm-yyyy HH:MM:SS'));

%% Fixation

Screen('FillOval', w, White, FixRct);
Screen('FrameOval', w, Black, FixRct);
Screen('Flip', w);

fprintf('Waiting for the task to begin in %g seconds...\n', round(T.Delay));
WaitSecs(T.Delay);

% Check loading time to subtract from ITI
tic;

%% COUNTDOWN?

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRIAL LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FirstTrial = (ThisRunNo - 1) * RunTrials + 1;

for trial = FirstTrial:FirstTrial + RunTrials - 1
    
    fprintf('--- TRIAL %g/%g ---\n', trial, RunTrials);
    
    %% Prepare for stimulus presentation
    
    if RealRun
        BitsiBB.clearResponses();
    end
    
    % Start ITI
    Screen('FillOval', w, White, FixRct);
    Screen('FrameOval', w, Black, FixRct);
    TStamp.event(trial,1) = Screen('Flip', w);
    
    %% Load images, create textures
    
    ThisView = Views{AllTrials.InitView(trial)};
    ThisSeq = AllTrials.Sequence{trial};
    
    ImgFile = sprintf('Scene%g%c', AllTrials.Scene(trial), ThisPos);
    
    Txtrs = nan(1, length(ThisSeq));
    
    for i = 1:length(ThisSeq)
        ThisImg = double(imread(fullfile(StimDir, [ImgFile '_' sprintf('%04d.png',Sequence(i))])));
        ThisImg = ThisImg * 0.5 + 64;
        Txtrs(i) = Screen('MakeTexture', w, ThisImg);
        clear ThisImg
    end
    
    % Load bounding box
    
    BBoxes = dlmread(fullfile(StimDir,[ImgFile '_BBox.txt']));
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
        
        for frame = 1:Frames.Probe
            
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
    WaitSecs(T.PreRespDelay);
      
    %% Get response
    
    ButPres = 0;
    
    for frame = 1:Frames.Resp
    
        Screen('FillOval', w, White, FixRct);
        Screen('FrameOval', w, Black, FixRct);
        Screen('Flip', w);
        
        if RealRun
            
            [wkey, timeStamp] = BitsiBB.getResponse(timeout, true);
            
            if ismember(wkey, RespKeys) && ~ButPres
                ButPres = 1;
                AllTrials.
    
end

