function RSceneTask_RunTraining(SubjNo, RunInfo, DataDir)

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

RunTrials = 360;

% Load settings
RST_Settings_General;
RST_Settings_Training;

%% Set diary and other files

LogDir = fullfile(RunDir, 'Logs');
if ~exist(LogDir, 'dir')
    mkdir(LogDir);
end

RunStart = datestr(now, 'dd-mm-yyyy_HH-MM-SS');
diary(fullfile(LogDir, sprintf('Subj%02d_%s_%g_%s.txt', SubjNo, RunType, ThisRunNo, RunStart)));

DataFile = fullfile(RunDir, sprintf('Subj%02d_%s_%g.mat', SubjNo, RunType, ThisRunNo));
BUpDir = fullfile(RunDir, 'Backup');
if ~exist(BUpDir, 'dir')
    mkdir(BUpDir);
end
BUpFile = fullfile(BUpDir, sprintf('Subj%02d_%s_%g_%s.mat', SubjNo, RunType, ThisRunNo, RunStart));

%% Get and report flip interval, and compensate waiting times in settings

ifi = Screen('GetFlipInterval', w);
fprintf('\nFlip interval: %2.4g ms (should be %2.4g ms).\n\n', ifi*1000, (1000/RefRate));

if ifi > (1/RefRate) * 1.5 % if flip interval is way too slow do not adapt timings
    warning('This is a serious timing problem!\n');
end

%% Message
% while loading images and waiting for scanner trigger

InstrTxt = sprintf(['Run %g/%g\n\n',...
    'Press a button every time the same image appears twice in a row.'], RunNo, NRuns);

DrawFormattedText(w, InstrTxt, 'center', 'center', White);

Screen('Flip', w);

%% Conditions

NStimTypes = 8;

[TrialMatrix, TargetTrials] = RSceneTask_CreateTrainRunMat;

assert(size(TrialMatrix, 1) == NStimTypes && size(TargetTrials, 1) == NStimTypes);

%% Create trial table
% total = 360 trials

% Determine order of stimulus types:

StimTypes = zeros(1, 40); % what stimulus types in what mini-blocks
StillToDo = ones(1, NStimTypes); % stimulus types still to assign

StimTypes(1) = randi(4);
StillToDo(StimTypes(1)) = 0;

for i = 2:40
    
    ChooseFrom = find(StillToDo);
    ChooseFrom = ChooseFrom(randperm(length(ChooseFrom)));

    StimTypes(i) = ChooseFrom(1);
    StillToDo(StimTypes(i)) = 0;
    
    if rem(i, 8) == 0 % reset stimulus types
        StillToDo = ones(1, 8);
    end
    
end

% 'Type' is the stimulus category
% 'Stimulus' is the specific stimulus

Variables = {'Type', 'Stimulus', 'Target', 'Response'};
AllTrials = nan(RunTrials, numel(Variables));
AllTrials = array2table(AllTrials, 'VariableNames', Variables);

for i = 1:RunTrials
    
    if rem(i, 9) == 1 % new mini-block
        whichtype = ceil(i/9);
        AllTrials.Type(i:i+8) = StimTypes(whichtype);
        TypeRep = sum(StimTypes(1:whichtype)==StimTypes(whichtype));
        AllTrials.Stimulus(i:i+8) = TrialMatrix(StimTypes(whichtype), (TypeRep-1)*9+1:(TypeRep-1)*9+9);
        AllTrials.Target(i:i+8) = TargetTrials(StimTypes(whichtype), (TypeRep-1)*9+1:(TypeRep-1)*9+9);
    end
    
end

%% Select images (load them later)    
%
% SHOULD WE SEPARATE BEDS AND COUCHES?
% 
% Stimulus types:
% 1. View A - Scene 30° - Object 30° 
% 2. View A - Scene 30° - Object 90°
% 3. View B - Scene 30° - Object 30° 
% 4. View B - Scene 30° - Object 90°
% 5. View A - Scene 90° - Object 30°
% 6. View A - Scene 90° - Object 90°
% 7. View B - Scene 90° - Object 30°
% 8. View B - Scene 90° - Object 90°
%

StimTypes_list = {'A', 30, 30; ...
                  'A', 30, 90; ...
                  'B', 30, 30; ...
                  'B', 30, 90; ...
                  'A', 90, 30; ...
                  'A', 90, 90; ...
                  'B', 90, 30; ...
                  'B', 90, 90};

assert(size(StimTypes_list, 1) == NStimTypes);
              
AllObjTxtrs = nan(RunTrials, 1);
AllSceneTxtrs = AllObjTxtrs;
ObjFiles = cell(RunTrials, 1);
SceneFiles = ObjFiles;

ObjFileName = 'Object%g%c_%04d.png';
SceneFileName = 'Background%g_%04d.png';

for i = 1:RunTrials

    StartPoint = StimTypes_list{AllTrials.Type(i), 1};
    SceneFrame = StimTypes_list{AllTrials.Type(i), 2};
    ObjFrame = StimTypes_list{AllTrials.Type(i), 3};
        
    ObjFileName_this = sprintf(ObjFileName, AllTrials.Stimulus(i), StartPoint, ObjFrame);
    SceneFileName_this = sprintf(SceneFileName, AllTrials.Stimulus(i), SceneFrame);
    
    ObjFiles{i} = fullfile(StimDir, ObjFileName_this);
    SceneFiles{i} = fullfile(StimDir, SceneFileName_this);
    
end

%% Load images for first group of blocks

thisblock = 1:72; % 8 mini-blocks of 9 stimuli each

for i = thisblock
    
    [ThisObj, ~, Alpha] = imread(ObjFiles{i});
    ThisObj = double(ThisObj) * 0.5 + 64;
    ThisObj = cat(3, ThisObj, Alpha);
    
    ThisScene = imread(SceneFiles{i});
    ThisScene = double(ThisScene) * 0.5 + 64;
    
    AllObjTxtrs(i) = Screen('MakeTexture', w, ThisObj);
    AllSceneTxtrs(i) = Screen('MakeTexture', w, ThisScene);
    
end

%% Create structure to save timestamps

TStamp = struct;

TStamp.trigger = NaN;
TStamp.event = nan(RunTrials, 1);
TStamp.offset = nan(RunTrials, 1);
TStamp.response = nan(RunTrials, 1);
TStamp.empty = nan(5, 1); % for each miniblock

%% Scanner trigger

if RealRun % dummy scanner or lab
        
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

% Start (sham) scanning
TStamp.trigger = GetSecs;

fprintf('\n SCANNER TRIGGER (training run %g): %s \n', ThisRunNo, datestr(now,'dd-mm-yyyy HH:MM:SS'));

%% Fixation
% Should I put a message to subjects that they have to wait?    

Screen('FillOval', w, White, FixRct);
Screen('FrameOval', w, Black, FixRct);
Screen('Flip', w);

fprintf('Waiting for the task to begin in %g seconds...\n', round(T.Delay));
WaitSecs(T.Delay);

%% COUNTDOWN?

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRIAL LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Block = 1; 

for trial = 1:RunTrials
    
    if rem(trial, 72) == 1, fprintf('--- BLOCK %g/%g ---\n', Block, 5); end
    
    %% Prepare for stimulus presentation
    
    if RealRun
        BitsiBB.clearResponses();
    end
    
    %% SHOW STUFF
    
    ButPres = 0;
    
    for frame = 1:(Frames.On + Frames.Off)
        if frame <= Frames.On
            Screen('DrawTexture', w, AllSceneTxtrs(trial), [], ImRect);
            Screen('DrawTexture', w, AllObjTxtrs(trial), [], ImRect);
        end
        Screen('FillOval', w, White, FixRct);
        Screen('FrameOval', w, Black, FixRct);
        tnow = Screen('Flip', w);
        
        % Timestamp
        
        if frame == 1
            TStamp.event(trial) = tnow;
        elseif frame == Frames.On + 1
            TStamp.offset(trial) = tnow;
        end
        
        % Get response:
        
        if RealRun
            
            [wkey, timeStamp] = BitsiBB.getResponse(timeout, true);
            
            if ismember(wkey, RespKeys) && ~ButPres
                ButPres = 1;
                AllTrials.Response(trial) = find(wkey==RespKeys);
                TStamp.response(trial) = timeStamp - TStamp.trigger;
            end
            
        end
        
        % Retrieve responses from keyboard (to abort run):
        if IsOSX
            [keydown, ~, keyCode] = KbCheck(-1);
        else
            [keydown, ~, keyCode] = KbCheck;
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
            AllTrials.Response(trial) = find(find(keyCode)==RespKeys);
        end % end of kbcheck
        
    end
    
    %% End of mini-block
    
    if rem(trial, 72) == 0
        
        %% Fixation (end of miniblock)
        
        Screen('FillOval', w, White, FixRct);
        Screen('FrameOval', w, Black, FixRct);
        TStamp.empty(Block) = Screen('Flip', w);
        tic; % to compensate
        Block = Block + 1;

        %% Monitor participant performance

        fprintf(' -> Hits: %g/%g\n', sum(~isnan(AllTrials.Response(trial-71:trial)) & AllTrials.Target(trial-71:trial)), sum(AllTrials.Target(trial-71:trial)));
        fprintf(' -> FAs: %g/%g\n',  sum(~isnan(AllTrials.Response(trial-71:trial)) & ~AllTrials.Target(trial-71:trial)), sum(~AllTrials.Target(trial-71:trial)));

        %% Monitor PTB performance
        
        EvDur = diff(TStamp.event(trial-71:trial));
        TrialDur = round(1000 * nanmean(EvDur));
        TrialStd = round(1000 * nanstd(EvDur));
        fprintf('((( Trial Duration %g/%g ms, SD = %g ms )))\n', TrialDur, TrialSet, TrialStd);
        
        %% Save
        
        save(DataFile, 'AllTrials', 'TStamp');
        save(BUpFile); % save everything
        
        %% Load images for next block
        
        Screen('Close', AllObjTxtrs(thisblock));
        Screen('Close', AllSceneTxtrs(thisblock));
        
        if trial ~= RunTrials
            
            thisblock = thisblock + 72;
            
            for i = thisblock
                
                [ThisObj, ~, Alpha] = imread(ObjFiles{i});
                ThisObj = double(ThisObj) * 0.5 + 64;
                ThisObj = cat(3, ThisObj, Alpha);
                
                ThisScene = imread(SceneFiles{i});
                ThisScene = double(ThisScene) * 0.5 + 64;
                
                AllObjTxtrs(i) = Screen('MakeTexture', w, ThisObj);
                AllSceneTxtrs(i) = Screen('MakeTexture', w, ThisScene);
                
            end
            
        end
        
        %% 
        
        if toc < T.Fix
            WaitSecs(T.Fix - toc);
        end
        
    end
end
           
%% Fixation time at end of run

Screen('FillOval',w, White, FixRct);
Screen('FrameOval', w, Black, FixRct);
Screen('Flip', w);

fprintf('\n\nWaiting for the run to end in %g seconds...', round(T.Delay));
WaitSecs(T.Delay - T.Fix);

%% Save stuff

save(DataFile, 'AllTrials', 'TStamp');
save(BUpFile); % save everything

%% EXIT MESSAGE

EndTxt = sprintf(['Well done! You completed run %g/%g.\n\n', ...
    'You can take a break while we prepare the next run.'], RunNo, NRuns);
DrawFormattedText(w, EndTxt, 'center', 'center', White);
Screen('Flip', w);

%% MONITOR STIMULUS PRESENTATION

% Run duration
RealDur = round(GetSecs - TStamp.trigger);
fprintf('\nRUN DURATION: %g seconds (expected duration: %g seconds). \n\n', RealDur, CompDur);

%% CLEAN UP

% Screen('Close', AllTxtrs);

close(BitsiBB);
if RealRun
    close(BitsiScanner);
    wmsg = 'scanner and response box';
else
    wmsg = 'virtual response box';
end
delete(instrfind);
fprintf('\nSerial ports (%s) CLOSED at %s. \n\n', wmsg, datestr(now));

% clean log
fprintf('Diary closed (%s)\n\n', datestr(now));
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

end
    