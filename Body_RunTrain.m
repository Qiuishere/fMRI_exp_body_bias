function Body_RunTrain(SubNo, RunInfo, DataDir)

%{
training session: isolated arms
%}
tic;

if ~exist( 'SubNo', 'var') % this will run if the "function" line is commented

    addpath(genpath('Functions'));
    StartUp;
    SubNo     = 00;
    RunNo     = 1;
    NRuns     = 13;
    ThisRunNo = 1;
    DataDir = 'Data_Demo';
    Environment = 1;
    RealRun = 0;
else
    RunNo = RunInfo(1);
    NRuns = RunInfo(2);
    ThisRunNo = RunInfo(3);
    global Environment
    global RealRun
end

RunType = 'Train';


try % open screen from here
    
    Settings_General;
    
    %% Instruction
    
    % message (while loading images and waiting for scanner trigger)
    Inst = 'Fixate at the center, \n Press any key as soon as possible if you see a blue arm';
    DrawFormattedText(w.Number, Inst, 'center', 'center', White);
    Screen('Flip', w.Number);
    
    
    %% Set diary and record the current script
    
    LogDir = fullfile(RunDir, 'Logs');
    if ~exist(LogDir, 'dir')
        mkdir(LogDir);
    end
    
    RunStart = datestr(now, 'dd-mm-yyyy_HH-MM-SS');
    diary(fullfile(LogDir, sprintf('Sub%02d_%s_%g_%s.txt', SubNo, RunType, ThisRunNo, RunStart)));
    
    
    ScriptName = {'Body_RunTrain','Settings_General', 'RUN'};
    myCode = cell(length(ScriptName),1);
    for i = 1:length(ScriptName)
        fid = fopen([ScriptName{i} '.m'],'r'); % 'r': open file for reading
        myCode{i} = fscanf(fid,'%300c');
        fclose('all');
    end
    
    %% Constructing all the trials=========================================
    
    changeAng = 6; % angle will change by 5 degrees. Only one step size as opposed to 3,6,9 in the behavior
    spatialAngs = [-45-changeAng; -45+changeAng; 45-changeAng; 45+changeAng]; % training angles, spatial angle on screen, neg:left, pos:right
    views  = {'L', 'R'}';
    directions = {'Front', 'Back'}';
    figures = {'Fe', 'Ma'}';
    
    NTrialPerBlock = 2 * 6; % should be a multiple of 2 (2 combinations)
    NBlock = length(spatialAngs) * 7;
    BlockOrder = randSamp(1: length(spatialAngs), NBlock, 'n'); % each miniblock only show one angle
    
    % a list for all the possible combinations for a particular angle
    list{1} = {'L', 'Front', 45+changeAng; 'R', 'Back', 45+changeAng;}; % for each spatial angle in a, there are two combinitions of view and direction
    list{2} = {'L', 'Front', 45-changeAng; 'R', 'Back', 45-changeAng;};
    list{3} = {'R', 'Front', 45+changeAng; 'L', 'Back', 45+changeAng;};
    list{4} = {'R', 'Front', 45-changeAng; 'L', 'Back', 45-changeAng;};
    
    
    % construct trial table
    T = table();
    for theblock = 1:NBlock
        theSpatialAng = BlockOrder(theblock);
        
        % assign the two combinations to a mini table of each spatial angle
        combinations = randSamp(1:2, NTrialPerBlock, 'n'); % within one mini block, randomize two possible combinations
        for thetrial = 1:NTrialPerBlock
            block.retinalAng(thetrial,1) = spatialAngs(theSpatialAng);
            
            thecombi = combinations(thetrial);
            block.view(thetrial,1)       = list{theSpatialAng}(thecombi,1);
            block.direction(thetrial,1)  = list{theSpatialAng}(thecombi,2);
            block.angle(thetrial,1)      = list{theSpatialAng}{thecombi, 3};
        end
        
        % randomly assign two figures to each combination
        for thecombi = 1: 2
            id = combinations==thecombi;
            thefigure = randSamp((1:length(figures))', length(find(id)), 'n');
            block.figure(id,1) = figures(thefigure);
        end
        
        % add target trial
        NTargetPerBlock = randperm(2,1); % each block has either one or two targets
        % make sure no two consequtive targets
        TargetRange = 2: NTrialPerBlock;
        targetId = [];
        for j = 1: NTargetPerBlock
            targetId(j) = randSamp(TargetRange, 1, 'n'); % so that the first and last trial won't be target
            toRemove = targetId(j) - 1: targetId(j) + 1; % the one before or follow the target can't be target
            [~, toRemoveId] = ismember(toRemove, TargetRange);
            TargetRange(toRemoveId((toRemoveId~=0))) = [];
        end
        
        % add target trial
        block.iftarget = zeros(NTrialPerBlock,1);
        block.iftarget(targetId) = 1;
        
        
        T = [T; struct2table(block)];
    end
    
    
    % make image textures for all images. This takes some time.
    T.image = strcat('Stimuli/Arms/Above-', T.direction, '/', T.figure, '_', T.direction, '_', num2str(T.angle), T.view, '.png');
    for i = 1:height(T)
        [X1, ~, alpha1] = imread(T.image{i});
        if T.iftarget(i)==1
            X1(:,:,1:2) = X1(:,:,1:2)/3*2; % reduce R G channel to make it blue
        end
        X1(:,:,4)  = alpha1;
        T.img1Txt(i)  =   Screen( 'MakeTexture', w.Number, X1);
    end
    
    T.correctKey(T.iftarget==1) = 1;
    T.resp = zeros(height(T),1);
    T.rt   = NaN(height(T),1);
    T.subjectid = SubNo * ones(height(T),1);
    T = movevars(T, 'subjectid','Before', 1);
    
    %% experimental parameters ============================================
    % temporal parameters
    dur.image  = 450;
    dur.wait   = 200;
    
    Nfr = structfun(@(x) round(x/1000 * w.RefreshRate), dur,'UniformOutput',0); % need to
    SwitchFr = cumsum(cell2mat(struct2cell(Nfr)));
    dur.trial = sum(cell2mat(struct2cell(dur)));
    dur.minimumRT = 0.25;
    
    time.task = height(T)*dur.trial/1000/60;
    time.rest = 2500;
    Nfr.rest = round(time.rest/1000 * w.RefreshRate);
    time.theory = (height(T)*dur.trial + time.before *2 + NBlock * time. rest)/1000/60;
    
    if ~RealRun
        T = T(1: 4 * NTrialPerBlock,:);
    end
    NTrial = height(T);
    
    % spatial parameters
    img.W     = size(X1,2);
    img.H     = size(X1,1);
    img.WAng  = 4;
    img.Scale = visAng2pix(img.WAng, 0, monitor)/img.W;
    img.Pos   = CenterRectOnPoint([0, 0, img.W * img.Scale, img.H * img.Scale], w.Center(1) , w.Center(2));
    
    %% Scanner trigger ====================================================
    toc;
    if RealRun % if dummy scanner or lab, get trigger to start the run
        warning('Waiting for scanner trigger...');
        BitsiScanner.clearResponses();
        firstScan = 0;
        
        while firstScan == 0
            while BitsiScanner.numberOfResponses() == 0  % what does this do??
                WaitSecs(0.001);
            end
            [resp] = BitsiScanner.getResponse(0.001, true);
            if resp == 97
                firstScan = 1;
            end
        end
    else
        warning('Sham run: press Spacebar to continue.');
        waitforspace; waitfornokey;
    end
    
    % !!  GET THE TRIGGER TIME
    time.trigger  = GetSecs;
    Screen('Flip', w.Number);
    
    % Start (sham) scanning
    fprintf('\n SCANNER TRIGGER (run %g): %s \n', RunNo, datestr(now,'dd-mm-yyyy HH:MM:SS'));
    
    
    %% Delay period at the beginning
    
    WaitSecs(time.before/1000 - time.countdown); % during this period the text is still shown
    fprintf('Waiting for the task to begin in %g seconds...\n', time.countdown);
    
    for n = 1:time.countdown
        CountDownTxt = sprintf('Starting in %g seconds', time.countdown - n);
        if n < time.countdown
            DrawFormattedText(w.Number, CountDownTxt, 'center', w.Center(2) + 30, White);
        else
            % show fixation at the last second
            Screen('DrawDots', w.Number, [0; 0], fix.Size, fix.Col, w.Center, 1);
        end
        Screen('Flip', w.Number);
        WaitSecs(1);
    end
    
    %% real task ==========================================================
    
    theblock = 1;
    vbl = Screen('Flip', w.Number);
    
    for thetrial = 1: NTrial
        ButPres     =   0;
        if RealRun BitsiBB.clearResponses(); end
        if rem(thetrial, NTrialPerBlock) == 1, fprintf('--- BLOCK %g/%g ---\n', theblock, NBlock); end
        
        for nf = 1: SwitchFr(end)
            if nf <= SwitchFr(1)
                % draw fixation all the time
                Screen('DrawDots', w.Number, [0; 0], fix.Size, fix.Col, w.Center, 1);
                Screen( 'DrawTexture', w.Number, T.img1Txt(thetrial), [], img.Pos);
            elseif nf > SwitchFr(1) && nf <= SwitchFr(2)
                if ButPres
                    Screen('DrawDots', w.Number, [0; 0], fix.Size, Black, w.Center, 1);
                else
                    Screen('DrawDots', w.Number, [0; 0], fix.Size, fix.Col, w.Center, 1);
                end
            end
            
            Screen('DrawingFinished', w.Number);
            vbl = Screen('Flip', w.Number, vbl + .5 * w.ifi);
            
            if nf == 1
                T.img1Onset(thetrial) = vbl - time.trigger;
            end
            
            %% collect response
            if RealRun
                % retreive responses from button box
                [wkey, timeStamp] = BitsiBB.getResponse(timeout, true);
                if ismember( wkey, RespKeys) && ~ButPres % only allow to response with the first key (index figure)
                    ButPres = 1;
                    % if the response is made within 300ms, it can't be for
                    % this trial
                    if timeStamp - T.img1Onset(thetrial) - time.trigger> dur.minimumRT
                        T.resp(thetrial)     =   1;          % 1 no matter which key it is
                        T.rt(thetrial)       =   timeStamp - T.img1Onset(thetrial) - time.trigger;  % record time relative to the previous image onset to avoid confusion
                    elseif isnan(T.rt(thetrial-1))
                        % when the response is made too fast, it probably belongs to the previous trial.
                        % But another case is the press from the last trial
                        % is not released (rt in thetrial - 1 was already recorded). Then we need to avoid override
                        % the previous rt
                        T.resp(thetrial-1)   =   1;          % 1
                        T.rt(thetrial-1)     =   timeStamp - T.img1Onset(thetrial - 1) - time.trigger;
                    end
                end
            end
            
            % check quit
            [keydown,respT,keyCode] = KbCheck;
            if keydown && keyCode(KbName('ESCAPE'))
                fprintf('\n\nExperiment terminated at %s, diary closed...\n', datestr(now));
                if RealRun
                    close(BitsiScanner);
                    fprintf('All serial ports closed.\n');
                end
                close(BitsiBB);
                ShowCursor;
                diary off
                clear temprun
                sca; return
            elseif keydown && any(keyCode(RespKeys)) && ~ButPres && ~RealRun
                ButPres = 1;
                % if the response is made within 300ms, it can't be for
                % this trial
                if respT - T.img1Onset(thetrial) - time.trigger> dur.minimumRT
                    T.resp(thetrial)     =   1;
                    T.rt(thetrial)       =   respT - T.img1Onset(thetrial) - time.trigger;
                elseif isnan(T.rt(thetrial-1)) % when the response is made too fast, it belongs to the previous trial
                    T.resp(thetrial-1)   =   1;
                    T.rt(thetrial-1)     =   respT - T.img1Onset(thetrial - 1) - time.trigger;
                end
            end % end of kbcheck
            
        end
        
        %% Fixation (end of miniblock)
        if rem(thetrial, NTrialPerBlock) == 0
            T.ifcorrect = T.resp == T.correctKey; % no reponse for non-target trials is also correct
            subT = T(thetrial - NTrialPerBlock + 1: thetrial,:);
            targetTrials = subT.ifcorrect(logical((subT.iftarget)));
            rt   = mean(subT.rt(subT.ifcorrect & subT.iftarget)); % omit the unresponded trials
            fprintf('((( Target detected; %g/%g. mean RT = %g s )))\n',length(find(targetTrials)),length(targetTrials), rt); % to print a % inside fprintf, you need %%
            
            for nfr = 1:Nfr.rest
                DrawFormattedText(w.Number, '+', 'center', 'center', Black);
                Screen('DrawingFinished', w.Number);
                vbl = Screen('Flip', w.Number, vbl + .5 * w.ifi);
            end
            
            theblock = theblock + 1;
            
        end
    end
    
    
    % Monitor PTB performance
    
    dur.trialActual = 1000 * diff(T.img1Onset); % get the actual duraction for each trial
    dur.trialActual = dur.trialActual(dur.trialActual> 0 & dur.trialActual< time.rest);
    TrialDur = round(nanmean(dur.trialActual));
    TrialStd = round(nanstd( dur.trialActual));
    fprintf('((( Trial Duration %g/%g ms, SD = %g ms )))\n', TrialDur, dur.trial, TrialStd);
    %     scatter(ones(size(dur.trialActual)), dur.trialActual)
    time.total = toc;
    

    %% EXIT and clean up===================================================
    targetTrials = T.ifcorrect(logical((T.iftarget)));
    EndTxt = sprintf(['Well done! You completed run %g/%g.\n\n', ...
        'You detected %g targets out of %g. \n\n', ...
        'You can take a break while we prepare the next run.'],...
        RunNo, NRuns, length(find(targetTrials)),length(targetTrials));
    DrawFormattedText(w.Number, EndTxt, 'center', 'center', text.Color);
    Screen('Flip', w.Number);
    
    fprintf('Waiting for the task to end in %g seconds...\n', time.before/1000);
    % only after the signal returned to baseline, close recording
    WaitSecs(time.before/1000);
    close(BitsiBB);
    
    
    if RealRun
        close(BitsiScanner);
        wmsg = 'scanner and response box';
    else
        wmsg = 'virtual response box';
    end
    fprintf('\nSerial ports (%s) CLOSED at %s. \n\n', wmsg, datestr(now));
    time.actual  =  round( GetSecs - time.trigger)/60;
    fprintf('RUN DURATION: %g mins (expected duration: %g mins). \n\n', time.actual, time.theory);
    
    
    %% Save 1. all data 2. resulttable=====================================
    
    fileName1 = fullfile(RunDir, sprintf('Sub%02d_%s_%g.mat', SubNo, RunType, ThisRunNo));
    save(fileName1);
    fileName2 = fullfile(RunDir, sprintf('ResultTable_Sub%02d_%s_%g.csv', SubNo, RunType, ThisRunNo));
    writetable(T, fileName2);
    
    fprintf('Diary closed (%s)\n\n', datestr(now));
    diary off
    
    % clean memory
    %  cleans up the Java heap from memory leaks, preventing the infamous
    %  Java OutOfMemory exception. https://nl.mathworks.com/matlabcentral/fileexchange/36757-java-heap-cleaner
    try  % use this to avoid a cessation from error
        jheapcl;
    end
    
    
    sca
    ShowCursor;
    Priority(0);
    toc;
    
catch ME
    sca
    ShowCursor;
    Priority(0);
    toc;
    rethrow(ME);
end