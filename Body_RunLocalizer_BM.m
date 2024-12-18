function Body_RunLocalizer_BM(SubNo, RunInfo, DataDir)

%{
localizer
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

RunType = 'Localizer_BM';


try
    
    Settings_General; % open screen from here
    
    %% Instruction
    
    % message (while loading images and waiting for scanner trigger)
    Inst = 'Fixate at the center, \n Press any key as soon as possible if you see a repeated motion pattern';
    DrawFormattedText(w.Number, Inst, 'center', 'center', White);
    Screen('Flip', w.Number);
    
    
    %% Set diary and record the current script
    
    LogDir = fullfile(RunDir, 'Logs');
    if ~exist(LogDir, 'dir')
        mkdir(LogDir);
    end
    
    RunStart = datestr(now, 'dd-mm-yyyy_HH-MM-SS');
    diary(fullfile(LogDir, sprintf('Sub%02d_%s_%g_%s.txt', SubNo, RunType, ThisRunNo, RunStart)));
    
    
    ScriptName = {'Body_RunLocalizer_BM','Settings_General','RUN'};
    myCode = cell(length(ScriptName),1);
    for i = 1:length(ScriptName)
        fid = fopen([ScriptName{i} '.m'],'r'); % 'r': open file for reading
        myCode{i} = fscanf(fid,'%300c');
        fclose('all');
    end
    
    %% Constructing all the trials=========================================
    
    fileDir = './Stimuli/Vanrie_mat/';
    files = listdir(fullfile(fileDir, '*.mat'));
    NPoint = 13;
    NAction = length(files);
    for themov = 1:NAction
        
        load(fullfile(fileDir,files{themov}))
        MOV =  RotationMatrix(0.5*pi,[0 0 1])* RotationMatrix(-pi,[0 1 0])* MOV ; 
        %this rotates the movie to face front: dim1: left-right; dim2:
        %depth; dim 3: up-down
        Motions{themov,1} = MOV;
        NFrame(themov) = size(MOV,2)/NPoint;
        
        % make scramble        
        ScramRange = range(MOV,2); % the range of all dots in each dimension
        randPos(:,:) = ScramRange.*(rand(3,NPoint)-0.5); % give a random position for each point within the range
        ScramMOV = MOV;
        for thepoint = 1:NPoint
            startPos(:,thepoint) = mean(MOV(:, thepoint:NPoint:end), 2); % the mean center location of each point
            ScramMOV(:, thepoint:NPoint:end) = MOV(:, thepoint:NPoint:end) - startPos(:,thepoint) + randPos(:,thepoint); % all move to the center and then to the random position
        end
        
        ScramMOV(3,:) = -ScramMOV(3,:); % 3rd dimension is height
        Motions{themov,2} = ScramMOV;
        
    end


    conditions = [1,2]; % {'Intact'; 'Scrambled'}; 

    NTrialPerBlock = 8;
    NRep           = 12; 
    NBlock         = length(conditions) * NRep;
    
    BlockOrder = randSamp(conditions, NBlock, 'n'); % each miniblock only show one angle

    % construct trial table
    T = table();
    for theblock = 1:NBlock
        theCate = BlockOrder(theblock);
        block.condition(1:NTrialPerBlock,:) = deal(conditions(theCate));
        
        % get the corresponding index of images of this block
     
        block.image       = randSamp(1:NAction, NTrialPerBlock, 'n'); 
        block.action      = files(block.image);
        % rotate
        block.rotate      = pi * (rand(NTrialPerBlock,1) - 0.5); % in rad, ranging from left 90 degree to right 90
        
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
        
        block.image(targetId)    = block.image(targetId - 1);
        block.rotate(targetId)   = block.rotate(targetId - 1);
        block.iftarget           = false(NTrialPerBlock,1);
        block.iftarget(targetId) = 1;
        
        T = [T; struct2table(block)];
    end
    
    T.correctKey(T.iftarget==1) = 1;
    T.resp = zeros(height(T),1);
    T.rt   = NaN(height(T),1);
    T.subjectid = SubNo * ones(height(T),1);
    T = movevars(T, 'subjectid','Before', 1);
    
    %% experimental parameters ============================================
    % temporal parameters
    dur.image  = 1000;
    dur.wait   = 200;
    
    Nfr = structfun(@(x) round(x/1000 * w.RefreshRate), dur,'UniformOutput',0); % need to
    SwitchFr = cumsum(cell2mat(struct2cell(Nfr)));
    dur.trial = sum(cell2mat(struct2cell(dur)));
    dur.minimumRT = 0.25;
    
    time.task   = height(T)*dur.trial/1000/60;
    time.rest   = 2500;
    Nfr.rest    = round(time.rest/1000 * w.RefreshRate);
    time.theory = (height(T)*dur.trial + time.before *2 + NBlock * time. rest)/1000/60;
    
    % spatial parameters
    
    img.HAng  = 7.125;
    img.Scale = visAng2pix(img.HAng, 0, monitor)/ScramRange(3);
    Motions   = cellfun(@(x) x*img.Scale, Motions,'UniformOutput', false);
    
    if ~RealRun
        T = T(1: 4 * NTrialPerBlock,:);
    end
    NTrial = height(T);
    

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
        movId = T.image(thetrial);
        theMOV = RotationMatrix(T.rotate(thetrial), [0,0,1]) * Motions{movId, T.condition(thetrial)};

        ButPres     =   0;
        
        if RealRun BitsiBB.clearResponses(); end
        if rem(thetrial, NTrialPerBlock) == 1, fprintf('--- BLOCK %g/%g ---\n', theblock, NBlock); end
        
        for nf = 1: SwitchFr(end)
            if nf <= SwitchFr(1)
                % draw fixation all the time
                thefr = ceil(nf/2); % the action data is at 30 hz
                while thefr > NFrame(movId); thefr = thefr - NFrame(movId); end % loop through the sequence
                
                if ButPres
                    DrawFormattedText(w.Number, '+', 'center', 'center', Black);
                else
                    DrawFormattedText(w.Number, '+', 'center', 'center', fix.Col);
                end
                Screen('DrawDots',w.Number, theMOV( [1,3], (thefr-1)*NPoint+1:thefr*NPoint), fix.Size, White, w.Center, 1);
            elseif nf > SwitchFr(1) && nf <= SwitchFr(2)
                if ButPres
                    DrawFormattedText(w.Number, '+', 'center', 'center', Black);
                else
                    DrawFormattedText(w.Number, '+', 'center', 'center', fix.Col);
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
                Screen('DrawDots', w.Number, [0; 0], fix.Size, Black, w.Center, 1);
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
        % clean log
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