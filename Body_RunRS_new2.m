% function Body_RunRS_new(SubNo, RunInfo, DataDir)

%{
1. set fixation time: whether to remove the period before? or restrict it
as 500ms?
2. check trial number and correctkeyS
%}
tic;

if ~exist( 'SubNo', 'var') % this will run if the "function" line is commented
    
    addpath(genpath('Functions'));
    StartUp;
    SubNo     = 00;
    RunNo     = 1;
    NRuns     = 13;
    ThisRunNo = 1;
    DataDir   = 'Data_Demo';
    Environment = 1;
    RealRun = 0;
    
else
    RunNo = RunInfo(1);
    NRuns = RunInfo(2);
    ThisRunNo = RunInfo(3);
    global Environment
    global RealRun
end
RunType = 'Test';


try % open screen from here
    
    Settings_General;
    upKey   = Keys.name{1}; % order already randomized.
    downKey = Keys.name{2};
    
    %% Instruction
    
    % message (while loading images and waiting for scanner trigger)
    Inst = sprintf('Fixate at the center, \n judge whether the arm in the second posture moves upward or downward.\n Up:     %s \n Down: %s', upKey, downKey);
    DrawFormattedText(w.Number, Inst, 'center', 'center', White);
    Screen('Flip', w.Number);
    
    
    %% Set diary and record the current script
    
    LogDir = fullfile(RunDir, 'Logs');
    if ~exist(LogDir, 'dir')
        mkdir(LogDir);
    end
    
    RunStart = datestr(now, 'dd-mm-yyyy_HH-MM-SS');
    diary(fullfile(LogDir, sprintf('sub-%03d_%s_%g_%s.txt', SubNo, RunType, ThisRunNo, RunStart)));
    
    
    ScriptName = {'Body_RunRS','Settings_General', 'RUN'};
    myCode = cell(length(ScriptName),1);
    for i = 1:length(ScriptName)
        fid = fopen([ScriptName{i} '.m'],'r'); % 'r': open file for reading
        myCode{i} = fscanf(fid,'%300c');
        fclose('all');
    end
    
    %% Constructing all the trials=========================================
    
    changeAng = 6; % angle will change by 6 degrees. Only one step size as opposed to 3,6,9 in the behavior
    
    % startAngs combined with movingDirs, this create 10 pairs of angle combination.
    % When analyzing, should compare 35->40 with 40->35 (not 35->40 and
    % 35->30), to ensure stimuli matches exactly in moving up and moving down
    % conditions
    movingDirs = [-changeAng, changeAng, 90]';
    correctKeys = [1, 2, nan];
    startAngles = [45,45+changeAng; 45,45-changeAng; 45,45]; % row1 for moving up, row 2 for moving down
    views = {'L', 'R'}';
    directions = {'Front', 'Back'}';
    figures = {'Fe', 'Ma'}';
    
    
    minNTrial = length(directions) * length(movingDirs) * length(views) * size(startAngles,2);
    NTrial.RS = minNTrial * 3;
    
    % create a list with all the possible combination of conditions
    clear list
    startAngid     = 1:size(startAngles,1);
    figure2use    = randSamp(1:length(figures), NTrial.RS/length(directions)/length(movingDirs), 'n'); % keep this list identical for key conditions
    list.movingDir= repelem(movingDirs, NTrial.RS/length(movingDirs),1);
    
    [list.direction, list.angle1, list.angle2, list.view] = deal(zeros(size(list.movingDir)));
    %first condi: movingdir
    for themovingdir = 1:length(movingDirs)
        id = list.movingDir == movingDirs(themovingdir);
        [list.correctKey(id,1)] = deal(correctKeys(themovingdir));
        % second condi: directions
        list.direction(id,1) = randSamp(1:length(directions), length(find(id)), 'n');
        
        for thedirection = 1:length(directions)
            id2 = id & list.direction==thedirection;
            list.figure(id2) = figure2use;
            
            % third level: start angle            
            list.angle1(id2) = randSamp(startAngles(themovingdir,:), length(find(id2)), 'n');
            list.angle2(id2) = list.angle1(id2) + movingDirs(themovingdir);
            if themovingdir==3 % large change, move 90 degrees
                list.angle2(id2) = 90 - list.angle1(id2);
            end
            
            for theangle = 1:length(startAngles(themovingdir,:))
                % fourth: views
                id3 = id2 & list.angle1==startAngles(themovingdir,theangle);
                list.view(id3) = randSamp(1:length(views), length(find(id3)), 'n');
            end
        end


    end
    
    % reverse the direction of large posture change
    list.direction2 = list.direction;
    list.direction2(list.movingDir==movingDirs(3)) = 3 - list.direction(list.movingDir==movingDirs(3));
   
    % convert to condition names
    list.view       = views(list.view);
    list.direction  = directions(list.direction);
    list.direction2 = directions(list.direction2);
    list.figure     = figures(list.figure);
    
%% convert to a table
    
    % construct the table for Repetition Suppresion
    T = struct2table(list); 

    % make image and loop through the list and create texture from image
    
    T.image1 = strcat('Stimuli/', T.direction,  '-Upper-shoulder-middle/', T.figure, '_', T.direction, '_', num2str(T.angle1), T.view, '.png');
    T.image2 = strcat('Stimuli/', T.direction2, '-Upper-shoulder-middle/', T.figure, '_', T.direction2, '_', num2str(T.angle2), T.view, '.png');
    
    for i = 1:height(T)
        [X1, ~, alpha1] = imread(T.image1{i});
        [X2, ~, alpha2] = imread(T.image2{i});
        X1(:,:,4)  = alpha1;
        X2(:,:,4)  = alpha2;
        T.img1Txt(i)  =   Screen( 'MakeTexture', w.Number, X1);
        T.img2Txt(i)  =   Screen( 'MakeTexture', w.Number, X2);
    end
    
    
    % also some none trials
    Tnone = T(1,:); % only one row
    [Tnone.angle1, Tnone.angle2, Tnone.movingDir, Tnone.img1Txt, Tnone.img2Txt]  = deal(NaN);
    
    % combine these three trial types together and randomize
    NTrial.none = 24;
    T = [T; repmat(Tnone, NTrial.none/height(Tnone), 1)];
    T = T(randperm(height(T)),:);
    
    T.resp = zeros(height(T),1);
    T.rt   = NaN(height(T),1);
    T.subjectid = SubNo * ones(height(T),1);
    T = movevars(T, 'subjectid','Before', 1);
    
    %% Experimental precedures ============================================
    
    % temporal parameters
    iti   = 2400 + [-600, 0, 600];  % make sure the length is a divisor of Ntrials to ensure identical time length
    if ~mod(height(T),length(iti))
        warning('ITI cannot achieve qual distribution')
    end
        
    T.iti = randSamp(iti, height(T), 'n');
    dur.fix    = 000;
    dur.image1 = 200;
    dur.blank  = 1000;
    dur.image2 = 200;
    dur.wait   = mean(T.iti);
    
    Nfr = structfun(@(x) round(x/1000 * w.RefreshRate), dur,'UniformOutput',0); % need to
    SwitchFr = cumsum(cell2mat(struct2cell(Nfr)));
    
    dur.trial = sum(cell2mat(struct2cell(dur)));
    time.task = height(T)*dur.trial/1000/60;
    time.rest = 12 * 1000; % try to use ms for all time variables
    
    Nfr.rest = round(time.rest/1000 * w.RefreshRate);
    time.theory = (height(T)*dur.trial + time.before *2 + time. rest *2)/1000/60;
    
    if ~RealRun
        T = T(1: 30,:);
    end
    NTrial.total = height(T);
    
    % spatial parameters
    img.W     = size(X1,2);
    img.H     = size(X1,1);
    img.WAng  = 4;
    img.Scale = visAng2pix(img.WAng, 0, monitor)/img.W;
    img.Pos   = CenterRectOnPoint([0, 0, img.W * img.Scale, img.H * img.Scale], w.Center(1) , w.Center(2));
    
    %% Scanner trigger ====================================================
    toc % get the time cost for code initiation before trigger
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
    
    vbl = Screen('Flip', w.Number);
    
    for thetrial = 1: NTrial.total
        
        SwitchFr(end) = SwitchFr(end-1) + T.iti(thetrial)/1000 * w.RefreshRate;
        ButPres     =   0;
        if RealRun BitsiBB.clearResponses(); end
        
        for nf = 1: SwitchFr(end)
            if nf <= SwitchFr(1)
                % draw dixation
                Screen('DrawDots', w.Number, [0; 0], fix.Size, fix.Col, w.Center, 1);
            elseif nf > SwitchFr(1) && nf <= SwitchFr(2)
                if ~isnan(T.angle1(thetrial))  % RS trial, draw image 1
                    Screen( 'DrawTexture', w.Number, T.img1Txt(thetrial), [], img.Pos);
                end
                Screen('DrawDots', w.Number, [0; 0], fix.Size, fix.Col, w.Center, 1);
            elseif nf > SwitchFr(2) && nf <= SwitchFr(3)
                % fixation
                Screen('DrawDots', w.Number, [0; 0], fix.Size, fix.Col, w.Center, 1);
            elseif nf > SwitchFr(3) && nf <= SwitchFr(4)
                if ~isnan(T.angle2(thetrial)) % RS trial
                    Screen('DrawTexture', w.Number, T.img2Txt(thetrial), [], img.Pos);
                end
                Screen('DrawDots', w.Number, [0; 0], fix.Size, fix.Col, w.Center, 1);
            else % the response stage
                if ButPres
                    Screen('DrawDots', w.Number, [0; 0], fix.Size, fix.Col, w.Center, 1);
                elseif ~isnan(T.angle2(thetrial)) % only show the prompt for RS trials, otherwise blank
                    prompt = sprintf('Please select:\n Up:     %s \n Down: %s', upKey, downKey);
                    DrawFormattedText(w.Number, prompt, 'center', 'center', text.Color);
                elseif nf -  SwitchFr(end - 1) > T.iti(thetrial)/2/1000 * w.RefreshRate
                    % for trials where no response is needed, show fixation for the last half of the iti
                    Screen('DrawDots', w.Number, [0; 0], fix.Size, fix.Col, w.Center, 1);
                end
            end
            
            Screen('DrawingFinished', w.Number);
            vbl = Screen('Flip', w.Number, vbl + .5 * w.ifi);
            
            if nf == SwitchFr(1) + 1
                T.img1Onset(thetrial) = vbl - time.trigger;
            elseif nf == SwitchFr(3) + 1
                T.img2Onset(thetrial) = vbl - time.trigger;
            end
            
            %% collect response
            if RealRun
                % retreive responses from button box
                [wkey, timeStamp] = BitsiBB.getResponse( timeout, true);
                if ismember( wkey, RespKeys) && ~ButPres && nf > SwitchFr(3)  % only response after image 2 is presented is valid
                    ButPres = 1;
                    T.resp(thetrial)   =   find(wkey==RespKeys);          % 1 or 2
                    T.rt(thetrial)     =   timeStamp - T.img2Onset(thetrial) - time.trigger;    % time from scanner trigger
                    fprintf('Responded %g, should be %g. \n', T.resp(thetrial), T.correctKey(thetrial));
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
            elseif keydown && ~RealRun && any(keyCode(RespKeys)) && ~ButPres && nf > SwitchFr(3)
                ButPres = 1;
                T.resp(thetrial)   =   find(find(keyCode)==RespKeys);
                T.rt(thetrial)     =   respT - T.img2Onset(thetrial) - time.trigger;
                fprintf('Responded %g, should be %g.\n', T.resp(thetrial), T.correctKey(thetrial));
            end % end of kbcheck
            
        end
        
        % at 1/3 and 2/3, have a long rest to release activation
        if ismember(thetrial,ceil(NTrial.total/3), ceil(Ntrial.total/3*2))
            msg = sprintf('Take a rest for %g seconds', time.rest/1000);
            for nfr = 1:Nfr.rest
                if nfr < Nfr.rest - w.RefreshRate*2
                    DrawFormattedText(w.Number, msg, 'center', 'center', Black);
                else
                    DrawFormattedText(w.Number, 'Be ready!', 'center', w.Center(2) - 30, Black);
                    Screen('DrawDots', w.Number, [0; 0], fix.Size, fix.Col, w.Center, 1);
                end
                
                Screen('DrawingFinished', w.Number);
                vbl = Screen('Flip', w.Number, vbl + .5 * w.ifi);
            end
        end
    end
    
    %% report participants'results
    RStrialId = ismember(T.movingDir,movingDirs(1:2)); % get the RS trials
    T.ifcorrect(RStrialId) = T.resp(RStrialId)== T.correctKey(RStrialId); % no reponse for MVPA trials are also correct
    figure;
    for thedir = 1:2
        for themovingdir = 1:2
            trialid = RStrialId & strcmp(T.direction,directions(thedir)) & T.movingDir==movingDirs(themovingdir);
            accRS(thedir,themovingdir) = mean(T.ifcorrect(trialid)); % only calculate acc for RS trials
            bar(thedir+(themovingdir-1.5)*0.4, accRS(thedir,themovingdir), 0.3); hold on
        end
    end
    
    rtRS   = mean(T.rt(RStrialId & T.rt ~=0)); % omit the unresponded trials
    fprintf('((( Accuracy for RS trials: %g %%, RT = %g s )))\n', 100 * mean(mean(accRS)), rtRS); % to print a % inside fprintf, you need %%
    
    
    % Monitor PTB performance
    
    dur.trialActual = 1000 * diff(T.img1Onset); % get the actual duraction for each trial
    dur.trialActual = dur.trialActual(dur.trialActual> 0 & dur.trialActual< time.rest);
    TrialDur = round(nanmean(dur.trialActual));
    TrialStd = round(nanstd( dur.trialActual));
    fprintf('((( Trial Duration %g/%g ms, SD = %g ms )))\n', TrialDur, dur.trial, TrialStd);
    %   figure;  scatter(ones(size(dur.trialActual)), dur.trialActual)
    time.total = toc/60;
    
    
    %% EXIT and clean up===================================================
    respondedRS = find((RStrialId & T.resp~=0));
    EndTxt = sprintf(['Well done! You completed run %g/%g.\n\n', ...
        'You responded %g trials out of %g trials. \n\n', ...
        'You can take a break while we prepare the next run.'],...
        RunNo, NRuns, length(respondedRS),length(find(RStrialId)));
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
    
    fileName1 = fullfile(RunDir, sprintf('sub-%03d_%s_%g.mat', SubNo, RunType, ThisRunNo));
    save(fileName1);
    fileName2 = fullfile(RunDir, sprintf('ResultTable_sub-%03d_%s_%g.csv', SubNo, RunType, ThisRunNo));
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