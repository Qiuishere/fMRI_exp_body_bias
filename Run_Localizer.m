function Run_Localizer(SubjNo, RunInfo, DataDir)

% Based as closely as possible on Surya Gayet's code (for consistency)
%
%

%% startup

global Environment
global RealRun

RunType = 'Localizer';

RunNo = RunInfo(1);
NRuns = RunInfo(2);
ThisRunNo = RunInfo(3);

% load settings
Settings_General;
Settings_Loc;

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

%% Get (and report) flip interval, and issue warning if huge deviation

ifi = Screen('GetFlipInterval',w);
fprintf('\nFlip interval: %2.4g ms (should be %2.4g ms).\n\n', ifi*1000, (1000/RefRate));

if ifi > (1/RefRate)*1.5
    warning( 'This is a serious timing problem! \n');
end

%% Message
% while loading images and waiting for scanner trigger

InstrTxt = sprintf(['Run %g/%g\n\n',...
    'Press a button every time\nthe same image appears twice in a row.'], RunNo, NRuns);

DrawFormattedText(w, InstrTxt, 'center', 'center', White);

Screen('Flip', w);

%% Conditions

% determine block order
blockorder = blockorders(:, :, ThisRunNo);

% add target trials (image repetitions)
ts                  =   zeros(imcount-1, 1); % first trial cannot be a repetition
ts(1:ntargets, 1)   =   1;
behav.blocks        =   [reshape( blockorder', 1, numel(blockorder)) 0];

% create target occurence matrix
behav.trials        =   zeros(imcount, length(behav.blocks));
for bi = 1: length( behav.blocks)
    if mod( bi, size(blockorder,2))~= 1
        behav.trials(2:end, bi) = Shuffle(ts);
    end
end

%% load images

% directories
cond_names  = listdir(StimDir);
nconditions = length(cond_names);

% determine size
myimsize = [cawidth(2) cawidth(1)]; %imresize uses [y x]

% empty array
filecount = imcount;

imTex = NaN( length(cond_names), filecount);

% loop through categories
for ncat = 1 : length(cond_names)

    % list files
    filelist = listdir(fullfile(StimDir, cond_names{ncat}, '*.jpg'));

    if length( filelist) > imcount
        filelist = filelist{1:imcount};
    end

    % loop through files and create texture from image
    for nim = 1:filecount

        thisim = imread(fullfile(StimDir, cond_names{ncat}, filelist{nim}));
        thisim = imresize( thisim(:,:,1), myimsize);

        imTex(ncat, nim) = Screen('MakeTexture', w, thisim);

    end % end of image loop
end % end of category loop


%% assign images to trials

behav.imag = NaN( size(behav.trials));
for nblock = 1 : size(behav.trials,2)

    % if NOT a baseline block
    if behav.blocks(nblock)

            % randomly pick 'imcount-ntargets' images (out of all cat. images)
            behav.imag( behav.trials(:, nblock)==0, nblock) = shufflemasks(imTex(behav.blocks(nblock), :), imcount-ntargets)';

            % and repeat image if target trial (and localizer run)
            behav.imag( logical( behav.trials(:, nblock)), nblock)    =  behav.imag( logical( [diff( behav.trials(:, nblock))==1; 0]), nblock);
    end
end


%% stimulus locations

sceRect         =   CenterRect( [0 0 fliplr( size( thisim)) ], windowRect);
fixRect1        =   CenterRect( [0 0 fixsize(1,:)           ], windowRect);
fixRect2        =   CenterRect( [0 0 fixsize(2,:)           ], windowRect);

fixrects        =   [fixRect1; fixRect2]';

% target rects: source rects
sourceRect(1,:) =   [0 0 fliplr( size( thisim)) ];                                               % LEFT
sourceRect(1,3) =   round( size( thisim, 2)/2);
sourceRect(2,:) =   [0 0 fliplr( size( thisim)) ];                                               % RIGHT
sourceRect(2,1) =   round( size( thisim, 2)/2);

% target rects: destination rects
devi.small      =   abs( (f.targs * size( thisim, 2)) - (f.targs * size( thisim, 2) * f.small)); % deviation due to rescaling
devi.large      =   abs( (f.targs * size( thisim, 2)) - (f.targs * size( thisim, 2) * f.large)); % deviation due to rescaling

smasize         =  f.small *  [0.5 * size( thisim, 2) size( thisim, 1)];                         % size of small rect
bigsize         =  f.large *  [0.5 * size( thisim, 2) size( thisim, 1)];                         % size of large rect

targRect(1,:,1) =  CenterRect( [0 0 smasize ], windowRect) - 0.5*[smasize(1) 0 smasize(1) 0];    % LEFT rect  (small)
targRect(1,:,1) =   targRect(1,:,1) - [devi.small 0 devi.small 0];

targRect(2,:,1) =  CenterRect( [0 0 smasize ], windowRect) + 0.5*[smasize(1) 0 smasize(1) 0];    % RIGHT rect (small)
targRect(2,:,1) =   targRect(2,:,1) + [devi.small 0 devi.small 0];

targRect(1,:,2) =  CenterRect( [0 0 bigsize ], windowRect) - 0.5*[bigsize(1) 0 bigsize(1) 0];    % LEFT rect  (large)
targRect(1,:,2) =   targRect(1,:,2) + [devi.large 0 devi.large 0];

targRect(2,:,2) =  CenterRect( [0 0 bigsize ], windowRect) + 0.5*[bigsize(1) 0 bigsize(1) 0];    % RIGHT rect (large)
targRect(2,:,2) =   targRect(2,:,2) - [devi.large 0 devi.large 0];

targRect        =   round( targRect);


%% scanner trigger

if RealRun

    warning('Waiting for scanner trigger');

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
tstamp.trigger  =   GetSecs;

fprintf('\n SCANNER TRIGGER (test run %g): %s \n', ThisRunNo, datestr(now,'dd-mm-yyyy HH:MM:SS'));

%% start message
for n = 1:ncount

    CountDownTxt = sprintf('\n\n\nStarting in %g', 3 - n + 1);
    DrawFormattedText(w, CountDownTxt, 'center', 'center', White);
    Screen('FillOval', w, fixcols, fixrects);
    Screen('Flip', w);
    WaitSecs(t.count);

end


%% pre-alocate vectors

% create matrices to save responses
behav.name       =   cond_names;
behav.resp       =   NaN( size( behav.trials));
behav.time       =   behav.resp;
behav.check      =   behav.resp;

% create matrices to save timestamps
tstamp.blocks    =   behav.blocks;
tstamp.event     =   NaN( size( behav.imag));
tstamp.offset    =   tstamp.event;

% set relevant variables
nblock      =   0;
oldkey      =   0;
butpres     =   0;
last_targ   =   tstamp.trigger;
if RealRun
    BitsiBB.clearResponses();
end

%% loop through stimulus presentation
% loop through mini-blocks

for nblock = 1 : length(tstamp.blocks)

    blockims = behav.imag(:, nblock);

    % print current block
    fprintf( '  - mini-block %g/%g', nblock, length(tstamp.blocks));

    % loop through images within block
    for nimage = 1: imcount

        % draw stimuli
        for nframe = 1 : nframes.on + nframes.off

            % draw image if block > 0
            if nframe <= nframes.on && tstamp.blocks(nblock)
                Screen( 'DrawTexture', w, blockims(nimage), [], sceRect);
            end

            % draw fixation
            Screen('FillOval', w, fixcols, fixrects);
            tnow = Screen('Flip', w);

            % log timestamp
            if nframe == 1
                tstamp.event(  nimage, nblock) = tnow;
                if logical( behav.trials( nimage, nblock))
                    last_targ   =   tnow;
                end
            elseif nframe == nframes.on + 1
                tstamp.offset( nimage, nblock) = tnow;
            end

            if RealRun

                % retreive responses from button box
                [wkey, timeStamp] = BitsiBB.getResponse( timeout, true);

                if ismember( wkey, RespKeys) && ~butpres
                    butpres = 1;
                    behav.resp(nimage, nblock)   =   find(wkey==RespKeys);          % 1 or 2
                    behav.time(nimage, nblock)   =   timeStamp - tstamp.trigger;    % time from scanner trigger
                    behav.check(nimage, nblock)  =   timeStamp - last_targ;
                elseif wkey == 0
                    butpres = 0;
                end
            end

            % retreive responses from keyboard (to abort run)
            if IsOSX
                [keydown, ~, keyCode] = KbCheck(-1);
            else
                [keydown, ~, keyCode] = KbCheck;
            end

            if keydown  && keyCode(EscKey)
                fprintf('\n\nExperiment terminated at %s, diary closed...\n', datestr(now));
                if RealRun
                    close(BitsiScanner);
                    fprintf('All serial ports closed.\n');
                end
                close(BitsiBB);
                delete(instrfind);
                ShowCursor;
                diary off
                clear temprun
                sca; return
            elseif keydown && ~RealRun && any(keyCode(RespKeys)) && ~butpres
                butpres = 1;
                behav.resp(nimage, nblock)   =   find(find(keyCode)==RespKeys);
            end % end of kbcheck

        end % end of current frame
    end % end of current image/trial

    % monitor participant performance
    fprintf(' -> Targets reported: %g/%g',  sum( ~isnan( behav.resp(:, nblock))), nansum( logical( behav.trials(:, nblock))));

    % monitor PTB performance
    evdur    =  diff( tstamp.event(:, nblock), 1);
    trialdur =  round( 1000 * nanmean( evdur(:)));
    trialstd =  round( 1000 * nanstd( evdur(:)));
    fprintf('(trial dur. %g/%g ms, SD = %g)\n', trialdur, trialset, trialstd);

    if ~tstamp.blocks(nblock)
        fprintf('\n');
    end

end % end of presentation loop (i.e., current block)

%% save stuff
tstamp.diff = tstamp.event - tstamp.trigger;

tstamp.names = cond_names;

save(DataFile, 'tstamp', 'behav');
save(BUpFile); % save everything


%% evaluate performance

if RealRun

    x1 = behav.check( ~isnan( behav.check));
    x1 = x1(:);
    perf = sum( x1 > 0.35 & x1 < 1.5); % random RT thresholds
    ntargs = length( behav.blocks(logical(behav.blocks))) .* ntargets;

    percorr = 100 * perf/ntargs;

    %% exit message

    % select message, depending on current run (w.r.t. planned runs)
    m1 = sprintf('Well done, you finished the run!\nTargets detected: %2.1f%%.', percorr);
    m2 = 'You can take a little break,\nwhile the next run is loading.';


    % display message
    DrawFormattedText(w, [m1, '\n', m2], 'center', 'center', White);
    Screen('Flip', w);

end

%% clean up

% close textures
Screen('Close',imTex(:));

% close trigger and response ports
close(BitsiBB);
if RealRun
    close(BitsiScanner);
    wmsg = 'scanner and response box';
else
    wmsg = 'virtual response box';
end
delete(instrfind);
fprintf('\nSerial ports (%s) CLOSED at %s. \n\n', wmsg, datestr(now));

% report run duration
compdur  =   round( size( tstamp.blocks, 2).*t.block + t.count.*ncount);
realdur  =   round( GetSecs - tstamp.trigger);
fprintf('RUN DURATION: %g seconds (expected duration: %g seconds). \n\n', realdur, compdur);

% report trial duration
evdur    =  diff( tstamp.event, 1);
trialdur =  round( 1000 * nanmean( evdur(:)));
trialstd =  round( 1000 * nanstd( evdur(:)));
fprintf('Single trial duration: %g ms, SD = %g (expected: %g ms). \n\n', trialdur, trialstd, trialset);

% clean log
fprintf('Diary closed (%s)\n\n', datestr(now));
diary off

% clean memory
try %#ok<TRYNC>
    jheapcl;
end

% wait for space press by experimenter
warning('Press space to continue...')
waitforspace; waitfornokey;

sca; ShowCursor;
