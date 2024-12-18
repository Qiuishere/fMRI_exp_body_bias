function [TrialMatrix, TargetTrials] = RSceneTask_CreateTrainRunMat
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 

% 4 repetitions

NStimTypes = 4; % n. of stimulus types

TotNTrials = 360;

NStimuli = 20;
NTargets = 10; % how many targets per stimulus type

BlockSize = 18; % how many trials in a mini-block 

TrialMatrix = zeros(NStimTypes, TotNTrials/NStimTypes);
TargetTrials = TrialMatrix; % where repetition occurs (and subject should respond)

% Generate each row of the matrices:

for i = 1:NStimTypes
    
    % 
    TargetsPerBlock = zeros(1, size(TrialMatrix, 2));
    TargetsPerBlock(1:NTargets) = 1;
    TargetsPerBlock = TargetsPerBlock(randperm(size(TrialMatrix, 2)));
    TargetsPerBlock = reshape(TargetsPerBlock, [BlockSize, size(TrialMatrix, 2)/BlockSize]);
    TargetsPerBlock = sum(TargetsPerBlock); 
    % 1 x NBlocks (5): target trials assigned to each block
    
    Stimuli = randperm(NStimuli);
    Indx = 1;
    
    for j = 1:BlockSize:size(TrialMatrix, 2)
        
        BlockNo = ceil(j/BlockSize); % which block we are in
        
        % get portion (row and columns) of matrices pertaining to this
        % block:
        TheseTargets = TargetTrials(i, j:(j+BlockSize-1));
        TheseTrials = TrialMatrix(i, j:(j+BlockSize-1));
        
        TargetsChosen = 0;
        
        while ~TargetsChosen % to make sure there aren't two targets in a row
            
            % Where target will be in a miniblock
            Targets = randperm(BlockSize - 1) + 1; 
            % (numbers from 2 to BlockSize, because of course the first
            % trial cannot be a target)
            
            Targets = Targets(1:TargetsPerBlock(BlockNo));
            % get the correct number of targets for this miniblock
            
            if ~any(diff(Targets)==1) % no consecutive numbers (2 targets in a row)
                TargetsChosen = 1; % we're good
            end
            
        end
        
        TheseTargets(Targets) = 1;
        
        for col = 1:size(TheseTrials, 2)
            
            if ismember(col, Targets) % if target
                
                TheseTrials(col) = TheseTrials(col-1); % repeat previous stimulus
                
            else
                
                TheseTrials(col) = Stimuli(Indx);
                Indx = Indx + 1;
                
            end
            
            if Indx == NStimuli
                
                Stimuli = randperm(NStimuli);
                Indx = 1;
                
            end
            
        end
        
        TargetTrials(i, j:(j+BlockSize-1)) = TheseTargets;
        TrialMatrix(i, j:(j+BlockSize-1)) = TheseTrials;
        
        clear TheseTargets TheseTrials
        
    end
end
