function [NextRunNo, NOccur] = EstiRun(SubjNo,Runs,DataDir)

% Estimate next run (number) based on previously saved files. 

RunFiles = cell(length(Runs),1);
RunsDone = nan(length(Runs),1);
NOccur = nan(length(Runs),1);

for run = 1:length(Runs)
    
    NOccur(run) = nnz(strcmp(Runs(1:run), Runs{run}));
    
    if ~strcmp(Runs{run}, 'Demo') && ~strcmp(Runs{run}, 'Resting state')
        
        RunFiles{run} = fullfile(DataDir,sprintf('Subj%02d', SubjNo), Runs{run}, ...
            sprintf('Subj%02d_%s_%g.mat', SubjNo, Runs{run}, NOccur(run)));
        
    else
        
        RunFiles{run} = fullfile(DataDir,sprintf('Subj%02d', SubjNo), Runs{run}, ...
            sprintf('Subj%02d_%s.mat', SubjNo, Runs{run}));
    
    end
    
    RunsDone(run) = exist(RunFiles{run}, 'file');
    
end

NextRunNo = find(RunsDone==0, 1, 'first');

end