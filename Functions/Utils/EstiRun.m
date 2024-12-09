function [NextRunNo, NOccur] = EstiRun(SubNo,Runs,DataDir)

% Estimate next run (number) based on previously saved files. 

RunFiles = cell(length(Runs),1);
RunsDone = nan(length(Runs),1);
NOccur = nan(length(Runs),1);

for run = 1:length(Runs)
    
    NOccur(run) = nnz(strcmp(Runs(1:run), Runs{run}));
    
    if ~strcmp(Runs{run}, 'Resting state')
        
        RunFiles{run} = fullfile(DataDir,sprintf('Sub%02d', SubNo), Runs{run}, ...
            sprintf('Sub%02d_%s_%g.mat', SubNo, Runs{run}, NOccur(run)));
        
    else
        
        RunFiles{run} = fullfile(DataDir,sprintf('Sub%02d', SubNo), Runs{run}, ...
            sprintf('Sub%02d_%s.mat', SubNo, Runs{run}));
    
    end
    
    RunsDone(run) = exist(RunFiles{run}, 'file');
    
end

NextRunNo = find(RunsDone==0, 1, 'first');

end