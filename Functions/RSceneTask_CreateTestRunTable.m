function AllTrials = RSceneTask_CreateTestRunTable
% 
% 
% 
% 
% 
% 
% 
% 

TotNTrials = 336;
NRuns = 7;
TrialsPerRun = TotNTrials/NRuns;

NScenes = 20;

Prob_cons = 0.75;

NCons = Prob_cons * TrialsPerRun; % consistent trials per run
assert(round(NCons)==NCons);

%% AllTrials table
% Create each run separately - we want everything balanced *within* run

Views = {'A', 'B'};
Rotations = [30, 90];

NConditions = numel(Views) * numel(Rotations);

% add option for different tasks? (e.g. 'sizes' instead of 'orients')
Variables = {'Scene', 'InitView', 'FinalView', 'Diff', ...
    'Orients', 'Consistent', 'Hit', 'RT'};

% Automatically if trials are evenly split between initial views, and then
% between final views, half will be wide and half narrow.
% Do that separately for consistent and inconsistent to avoid imbalances.

AllTrials = nan(TotNTrials, numel(Variables));

for i = 1:TrialsPerRun:TotNTrials
    
    ThisRun = AllTrials(i:i + TrialsPerRun - 1, :);
    ThisRun(1:NCons, 2:3) = repmat(allcomb(1:2, Rotations), NCons/NConditions, 1);
    ThisRun(NCons + 1:end, 2:3) = repmat(allcomb(1:2, Rotations), (TrialsPerRun - NCons)/4, 1);
    ThisRun(1:NCons, 6) = 1;
    ThisRun(NCons + 1:end, 6) = 0;
    ThisRun = ThisRun(randperm(TrialsPerRun), :);
    AllTrials(i:i + TrialsPerRun - 1, :) = ThisRun;
    
end

AllTrials = array2table(AllTrials, 'VariableNames', Variables);

%% Scenes' order:

SceneOrder = repmat((1:NScenes)', ceil(TotNTrials/NScenes), 1);
SceneOrder = SceneOrder(randperm(length(SceneOrder)));

AllTrials.Scene = SceneOrder(1:TotNTrials);

% Add sequence:

Seq = nan(TotNTrials, 5);

Seq(:, 1) = 0;
for i = 1:TotNTrials
    if AllTrials.FinalView(i)==30, Limit = 25; else, Limit = 60; end
    Seq(i, 2:4) = sort(randsample(15:5:Limit, 3));
end
Seq(:, 5) = AllTrials.FinalView;

AllTrials.Sequence = num2cell(Seq, 2);


