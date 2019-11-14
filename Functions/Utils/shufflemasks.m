function maskersnew = shufflemasks(maskersold,lengthnew)
%SHUFFLEMASKS: shuffle and repeat array without repetions
%
% INPUT
% maskersold: array of numbers
% lengthnew: length of new array
%
% OUTPUT
% maskersnew: new array of shuffled numbers of length n
%
%
% SG 24-04-2018


%% preparation
% count number of elements (once)
elems       =   numel(maskersold);

% count number of iterations needed to reach requested number of items
nitis       =   ceil( lengthnew/elems);

% make empty array
allmasks    =   NaN( 1, nitis*elems);

%% loop through iterations
for niti    =   1 : nitis
    
    % shuffle items
    maskershuf  =   maskersold( randperm( elems));
    
    % if repetition, swap first item with random item from array
    if maskershuf(1) == maskersold(end)
        newpos  =   (randperm(elems-1,1)+1);                % determine random item
        maskershuf([1 newpos]) = maskershuf([newpos 1]);    % swap
    end
    
    % append shuffled array
    allmasks( (niti-1)*elems+1 : niti*elems) = maskershuf;
    
    % replace old array with (shuffled) new array
    maskersold  =   maskershuf;
end

% keep only requested number of items
maskersnew  =   allmasks(1, 1 : lengthnew);


