%% Slightly different things in Surya's code

%% Screen settings

screensize = [369 277];
dispdist = 955;

%%
% text settings
fontsize    =   0.4;
textcol     =   [255 255 255];
textcol2    =   [  0   0   0];
ydev        =   0.5;

% background
backgr = ones(1,3) .* 104;

fixsize(1,:)    =   [0.2 0.2];
fixsize(2,:)    =   [0.1 0.1];
fixcols(1,:)    =   [255 255 255];
fixcols(2,:)    =   [  0   0   0];

ncount      =   3;  % countdown from ncount to 1 before run starts (after trigger)
t.count     =   0.8;  % each count takes t.count seconds

%% hard coded experiment settings (LOCALIZER & TRAINING runs)

% general settings
blockreps =     4;  % blockorders are fixed manually (adjust in section below)
imcount   =     20; % number of images per mini-block
ntargets  =     2;  % number of 1-back targets per block

% timing
t.off   =   0.300;
t.on    =   0.450;

% sizes (dva) and colors (rgb)
cawidth         =   [12 12];     % Based on Peelen, Fei Fei, & Kastner, 2009

% training task 1
mask.size       =   3;
mask.ecc        =   4;
mask.col        =   backgr(1);
mask.opac       =   0.7;

% training task 2
f.large         =   1.30;
f.small         =   0.70;
f.targs         =   0.1634; % target eccentricity (fraction of image width)


%% Block orders for localizers

% Block orders for 2 x 2 localizers (run 1/2, session 1/2)
% ONLY ONE SESSION

blockorders(:,:,1)  =   [...
    0 1 2 3 4; ...
    0 3 1 4 2; ...
    0 4 3 2 1; ...
    0 2 4 1 3];

blockorders(:,:,2)  =   [...
    0 4 2 1 3; ...
    0 2 3 4 1; ...
    0 1 4 3 2; ...
    0 3 1 2 4];

% blockorders(:,:,1,2)  =   [...
%     0 2 1 4 3; ...
%     0 4 2 3 1; ...
%     0 3 1 4 2; ...
%     0 1 3 2 4];
% 
% blockorders(:,:,2,2)  =   [...
%     0 3 2 4 1; ...
%     0 1 4 2 3; ...
%     0 2 3 1 4; ...
%     0 4 1 3 2];


%% conversions
   
% convert dva to pix
ydev            =   fliplr(round( ang2pix3( ydev,         dispdist, windowRect(3:4), screensize)));
fontsize        =   fliplr(round( ang2pix3( fontsize,     dispdist, windowRect(3:4), screensize)));
fixsize(1,:)    =   round(ang2pix3( fixsize(1,:), dispdist, windowRect(3:4), screensize));
fixsize(2,:)    =   round(ang2pix3( fixsize(2,:), dispdist, windowRect(3:4), screensize));
cawidth         =   round(ang2pix3( cawidth,   dispdist, windowRect(3:4), screensize));
mask.size       =   round(ang2pix3( mask.size, dispdist, windowRect(3:4), screensize));
mask.ecc        =   round(ang2pix3( mask.ecc,  dispdist, windowRect(3:4), screensize));

mask.ecc        =   mask.ecc(1);

% select y component
ydev        =   ydev(2);
fontsize    =   fontsize(2);

% transpose matrices
fixcols     =   fixcols';

% convert seconds to frames
nframes.on  =   round(t.on  .* RefRate);
nframes.off =   round(t.off .* RefRate);

% timing
t.block =   (t.off + t.on).*imcount;
trialset =  round( (nframes.on + nframes.off)/RefRate*1000);

%% Stimulus directory

if IsOSX
    StimDir = '/Users/u010155/Desktop/localizer';
else
    StimDir = fullfile('Stimuli', RunType);
end
