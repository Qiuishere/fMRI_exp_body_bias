function [pix, dva, x] = ang2pix3( dvaorpix, dispdist, sspix, scrmm)


%%
% Convert degrees of visual angle into number of pixels or vice versa [x y]
%
%
% INPUT:
%
% - [x y] in degrees of visual angle OR pixels
% - Distance between the observer and the display in mm
% - Screen size [x y] in pixels
% - Screen size [x y] in mm
%
%
% OUTPUT:
%
% - Stimulus size in pixels,  assuming input in dva     [x y]
% - Stimulus size in dva,     assuming input in pixels  [x y]
%
% - Extra output in struct:
%                             - x.ppmm: number of pixels per mm
%                             - x.mmpp: pixel size in mm
%                             - x.smm:  screensize in mm
%                             - x.spx:  screensize in pixels
%
% -------------------------------------------------------------------------


%% retreive screen properties
x(1).smm(1,1:2)             =   scrmm;                      % Screen size in mm
x.spx                       =   sspix;                      % Screen size in pix


%% convert dva to pixels
osmm      =   2 .* dispdist .* tand(0.5 .* dvaorpix);       
x.ppmm    =   x.spx ./ x.smm;                               % Number of pixels per mm
pix       =   osmm .* x.ppmm;


%% convert pixels to dva
x.mmpp    =   1 ./ x.ppmm;                                  % pixel size in mm
hangle    =   (0.5 .* dvaorpix .* x.mmpp) ./ dispdist;
dva       =   2 .* atand(hangle);


