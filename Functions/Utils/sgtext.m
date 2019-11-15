function sgtext( w, textstr, screenrect, fontsize, textfont, textcol, xydeviation)
% Draw text (relative to center of screen)

% text settings
Screen( 'TextFont', w, textfont);
Screen( 'TextSize', w, fontsize);
%Screen( 'TextStyle', w, 1+2);

% get text box
bbox        =   Screen( 'TextBounds', w, textstr); %bbox = bbox + [0 0 5 5];
textrect    =   CenterRect( bbox, screenrect);

if exist('xydeviation','var')
   textrect(1,1:2) = textrect(1,1:2) + xydeviation; 
end

% draw text
Screen( 'DrawText', w, textstr, textrect(1), textrect(2), textcol);

end

