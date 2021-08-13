
function newFrame = OccBuildFrame(rawImg, bwImg, occIm, txt)
% newFrame = OccBuildFrame(rawImg, bwImg, occIm, textStructure)
% 
% 
% ES Dickinson, Yale University, Aug 2021


newFrame = ([txt.vert, rawImg, txt.vert, bwImg, txt.vert, occIm, txt.vert]);
newFrame = [txt.hori; newFrame; txt.hori];    
newFrame = [newFrame,txt.buffZone,txt.cbar];
newFrame = [zeros(txt.thickBuffer*txt.offsets(1),size(newFrame,2)); newFrame]; % top buffer
newFrame = [zeros(size(newFrame,1), txt.thickBuffer),newFrame];%right side buffer
newFrame = [newFrame, zeros(size(newFrame,1), txt.offsets(2)*txt.thickBuffer)]; % left side buffer
newFrame = [newFrame; zeros(txt.thickBuffer, size(newFrame,2))]; % bottom buffer

% insert title and text...
newFrame = (insertText(newFrame, txt.loc, txt.str,'FontSize',40,'BoxColor',...
    txt.colors,'BoxOpacity',1,'TextColor','white', 'anchorpoint', 'center'));

title_loc = [size(newFrame,2)-txt.thickBuffer, txt.thickBuffer+round(txt.height/2)];
newFrame = (insertText(newFrame, title_loc, 'Prob.','FontSize',40,'BoxColor',...
    'black','BoxOpacity',1,'TextColor','white', 'anchorpoint', 'center'));

end