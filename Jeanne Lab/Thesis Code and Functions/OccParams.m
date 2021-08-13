
function textStructure = OccParams(vidPath, thinBuffer, thickBuffer)
% textStucture = OccParams(vidPath, thinBuffer, thickBuffer)
% 
% 
% ES Dickinson, Yale University, Aug 2021


textStructure.colors = repmat([0,0,0],[5,1]);
textStructure.str = {'raw', 'binary', 'occupation probabilty', '1', '0'};
textStructure.thinBuffer = thinBuffer;
textStructure.thickBuffer = thickBuffer;

movieInfo = VideoReader(vidPath);
height = movieInfo.Height;
width = movieInfo.Width;

% space buffers
MT = zeros(height, width);
line_size = thinBuffer;
vert = ones(height,line_size);
hori = ones(line_size,(width*3)+(line_size*4));

% build a colorbar
barWidth = thickBuffer;
horimini = ones(line_size,barWidth+(line_size*2));
buffZone = zeros(height+(2*line_size), barWidth);
clims = [1,0];
cmap = repmat(linspace(clims(1),clims(2), height)',1,barWidth);
cbar = [horimini; vert,cmap,vert; horimini];

% set up text locations:
offsets = [2,4]; %top&bottom, right&left
blankFrame = ([vert, MT, vert, MT,vert, MT, vert]);
blankFrame = [hori;blankFrame; hori];    
blankFrame = [blankFrame,buffZone,cbar];
blankFrame = [zeros(barWidth*offsets(1),size(blankFrame,2)); blankFrame]; % top buffer
blankFrame = [zeros(size(blankFrame,1), barWidth),blankFrame];%right side buffer
blankFrame = [blankFrame, zeros(size(blankFrame,1), offsets(2)*barWidth)]; % left side buffer
blankFrame = [blankFrame; zeros(barWidth,size(blankFrame,2))]; % bottom buffer
frameSize = size(blankFrame);

endloc = frameSize(2)-barWidth*(offsets(2)-1);
pos = [barWidth, barWidth+line_size+round(width/2);...
        barWidth, barWidth+2*line_size+round(width*1.5);...
        barWidth, barWidth+3*line_size+round(width*2.5);...
        2*barWidth, endloc;...
        frameSize(1)-barWidth, endloc];
textStructure.loc = [pos(:,2), pos(:,1)];
textStructure.frameSize = frameSize;
textStructure.vert = vert;
textStructure.hori = hori;
textStructure.cbar = cbar;
textStructure.buffZone = buffZone;
textStructure.offsets = offsets;
textStructure.height = height;
textStructure.width = width;

end