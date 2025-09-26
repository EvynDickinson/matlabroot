
function getArenaOutlines(img,savePath)

baseFolder = getCloudPath;

% Find the arena centers:
try load([savePath ' parameters.mat']);
catch 
    disp('Save parameters and then retry')
    return
end

demo_well_loc = load([baseFolder(1:end-6) '/demo_well_loc.mat']);
demo_well_label = [];
for i = 1:16
    demo_well_label{i} = num2str(i);
end
RGB = insertText(img,demo_well_loc.demo_well_loc,demo_well_label,'FontSize', 24); % previously 15

% get numbers for arena sizes / labels
radii = 165; %well surround regions
r = 435; % radius of the arena
arenaIdx = {'A', 'B', 'C', 'D'};

% Select the well centers from top right corner and go clockwise:
well_loc_file = [savePath ' well_locations.mat'];
if ~exist(well_loc_file,'file')
    h = warndlg('Start at top left arena and click wells in clockwise order starting at 12. Proceed clockwise through the arenas.'); 
    uiwait(h)
    % add text hints for outline order:
    well_loc = readPoints(RGB,16); % get locations of the wells from the image

    % Allocate well locations from manually selected locations
    decoder_new = [11, 12, 9, 10, 15, 16, 13, 14, 7, 8, 5, 6, 3, 4, 1, 2];
    wellcenters = zeros(2,16); %blank locations to fill with coordinates
    for arena = 1:4
        % Divide well centers by arena
        roi = (arena-1)*4 + 1 : (arena-1)*4 + 4;
        % check if the experiment is from the previous version for arena alignment
        wellcenters(:,roi) =  well_loc(:,decoder_new(roi));
    end
        
    % save the arena well locations:
    save(well_loc_file, 'wellcenters','arenaIdx','r');

else
    load(well_loc_file);
end

arenaData = struct;
CList = {'DeepPink','Orange', 'Lime', 'DodgerBlue'};
for arena = 1:4
    kolor = Color(CList{arena});
    arenaData(arena).r = r;
    arenaData(arena).name = ['Arena ' arenaIdx{arena}];
    arenaData(arena).color = kolor;
    % Divide well centers by arena
    roi = (arena-1)*4 + 1 : (arena-1)*4 + 4;
    arenaData(arena).wellcenters = wellcenters(:,roi);

    % Find the center of each arena
    x1 = arenaData(arena).wellcenters(1,1:2:4);
    y1 = arenaData(arena).wellcenters(2,1:2:4);
    x2 = arenaData(arena).wellcenters(1,2:2:4);
    y2 = arenaData(arena).wellcenters(2,2:2:4);
    [xi,yi] = polyxpoly(x1,y1,x2,y2);
    arenaData(arena).centre = [xi;yi];

    % Well labels
    AI = ['Arena' arenaIdx{arena}];
    wellLabels = {parameters.(AI).well_1;...
                  parameters.(AI).well_2;...
                  parameters.(AI).well_3;...
                  parameters.(AI).well_4};   
    arenaData(arena).wellLabels = wellLabels;
end


% Visual check of alignment, ID, and well contents
fig = figure; 
    imshow(img); set(fig, 'color', 'k', 'pos', [341 50 1101 946]);
    hold on
    for arena = 1:4 
        kolor = arenaData(arena).color;
        centre = arenaData(arena).centre;
        WC = arenaData(arena).wellcenters;
        WL = arenaData(arena).wellLabels;
        % draw arena circles
        viscircles(centre', r, 'color', kolor);
        % label arenas to show that it's what we think:
        text(centre(1),centre(2),arenaData(arena).name,'color',...
             kolor,'fontsize', 11,'horizontalAlignment', 'center')
        for well = 1:4
            % label wells
            text(WC(1,well),WC(2,well)-75,strrep(WL{well},'_',' '),'color',...
             'w','fontsize',10,'horizontalAlignment','center')
            % draw well ROIs
            viscircles(WC(:,well)', radii, 'color', 'w','linestyle','--','linewidth', 0.1);
        end
    end
    tltstr = strsplit(savePath,'\');
    titlestr = tltstr{end};
    titlestr = strrep(tltstr,'_', ' ');
    title(titlestr,'color','w')
figure(fig)


% Does the arena fit work???
answer = questdlg('Okay arena alignment?');
switch answer
    case 'Yes'
        export_fig(fig, [savePath ' arena outlines.png'],'-png' , '-nocrop', '-r80' , '-painters', '-rgb');
        close(fig)
        % Save well labels and alignments
        save([savePath ' arena coordinates.mat'], 'arenaData')
        % save_figure(fig, [analysisDir expName ' arena outlines'],'-png',true);
    case 'Cancel'
        % don't save figure but keep going
    case 'No'
        return

end



end