

% plot temperature protocols that are used in the Behavior paper 2024

%% select & load example protocols:

 % Import and Load Temperature Protocol Options
[excelfile, Excel, xlFile] = load_QuadBowlExperiments;

[~,~,temp_protocol_XL] = xlsread(xlFile, 'Temp Protocols');
protocolList = temp_protocol_XL(2:end,1);
idx = listdlg("PromptString",'Select the temp protocols to plot',...
                    'SelectionMode', 'multiple', 'ListString',protocolList,...
                    'ListSize',[250,500]);
temp_protocols = protocolList(idx);
num = length(temp_protocols);

folder = getDataPath(1,0);

data = [];

% find experiments with the given temperature protocols
for i = 1:num

    % find a trial with the appropriate temp protocol
    TP = temp_protocols{i};
    all_protocols = excelfile(:,Excel.protocol);
    rows = find(strcmp(all_protocols,TP));
    row = rows(end);
    trial_ID = excelfile{row, Excel.trialID};
    exp_ID = excelfile{row, Excel.expID};
    arena = excelfile{row, Excel.arena};


    % load the file
    dummy = load([folder trial_ID '/' exp_ID arena ' timecourse data.mat']);
    data(i).temp = dummy.data.occupancy.temp;
    data(i).time = dummy.data.occupancy.time;
    
    clear dummy
end

%% Plot the temperature protocols

cMap = colormap(turbo(num+2));
fig = gcf;
close(fig)
cMap(1,:) = []; cMap(num+1,:) = []; %remove ends since they are too dark for black background
sSpan = 180;

fig = getfig('',1);
    hold on
    for i = 1:num
        plot(data(i).time, smooth(data(i).temp,sSpan),'linewidth', 1,'Color',cMap(i,:))
    end
    
    xlabel('time (min)')
    ylabel('temp (\circC)')
    
    fig = formatFig(fig, false);
    set(gca, 'TickDir','out')
    legend(strrep(temp_protocols,'_',' '),'fontsize', 8,'box', 'off');

% Save figure
figDir = getCloudPath;
saveDir = figDir(1:end-5);
save_figure(fig,[saveDir 'Manuscripts\2024 Behavior Paper\Figures/Temp Protocols'],'-pdf');

fig_path = 'D:\Evyn Lab Data\Manuscript Figures\Temp Protocols';
save_figure(fig,fig_path);
































