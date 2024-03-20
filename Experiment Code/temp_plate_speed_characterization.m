
file_path = 'G:\My Drive\Jeanne Lab\DATA\02.29.2024\Speed_test_4.csv';
tempLog = readmatrix(file_path);

time = tempLog(:,1);
target = tempLog(:,2);
temp = tempLog(:,3);
work = tempLog(:,4);


%% Graph timecourse
LW = 2; %linewidth

fig = getfig('',true);
hold on
plot(time, target,'Color','white','LineWidth',LW)
plot(time, temp, 'Color',Color('teal'),'LineWidth',LW)
xlabel('time (sec)')
ylabel('temp (\circC)')
formatFig(fig,true);

xlim([0,3000])

save_figure(fig,[file_path(1:end-4) ' timecourse.png'],'-png')

