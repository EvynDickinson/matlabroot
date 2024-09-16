[vidname,vidpath] = uigetfile('*_1.avi');
videoList = dir([vidpath,vidname]);
videoStartTime = videoList(1).date(end-7:end);
disp(videoStartTime)

% % write the experiment details into the excel sheet
% isExcelFileOpen(xlFile); % check that file details can be written to spreadsheet
% for arena = 1:4
%     writecell({videoStartTime},xlFile,'Sheet','Exp List','Range', [Alphabet(Excel.starttime) num2str(XLrow(arena))])
%     writecell({facility},xlFile,'Sheet','Exp List','Range', [Alphabet(Excel.facility) num2str(XLrow(arena))])
%     writecell({'Y'}, xlFile, 'Sheet','Exp List','Range',[Alphabet(Excel.step1) num2str(XLrow(arena))]);
% end