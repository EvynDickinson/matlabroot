
function tempPoints = getTempTurnPoints(TempProtocolString)
% tPoints = getTempTurnPoints(tempList{ii})


% pull the appropriate start|end points:
switch TempProtocolString
    case 'double_eaton_cooling_warming_ramp'
        tempPoints.down = [5720,34628;...
                         64819,93739];
        tempPoints.up = [34629,64818;...
                           93740,141259];
    case 'double_eaton_warming_cooling_ramp'
        tempPoints.down = [5720,34628;...
                         64819,93739];
        tempPoints.up = [34629,64818;...
                           93740,141259];
    case 'eaton_cooling_warming_ramp'
        tempPoints.down = [5700,34799];
        tempPoints.up = [34800,64000];
    case 'eaton_warming_cooling_ramp'
        tempPoints.down = [34800,64000];
        tempPoints.up = [5700,34799];
end

% increasing|decreasing ROI:
nUp = size(tempPoints.up,1);
nDown = size(tempPoints.down,1);
tempPoints.nUp = nUp;
tempPoints.nDown = nDown;
UpROI = [];
for ii = 1:nUp
    roi = tempPoints.up(ii,1):tempPoints.up(ii,2);
    UpROI = [UpROI,roi];
end
DownROI = [];
for ii = 1:nDown
    roi = tempPoints.down(ii,1):tempPoints.down(ii,2);
    DownROI = [DownROI,roi];
end
tempPoints.UpROI = UpROI;
tempPoints.DownROI = DownROI;