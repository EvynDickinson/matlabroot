
function tempPoints = getTempTurnPoints(TempProtocolString)
% tPoints = getTempTurnPoints(tempList{ii})



switch TempProtocolString
    case 'double_eaton_cooling_warming_ramp'
        tempPoints.up = [5720,34628;...
                         64819,93739];
        tempPoints.down = [34629,64818;...
                           93740,141259];
end