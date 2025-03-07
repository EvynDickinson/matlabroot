

function param = getTempParamsHighRes(protocol)

switch protocol
    case 'high_res_LTS_35-15'
        param.ntrans = 4;
        param.labelstr = 'click (1) the start of the ramp, (2) the peak, (3) the trough, and (4) the peak';
        % find the x-time value for each time period
        tRate = struct;
        tRate(1).name = 'start hold';
        tRate(1).color = Color('grey');
        tRate(2).name = 'warming';
        tRate(2).color = Color('red');
        tRate(3).name = 'cooling';
        tRate(3).color = Color('dodgerblue');
        tRate(4).name = 'warming';  
        tRate(4).color = Color('red');
        tRate(5).name = 'cooling';
        tRate(5).color = Color('dodgerblue');
        param.tRate = tRate;


    case 'courtship_F_LRR_25-17'
        param.ntrans = 3;
        param.labelstr = 'click (1) the start of the ramp, (2) the bottom, and (3) the end of the ramp';
        % find the x-time value for each time period
        tRate = struct;
        tRate(1).name = 'start hold';
        tRate(1).color = Color('grey');
        tRate(2).name = 'cooling';
        tRate(2).color = Color('dodgerblue');
        tRate(3).name = 'warming';
        tRate(3).color = Color('red');
        tRate(4).name = 'end hold';
        tRate(4).color = Color('grey');
        param.tRate = tRate;

end