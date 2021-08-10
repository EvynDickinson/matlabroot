   
clear all; close all; clc 
strtpoint = 470; %line in excel file to start on

% Find files from Excel:     
[excelfile, Excel] = load_flysummary;

% Pull out the structure names from the file:
structure_names.excelfile = excelfile(strtpoint:end, Excel.structure);
ind = 1;
for ii = 1:length(structure_names.excelfile)
    if ischar(structure_names.excelfile{ii})
        structure_names.text{ind} = (structure_names.excelfile{ii});
        ind = ind+1;
    end
end; clear ind

%find the unique structure names:
structure_names.unique = unique(structure_names.text);
indx = listdlg('ListString', structure_names.unique, 'SelectionMode', 'Multiple', 'ListSize', [300, 400]);

%for each structure selected, load the data corresponding to that fly
for ii = 1:length(indx)
   %FIND STRUCTURE DETAILS
    structure_name = structure_names.unique{indx(ii)};
    location = find(strcmpi(structure_name, excelfile(strtpoint:end,Excel.new_struct_name)))+strtpoint-1;
    bb = 1;
    for tt = 1:length(location)
        a = cell2mat(excelfile(location(tt), Excel.structurenum));
        if ischar(a) %no value for the structure
        else %has a number in the structure
            flynum(bb) = a; %find total number of flies
            Location.rownum(bb) = location(tt); %get the rownumber for each fly
            bb = bb+1; %increase the index for the number of flies
        end
    end
    num.flies = max(flynum); clear flynum bb a tt
    % Error Check
    if sum(num.flies) == 0
        warndlg('Check name of flies in Excel File, no matching files found')
        return
    end
    
    %LOAD DATA INTO MATLAB
    % begin loop to load data from each fly in 'fly.num'
    for kk = 1:num.flies

        % Find Specific Fly & Load Data:
        fprintf('\n Fly: \n')
        disp([excelfile(Location.rownum(kk), Excel.date) excelfile(Location.rownum(kk), Excel.flynum)])
        
        % Date and number of the fly that matches the desired one
        %get row number of the fly(kk) 1_0 fly:
        %gives a logical index of NaN = 1
        NaN_index = cellfun(@(excelfile) any(isnan(excelfile(:,Excel.date))), excelfile(:,Excel.date));
        date_rownum = Location.rownum(kk);
%         while NaN_index(date_rownum) == 1
%             date_rownum = date_rownum-1;
%         end
        fly_date = string(excelfile(date_rownum, Excel.date));
        fly_num = string(excelfile(Location.rownum(kk), Excel.flynum));
        genetic_cross = string(excelfile(Location.rownum(kk), Excel.genetic_cross));
        %convert the number into fly ID date format, e.g 03282018 from 3.28.18
        c = strsplit(fly_date,'.');
        %month conversion:
        switch str2double(c(1))
            case num2cell(1:9)
                C{1} = ['0' c(1)];
            case num2cell(10:12)
                C{1} = c(1);
            otherwise
        end
        %day conversion:
        switch str2double(c(2))
            case num2cell(1:9)
                C{2} = ['0' c(2)];
            case num2cell(10:31)
                C{2} = c(2);
            otherwise
        end

        datestring = cell2mat([cell2mat(C{1}) cell2mat(C{2}) '20' c(3)]);
        fly_ID = [datestring '_fly' cell2mat(fly_num)];

        % Load analyzed data *can take a sec or two:
        folder_name = cell2mat(fly_date);
        try
            directory = ['D:\Evyn Data Files\' folder_name '\Fly ' char(fly_num) '\Analysis\'];
            load([directory fly_ID '.mat']);
        catch
            try
                directory = ['F:\Evyn Data Files\' folder_name '\Fly ' char(fly_num) '\Analysis\'];
                load([directory fly_ID '.mat']);
            catch
                try
                    directory = ['E:\FicTrac Raw Data\' folder_name '\Analysis\'];
                    load([directory fly_ID '.mat']);
                catch
                    directory = ['C:\matlabroot\FicTrac Raw Data\' folder_name '\Analysis\'];
                    load([directory fly_ID '.mat']);
                end
            end
        end
        fprintf(['\n Loaded fly data: ' fly_ID '\n'])
        % Add a video category to the data
        if ~isfield(fly,'Video')
            fprintf('\n Adding Video field\n')
            Sstrt = 1;Sstp = 44;Cstrt = 46;Cstp = 61;
            for cond = 1:28
                for rep = 1:3
                    fly.Video.data(cond, rep).raw = ...
                        [fly.Control.data(cond,rep).raw(Cstrt:Cstp,:); fly.Stim.data(cond,rep).raw(Sstrt:Sstp,:)];

                end
                fly.Video.speed(cond).data = ...
                    [fly.Control.speed(cond).data(Cstrt:Cstp,:); fly.Stim.speed(cond).data(Sstrt:Sstp,:)];
                fly.Video.rotvelocity(cond).data = ...
                    [fly.Control.rotvelocity(cond).data(Cstrt:Cstp,:); fly.Stim.rotvelocity(cond).data(Sstrt:Sstp,:)];
            end
            % Save the updated fly video:
            newfile = [directory, fly_ID];
            save(newfile, 'fly');
        end
         
        % Copy the analysis folder to google drive:
        % Destination:
        dest = ['G:\My Drive\Data\FicTrac Raw Data\', fly.param.folder_date, '\Fly ', fly_ID(end-2:end), '\Analysis'];
        copyfile(directory, dest)
        
        %ADD fly to FLY structure:
        try
        a = length(FLY);
        FLY(a+1) = fly;
        catch %first fly added to structure
            FLY(1) = fly;
        end
        clear a C c fly
        fprintf([' Number: ' num2str(kk) '/' num2str(num.flies) '\n '])
    end
    
end
beep
% Save Structure:
switch questdlg(['Save flies to structure: ' structure_name '?'], 'Save Data', 'Yes', 'Different Name', 'No', 'Yes')
    case 'Yes'
        save(structure_name,'FLY', 'structure_name', '-v7.3')
        fprintf('\n Structure saved \n')
    case 'Different Name'
        structure_name = cell2mat(inputdlg('Name to save?'));
        save(structure_name, 'FLY', 'structure_name')
        fprintf('\n Structure saved \n')
    case 'No'
        return
end
beep 
    
    
    

