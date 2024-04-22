
filelist = dir('G:\My Drive\Jeanne Lab\Data');
names = {filelist(:).name};
checklist = [];
for i = 1:length(names)
    if strfind(names{i}, '2021')==7
        checklist(i) = true;
    elseif strfind(names{i}, '2022')==7
        checklist(i) = true;
    else
        checklist(i) = false;
    end
end

sum(checklist)



% do the cross comparision

listdirs