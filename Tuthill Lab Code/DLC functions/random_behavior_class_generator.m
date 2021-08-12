
% num.flies = 8;
% num.conds = 20;
% num.reps = 6;

structure_name = '';

behavior_list = {'walking', 'stationary', 'grooming', 'other'};
phase_list = {'stance', 'swing'};

for ifly = 1:num.flies
    for cond = 1:num.conds
        for rep = 1:num.reps
           ii = randi([1,4],1,1);
           II = randi([1,2],1,1);
           group(ifly).behavior{cond,rep} = behavior_list{ii};
           group(ifly).phase{cond,rep} = phase_list{II};
        end
    end
end

save(['C:\matlabroot\behavior class\' structure_name ' behavior class'], 'group')
    



