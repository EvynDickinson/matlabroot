
clear
clc

%% Display well contents (no changes, just display)
arenaIdx = {'A', 'B', 'C', 'D'};
A = [];
for arena = 1:4
    wellLabels = {parameters.(['Arena' arenaIdx{arena}]).well_1;...
                  parameters.(['Arena' arenaIdx{arena}]).well_2;...
                  parameters.(['Arena' arenaIdx{arena}]).well_3;...
                  parameters.(['Arena' arenaIdx{arena}]).well_4};   
    A = [A, wellLabels];
end
wellLabels = array2table(A,"VariableNames",["Arena A","Arena B","Arena C","Arena D"]);
disp(wellLabels) 

% Clear extra variables:
clear i ii ans wellLabels A arenaIdx arena

%% UPDATE THE NAMES IN FOOD WELLS
%Blank all of them
for i = 1:4
    for ii = 1:4
        parameters.(['Arena' Alphabet(i)]).(['well_' num2str(ii)]) = 'Empty';
    end
end

%Readd the appropriate food labels:
parameters.ArenaA.well_3 = 'caviar';
parameters.ArenaB.well_3 = 'caviar';
parameters.ArenaC.well_1 = 'caviar';
parameters.ArenaD.well_1 = 'caviar';

% disp well contents:
arenaIdx = {'A', 'B', 'C', 'D'};
A = [];
for arena = 1:4
    wellLabels = {parameters.(['Arena' arenaIdx{arena}]).well_1;...
                  parameters.(['Arena' arenaIdx{arena}]).well_2;...
                  parameters.(['Arena' arenaIdx{arena}]).well_3;...
                  parameters.(['Arena' arenaIdx{arena}]).well_4};   
    A = [A, wellLabels];
end
wellLabels = array2table(A,"VariableNames",["Arena A","Arena B","Arena C","Arena D"]);
disp(wellLabels) 

% Clear extra variables:
clear i ii ans wellLabels A arenaIdx arena

%% UPDATE THE GENOTYPE

%Blank all of them
for i = 1:4
        parameters.(['Arena' Alphabet(i)]).genotype = 'UAS-Kir2.1_B_F8';
        disp(['Arena ' Alphabet(i) ' genotype = ' parameters.(['Arena' Alphabet(i)]).genotype])
end

% %Custom (not all) updates to genotype the appropriate food labels:
% parameters.ArenaA.genotype = '0.5M_glucose_15%ACV';
% parameters.ArenaB.genotype = '0.5M_glucose_15%ACV';
% parameters.ArenaC.genotype = '0.5M_glucose_5%ACV';
% parameters.ArenaD.genotype = '0.5M_glucose_5%ACV';

% Clear extra variables:
clear i ii ans

%%

parameters.protocol = 'Large_temp_sweep_15_35_FPS6';

%% 

expData.parameters.protocol = 'Large_temp_sweep_15_35_FPS6';
