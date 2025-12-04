
function output = pullGroupInfo(expGroup)
% output = pullGroupInfo(expGroup)
%
% includes general features: 
% --------------------------------------
% *** expOrder -- the preferred plotting order of the subgroups
% *** colors -- the colors for the groups based on their expOrder 
%
% includes stats comparison for temp tuning curves: 
% ------------------------------------------------------------------------
% *** comp_pairs -- what pairs of data should be compared 
%       for temp tuning curve stats
% *** required_comp_pairs -- which groups must all be 
%       significant for test data significance
% *** nComp -- how many pairs of stats paired comparisions there are
% 


% Color selections
switch expGroup
   case 'Berlin LTS 15-35 caviar vs empty'
        expOrder = 1:2;
        colors = {'DodgerBlue', 'Gray'};

     case 'Berlin LTS 15-35 plate 1 vs plate 2'
        expOrder = 1:4;
        colors = {'Tomato', 'Dodgerblue', 'Peachpuff', 'Powderblue'};  

    case 'Berlin LTS 15-35 intact vs no antenna no food'
        expOrder = 1:2;
        colors = {'dodgerblue', 'peachpuff'};
        comp_pairs = [1 2 2]; % n1 & n2 : exp groups to compare, n3: color to plot in
        required_comp_pairs = [1, 1, 2]; % comp pairs (n1 n2; above) that must both be true for the test significance in n3

    case 'Berlin F LRR 25-17 caviar intact vs wax vs hold'
        expOrder = 1:4;
        colors = {'Dodgerblue', 'Tomato', 'Grey', 'Grey'};

    case 'Berlin F LRR caviar mechanical removal comparisons'
        expOrder = [1,2,3]; % berlin, no antenna, no arista
        colors = {'Grey', 'MediumSlateBlue', 'Gold'};
        comp_pairs = [1 2 2; 1 3 3]; % n1 & n2 : exp groups to compare, n3: color to plot in
        required_comp_pairs = [1, 1, 2; 2, 2, 3]; % comp pairs (n1 n2; above) that must both be true for the test significance in n3

    case 'Berlin temp rate caviar'
        expOrder = [5, 3, 2, 1, 4]; % slow to fast
        colors = {'Deeppink','Gold','MediumSpringGreen','mediumslateblue', 'dodgerblue'};

    case 'Berlin LTS 15-35 no food mechanical removal comparisons'
        expOrder = [1,2,3]; % berlin, no antenna, no arista
        colors = {'Grey', 'MediumSlateBlue', 'Gold'};
        comp_pairs = [1 2 2; 1 3 3]; % n1 & n2 : exp groups to compare, n3: color to plot in
        required_comp_pairs = [1, 1, 2; 2, 2, 3]; % comp pairs (n1 n2; above) that must both be true for the test significance in n3

    case 'Berlin temperature hold caviar'
        expOrder = [1 2 3 8 4 5 9 6 7]; %out of order bc 23C and 30C were added later
        % colors = {'Dodgerblue', 'Tomato', 'Grey', 'Grey'};

    case 'Berlin temperature holds' % includes food and empty trials
        expOrder = 1:num.exp;
        kolor1 = Color('Blue', 'LightGrey', 5); % ultimately 5 (with all temp trials added)
        kolor2 = Color('LightGrey', 'Red', 5);  % ultimately 5 (with all temp trials added)
        colors = nan([18, 3]); % ultimately 18x3 (with all temp trials added)
        colors(1:2:end,:) = [kolor1(1:end-1,:); kolor2];
        colors(2:2:end,:) = [kolor1(1:end-1,:); kolor2];
        % figure; scatter(1:num.exp, 1:num.exp, 50, colors(1:num.exp,:), 'filled'); % color test

    case 'Berlin Temp Holds Caviar' % includes food and empty trials
        expOrder = [1,2,3,8,4,5,9,6,7];
        kolor1 = Color('Blue', 'LightGrey', 5); % ultimately 5 (with all temp trials added)
        kolor2 = Color('LightGrey', 'Red', 5);  % ultimately 5 (with all temp trials added)
        colors = nan([18, 3]); % ultimately 18x3 (with all temp trials added)
        colors(1:2:end,:) = [kolor1(1:end-1,:); kolor2];
        colors(2:2:end,:) = [kolor1(1:end-1,:); kolor2];
        % figure; scatter(1:num.exp, 1:num.exp, 50, colors(1:num.exp,:), 'filled'); % color test
        
    case {'TrpA1-Gal4 x UAS-Kir2.1 LTS 15-35 no food','TrpA1-gal4 LTS 15-35 no food'}
        expOrder = 1:3; % UAS control, GAL4 control, GAL4>UAS
        colors = {'LightPink', 'HotPink', 'DodgerBlue'};
        comp_pairs = [3 1 1; 3 2 2]; % n1 & n2 : exp groups to compare, n3: color to plot in
        required_comp_pairs = [1, 2, 3]; % comp pairs (n1 n2; above) that must both be true for the test significance in n3

    case {'TrpA1-Gal4 x UAS-Kir2.1_A1 LTS 15-35 caviar','TrpA1-Gal4 x UAS-Kir2.1_A1 LTS 15-35 caviar plate 1'}
        expOrder = 1:3; % UAS control, GAL4 control, GAL4>UAS
        colors = {'LightPink', 'HotPink', 'DodgerBlue'};
        comp_pairs = [3 1 1; 3 2 2]; % n1 & n2 : exp groups to compare, n3: color to plot in
        required_comp_pairs = [1, 2, 3]; % comp pairs (n1 n2; above) that must both be true for the test significance in n3 
    
    case 'TrpA1-gal4 x Kir2.1 no antenna LTS 15-35 no food' 
        expOrder = [1,2]; % intact, antenna removed
        comp_pairs = [2 1 1]; % n1 & n2 : exp groups to compare, n3: color to plot in
        required_comp_pairs = [1, 1, 2]; % comp pairs (n1 n2; above) that must both be true for the test significance in n3 
        colors = {'Grey', 'DodgerBlue'};

    case 'TrpA1-Gal4 x UAS-Kir2.1_A1  LTS 15-35 no food mechanical removal comparisons' 
        expOrder = 1:2; % intact vs antenna-less
        colors = {'Grey', 'DodgerBlue'};
        comp_pairs = [2, 1, 1]; % n1 & n2 : exp groups to compare, n3: color to plot in
        required_comp_pairs = [1, 1, 2]; % comp pairs (n1 n2; above) that must both be true for the test significance in n3 

   case 'TrpA1-Gal4 x UAS-Kir2.1_A1  LTS 15-35 food vs no food' 
        expOrder = 1:2; % caviar vs no food
        colors = {'DodgerBlue', 'Grey'};
         comp_pairs = [2 1 1]; % n1 & n2 : exp groups to compare, n3: color to plot in
        required_comp_pairs = [1, 1, 1]; % comp pairs (n1 n2; above) that must both be true for the test significance in n3 

    case 'Berlin vs UAS-Kir2.1 caviar background comparison'
        expOrder = 1:num.exp;
        colors = {'turquoise', 'DarkOrchid'};
    
    case 'IR25a-gal4 x Kir2.1 and controls LTS 15-25 caviar'
        expOrder = 1:3; % UAS, GAL, UAS-GAL
        colors = {'LightPink', 'HotPink', 'DodgerBlue'};
        comp_pairs = [3 1 1; 3 2 2]; % n1 & n2 : exp groups to compare, n3: color to plot in
        required_comp_pairs = [1, 2, 3]; % comp pairs (n1 n2; above) that must both be true for the test significance in n3

    case 'IR25a-gal4 x Kir2.1 and controls F LRR caviar'
        expOrder = 1:3; % UAS, GAL, UAS-GAL
        colors = {'LightPink', 'HotPink', 'DodgerBlue'};
        comp_pairs = [3 1 1; 3 2 2]; % n1 & n2 : exp groups to compare, n3: color to plot in
        required_comp_pairs = [1, 2, 3]; % comp pairs (n1 n2; above) that must both be true for the test significance in n3
    
    case 'R77C10-gal4 x Kir2.1 and controls F LRR 17-25 caviar'
        expOrder = 1:3; % gal4 control, UAS control, 
        colors = {'LightPink', 'HotPink', 'DodgerBlue'};
    
    case 'IR40a_TM2-gal4 x Kir2.1 and controls F LRR 17-25 caviar'
        expOrder = 1:3; % gal4 control, UAS control, 
        colors = {'LightPink', 'HotPink', 'DodgerBlue'};
    
    case 'IR40a_TM2-gal4 x Kir2.1 no arista LTS 15-35 empty'
        expOrder = 1:3;
        colors = {'magenta', 'grey', 'dodgerblue'};
        
    case {'Orco-Gal4 x UAS-Kir2.1_A1 LTS 15-35 caviar','Orco-Gal4 x UAS-Kir2.1_A1 LTS 15-35 caviar plate 1'} %no uas control yet
        % expOrder = 1:3; % UAS control, GAL4 control, GAL4>UAS
        % colors = {'LightPink', 'HotPink', 'DodgerBlue'};
        % comp_pairs = [3 1 1; 3 2 2]; % n1 & n2 : exp groups to compare, n3: color to plot in
        % required_comp_pairs = [1, 2, 3]; % comp pairs (n1 n2; above) that must both be true for the test significance in n3 
        expOrder = 1:2; % GAL4 control, GAL4>UAS
        colors = {'LightPink', 'DodgerBlue'};
        comp_pairs = [2 1 1]; 
        required_comp_pairs = [1 3]; 
end

% set default experiment pairs to be 1 and 2 for statistical significance
if exist('comp_pairs','var') 
    nComps = size(comp_pairs,1);
    if ~exist('required_comp_pairs','var') 
        required_comp_pairs = [];
    end
    %  save statistical comparison pairs to the output function
    output.required_comp_pairs = required_comp_pairs;
    output.comp_pairs = comp_pairs;
    output.nComps = nComps;
end

% test that the group struct exists in the switch: 
if exist('expOrder', 'var')
    output.expOrder = expOrder;
    output.color = colors;
else 
    disp('Group structure not yet in ''pullGroupInfo.m''')
    output = [];
end


%     case 'Wild Types S LRR 23-15 caviar'
%         expOrder = 1:3; % UAS, GAL, UAS-GAL
%         colors = {'LightPink', 'HotPink', 'DodgerBlue'};
%         comp_pairs = [3 1 1; 3 2 2]; % n1 & n2 : exp groups to compare, n3: color to plot in
%         required_comp_pairs = [1, 2, 3]; % comp pairs (n1 n2; above) that must both be true for the test significance in n3
% 
% % lat_list = {'Swedish', 'Berlin', 'Oregon','Canton','Malawi', 'Zimbabwe'};
% % latitudes = [60.1282,  52.5200,    43.8041,   40.7989, -13.2543, -19.0154];..
%                    
