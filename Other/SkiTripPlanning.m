

%% Parameters (fixed costs)
clear

%  Ticket prices
adult_full_ticket = 129;
adult_half_ticket = 106;
senior_full_ticket = 111;
senior_half_ticket = 93;

groups = {'AF','AH','SF','SH'};
max_discount = 0.8;

% Pass price
ikon_pass = 829;
buddy_passes_per_pass= 8;
buddy_pass_discount = 0.75;

tax  = 0.89;
cost = struct;

%% Conditions : uniform ski days between all people
skier = struct;
cost = struct;
skier(1).name = 'Patsy';
skier(1).g = 'senior';
skier(2).name = 'Clarissa';
skier(2).g = 'adult';
skier(3).name = 'Evyn';
skier(3).g = 'adult';
skier(4).name = 'Ravella';
skier(4).g = 'adult';
skier(5).name = 'Matt';
skier(5).g = 'adult';


for cond = 1:9
    
    switch cond
        case 1 % no passes & 6 full days skiing all
            passholders = {' '};
            full_days = 6;
            half_days = 0;
            for ii = 1:length(skier)
                skier(ii).full_days = full_days;
                skier(ii).half_days = half_days;
            end
    
        case 2 % 1 pass & 6 full days skiing all
            passholders = {'Evyn'};
            full_days = 6;
            half_days = 0;
            for ii = 1:length(skier)
                skier(ii).full_days = full_days;
                skier(ii).half_days = half_days;
            end
        case 3 % 2 passes & 6 full days skiing all
            passholders = {'Evyn','Patsy'};
            full_days = 6;
            half_days = 0;
            for ii = 1:length(skier)
                skier(ii).full_days = full_days;
                skier(ii).half_days = half_days;
            end
        case 4 % 3 passes & 6 full days skiing all
            passholders = {'Evyn','Patsy','Ravella'};
            full_days = 6;
            half_days = 0;
            for ii = 1:length(skier)
                skier(ii).full_days = full_days;
                skier(ii).half_days = half_days;
            end
        case 5 % no passes & 5 full days skiing all
            passholders = {' '};
            full_days = 5;
            half_days = 1;
            for ii = 1:length(skier)
                skier(ii).full_days = full_days;
                skier(ii).half_days = half_days;
            end
    
        case 6 % 1 pass & 5 full days skiing all
            passholders = {'Evyn'};
            full_days = 5;
            half_days = 1;
            for ii = 1:length(skier)
                skier(ii).full_days = full_days;
                skier(ii).half_days = half_days;
            end
        case 7 % 2 passes & 5 full days skiing all
            passholders = {'Evyn','Patsy'};
            full_days = 5;
            half_days =1;
            for ii = 1:length(skier)
                skier(ii).full_days = full_days;
                skier(ii).half_days = half_days;
            end
        case 8 % 3 passes & 5 full days skiing all
            passholders = {'Evyn','Patsy','Ravella'};
            full_days = 5;
            half_days = 1;
            for ii = 1:length(skier)
                skier(ii).full_days = full_days;
                skier(ii).half_days = half_days;
            end
%         case 9 % 3 passes & 5 full days skiing all
%             passholders = {'Evyn','Patsy','Ravella'};
%             full_days = 5;
%             half_days = 1;
%             for ii = 1:length(skier)
%                 skier(ii).full_days = full_days;
%                 skier(ii).half_days = half_days;
%             end
    end

    nskiers = size(skier,2);
    npasses = size(passholders,2);

    % Ticket counts calculation
    [adult_full, adult_half, senior_full, senior_half, passes] = deal(0);
    for i = 1:nskiers
        % ski pass holder ID
        if any(strcmpi(skier(i).name,passholders))
            skier(i).full_days = 0;
            skier(i).half_days = 0;
            skier(i).pass = 1;
            passes = passes + 1;
        else
            skier(i).pass = 0;
        end
        % ticket count
        switch skier(i).g
            case 'adult'
                adult_full = adult_full + skier(i).full_days;
                adult_half = adult_half + skier(i).half_days;
            case  'senior'
                senior_full = senior_full + skier(i).full_days;
                senior_half = senior_half + skier(i).half_days;
        end
    end

    % Ticket distribution
    nBuddyTix = passes*buddy_passes_per_pass;

    % distribute tickets to each group according to savings
    [n_AF_buddy, n_SF_buddy, n_AH_buddy, n_SH_buddy] = deal(0);
    while nBuddyTix>0
        % Adult full 
        if adult_full>0
            n_AF_buddy = n_AF_buddy+1;
            adult_full = adult_full-1;
            nBuddyTix = nBuddyTix-1;
            continue
        end
        % Senior full
        if senior_full>0
            n_SF_buddy = n_SF_buddy+1;
            senior_full = senior_full-1;
            nBuddyTix = nBuddyTix-1;
            continue
        end
        % Adult half
         if adult_half>0
            n_AH_buddy = n_AH_buddy+1;
            adult_half = adult_half-1;
            nBuddyTix = nBuddyTix-1;
            continue
         end
        % Senior half
        if senior_half>0
            n_SH_buddy = n_SH_buddy+1;
            senior_half = senior_half-1;
            nBuddyTix = nBuddyTix-1;
            continue
        end

        if sum(senior_half+adult_half+senior_full+adult_full)==0
            break
        end
    end

    % -- Cost calculations --- :
    %buddy tickets
    cost(cond).n_AF_buddy = n_AF_buddy;
    cost(cond).AF_buddy = n_AF_buddy*(adult_full_ticket*buddy_pass_discount)*tax;
    cost(cond).n_AH_buddy = n_AH_buddy;
    cost(cond).AH_buddy = n_AH_buddy*(adult_half_ticket*buddy_pass_discount)*tax;
    cost(cond).n_SF_buddy = n_SF_buddy;
    cost(cond).SF_buddy = n_SF_buddy*(senior_full_ticket*buddy_pass_discount)*tax;
    cost(cond).n_SH_buddy = n_SH_buddy;
    cost(cond).SH_buddy = n_SH_buddy*(senior_half_ticket*buddy_pass_discount)*tax;

    %discount tickets
    cost(cond).n_AF_reduced = adult_full;
    cost(cond).AF_reduced = adult_full*(adult_full_ticket*max_discount)*tax;
    cost(cond).n_AH_reduced = adult_half;
    cost(cond).AH_reduced = adult_half*(adult_half_ticket*max_discount)*tax;
    cost(cond).n_SF_reduced = senior_full;
    cost(cond).SF_reduced = senior_full*(senior_full_ticket*max_discount)*tax;
    cost(cond).n_SH_reduced = senior_half;
    cost(cond).SH_reduced = senior_half*(senior_half_ticket*max_discount)*tax;

    %full price tickets
    cost(cond).n_AF_full = adult_full;
    cost(cond).AF_full = adult_full*(adult_full_ticket)*tax;
    cost(cond).n_AH_full = adult_half;
    cost(cond).AH_full = adult_half*(adult_half_ticket)*tax;
    cost(cond).n_SF_full = senior_full;
    cost(cond).SF_full = senior_full*(senior_full_ticket)*tax;
    cost(cond).n_SH_full = senior_half;
    cost(cond).SH_full = senior_half*(senior_half_ticket)*tax;

    %passes
    cost(cond).n_passes = passes;
    cost(cond).passes = passes*ikon_pass;

    %total costs:
    n = 0;
    for g = 1:length(groups)
        n = n+cost(cond).([groups{g} '_buddy']);
    end
    cost(cond).buddy_tix = n;

    n = 0;
    for g = 1:length(groups)
        n = n+cost(cond).([groups{g} '_reduced']);
    end
    cost(cond).reduced_tix = n;

    n = 0;
    for g = 1:length(groups)
        n = n+cost(cond).([groups{g} '_full']);
    end
    cost(cond).full_tix = n;

    cost(cond).min_total = cost(cond).passes + cost(cond).buddy_tix + cost(cond).reduced_tix;
    cost(cond).max_total = cost(cond).passes + cost(cond).buddy_tix + cost(cond).full_tix;

end


%% Visualize the cost scenarios:


fig = getfig('',1);
[data,x] = deal([]);
buff = 0.15;
idx = 0;
for cond = 1:length(cost)
    
    idx = idx + 1;
    x_max(cond) = cond+buff;
    data(idx,:) = [cost(cond).reduced_tix, cost(cond).buddy_tix, cost(cond).passes];
    x(idx) = cond-buff;
    
    idx = idx + 1;
    data(idx,:) = [cost(cond).full_tix, cost(cond).buddy_tix, cost(cond).passes];
    x(idx) = cond+buff;
end
h = bar(x, data,'stacked'); 

formatFig(fig,true);
ylabel('$      ','fontsize',30,'Rotation',0)
set(gca,'xticklabel',[],'xtick',[])

h(1).FaceColor = Color('teal');%tickets
h(2).FaceColor = Color('orange'); %buddy tix
h(3).FaceColor = Color('lavender'); %passes











    



















