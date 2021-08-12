

function [fig_handle, stance] = DLC_plot_swingNstance(data,Center,Radius,threshold,SZ)
% 
% [fig_handle, stance] = DLC_plot_swingNstance(data,Center,Radius,threshold,SZ)
% OR
% [fig_handle, stance] = DLC_plot_swingNstance(stance)
% already generated logical of stance|swing
%
% Plot swing/stance from coordinate data for given legs
% 
% ES Dickinson,
% University of Washington, 2020

if nargin > 1
    % check for stance height thresholds
    try a = threshold<1 & threshold>=0;
        if a == false
            threshold = 0.02; %2percent threshold
        end
    catch
        threshold = 0.02; %2percent threshold
    end
    % check for user plot size:
    try 
        h = SZ(1);
        w = SZ(2);
    catch    
        h = 300;    % default base height in pixels
        w = 500;    % default base width in pixels
    end

    % find the euclidian distance for each leg in each frame:
    for ii = 1:length(data)
        dist(ii,:) = pdist2(Center,data(ii).raw); 
    end
    % find swing vs stance
    stance = dist<Radius+(Radius*threshold); %1=in stance
% previously generated logical index of swing and stance:
elseif nargin == 1
    stance = data;
    h = 300;    % default base height in pixels
    w = 500;    % default base width in pixels
end


% image info:
ncol = size(stance,2);              % nunber of frames
nrows = size(stance,1);             % number of legs
uh = round(h/nrows);                % unit height
uw = round(w/ncol);                 % unit width
swing_unit = zeros(uh,uw);          % swing black block
stance_unit = ones(uh,uw);          % stance white block

% build image:
blankIM = [];    
for jj = 1:nrows % for each leg
    temp = [];
    for ii = 1:ncol % for each frame
        switch stance(jj,ii)
            case true
                a = stance_unit;
            case false
                a = swing_unit;
        end
        temp = [temp, a];
    end
    blankIM = [blankIM; temp];
end
IM = repmat(blankIM,[1,1,3]); % B&W to color

% Display
fig_handle = getfig('',1);
image(IM);


end

