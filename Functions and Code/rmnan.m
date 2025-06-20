function [yA, yB] = rmnan(A, B, dim, bidirectional)
% [yA, yB] = rmnan(A, B, dim, bidirectional)
% Remove NaNs from vector A and corresponding entries in B.
% If bidirectional is true, also remove entries where B has NaNs.
% If A is a vector and B is a matrix, 'dim' specifies which dimension to filter.
%
% Inputs:
%   A - vector
%   B - (optional) vector or matrix
%   dim - (optional) dimension for filtering B if B is matrix (default = 2)
%   bidirectional - (optional) if true, remove where A or B has NaNs (default = false)
%
% Outputs:
%   yA - cleaned A
%   yB - cleaned B (only if B is provided)
%
% ES Dickinson (and ChatGPT)

% set across rows are the default direction if none selected and B is a
% matrix while A is a vector
if nargin < 3 || isempty(dim)
    dim = 2;
end
if nargin < 4
    bidirectional = false; % only match B to A
end

Aloc = isnan(A);
Bloc = false(size(Aloc)); 

if nargin > 1 && bidirectional
    if isvector(B) && isvector(A) % both A and B are vectors
        Bloc = isnan(B);
    elseif ~isvector(A) && ~isvector(B)
        Bloc = isnan(B);
    elseif isvector(A) && ~isvector(B)
        % A is vector, B is matrix
        if dim == 1
            Bloc = any(isnan(B), 2);
        else
            Bloc = any(isnan(B), 1);
        end
        Bloc = Bloc(:); % ensure shape matches A
    end
end

try nanLoc = Aloc | Bloc;
catch
   warndlg('Mismatching vectors A and B for removing nans')
end

yA = A(~nanLoc);

if nargin > 1 
    if isvector(B)
        yB = B(~nanLoc);
    elseif isvector(A)
        if dim == 1
            yB = B(~nanLoc, :);
        else
            yB = B(:, ~nanLoc);
        end
    else
        error('If A is not a vector, B must be a vector or dim must be specified properly.');
    end
else
    yB = [];
end


end
