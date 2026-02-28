
function letter = Alphabet(number)
% letter = Alphabet(number)
% Get the letter of the alphabet that corresponds to the input number
% Useful for writing into Excel sheets
%
% INPUT
%   'number' : number(s) in the alphabet to find the corresponding letter(s)
%           Can be a single number or array of numbers (double)
%           Can also be a character to get the number (reverse lookup)
%                       e.g., Alphabet('C') --> 3
%
% OUTPUT
%   'letter' : letter(s) corresponding to the input number
%           Single number input --> character (e.g., Alphabet(1) --> 'A')
%           Array input --> cell array (e.g., Alphabet([1,3,5]) --> {'A','C','E'})
%           Character input --> number (e.g., Alphabet('A') --> 1)
%           Will double up letters as n increases beyond length of alphabet
%           e.g., Alphabet(30) --> 'AD'
%
% ES DICKINSON, Dec 2018
% Updated 2026 - added array input support

%%
alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
nLetters = length(alphabet);

% Handle character input (reverse lookup)
if ischar(number)
    letter = find(alphabet == number);
    return;
end

% Handle array input
if numel(number) > 1
    letter = cell(size(number));
    for i = 1:numel(number)
        letter{i} = Alphabet(number(i));  % Recursive call for each element
    end
    return;
end

% Handle single number input
if number > nLetters
    start = fix(number / nLetters);
    last = rem(number, nLetters);
    if last == 0 
        start = start - 1;
        last = nLetters;
    end
    letter = [alphabet(start) alphabet(last)];
else
    letter = alphabet(number);
end

%% OLD CODE pre 2026
% %%
% alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
% nLetters = length(alphabet);
% % determine if letter or number:
% if ischar(number)
%     letter = find(alphabet==number);
% else
%     if number > nLetters
%         start = fix(number/nLetters);
%         last = rem(number,nLetters);
%         letter = [alphabet(start) alphabet(last)];
%     else
%         letter = alphabet(number);
%     end
% end



