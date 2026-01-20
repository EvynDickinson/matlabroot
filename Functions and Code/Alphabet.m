
function letter = Alphabet(number)
% letter = Alphabet(number)
% Get the letter of the alphabet that corresponds to the input number
% Useful for writting into Excel sheets
% 
% INPUT
%   'number' : number in the alphabet to find the corresponding letter in
%           the alphabet. Input type is double
%
% OUTPUT
%   'output' : letter corresponding to the input number (in character)
%           e.g. : Alphabet(1) --> 'A'
%           will double up numbers as the n increases beyond length of alphabet
%           e.g. : Alphabet(30) --> 'AD'
% 
% ES DICKINSON, Dec 2018

%%
alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
nLetters = length(alphabet);
% determine if letter or number:
if ischar(number)
    letter = find(alphabet==number);
else
    if number > nLetters
        start = fix(number/nLetters);
        last = rem(number,nLetters);
        letter = [alphabet(start) alphabet(last)];
    else
        letter = alphabet(number);
    end
end



