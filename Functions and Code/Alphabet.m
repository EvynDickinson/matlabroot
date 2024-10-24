function output = Alphabet(input)
% output = Alphabet(input)
% Get the letter of the alphabet that corresponds to the input number
% Useful for writting into Excel sheets
% Input:
% 'input' [number of the letter in the alphabet, as double]
% 'input' [letter in the alphabet, as character]
% Output:
% 'output' [letter, in char | number, in double]
% 
% 
% ES Dickinson, University of Washington, Dec 2018

alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
% determine if letter or number:
if ischar(input)
    output = find(alphabet==input);
else
    if input > 26
        start = fix(input/26);
        last = rem(input,26);
        output = [alphabet(start) alphabet(last)];
    else
        output = alphabet(input);
    end
end



