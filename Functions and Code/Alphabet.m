function letter = Alphabet(input_num)
% letter = Alphabet(input_num)
% Get the letter of the alphabet that corresponds to the input number
% Useful for writting into Excel sheets
% Input:
% 'input_num' [number of the letter in the alphabet, as double]
% Output:
% 'letter' [letter, in char]
% 
% 
% ES Dickinson, University of Washington, Dec 2018
alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
if input_num > 22
    start = fix(input_num/22);
    last = rem(input_num,22);
    letter = [alphabet(start) alphabet(last)];
else
    letter = alphabet(input_num);
end

end