

function newData = autoCat(data,newMat,dim,loc)
% newData = autoCat(data,newMat,dim)
% automatically adds nans to the matrix to expand it to account for larger
% data set sizes 
%
% dim selects how the matrix should be concatenated:
% dim = 1 --> add rows to matrix
% dim = 2 --> add columns to matrix
% loc selects whether to add to the beginning or end of the matrix to make
% it all fit
% loc = 1 --> add to the beginning (top or left) of the matrix
% loc = 2 --> add to the end (bottom or right) of the matrix
%
% data is the current matrix
% newMat is the new data to be concatenated with 'data'
%
%
% Evyn Dickinson, Yale University 2022


