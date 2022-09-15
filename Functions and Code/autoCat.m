

function newData = autoCat(data,newMat,rows,bottom,DimFlip)
%
% newData = autoCat(data,newMat,rows,bottom,DimFlip)
%
% automatically adds nans to the matrix to expand it to account for larger
% data set sizes 
%
% rows selects how the matrix should be concatenated:
% rows = true --> add rows to matrix --> [x ; y]
% rows = false --> add columns to matrix --> [x , y]
% bottom selects whether to add to the beginning or end of the matrix to make
% it all fit
% bottom = true --> add to the end (bottom or right) of the matrix
% bottom = false --> add to the beginning (top or left) of the matrix
% 
% DimFlip = true --> able to change the dimensions of the new matrix to
% match the old matrix (e.g., a column vector can become a row vector)
% 
% data is the current matrix
% newMat is the new data to be concatenated with 'data'
%
%
% Evyn Dickinson, Yale University 2022

% determine the dimensions:
dataDim = size(data);
addDim = size(newMat);

% default dimensions and locations if none provided
if ~exist('rows','var')
    rows = true;
end
if ~exist('bottom','var')
    bottom = true;
end

if ~exist('DimFlip','var')
    DimFlip = true; % vector orientation can be flipped
elseif ~any(addDim==1)
    DimFlip = false; % don't auto-flip matrices, only vectors
end

if all(dataDim==0)
    newData = newMat;
else
    switch rows
        case true % ADD COLUMNS
            % Size matches
            if dataDim(2)==addDim(2)
            elseif DimFlip && dataDim(2)==addDim(1)
               newMat = newMat';
            % case: size does not match
            elseif dataDim(2)>addDim(2) %new data is smaller than old data
               newMat = [newMat, nan(addDim(1),dataDim(2)-addDim(2))];
            elseif dataDim(2)<addDim(2) %existing data is smaller than new data
               data = [data, nan(dataDim(1),addDim(2)-dataDim(2))];
            else
                warndlg('Error in autoconcatenation')
                return
            end
            if bottom
                newData = [data; newMat];
            else
                newData = [newMat; data];
            end
        case false % ADD ROWS
            % Size matches
            if dataDim(1)==addDim(1)
            elseif DimFlip && dataDim(1)==addDim(2)
               newMat = newMat';
            % case: size does not match
            elseif dataDim(1)>addDim(1) %new data is smaller than old data
               newMat = [newMat; nan(dataDim(1)-addDim(1),addDim(2))];
            elseif dataDim(1)<addDim(1) %existing data is smaller than new data
               data = [data; nan(addDim(1)-dataDim(1),dataDim(2))];
            else
                warndlg('Error in autoconcatenation')
                return
            end
            if bottom
                newData = [data, newMat];
            else
                newData = [newMat, data];
            end

    end
end


























        