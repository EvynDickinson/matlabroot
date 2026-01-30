

 function v = shuffle_data(v)
 % v = shuffle_data(v);
 % 
 % PURPOSE
 % shuffle the elements in the vector V 
 % (useful for shuffling x-values in a scatter plot)
 %
 % INPUTS
 %   'v' : vector of values to be shuffled
 %
 % OUTPUTS
 % 'v' : shuffled data vector
 %
 % EXAMPLE
 %    v = shuffle_data(1:10)
 %    --> v = [ 6  3  7 8  5  1  2  4  9  10]
 %
 % ES DICKINSON, 2022

%%
 
 v = v(randperm(length(v)));


 