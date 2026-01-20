

 function v = shuffle_data(v)
 % v = shuffle_data(v);
 % shuffle the elements in the vector V 
 % (useful for shuffling x-values in a scatter plot)
 %
 % INPUTS
 %   'v' : vector of values to be shuffled
 %
 % OUTPUTS
 % 'v' : shuffled data vector
 %
 % ES DICKINSON, 2022

%%
 
 v = v(randperm(length(v)));


 end