

 function v = shuffle_data(v)
 % v = shuffle_data(v);
 % shuffle the elements in the vector V
     v = v(randperm(length(v)));
 end