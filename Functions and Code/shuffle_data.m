

 function v = shuffle_data(v)
 % v = shuffle(v);
 % shuffle the elements in the vector V
     v = v(randperm(length(v)));
 end