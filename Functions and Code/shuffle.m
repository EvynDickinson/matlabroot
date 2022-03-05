

 function v = shuffle(v)
 % v = shuffle(v);
 % shuffle the elements in the vector V
     v = v(randperm(length(v)));
 end