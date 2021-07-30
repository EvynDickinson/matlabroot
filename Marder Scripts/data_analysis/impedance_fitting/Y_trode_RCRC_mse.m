function mse=f(theta,f,Y_true)

Y=Y_trode_RCRC(theta,f);
res=abs(Y-Y_true);
mse=mean(res.^2);
