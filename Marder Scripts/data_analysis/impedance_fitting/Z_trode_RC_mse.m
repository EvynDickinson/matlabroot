function mse=f(theta,f,Z_true)

Z=Z_trode_RC(theta,f);
res=abs(Z-Z_true);
mse=mean(res.^2);