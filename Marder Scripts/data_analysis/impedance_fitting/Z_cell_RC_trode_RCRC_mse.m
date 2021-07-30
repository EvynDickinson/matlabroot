function mse=f(theta,f,Z_true)

Z=Z_cell_RC_trode_RCRC(theta,f);
res=abs(Z-Z_true);
mse=mean(res.^2);
