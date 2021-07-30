function mse=f(theta,f,Z_true)

Z=Z_cell_RC_two_trodes_RC(theta,f);
res=abs(log(Z)-log(Z_true));
mse=mean(res.^2);
