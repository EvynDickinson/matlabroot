function mse=f(theta,f,Z_true)

Z=Z_cell_three_taus_two_trodes_RC_Cx(theta,f);
res=abs(log(Z)-log(Z_true));
mse=mean(res.^2);
