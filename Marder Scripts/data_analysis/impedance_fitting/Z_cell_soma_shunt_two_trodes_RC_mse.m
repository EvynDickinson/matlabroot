function mse=f(theta,f,Z_true)

Z=Z_cell_soma_shunt_two_trodes_RC(theta,f);
res=abs(Z-Z_true);
mse=mean(res.^2);
