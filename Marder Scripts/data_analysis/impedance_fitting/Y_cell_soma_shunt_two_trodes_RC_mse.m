function mse=f(theta,f,Y_true)

Y=Y_cell_soma_shunt_two_trodes_RC(theta,f);
res=abs(Y-Y_true);
mse=mean(res.^2);
