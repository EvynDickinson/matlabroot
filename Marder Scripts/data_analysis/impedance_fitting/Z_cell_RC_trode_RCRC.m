function Z=f(theta,f)

R_m=theta(1);  % MOhm
C_m=theta(2);  % pF
R1_e=theta(3);  % MOhm
C1_e=theta(4);  % pF
R2_e=theta(5);  % MOhm
C2_e=theta(6);  % pF

[Z_x,Z_in,Z_out]=Z_x_trode_RCRC(R1_e,C1_e,R2_e,C2_e,f);
Y_C_m=1e-6*j*2*pi*C_m*f;  % uS
Z_m=R_m./(1+R_m*Y_C_m);  % MOhm
Z=Z_in-Z_x.^2./(Z_m+Z_out);
