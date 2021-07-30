function Z=f(theta,f)

R1=theta(1);  % MOhm
C1=theta(2);  % pF
R2=theta(3);  % MOhm
C2=theta(4);  % pF
Y_C1=1e-6*j*2*pi*f*C1;  % uS
Y_C2=1e-6*j*2*pi*f*C2;  % uS
Z=(R1+R2+R1*R2*Y_C2)./(1+R2*Y_C2+R2*Y_C1+R1*Y_C1+R1*R2*(Y_C1.*Y_C2));

