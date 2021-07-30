function Z=f(theta,f)

Rm=theta(1);  % MOhm
Cm=theta(2);  % pF
R1=theta(3);  % MOhm
C1=theta(4);  % pF
R2=theta(5);  % MOhm
C2=theta(6);  % pF

s=j*2*pi*f;  % rad/s
Y_m=1/Rm+1e-6*Cm*s;  % uS
Y_C1=1e-6*C1*s;  % uS
Y_C2=1e-6*C2*s;  % uS
U_1=1+R1*Y_C1;  % unitless
U_2=1+R2*Y_C2;  % unitless
Y=U_2.*Y_C1+U_1.*Y_C2+U_1.*U_2.*Y_m;  % uS
Z=1./Y;  % MOhms
