function Z=f(theta,f)

Rin1=theta(1);  % input resistance, MOhm
tau1=1e-3*theta(2);  % one time constant, ms->s
Rin2=theta(3);  % input resistance, MOhm
tau2=1e-3*theta(4);  % one time constant, ms->s
R1=theta(5);  % MOhm
C1=1e-6*theta(6);  % pF->uF
R2=theta(7);  % MOhm
C2=1e-6*theta(8);  % pF->uF
Cx=1e-6*theta(9);  % pF->uF

s=j*2*pi*f;  % rad/s
Z_m=Rin1./(1+tau1*s)+Rin2./(1+tau2*s);  % MOhm
Y_m=1./Z_m;  % uS

Y_C1=C1*s;  % uS
Y_C2=C2*s;  % uS
Y_Cx=Cx*s;  % uS

U_1=1+R1*Y_C1;  % unitless
U_2=1+R2*Y_C2;  % unitless

H11=R1.*Y_C2+(1+R1*Y_m).*U_2;
H12=-(R1+R2+R1*R2*Y_m);
H21=U_2.*Y_C1+U_1.*Y_C2+U_1.*U_2.*Y_m;  % uS
H22=-(R2*Y_C1+U_1.*(1+R2*Y_m));

Y=(H21+(H11-H22-2).*Y_Cx)./(1-H12.*Y_Cx);  % uS
Z=1./Y;  % MOhms
