function Y=f(theta,f)

Rin=theta(1);  % input resistance, MOhm
rho=theta(2);  % dendritic dominace == R_s/R_d, unitless
L=theta(3);  % electrotonic length, unitless
tau_d=theta(4);  % dendritic time constant, ms
epsilon=theta(5);  % == tau_soma/tau_d, unitless
R1=theta(6);  % MOhm
C1=theta(7);  % pF
R2=theta(8);  % MOhm
C2=theta(9);  % pF

s=j*2*pi*f;  % rad/s
s_term=sqrt(1+1e-3*tau_d*s);
Y_m=(1+epsilon*1e-3*tau_d.*s+rho*coth(L)*s_term.*tanh(L*s_term)) / ...
                                                     (Rin*(rho+1)) ;  % uS
Y_C1=1e-6*C1*s;  % uS
Y_C2=1e-6*C2*s;  % uS
U_1=1+R1*Y_C1;  % unitless
U_2=1+R2*Y_C2;  % unitless
Y=U_2.*Y_C1+U_1.*Y_C2+U_1.*U_2.*Y_m;  % uS
%Z=1./Y;  % MOhms
