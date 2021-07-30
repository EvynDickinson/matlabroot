function Z=f(theta,f)

R=theta(1);  % MOhm
C=theta(2);  % pF
tau=1e-6*R*C;  % s
Z=R./(1+j*2*pi*tau*f);
