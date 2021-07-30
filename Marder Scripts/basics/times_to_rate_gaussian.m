function [t,r]=f(dt,T,t0,ts,sigma_t)

N=round(T/dt);
T=N*dt;  % make sure dt and T are consistent
t=(t0:dt:t0+T-dt)';  % timeline
is=floor((ts-t0)/dt)+1;  % index into t of each spike
is(is>N)=N;  % in case there's a ts(i)==T
r_proto=zeros([N 1]);
r_proto(is)=1;
sigma_i=sigma_t/dt;
kernel=1/dt*gaussian_kernel_1d(sigma_i,2*ceil(4*sigma_i)+1);
r=convn(r_proto,kernel,'same');
