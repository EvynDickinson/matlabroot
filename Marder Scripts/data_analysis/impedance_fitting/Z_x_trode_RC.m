function [Z_x,Z_in,Z_out]=f(R,C,f)

% R's in Mohms, C's in pF, f in Hz
% Y's are in uS

Z_C=1./(1e-6*j*2*pi*f*C);  % MOhm
Z_x=Z_C;
Z_in=Z_C;
Z_out=Z_C+R;
