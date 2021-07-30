function [Z_x,Z_in,Z_out]=f(R1,C1,R2,C2,f)

% R's in Mohms, C's in pF, f in Hz
% Y's are in uS

% There's surely a better way to compute this, but oh well

Y_C1=1e-6*j*2*pi*C1*f;  % uS
Y_C2=1e-6*j*2*pi*C2*f;  % uS
F_in=1+R2*Y_C2+R2*Y_C1+R1*Y_C1+R1*R2*Y_C1.*Y_C2;
F_out=(1+R1*Y_C2);
F_x=R1+R2+R1*R2*Y_C2;
Z_x=F_x./(1-F_in.*F_out);
Z_in=-F_out.*Z_x;
Z_out=-F_in.*Z_x;
