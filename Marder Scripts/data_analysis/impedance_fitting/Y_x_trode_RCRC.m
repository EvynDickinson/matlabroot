function [Y_x,Y_in,Y_out]=f(R1,C1,R2,C2,f)

% R's in Mohms, C's in pF, f in Hz
% Y's are in uS

Y_C1=1e-6*j*2*pi*C1*f;  % uS
Y_C2=1e-6*j*2*pi*C2*f;  % uS
Y_x=1./(R1+R2+R1*R2*Y_C2);
Y_in=(1+R2*Y_C2+R2*Y_C1+R1*Y_C1+R1*R2*Y_C1.*Y_C2).*Y_x;
Y_out=(1+R1*Y_C2).*Y_x;
