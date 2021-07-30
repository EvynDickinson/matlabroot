function [Y_x,Y_in,Y_out]=f(R,C,f)

% R's in Mohms, C's in pF, f in Hz
% Y's are in uS

Y_C=1e-6*j*2*pi*f*C;  % uS
Y_x=repmat(-1/R,size(f));
Y_out=repmat(+1/R,size(f));
Y_in=1/R+Y_C;
