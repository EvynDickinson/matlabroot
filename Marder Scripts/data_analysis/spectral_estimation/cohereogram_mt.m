function [f,t,...
          x_mean,...
          P,...
          C_mag,C_phase,...
          x_mean_ci,...
          P_ci,...
          C_mag_ci,C_phase_ci,...
          P_taos]=f(dt,x,T_window,NW,K,F,p_FFT_extra,conf_level)

% On return:        
%   P is of shape [N_f n_windows Q]
%   C_mag is of shape [N_f n_windows Q Q]
%   C_phase is of shape [N_f n_windows Q Q]
%   f is of shape [N_f 1]
%   t is of shape [n_windows 1]
%   x_mean is of shape [n_windows Q]
%   x_mean_se is of shape [n_windows Q]
%   P_taos is of shape [N_f n_windows Q K]

% get dimensions
N_total=size(x,1);
Q=size(x,2);

% convert window size from time to elements 
N_window=round(T_window/dt);
if mod(N_window,2)~=0
  N_window=N_window-1;  % make even
end
n_windows_no_overlap=floor(N_total/N_window);
n_windows=2*n_windows_no_overlap-1;  % we do 50% overlaps

% truncate data so we have integer number of windows
N_total=n_windows_no_overlap*N_window;
x=x(1:N_total,:);

% figure out the dimensions of the output
% will be N_f x n_windows
N_fft=2^(ceil(log2(N_window))+p_FFT_extra);
df=1/(N_fft*dt);
N_f=ceil(F/df)+1;
F=(N_f-1)*df;  % make consistent

% set the time base (these are the centers of each window)
t=dt*(N_window/2)*(1:n_windows)';

% do it
x_mean=zeros(n_windows,Q);
P=zeros(N_f,n_windows,Q);
C_mag=zeros(N_f,n_windows,Q,Q);
C_phase=zeros(N_f,n_windows,Q,Q);
x_mean_ci=zeros(n_windows,Q,2);
P_ci=zeros(N_f,n_windows,Q,2);
C_mag_ci=zeros(N_f,n_windows,Q,Q,2);
C_phase_ci=zeros(N_f,n_windows,Q,Q,2);
P_taos=zeros(N_f,n_windows,Q,K);
for j=1:n_windows
  i_start=(j-1)*(N_window/2)+1;
  i_end=i_start+N_window-1;
  x_this=x(i_start:i_end,:);
  x_this_mean=mean(x_this,1);
  x_this_cent=x_this-repmat(x_this_mean,[N_window 1]);
  [f,...
   P_this,...
   C_mag_this,C_phase_this,...
   N_fft,f_res_diam,...
   P_ci_this, ...
   C_mag_ci_this, ...
   C_phase_ci_this, ...
   dummy, dummy, ...
   dummy, dummy, ...
   dummy, ...
   P_taos_this]=...
    coh_mt(dt,x_this_cent,...
           NW,K,F,conf_level,p_FFT_extra);
  clear dummy;
  
%   if j==round(n_windows/2)
%     % plot spectra
%     for q=1:Q
%       figure;
%       hold on;
%       plot(f,log10(P_ci_this(:,q,1)),'color',[0.5 0.5 1]);
%       plot(f,log10(P_ci_this(:,q,2)),'color',[0.5 0.5 1]);
%       plot(f,log10(P_this(:,q)));
%       hold off;
%       box on;
%       title(sprintf('q = %d ,t = %f',q,t(j)));
%       figure;
%       hold on;
%       plot(f,squeeze(log10(P_taos_this(:,q,:))),'color',[0.5 0.5 1]);
%       plot(f,log10(P_this(:,q)));
%       hold off;
%       box on;
%       title(sprintf('q = %d ,t = %f',q,t(j)));
%     end
%   end

  x_mean(j,:)=x_this_mean;
  P(:,j,:)=reshape(P_this,[N_f 1 Q]);
  C_mag(:,j,:,:)=reshape(C_mag_this,[N_f 1 Q Q]);
  C_phase(:,j,:,:)=reshape(C_phase_this,[N_f 1 Q Q]);
  x_this_mean_se=sqrt(0.5*P_this(1,:)/T_window);
    % P&W p189.  This is only an approx, but I don't know what else to 
    % use...
  ci_factor=norminv((1+conf_level)/2);
    % Seems like I should be using a T, not a normal, but I'm not clear
    % on how many degrees of freedom I'd use...
  x_mean_ci(j,:,1)=x_this_mean-ci_factor*x_this_mean_se;
  x_mean_ci(j,:,2)=x_this_mean+ci_factor*x_this_mean_se;
  P_ci(:,j,:,:)=P_ci_this;
  C_mag_ci(:,j,:,:,:)=C_mag_ci_this;
  C_phase_ci(:,j,:,:,:)=C_phase_ci_this;
  P_taos(:,j,:,:)=reshape(P_taos_this,[N_f 1 Q K]);  
end
f_res_diam

