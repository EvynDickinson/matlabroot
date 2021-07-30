function [gastricity,pyloricity,...
          A_gastric,A_pyloric,A_total,...
          gastricity_ci,pyloricity_ci,...
          A_gastric_ci,A_pyloric_ci,A_total_ci]=...
  gpicity(f,Prr,Prr_taos,...
          f_gastric_cutoff,...
          f0_pyloric,...
          n_harmonics_pyloric,...
          W_pyloric,...
          conf_level)

% Prr is n_f x n_t x n_signals
% Prr_taos is n_f x n_t x n_signals x n_tapers
%
% On return:
%   gastricity is n_t x n_signals
%   pyloricity is n_t x n_signals
%   A_gastric is n_t x n_signals
%   A_pyloric is n_t x n_signals
%   A_total is n_t x n_signals
%

% get dims
df=f(2);
n_f=size(Prr,1);
n_t=size(Prr,2);
n_signals=size(Prr,3);
n_tapers=size(Prr_taos,4);

% get total power, RMS amplitude
P_total=reshape(df*sum(Prr,1),[n_t n_signals]);
A_total=sqrt(P_total);

% estimate gastricity
f_range_gastric=(f<f_gastric_cutoff);
P_gastric=reshape(df*sum(Prr(f_range_gastric,:,:),1),...
                  [n_t n_signals]) ; % Hz^2
A_gastric=sqrt(P_gastric);  %  Hz
gastricity=A_gastric./A_total;  % pure

% estimate pyloricity
P_pyloric=zeros(n_t,n_signals);
for i=1:(n_harmonics_pyloric+1)
  f_center=f0_pyloric*i;
  f_range=(f_center-W_pyloric<=f)&(f<f_center+W_pyloric);
  P_pyloric=P_pyloric+...
            reshape(df*sum(Prr(f_range,:,:),1),[n_t n_signals]);  % Hz^2
end
A_pyloric=sqrt(P_pyloric);  %  Hz
pyloricity=A_pyloric./A_total;  % pure


%
% now calc all this stuff for the take-away-one estimates
%

% get total power, RMS amplitude
P_total_taos=reshape(df*sum(Prr_taos,1),[n_t n_signals n_tapers]);
A_total_taos=sqrt(P_total_taos);

% estimate gastricity
P_gastric_taos=reshape(df*sum(Prr_taos(f_range_gastric,:,:,:),1),...
                       [n_t n_signals n_tapers]) ; % Hz^2
A_gastric_taos=sqrt(P_gastric_taos);  %  Hz
gastricity_taos=A_gastric_taos./A_total_taos;  % pure

% estimate pyloricity
P_pyloric_taos=zeros(n_t,n_signals,n_tapers);
for i=1:(n_harmonics_pyloric+1)
  f_center=f0_pyloric*i;
  f_range=(f_center-W_pyloric<=f)&(f<f_center+W_pyloric);
  P_pyloric_taos=P_pyloric_taos+...
                 reshape(df*sum(Prr_taos(f_range,:,:,:),1),...
                         [n_t n_signals n_tapers]);  % Hz^2
end
A_pyloric_taos=sqrt(P_pyloric_taos);  %  Hz
pyloricity_taos=A_pyloric_taos./A_total_taos;  % pure


%
% transform things, and calc the SEs
%

% for A_'s, use log transform
[A_total_log10,A_total_log10_se]=...
  log10_se_jackknife(A_total,A_total_taos);
[A_gastric_log10,A_gastric_log10_se]=...
  log10_se_jackknife(A_gastric,A_gastric_taos);
[A_pyloric_log10,A_pyloric_log10_se]=...
  log10_se_jackknife(A_pyloric,A_pyloric_taos);

% for icity's use antilogisitic-of-square transform
[gastricity_algcsqr,gastricity_algcsqr_se]=...
  algcsqr_se_jackknife(gastricity,gastricity_taos);
[pyloricity_algcsqr,pyloricity_algcsqr_se]=...
  algcsqr_se_jackknife(pyloricity,pyloricity_taos);


%
% translate SE's into CI's
%

ci_factor=tinv((1+conf_level)/2,n_tapers-1);
A_total_ci(:,:,1)=10.^(A_total_log10-ci_factor*A_total_log10_se);
A_total_ci(:,:,2)=10.^(A_total_log10+ci_factor*A_total_log10_se);
A_gastric_ci(:,:,1)=10.^(A_gastric_log10-ci_factor*A_gastric_log10_se);
A_gastric_ci(:,:,2)=10.^(A_gastric_log10+ci_factor*A_gastric_log10_se);
A_pyloric_ci(:,:,1)=10.^(A_pyloric_log10-ci_factor*A_pyloric_log10_se);
A_pyloric_ci(:,:,2)=10.^(A_pyloric_log10+ci_factor*A_pyloric_log10_se);
gastricity_ci(:,:,1)=...
  sqrt(logistic(gastricity_algcsqr-ci_factor*gastricity_algcsqr_se));
gastricity_ci(:,:,2)=...
  sqrt(logistic(gastricity_algcsqr+ci_factor*gastricity_algcsqr_se));
pyloricity_ci(:,:,1)=...
  sqrt(logistic(pyloricity_algcsqr-ci_factor*pyloricity_algcsqr_se));
pyloricity_ci(:,:,2)=...
  sqrt(logistic(pyloricity_algcsqr+ci_factor*pyloricity_algcsqr_se));

