function [A_log10,A_log10_se]=log10_se_jackknife(A,A_taos)

% A is n_t x n_signals
% A_taos is n_t x n_signals x n_tapers

n_t=size(A,1);
n_signals=size(A,2);
n_tapers=size(A_taos,3);

A_log10=log10(A);
A_taos_log10=log10(A_taos);
A_taos_log10_mean=mean(A_taos_log10,3);  % mean across TAOs
A_log10_se=...
  sqrt((n_tapers-1)/n_tapers*...
         sum((A_taos_log10-...
              repmat(A_taos_log10_mean,[1 1 n_tapers])).^2,3));
