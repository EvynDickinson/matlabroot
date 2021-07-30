function C_mag_thresh=f(t,y,x,nw,f_star,alpha_thresh,method)

% process args
if nargin<7
  method='analytic';
end

% calc and print the alpha level for multiple comparisons
N_signals=size(y,2);
alpha_thresh_all=1-(1-alpha_thresh)^N_signals;
if alpha_thresh_all>0.0501
  warning('alpha_thresh_all is greater than 0.05!');
end

% depending on method, calc C_mag_thresh the easy way or the hard way
if strcmp(method,'analytic')
  K=2*nw-1;
  dof=2*K;  % degrees of freedom
  C_mag_thresh=sqrt(1-alpha_thresh^(1/(dof/2-1)));
else
  % get dims
  N=length(t);

  % figure out number of perms to do
  N_perms_approx=100/alpha_thresh;
  N_passes=ceil(N_perms_approx/N_signals)
  N_perms=N_passes*N_signals;

  % do the perms
  x=repmat(x,[1 N_signals]);
  C_mag_samples=zeros([1 N_signals N_passes]);
  %figure;
  for j=1:N_passes
    fprintf(1,'.');  if mod(j,50)==0 fprintf(1,'\n'); end
    perm=randperm(N);
    shuffled_y=y(perm,:);
    %plot(t,shuffled_y);
    C_mag_samples(:,:,j)=...
      coh_mtm(t,shuffled_y,x,nw,[],f_star);
  end
  if mod(j,50)~=0 fprintf(1,'\n'); end

  % collect the samples into a vector
  C_mag_samples=C_mag_samples(:);

  % the histo
  n_bins=100;
  dc=1/n_bins;
  C_mag_grid=0:dc:1;
  C_mag_grid_mod=[ -Inf C_mag_grid +Inf ];
  C_mag_hist=histc(C_mag_samples,C_mag_grid_mod);
  C_mag_hist=C_mag_hist(1:end-2);  % chop off garbage values
  C_mag_cdf_est=cumsum(C_mag_hist)/N_perms;  % add up the values
  C_mag_grid_centers=(C_mag_grid(1:end-1)+C_mag_grid(2:end))/2;
  C_mag_pdf_est=C_mag_hist(2:end)/N_perms/dc;

  % how is this different?
  [C_mag_hist_alt,C_mag_grid_alt]=hist(C_mag_samples,n_bins);
  dc_alt=(C_mag_grid_alt(end)-C_mag_grid_alt(1))/(length(C_mag_grid_alt)-1)
  C_mag_pdf_est_alt=C_mag_hist_alt/N_perms/dc_alt;

  % the theoretical distribution
  K=2*nw-1;
  dof=2*K;  % degrees of freedom
  C_mag_pdf=(dof-2)*C_mag_grid.*(1-C_mag_grid.^2).^(dof/2-2);
  C_mag_cdf=1-(1-C_mag_grid.^2).^(dof/2-1);

  % plot cdf & est
  figure;
  plot(C_mag_grid,C_mag_cdf,'r',...
       C_mag_grid,C_mag_cdf_est,'b');

  % plot pdf & est
  figure;
  plot(C_mag_grid,C_mag_pdf,'r',...
       C_mag_grid_centers,C_mag_pdf_est,'b',...
       C_mag_grid_alt,C_mag_pdf_est_alt,'g');

  % plot pdf & est as bar chart
  figure;
  h=bar(C_mag_grid_centers,C_mag_pdf_est);
  set(h,'EdgeColor','none');
  hold on;
  plot(C_mag_grid,C_mag_pdf,'r');
  hold off;

  % compute the significance threshold
  C_mag_samples_sorted=sort(C_mag_samples);
  C_mag_thresh=interp1(linspace(0,1,N_perms),...
                       C_mag_samples_sorted,...
                       1-alpha_thresh,...
                       'linear')
end