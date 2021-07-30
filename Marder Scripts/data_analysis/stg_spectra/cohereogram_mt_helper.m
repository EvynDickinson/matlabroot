% check for f_scale_log
if ~exist('f_scale_log','var')
  f_scale_log=true;
end
% check for W_show_low
if ~exist('W_show_low','var')
  W_show_low=0.04;
end

% get dimensions
n_files=i_max-i_min+1;
n_cells=length(cell_names);

% load the files
ts=cell(n_cells,1);
for j=1:n_cells
  ts{j}=zeros(0,1);
end
t_offset=0;
for i=i_min:i_max
  file_name=sprintf('%s %04d.smr',exp_name,i);
  for j=1:n_cells
    channel_name=sprintf('%sspikes',cell_names{j});
    [ts_this,T_this]=...
      load_named_event_channel_from_smr(file_name,channel_name);
    ts{j}=[ts{j};ts_this+t_offset];
  end
  t_offset=t_offset+T_this;
end

% set t0, T for all files
t0=0;
T=t_offset;

% filter the spikes according to the given refractory period
for j=1:n_cells
  ts_new=not_too_soon(ts{j},T_refract);
  ts{j}=ts_new;
end

% make a rate signal, timeline
N=round(T/dt);
T=dt*N;
r=zeros(N,n_cells);
for j=1:n_cells
  [t,r_this]=times_to_rate_simple(dt,T,t0,ts{j});
  r(:,j)=r_this;
end
t=double(t);
r=double(r);
r_mean=mean(r,1)
fs=1/dt;  % Hz

% determine window size
n_windows=round(T/T_window_desired)
N_window=floor(N/n_windows);
T_window=dt*N_window

% want N to be integer multiple of N_window
N=N_window*n_windows;
t=t(1:N);
r=r(1:N,:);
T=dt*N

% calc the spectrum for each segment, using multitaper routine
[f,t,...
 r_mean,...
 P,...
 C_mag,C_phase,...
 r_mean_ci,...
 P_ci,...
 C_mag_ci,C_phase_ci,...
 P_taos]=...
  cohereogram_mt(dt,r,T_window,NW,K,W_keep,p_FFT_extra,conf_level);
n_windows_overlapped=length(t);

% calculate gastricity, pyloricity
[gastricity,pyloricity,...
 A_gastric,A_pyloric,A_total,...
 gastricity_ci,pyloricity_ci,...
 A_gastric_ci,A_pyloric_ci,A_total_ci]=...
  gpicity(f,P,P_taos,...
          f_gastric_cutoff,...
          f0_pyloric,...
          n_harmonics_pyloric,...
          W_pyloric,...
          conf_level);

% % plot A_total vs. r_mean
% figure;
% plot(r_mean(1:2:end,:),A_total(1:2:end,:),'.');
% xlabel('Mean rate (Hz)');
% ylabel('Sigma (Hz)');
% xlim([0 10]);
% ylim([0 10]);
% axis square;

% calc means, SDs
r_mean_mean=mean(r_mean(1:2:end,:),1);
r_mean_sd=std(r_mean(1:2:end,:),[],1);
A_total_mean=mean(A_total(1:2:end,:),1);
A_total_sd=std(A_total(1:2:end,:),[],1);
A_gastric_mean=mean(A_gastric(1:2:end,:),1);
A_gastric_sd=std(A_gastric(1:2:end,:),[],1);
A_pyloric_mean=mean(A_pyloric(1:2:end,:),1);
A_pyloric_sd=std(A_pyloric(1:2:end,:),[],1);
pyloricity_mean=mean(pyloricity(1:2:end,:),1);
pyloricity_sd=std(pyloricity(1:2:end,:),[],1);
gastricity_mean=mean(gastricity(1:2:end,:),1);
gastricity_sd=std(gastricity(1:2:end,:),[],1);

% calc normalization for Prr
P_total_mean=mean(A_total(1:2:end,:).^2,1);

% plot spectrograms
for j=1:n_cells
  label=sprintf('%s %04d-%04d %s',exp_name,i_min,i_max,cell_names{j});
  figure;
  %plot_spectrogram(t,f,r_mean(:,j),P_dB(:,:,j),W_show,label,r_lim,P_dB_lim);
  plot_spectrogram_gp_log(t0,T,...
                          t,f,...
                          r_mean(:,j),...
                            reshape(r_mean_ci(:,j,:),...
                                    [n_windows_overlapped 2]),...
                          P(:,:,j),P_total_mean(j),...
                          W_show,label,...
                          r_lim,P_norm_dB_lim,A_lim,...
                          gastricity(:,j),...
                            reshape(gastricity_ci(:,j,:),...
                                    [n_windows_overlapped 2]),...
                          pyloricity(:,j),...
                            reshape(pyloricity_ci(:,j,:),...
                                    [n_windows_overlapped 2]),...
                          A_gastric(:,j),...
                            reshape(A_gastric_ci(:,j,:),...
                                    [n_windows_overlapped 2]),...
                          A_pyloric(:,j),...
                            reshape(A_pyloric_ci(:,j,:),...
                                    [n_windows_overlapped 2]),...
                          A_total(:,j),...
                            reshape(A_total_ci(:,j,:),...
                                    [n_windows_overlapped 2]),...
                          r_mean_mean(j),r_mean_sd(j),...
                          A_total_mean(j),A_total_sd(j),...
                          A_gastric_mean(j),A_gastric_sd(j),...
                          A_pyloric_mean(j),A_pyloric_sd(j),...
                          pyloricity_mean(j),pyloricity_sd(j),...
                          gastricity_mean(j),gastricity_sd(j),...
                          f_gastric_cutoff,...
                          f0_pyloric,...
                          n_harmonics_pyloric,...
                          W_pyloric,...
                          f_scale_log,...
                          W_show_low);
  if exist('make_pdfs','var') && make_pdfs
    print_pdf(gcf,label);
  end
end

% plot cohereograms
for j=1:n_cells
  for k=1:n_cells
    if j>k
      label=sprintf('%s %04d-%04d C_%s_%s',exp_name,...
                    i_min,i_max,cell_names{j},cell_names{k});
      figure;
      plot_cohereogram_gp_log(t0,T,...
                              t,f,C_mag(:,:,j,k).^2,C_phase(:,:,j,k),...
                              W_show,label,...
                              f_gastric_cutoff,...
                              f0_pyloric,...
                              n_harmonics_pyloric,...
                              W_pyloric,...
                              f_scale_log,...
                              W_show_low);
      if exist('make_pdfs','var') && make_pdfs
        print_pdf(gcf,label);
      end
    end
  end
end

% plot the legend for the cohereograms
alpha_thresh=0.05;
dof=2*K;  % degrees of freedom
C_mag_thresh=sqrt(1-alpha_thresh^(1/(dof/2-1)));
figure;
plot_coh2l75_border_legend(C_mag_thresh);



%
% now calculate spectra for the whole data set
%

% subtract off mean rate
r_mean=mean(r);
r_cent=r-repmat(r_mean,[N 1]);

% put windows into the third index
r_cent_windowed_proto=reshape(r_cent,[N_window n_windows n_cells]);
r_cent_windowed=permute(r_cent_windowed_proto,[1 3 2]);

% calc the spectra+coherences, using multitaper routine
[f,Prr,C_mag,C_phase,...
 N_fft,f_res_diam,...
 Prr_ci ...
 C_mag_ci ...
 C_phase_ci ...
 dummy dummy ...
 dummy dummy ...
 dummy ...
 Prr_taos]=...
  coh_mt(dt,r_cent_windowed,...
         NW,K,W_keep,conf_level,p_FFT_extra);
clear dummy;
N_fft
f_res_diam
n_f=length(f);

% calc amplitudes
Arr=sqrt(Prr);
Arr_ci=sqrt(Prr_ci);

% find the fundamental pyloric band
f_min_pyl_fund=f0_pyloric-W_pyloric;
f_max_pyl_fund=f0_pyloric+W_pyloric;
j_f=(1:n_f)';
in_pyl_fund_band=((f_min_pyl_fund<=f)&(f<f_max_pyl_fund));
j_in_pyl_fund_band=j_f(in_pyl_fund_band);

% get the mag of the peak amp. within the fundamental pyloric band
f_Arr_pyl_peak=nan(n_cells);
Arr_pyl_peak=nan(n_cells);
for k=1:n_cells
  [dummy,j_temp]=max(Arr(in_pyl_fund_band,k));
  j_pyl_peak=j_in_pyl_fund_band(j_temp);
  f_Arr_pyl_peak(k)=f(j_pyl_peak);
  Arr_pyl_peak(k)=Arr(j_pyl_peak,k);
end

% get the mag, angle of the peak coh within the fundamental pyloric band
f_C_mag_pyl_peak=nan(n_cells,n_cells);
C_mag_pyl_peak=nan(n_cells,n_cells);
C_phase_pyl_peak=nan(n_cells,n_cells);
for k=1:n_cells
  for l=1:n_cells
    [dummy,j_temp]=max(C_mag(in_pyl_fund_band,k,l));
    j_pyl_peak=j_in_pyl_fund_band(j_temp);
    f_C_mag_pyl_peak(k,l)=f(j_pyl_peak);
    C_mag_pyl_peak(k,l)=C_mag(j_pyl_peak,k,l);
    C_phase_pyl_peak(k,l)=C_phase(j_pyl_peak,k,l);
    if ((k==2)&&(l==1))
      drawnow;
    end
  end
end

% % calculate gastricity, pyloricity
% [gastricity,pyloricity,...
%  A_gastric,A_pyloric,A_total,...
%  gastricity_ci,pyloricity_ci,...
%  A_gastric_ci,A_pyloric_ci,A_total_ci]=...
%   gpicity(f,...
%           reshape(Prr,[n_f 1 n_cells]),...
%           reshape(Prr_taos,[n_f 1 n_cells n_windows*K]),...
%           f_gastric_cutoff,...
%           f0_pyloric,...
%           n_harmonics_pyloric,...
%           W_pyloric,...
%           conf_level);
% gastricity=reshape(gastricity,[n_cells 1]);
% gastricity_ci=reshape(gastricity_ci,[n_cells 2]);
% pyloricity=reshape(pyloricity,[n_cells 1]);
% pyloricity_ci=reshape(pyloricity_ci,[n_cells 2]);

% plot spectra
gray=[0.5 0.5 0.5];
for j=1:n_cells
  label=sprintf('%s %04d-%04d %s spectrum',exp_name,...
                i_min,i_max,cell_names{j});
  figure;
  line(f,Arr_ci(:,j,1),'color',gray);
  line(f,Arr_ci(:,j,2),'color',gray);
  line(f,Arr(:,j),'color','k');
  line(f_Arr_pyl_peak(j),Arr_pyl_peak(j),...
       'marker','o','linestyle','none',...
       'color','r');
  text(f_Arr_pyl_peak(j),0,...
       sprintf('%0.2f',f_Arr_pyl_peak(j)),...
       'horizontalalignment','center',...
       'verticalalignment','bottom');
  yl=ylim;
  ylim([0 yl(2)]);
  xlim([0 W_show]);
  text(W_show,1.01*yl(2),...
       sprintf('%0.2f',Arr_pyl_peak(j)),...
       'horizontalalignment','right',...
       'verticalalignment','bottom');
  ylabel('Amplitude (Hz/Hz^{0.5})');
  xlabel('Frequency (Hz)');
  title(label,'interpreter','none');
  box on;
  if exist('make_pdfs','var') && make_pdfs
    print_pdf(gcf,label);
  end
end

% % plot icities
% fig_h=figure;
% text(gastricity,pyloricity,cell_names);
% xlim([0 1.1]);
% ylim([0 1.1]);
% xlabel('Gastricity');
% ylabel('Pyloricity');
% axis square;
% title('Pyloricity vs. gastricity');
% if exist('make_pdfs','var') && make_pdfs
%   print_pdf(fig_h,sprintf('%s %04d icities',exp_name,trial_index));
% end

% calc the significance threshold for coherence
alpha_thresh=0.05;
dof=2*K*n_windows;  % degrees of freedom
C_mag_thresh=sqrt(1-alpha_thresh^(1/(dof/2-1)));

% plot coherences
purple=[0.5 0 0.5];
light_purple=0.9*[1 0.8 1];
for j=1:n_cells
  for k=1:n_cells
    if j>k
      label=sprintf('%s %04d-%04d C_%s_%s',exp_name,...
                    i_min,i_max,cell_names{j},cell_names{k});
      figure;
      subplot(2,1,1);
      line(f,repmat(C_mag_thresh,size(f)),...
           'linestyle',':','color',0.7*[1 1 1]);
      line(f,C_mag_ci(:,j,k,1),'color',light_purple);
      line(f,C_mag_ci(:,j,k,2),'color',light_purple);
      line(f,C_mag(:,j,k),'color',purple);
      line(f_C_mag_pyl_peak(j,k),C_mag_pyl_peak(j,k),...
           'marker','o','linestyle','none',...
           'color',purple);
      ylim([0 1.05]);
      xlim([0 W_show]);
      text(W_show,1.05,...
           sprintf('%0.2f',C_mag_pyl_peak(j,k)),...
           'horizontalalignment','right',...
           'verticalalignment','bottom');
      ylabel('C mag');
      title(label,'interpreter','none');
      box on;
      subplot(2,1,2);
      if all(isfinite(C_phase_ci(:,j,k,1)))
        line_wrap(f,-0.5/pi*C_phase_ci(:,j,k,1),[0 1],...
                  'color',light_purple);
      end
      if all(isfinite(C_phase_ci(:,j,k,2)))
        line_wrap(f,-0.5/pi*C_phase_ci(:,j,k,2),[0 1],...
                  'color',light_purple);
      end
      if all(isfinite(C_phase(:,j,k)))
        line_wrap(f,-0.5/pi*C_phase(:,j,k),[0 1],...
                  'color',purple);
      end
      line(f_C_mag_pyl_peak(j,k),wrap01(-C_phase_pyl_peak(j,k)),...
           'marker','o','linestyle','none',...
           'color',purple);
      ylim([0 1]);
      xlim([0 W_show]);
      text(W_show,1,...
           sprintf('%0.2f',wrap01(-C_phase_pyl_peak(j,k))),...
           'horizontalalignment','right',...
           'verticalalignment','bottom');
      ylabel('C phase (cycles)');
      xlabel('Frequency (Hz)');
      box on;
      if exist('make_pdfs','var') && make_pdfs
        print_pdf(gcf,label);
      end
    end
  end
end


