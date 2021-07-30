function plot_cohereogram_gp(t0,T,...
                             t,f,C_mag,C_phase,W_show,label,...
                             f_gastric_cutoff,...
                             f0_pyloric,...
                             n_harmonics_pyloric,...
                             W_pyloric,...
                             f_scale_log,...
                             W_show_low)

% designed to align with plot_spectrogram_gp_log

% deal w/ args
if nargin<13 || isempty(f_scale_log)
  f_scale_log=true;
end

C=C_mag.*exp(i*C_phase);  % want complex coherence
n_windows=length(t);
if n_windows>1
  dt=(t(end)-t(1))/(length(t)-1);  % need this later
else
  dt=T/2;
end

set_figure_size([10 7.5]);
set(gcf,'PaperOrientation','landscape');

% do the cohereogram itself
axes_sg=subplot(6,1,1:3);
cmap_phase=l75_border(256);  % to show abs(C)==1 colors
if f_scale_log
  % want to resample f, C onto a log scale
  f_log10_min=log10(f(2));
  f_log10_max=log10(W_show);
  df_log10_want=0.001;
  n_f_log10=ceil((f_log10_max-f_log10_min)/df_log10_want)+1;
  df_log10=(f_log10_max-f_log10_min)/(n_f_log10-1);
  f_log10=f_log10_min+df_log10*(0:(n_f_log10-1))';
  f_resamp=10.^f_log10;
  C_resamp=interp1(f,C,f_resamp,'linear');

  % now, want to duplicate first and last cols
  t_resamp=[t(1)-dt;t;t(end)+dt];
  C_resamp=[C_resamp(:,1) C_resamp C_resamp(:,end)];

  im=coh2l75_border(C_resamp);
  image('cdata',im,...
        'xdata',[t_resamp(1) t_resamp(end)],...
        'ydata',[f_resamp(1) f_resamp(end)]);
  %axis xy;
  set(gca,'yscale','log');
  set(gca,'layer','top');
  ylim([W_show_low W_show]);
else
  im=coh2l75_border(C);
  image(t,f,im);
  axis xy;
  ylim([0 W_show]);
end
xlim(t0+[0 T]);
ylabel('Frequency (Hz)');
xlabel('Time (s)');
title(label,'interpreter','none');

% draw gastricity, pyloricity cutoff lines
hold on;
line(xlim,[f_gastric_cutoff f_gastric_cutoff],'color',[1 0 0]);
for j=1:(n_harmonics_pyloric+1)
  f_center=f0_pyloric*j;
  f_lower=f_center-W_pyloric;
  f_upper=f_center+W_pyloric;
  line(xlim,[f_lower f_lower],'color',[0 0.8 0]);
  line(xlim,[f_upper f_upper],'color',[0 0.8 0]);
end
hold off;

% draw the colorbar
colorbar_axes_h=colorbar;
colorbar_image_h=findobj(colorbar_axes_h,'Tag','TMW_COLORBAR');
set(colorbar_image_h,'YData',[-180 +180]);
set(colorbar_axes_h,'YLim',[-180 +180]);
set(colorbar_image_h,'CData',reshape(cmap_phase,[256 1 3]));
set(colorbar_axes_h,'YTick',[-180 -90 0 +90 +180]);

drawnow;
