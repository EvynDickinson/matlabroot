function plot_spectrogram_gp_log(t0,T,...
                                 t,f,...
                                 r_mean,r_mean_ci,...
                                 Prr,Prr_total_mean,...
                                 W_show,label,...
                                 r_lim,Prr_norm_dB_lim,A_lim,...
                                 gastricity,gastricity_ci,...
                                 pyloricity,pyloricity_ci,...
                                 A_gastric,A_gastric_ci,...
                                 A_pyloric,A_pyloric_ci,...
                                 A_total,A_total_ci,...
                                 r_mean_mean,r_mean_sd,...
                                 A_total_mean,A_total_sd,...
                                 A_gastric_mean,A_gastric_sd,...
                                 A_pyloric_mean,A_pyloric_sd,...
                                 pyloricity_mean,pyloricity_sd,...
                                 gastricity_mean,gastricity_sd,...
                                 f_gastric_cutoff,...
                                 f0_pyloric,...
                                 n_harmonics_pyloric,...
                                 W_pyloric,...
                                 f_scale_log,...
                                 W_show_low)

% deal w/ args
if nargin<40 || isempty(f_scale_log)
  f_scale_log=true;
end

% includes a plot of the pyloric & gastric amplitude
% this uses a log scale for the frequency

Prr_norm=Prr./repmat(Prr_total_mean,size(Prr));
Prr_norm_dB=10*log10(Prr_norm);  % convert to dB
n_windows=length(t);
if n_windows>1
  dt=(t(end)-t(1))/(length(t)-1);  % need this later
else
  dt=T/2;
end
set_figure_size([10 7.5]);
set(gcf,'PaperOrientation','landscape');

axes_sg=subplot(6,1,1:3);
if ~f_scale_log
  if isempty(Prr_norm_dB_lim)
    imagesc(t,f,Prr_norm_dB);
  else
    imagesc(t,f,Prr_norm_dB,Prr_norm_dB_lim);
  end  
  axis xy;
  ylim([0 W_show]);
else
  % want to resample f, Prr_dB onto a log scale
  f_log10_min=log10(f(2));
  f_log10_max=log10(W_show);
  df_log10_want=0.001;
  n_f_log10=ceil((f_log10_max-f_log10_min)/df_log10_want)+1;
  df_log10=(f_log10_max-f_log10_min)/(n_f_log10-1);
  f_log10=f_log10_min+df_log10*(0:(n_f_log10-1))';
  f_resamp=10.^f_log10;
  Prr_norm_dB_resamp=interp1(f,Prr_norm_dB,f_resamp,'linear');

  % now, want to duplicate first and last cols
  t_resamp=[t(1)-dt;t;t(end)+dt];
  Prr_norm_dB_resamp=...
    [Prr_norm_dB_resamp(:,1) Prr_norm_dB_resamp Prr_norm_dB_resamp(:,end)];
  
  % okay, now do actual plotting
  if isempty(Prr_norm_dB_lim)
    image('cdata',Prr_norm_dB_resamp,...
          'xdata',[t_resamp(1) t_resamp(end)],...
          'ydata',[f_resamp(1) f_resamp(end)],...
          'cdatamapping','scaled');
  else
    image('cdata',Prr_norm_dB_resamp,...
          'xdata',[t_resamp(1) t_resamp(end)],...
          'ydata',[f_resamp(1) f_resamp(end)],...
          'cdatamapping','scaled');
    caxis(Prr_norm_dB_lim);
  end
  box on;
  set(gca,'yscale','log');
  set(gca,'layer','top');
  ylim([W_show_low W_show]);
end
xlim(t0+[0 T]);
ylabel('Frequency (Hz)');
title(label,'interpreter','none');
colormap(spectrum_smooth(256));
colorbar;

% draw gastricity, pyloricity cutoff lines
hold on;
line(xlim,[f_gastric_cutoff f_gastric_cutoff],'color',[1 0 0]);
for i=1:(n_harmonics_pyloric+1)
  f_center=f0_pyloric*i;
  f_lower=f_center-W_pyloric;
  f_upper=f_center+W_pyloric;
  line(xlim,[f_lower f_lower],'color',[0 0.8 0]);
  line(xlim,[f_upper f_upper],'color',[0 0.8 0]);
end
hold off;

axes_Agp=subplot(6,1,4);
hold on;
plot(t,A_gastric_ci(:,1),'color',[1 0.75 0.75]);
plot(t,A_gastric_ci(:,2),'color',[1 0.75 0.75]);
plot(t,A_pyloric_ci(:,1),'color',[0.75 1 0.75]);
plot(t,A_pyloric_ci(:,2),'color',[0.75 1 0.75]);
plot(t,A_total_ci(:,1),'color',[0.75 0.75 0.75]);
plot(t,A_total_ci(:,2),'color',[0.75 0.75 0.75]);
plot(t,A_gastric,'color',[1 0 0]);
%plot(t,A_gastric,'.','color',[1 0 0]);
plot(t,A_pyloric,'color',[0 0.8 0]);
%plot(t,A_pyloric,'.','color',[0 0.8 0]);
plot(t,A_total,'color',[0 0 0]);
hold off;
box on;
pos_sg=get(axes_sg,'position');
pos_Agp=get(gca,'position');
pos_Agp(1)=pos_sg(1);
pos_Agp(3)=pos_sg(3);
set(gca,'position',pos_Agp);
ylabel('Amp. (Hz)');
xlim(get(axes_sg,'xlim'));
if ~isempty(A_lim)
  ylim(A_lim);
end
xl=xlim;
yl=ylim;
x_text=xl(2)+0.01*diff(xl);
dy=-diff(yl)*0.25;
y_origin=yl(2)+dy;
text(x_text,y_origin,...
     sprintf('%0.2f \\pm %0.2f',A_total_mean,A_total_sd),...
     'color','k');
text(x_text,y_origin+dy,...
     sprintf('%0.2f \\pm %0.2f',A_pyloric_mean,A_total_sd),...
     'color',[0 0.8 0]);
text(x_text,y_origin+2*dy,...
     sprintf('%0.2f \\pm %0.2f',A_gastric_mean,A_gastric_sd),...
     'color','r');

axes_gp=subplot(6,1,5);
hold on;
plot(t,gastricity_ci(:,1),'color',[1 0.75 0.75]);
plot(t,gastricity_ci(:,2),'color',[1 0.75 0.75]);
plot(t,pyloricity_ci(:,1),'color',[0.75 0.9 0.75]);
plot(t,pyloricity_ci(:,2),'color',[0.75 0.9 0.75]);
plot(t,gastricity,'color',[1 0 0]);
plot(t,pyloricity,'color',[0 0.8 0]);
hold off;
box on;
pos_sg=get(axes_sg,'position');
pos_gp=get(gca,'position');
pos_gp(1)=pos_sg(1);
pos_gp(3)=pos_sg(3);
set(gca,'position',pos_gp);
ylabel('Normed Amp.');
ylim([0 1]);
xlim(get(axes_sg,'xlim'));
xl=xlim;
yl=ylim;
x_text=xl(2)+0.01*diff(xl);
dy=-diff(yl)*0.25;
y_origin=yl(2)+dy;
text(x_text,y_origin,...
     sprintf('%0.3f \\pm %0.3f',pyloricity_mean,pyloricity_sd),...
     'color',[0 0.8 0]);     
text(x_text,y_origin+dy,...
     sprintf('%0.3f \\pm %0.3f',gastricity_mean,gastricity_sd),...
     'color','r');     

axes_rate=subplot(6,1,6);
hold on;
plot(t,r_mean_ci(:,1),'color',0.75*[1 1 1]);
plot(t,r_mean_ci(:,2),'color',0.75*[1 1 1]);
plot(t,r_mean,'k');
hold off;
box on;
xlabel('Time (s)');
ylabel('Rate (Hz)');
pos_sg=get(axes_sg,'position');
pos_rate=get(gca,'position');
pos_rate(1)=pos_sg(1);
pos_rate(3)=pos_sg(3);
set(gca,'position',pos_rate);
xlim(get(axes_sg,'xlim'));
if isempty(r_lim)
  yl=get(gca,'ylim');
  set(gca,'ylim',[0 yl(2)]);
else
  set(gca,'ylim',r_lim);
end
xl=xlim;
yl=ylim;
x_text=xl(2)+0.01*diff(xl);
dy=-diff(yl)*0.25;
y_origin=yl(2)+dy;
text(x_text,y_origin,...
     sprintf('%0.2f \\pm %0.2f',r_mean_mean,r_mean_sd),...
     'color','k');

drawnow;
