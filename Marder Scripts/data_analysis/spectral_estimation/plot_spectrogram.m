function plot_spectrogram(t,f,x_mean,Pxx_dB,F,label,units,x_lim,P_dB_lim)

if nargin<7
  units='';
end
if nargin<8
  x_lim=[];
end
if nargin<9
  P_dB_lim=[];
end

set_figure_size([10 7.5]);
set(gcf,'PaperOrientation','landscape');

axes_sg=subplot(5,1,1:4);
if isempty(P_dB_lim)
  imagesc(t,f,Pxx_dB);
else
  imagesc(t,f,Pxx_dB,P_dB_lim);
end  
axis xy;
%colormap(gray(256));
colormap(spectrum_smooth(256));
ylabel('Frequency (Hz)');
ylim([0 F]);
colorbar;
title(label,'interpreter','none');

subplot(5,1,5);
plot(t,x_mean,'k');
xlabel('Time (s)');
if isempty(units)
  ylabel('Mean');
else
  ylabel(sprintf('Mean (%s)',units));
end
pos_sg=get(axes_sg,'position');
pos_rate=get(gca,'position');
pos_rate(1)=pos_sg(1);
pos_rate(3)=pos_sg(3);
set(gca,'position',pos_rate);
xlim(get(axes_sg,'xlim'));
if ~isempty(x_lim)
  set(gca,'ylim',x_lim);
end  
drawnow;
