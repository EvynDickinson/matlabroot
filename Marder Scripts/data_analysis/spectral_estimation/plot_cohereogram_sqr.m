function plot_cohereogram_sqr(t,f,C_mag,C_phase,F,label,C_mag_thresh)

% this version combines magnitude and phase info into one plot

C=C_mag.*exp(1i*C_phase);  % want complex coherence

set_figure_size([10 7.5]);
set(gcf,'PaperOrientation','landscape');

axes_sg=subplot(5,1,1:4);
cmap_phase=l75_border(256);  % to show abs(C)==1 colors
im=cohsqr2l75_border(C,C_mag_thresh);
image(t,f,im);
axis xy;
ylim([0 F]);
ylabel('Frequency (Hz)');
title(label,'interpreter','none');
colorbar_axes_h=colorbar;
colorbar_image_h=findobj(colorbar_axes_h,'Tag','TMW_COLORBAR');
set(colorbar_image_h,'YData',[-180 +180]);
set(colorbar_axes_h,'YLim',[-180 +180]);
set(colorbar_image_h,'CData',reshape(cmap_phase,[256 1 3]));
set(colorbar_axes_h,'YTick',[-180 -90 0 +90 +180]);

drawnow;
