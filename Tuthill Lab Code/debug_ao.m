
clear all; close all; daqreset;
s_vid_light = daq.createSession('ni');
s_vid_light.Rate = 100;
intensity = 1.5    ; 
light_length = 2  ;
on = [zeros(s_vid_light.Rate,1)' intensity*ones(s_vid_light.Rate,1)']'; % ends on
off = [intensity*ones(s_vid_light.Rate*light_length,1)' zeros(s_vid_light.Rate,1)']'; %stays off 


laser =  on;

laser =  off;  


addAnalogOutputChannel(s_vid_light,'Dev1', 'ao0', 'Voltage');
queueOutputData(s_vid_light, laser) 
startBackground(s_vid_light); 

     
t = datetime('now');
fprintf(['\n   Finished ' num2str(intensity) 'V: \n']) 
disp(t)  


%% ERROR CHECKING: test for missed AO signals
figure(2); hold all;
    plot(data(8,:), 'color', color.orange); %Basler camera trigger
    plot(data(5,:), 'k'); %condition signal
    
        vline(peaks)
        vline(peaks2)
 
    ylabel('volts')
    xlabel('time')
%     h = questdlg('Any slips in the AO signal?', 'Error', 'No');
%     fprintf('\n test')

a = num.starts == num.peaks ;
if a == 0
    h = warndlg('Missing a peak!');
    uiwait(h)
    hh = figure(2); hold all;
        plot(data(8,:), 'color', color.orange); %Basler camera trigger
        plot(data(5,:), 'k'); %condition signal
        vline(peaks)
        vline(peaks2)
        ylabel('volts')
        xlabel('time')
    uiwait(hh)
    answer = questdlg('Did the AO signal return to baseline?', 'Error', 'No');
    if answer == 'No'
        pd = diff(peaks);
        wrong_peak = find(pd > mean(pd)+std(pd));
        hh = figure(2); hold all;
            plot(data(8,:), 'color', color.orange); %Basler camera trigger
            plot(data(5,:), 'k'); %condition signal
            vline(peaks)
            vline(peaks2)
            vline(peaks2(wrong_peak), 'g')
            ylabel('volts')
            xlabel('time')
        uiwait(hh)
        answer2 = questdlg('Is this the correct AO condition signal?', 'Error', 'Yes');
        if answer2 == 'Yes'
           peaks = [peaks(1:wrong_peak), peaks2(wrong_peak), peaks(wrong_peak+1:end)]; 
           hhh = figure; hold all; plot(peaks); plot(peaks2)
           uiwait(hhh)
        end
    end
end
clear a answer answer2 h hh hhh

figure; hold all
plot(pd)
hline(mean(pd)-std(pd))
hline(mean(pd)+std(pd))

%% AO signal adjustments

clear peaks peaks2
peaks = find(diff(cs) > 2*std(diff(cs)));
peaks([diff(peaks) < 10]) = [];
peaks2 = find(diff(cs) < -2*std(diff(cs)));
peaks2([diff(peaks2) < 10]) = [];

figure; hold all;
    plot(data(8,:), 'color', color.orange); %Basler camera trigger
    plot(data(5,:), 'k'); %condition signal
    vline(peaks, 'r:') %red = start of signal
    vline(peaks2, 'k:') %black = end
    ylabel('volts')
    xlabel('time')
    
% -------------------------------------------------------------------------  
% REMOVE bad Condition AO peaks (create new extrapeak variable from graph)
clear a A
buffer = 2000;
a = extrapeak(1,1) + buffer;
b = extrapeak(1,1) - buffer;
% peaks2
A = find(peaks2<=a & peaks2>=b)
% A = A+1;
peaks2(A(1)) = [];
% peaks
A = find(peaks<=a & peaks>=b)

peaks(A(1)) = [];

peaks(A(1)+1) = [];
    
figure; hold all;
    plot(data(8,:), 'color', color.orange); %Basler camera trigger
    plot(data(5,:), 'k'); %condition signal
    vline(peaks, 'r:')
    vline(peaks2, 'k:')
    ylabel('volts')
    xlabel('time')
    
    vline(peaks2(A), 'b')
    
    vline(peaks(A), 'b')
    
    vline(extrapeak, 'b')
    
    
% -------------------------------------------------------------------------
% ADD in new peaks
%TYPE ONE: when the AO signal is overly long
extrapeak = extrapeak(1,1);

peaks = [peaks(1:A), peaks2(A), peaks(A+1:end)]; 
% figure; scatter(1:96, [peaks(1:A), peaks2(A), peaks(A+1:end)])


%TYPE TWO: blank adding in new position:
%check postion with graph
clear a A b
extrapeak = extrapeak(1,1);
figure; hold all;
    plot(data(8,:), 'color', color.orange); %Basler camera trigger
    plot(data(5,:), 'k'); %condition signal
    vline(peaks, 'r:')
    vline(peaks2, 'k:')
    ylabel('volts')
    xlabel('time')
vline(extrapeak, 'b')

% PEAKS
%find the peak positions within the broad spectrum of other peaks
for pp = 1:length(peaks)
    a(pp,1) = extrapeak>peaks(pp);
end
%missing peak is in position 1
if sum(a)==0
    peaks = [extrapeak, peaks];
end
%missing peak is in another positon

% PEAKS2
%find the peak positions within the broad spectrum of other peaks
for pp = 1:length(peaks2)
    a(pp,1) = extrapeak>peaks2(pp);
end
%missing peak is in position 1
if sum(a)==0
    peaks = [extrapeak, peaks2];
end
%missing peak is in another positon
A = find(a==0);
A = A(1);
figure; hold all;
    plot(data(8,:), 'color', color.orange); %Basler camera trigger
    plot(data(5,:), 'k'); %condition signal
    vline(peaks, 'r:')
    vline(peaks2, 'k:')
    ylabel('volts')
    xlabel('time')
    
    vline([peaks(1:A), extrapeak, peaks(A+1:end)], 'g')
    
    vline([peaks2(1:A-1), extrapeak, peaks2(A:end)], 'g')

peaks = [peaks(1:A), extrapeak, peaks(A+1:end)]; 
peaks2 = [peaks2(1:A-1), extrapeak, peaks2(A:end)];











