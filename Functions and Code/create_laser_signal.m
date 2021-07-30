
function [conds_matrix, LED] = create_laser_signal(conds_matrix, param)
% 
% [conds_matrix, LED] = create_laser_signal(conds_matrix, param)
% 
% Add the laser and basler camera signals into the conditions 
% matrix with the appropriate trigger timing. 
% 
% Inputs:
% 'conds_matrix' [structure with the condition information]
% 'param' [structure with parameter information]
% 
% Outputs:
% 'conds_matrix' [structure with signals added to each cond]
% 'LED' [structure with all the signal information per cond]
% 
% ES Dickinson, University of Washington, 2019

% % % LASER INTENSITY SCRIPT
% % Set up LED and Basler data queues -- ONE/STIM TYPE
% samplerate = param.basler_samplerate;
% signal_length = samplerate*param.basler_length;
% 
% for cond = 2:length(param.LED_intensity)
%  
%     % ---- LED trigger data ---- %
%     %ratio: ON:OFF
%     LED(cond).frequency = param.laser_freq; %hz of LED frequency
%     LED(cond).ratio_num = param.laser_ratio; % 1:num ratio of light off
%     LED(cond).rate = round(samplerate/LED(cond).frequency); %adjustment ratio to get LED signal speed correct 
%     LED(cond).DutyCycle = 1/(1+LED(cond).rate);
%     LED(cond).light_length = param.light_length;
%       
%     % creating the LED trigger signal
%     LED(cond).pulse = [(ones(LED(cond).rate,1)*param.LED_intensity(cond)); zeros(round(LED(cond).ratio_num*LED(cond).rate),1)];
%     LED(cond).pulse_length = length(LED(cond).pulse);
%     LED(cond).pulse_num = round(samplerate/LED(cond).pulse_length)*param.light_length; %should equal the desired light Hz  
%     LED(cond).sig = LED(cond).pulse;
%     %concatenate the individual pulses to reach the total length of light exposure
%     for ii = 1:(LED(cond).pulse_num-1)
%         LED(cond).sig = [LED(cond).sig; LED(cond).pulse];
%     end  
% 
%     LED(cond).sig(end-(2*LED(cond).pulse_length-1):end) = 0; %two units of off at the end of the signal
%     LED(cond).pre_sig = zeros(samplerate*param.basler_delay,1); %actual signal
%     LED(cond).post_sig = param.basler_length-param.basler_delay-param.light_length; %time (sec) of post LED w/basler on
% 
%     LED(cond).ON_outsig = [LED(cond).pre_sig; LED(cond).sig; zeros(round(samplerate*LED(cond).post_sig), 1)]; %; zeros((10*param.light_length(cond)),1)
%     xxlength = signal_length-size(LED(cond).ON_outsig,1);
%     LED(cond).ON_outsig = [LED(cond).ON_outsig; zeros(xxlength,1)];
%     
%     a = size(LED(cond).ON_outsig);
%     LED(cond).OFF_outsig = zeros(a(1), 1);
%     
% end
% 
% LED(1).light_length = param.light_length(1);
% LED(1).ON_outsig = LED(2).OFF_outsig;
% LED(1).OFF_outsig = LED(2).OFF_outsig;
% 
% % ---- Basler trigger data ---- %
% basler_volts = param.basler_volts;
% basler_outsig = zeros(a(1), 1);
% basler_rate = round(samplerate/param.Basler_fps);
% basler_outsig(1:basler_rate:end) = basler_volts;
% clear a xxlength
% 
% %error checking analogue output signal
% if size(LED(3).ON_outsig) == size(basler_outsig)
%     fprintf('\n Laser signal and camera signal aligned \n')
% else 
%     warndlg('Mismatched laser and camera signal')
%     return
% end
% 
% % ---- COMBINED DATA QUEUE ---- %
% for cond = 1:param.num_conds
%     switch conds_matrix(cond).opto
%         case param.LED_intensity(1)
%             conds_matrix(cond).LED = [LED(1).ON_outsig, basler_outsig];
%         case param.LED_intensity(2)
%             conds_matrix(cond).LED = [LED(2).ON_outsig, basler_outsig]; 
%         case param.LED_intensity(3)
%             conds_matrix(cond).LED = [LED(3).ON_outsig, basler_outsig]; 
%         case param.LED_intensity(4)
%             conds_matrix(cond).LED = [LED(4).ON_outsig, basler_outsig]; 
%         case param.LED_intensity(5)
%             conds_matrix(cond).LED = [LED(5).ON_outsig, basler_outsig]; 
%         case param.LED_intensity(6)
%             conds_matrix(cond).LED = [LED(6).ON_outsig, basler_outsig]; 
%         case param.LED_intensity(7)
%             conds_matrix(cond).LED = [LED(7).ON_outsig, basler_outsig];
%     end
% end

% figure;
% for ii = 1:7
%     subplot(2,4,ii)
%     plot(conds_matrix(ii).LED(:,1))
% end


% LASER DURATION SCRIPT:
% Set up LED and Basler data queues -- ONE/STIM TYPE
samplerate = param.basler_samplerate;
signal_length = samplerate*param.basler_length;

for cond = 2:length(param.light_length)
 
    % ---- LED trigger data ---- %
    %ratio: ON:OFF
    LED(cond).frequency = param.laser_freq; %hz of LED frequency
    LED(cond).ratio_num = param.laser_ratio; % 1:num ratio of light off
    LED(cond).rate = round(samplerate/LED(cond).frequency); %adjustment ratio to get LED signal speed correct 
    LED(cond).DutyCycle = 1/(1+LED(cond).rate);
    LED(cond).light_length = param.light_length(cond);
      
    % creating the LED trigger signal
    LED(cond).pulse = [(ones(LED(cond).rate,1)*param.LED_intensity); zeros(round(LED(cond).ratio_num*LED(cond).rate),1)];
    LED(cond).pulse_length = length(LED(cond).pulse);
    LED(cond).pulse_num = round(samplerate/LED(cond).pulse_length)*param.light_length(cond); %should equal the desired light Hz  
    LED(cond).sig = LED(cond).pulse;
    %concatenate the individual pulses to reach the total length of light exposure
    for ii = 1:(LED(cond).pulse_num-1)
        LED(cond).sig = [LED(cond).sig; LED(cond).pulse];
    end  

    LED(cond).sig(end-(2*LED(cond).pulse_length-1):end) = 0; %two units of off at the end of the signal
    LED(cond).pre_sig = zeros(samplerate*param.basler_delay,1); %actual signal
    LED(cond).post_sig = param.basler_length-param.basler_delay-param.light_length(cond); %time (sec) of post LED w/basler on

    LED(cond).ON_outsig = [LED(cond).pre_sig; LED(cond).sig; zeros(round(samplerate*LED(cond).post_sig), 1)]; %; zeros((10*param.light_length(cond)),1)
    xxlength = signal_length-size(LED(cond).ON_outsig,1);
    LED(cond).ON_outsig = [LED(cond).ON_outsig; zeros(xxlength,1)];
    
    a = size(LED(cond).ON_outsig);
    LED(cond).OFF_outsig = zeros(a(1), 1);
    
end

LED(1).light_length = param.light_length(1);
LED(1).ON_outsig = LED(2).OFF_outsig;
LED(1).OFF_outsig = LED(2).OFF_outsig;

% ---- Basler trigger data ---- %
basler_volts = param.basler_volts;
basler_outsig = zeros(a(1), 1);
basler_rate = round(samplerate/param.Basler_fps);
basler_outsig(1:basler_rate:end) = basler_volts;

basler_outsig(2:basler_rate:end) = basler_volts;


clear a xxlength

%error checking analogue output signal
if size(LED(3).ON_outsig) == size(basler_outsig)
    fprintf('\n Laser signal and camera signal aligned \n')
else 
    warndlg('Mismatched laser and camera signal')
    return
end

% ---- COMBINED DATA QUEUE ---- %
for cond = 1:param.num_conds
    switch conds_matrix(cond).opto
        case param.light_length(1)
            conds_matrix(cond).LED = [LED(1).ON_outsig, basler_outsig];
        case param.light_length(2)
            conds_matrix(cond).LED = [LED(2).ON_outsig, basler_outsig]; 
        case param.light_length(3)
            conds_matrix(cond).LED = [LED(3).ON_outsig, basler_outsig]; 
        case param.light_length(4)
            conds_matrix(cond).LED = [LED(4).ON_outsig, basler_outsig]; 
        case param.light_length(5)
            conds_matrix(cond).LED = [LED(5).ON_outsig, basler_outsig]; 
        case param.light_length(6)
            conds_matrix(cond).LED = [LED(6).ON_outsig, basler_outsig]; 
        case param.light_length(7)
            conds_matrix(cond).LED = [LED(7).ON_outsig, basler_outsig];
    end
end



end
% for cond = 1:param.num_conds
%     conds_matrix(cond).LED = [LED(3).ON_outsig, basler_outsig]; 
% end


% figure;
% for ii = 1:8
%     subplot(2,4,ii)
%     plot(param.conds_matrix(ii).LED(:,2))
% end


% signal_length = s_vid_light.Rate*param.basler_length;
% 
% for jj = 2:length(param.light_length)
%     
%     % ---- LED trigger data ---- %
%     %ratio: ON:OFF
%     LED(jj).frequency = 1200; %hz of LED frequency
%     LED(jj).ratio_num = 0.5; % 1:num ratio of light off
%     LED(jj).rate = round(s_vid_light.Rate/LED(jj).frequency); %adjustment ratio to get LED signal speed correct 
%     LED(jj).DutyCycle = 1/(1+LED(jj).rate);
%     LED(jj).light_length = param.light_length(jj);
%       
%     % creating the LED trigger signal
%     LED(jj).pulse = [(ones(LED(jj).rate,1)*param.LED_intensity); zeros(round(LED(jj).ratio_num*LED(jj).rate),1)];
%     LED(jj).pulse_length = length(LED(jj).pulse);
%     LED(jj).pulse_num = round(s_vid_light.Rate/LED(jj).pulse_length)*param.light_length(jj); %should equal the desired light Hz  
%     LED(jj).sig = LED(jj).pulse;
%     %concatenate the individual pulses to reach the total length of light exposure
%     for ii = 1:(LED(jj).pulse_num-1)
%         LED(jj).sig = [LED(jj).sig; LED(jj).pulse];
%     end  
%     
%     LED(jj).sig(end-(2*LED(jj).pulse_length-1):end) = 0; %two units of off at the end of the signal
%     LED(jj).pre_sig = zeros(s_vid_light.Rate*param.basler_delay,1); %actual signal
%     LED(jj).post_sig = param.basler_length-param.basler_delay-param.light_length(jj); %time (sec) of post LED w/basler on
% 
%     LED(jj).ON_outsig = [LED(jj).pre_sig; LED(jj).sig; zeros(round(s_vid_light.Rate*LED(jj).post_sig), 1)]; %; zeros((10*param.light_length(jj)),1)
%     xxlength = signal_length-size(LED(jj).ON_outsig,1);
%     LED(jj).ON_outsig = [LED(jj).ON_outsig; zeros(xxlength,1)];
%     
%     a = size(LED(jj).ON_outsig);
%     LED(jj).OFF_outsig = zeros(a(1), 1);
%     
% end
% 
% LED(1).light_length = param.light_length(1);
% LED(1).ON_outsig = LED(2).OFF_outsig;
% LED(1).OFF_outsig = LED(2).OFF_outsig;
% 
% % ---- Basler trigger data ---- %
% basler_volts = 9;
% basler_outsig = zeros(a(1), 1);
% basler_rate = round(s_vid_light.Rate/param.Basler_fps);
% basler_outsig(1:basler_rate:end) = basler_volts;
% clear a xxlength
% 
% %error checking analogue output signal
% if size(LED(3).ON_outsig) == size(basler_outsig)
%     fprintf('\n Laser signal and camera signal aligned \n')
% else 
%     warndlg('Mismatched laser and camera signal')
%     return
% end
% 
% % ---- COMBINED DATA QUEUE ---- %
% for jj = 1:param.num_conds
%     switch conds_matrix(jj).opto
%         case param.light_length(1)
%             conds_matrix(jj).LED = [LED(1).ON_outsig, basler_outsig];
%         case param.light_length(2)
%             conds_matrix(jj).LED = [LED(2).ON_outsig, basler_outsig]; 
%         case param.light_length(3)
%             conds_matrix(jj).LED = [LED(3).ON_outsig, basler_outsig]; 
%         case param.light_length(4)
%             conds_matrix(jj).LED = [LED(4).ON_outsig, basler_outsig]; 
%         case param.light_length(5)
%             conds_matrix(jj).LED = [LED(5).ON_outsig, basler_outsig]; 
%         case param.light_length(6)
%             conds_matrix(jj).LED = [LED(6).ON_outsig, basler_outsig]; 
%         case param.light_length(7)
%             conds_matrix(jj).LED = [LED(7).ON_outsig, basler_outsig];
%     end
% end
% clear ind



