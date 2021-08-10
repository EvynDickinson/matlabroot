

%% SET UP DAQ SESSIONS FOR VIDEOS AND LASER
% LED|Basler session
s_vid_light = daq.createSession('ni');     
s_vid_light.Rate = 10000;

% add analog output channels for LED|Basler
addAnalogOutputChannel(s_vid_light,'Dev1', 'ao0', 'Voltage'); %LASER output
addAnalogOutputChannel(s_vid_light,'Dev1', 'ao1', 'Voltage'); %Basler outputs


LED_intensity = 4;
light_length = 0.5;    
basler_length = 2;
basler_delay = 0.5;
Basler_fps = 300;

% ---- LASER trigger data ---- %
%ratio: ON:OFF
LED.frequency = 1200; %hz of LED frequency
LED.ratio_num = 0.5; % 1:num ratio of light on:off
LED.rate = round(s_vid_light.Rate/ LED.frequency); %adjustment ratio to get LED signal speed correct
% creating the LED trigger signal
LED.pulse = [(ones(LED.rate,1)*LED_intensity); zeros(round(LED.ratio_num*LED.rate),1)];
LED.pulse_length = length(LED.pulse);
LED.pulse_num = round(s_vid_light.Rate/LED.pulse_length)*light_length; %should equal the desired light Hz  
LED.sig = LED.pulse;
%concatenate the individual pulses to reach the total length of light exposure
for ii = 1:(LED.pulse_num-1)
    LED.sig = [LED.sig; LED.pulse];
end
LED.sig(end-(2*LED.pulse_length-1):end) = 0; %two units of off at the end of the signal
LED.pre_sig = zeros(s_vid_light.Rate*basler_delay,1);
LED.post_sig = basler_length-basler_delay-light_length; %sec post LED w/basler on
temp.A = zeros(s_vid_light.Rate*basler_delay,1);
temp.B = LED.sig;
temp.C = zeros(s_vid_light.Rate*LED.post_sig, 1);
temp.D = zeros((10*light_length),1);
LED.ON_outsig = [temp.A; temp.B; temp.C; temp.D];  clear temp;
a = size(LED.ON_outsig);
LED.OFF_outsig = zeros(a(1), 1);

% ---- Basler trigger data ---- %
basler_volts = 9;
basler_outsig = zeros(a(1), 1);
basler_rate = round(s_vid_light.Rate/Basler_fps);
basler_outsig(1:basler_rate:end) = basler_volts;
outsig = [LED.ON_outsig, basler_outsig];

 
    
%%
    
    
% Run the program:

queueOutputData(s_vid_light, outsig)     
startBackground(s_vid_light)
% start(A_vid); start(B_vid); start(C_vid); start(D_vid); start(E_vid); start(F_vid);   
    
pause(basler_length * 1.0)   

% stop(A_vid); stop(B_vid); stop(C_vid); stop(D_vid); stop(E_vid); stop(F_vid);
% 
s_vid_light.stop  
    

  















