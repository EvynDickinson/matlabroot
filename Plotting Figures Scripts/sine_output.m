s = daq.createSession('ni');
addAnalogOutputChannel(s,'Dev1',0,'Voltage');

s.Rate = 1000;

amplitude = 1


outputSignal1 = amplitude*sin(linspace(0,pi*2,s.Rate)');

outputSignal1 = 5*ones(1,s.Rate*5)';
plot(outputSignal1);
xlabel('Time');
ylabel('Voltage');
legend('Analog Output 0', 'Analog Output 1');
queueOutputData(s,[outputSignal1; outputSignal1; outputSignal1]);

s.startBackground;

    
    Panel_com('stop');%% start the pattern
  Panel_com('set_pattern_id', 2);   
    Panel_com('set_position', [48 1]);
    Panel_com('send_gain_bias', [10, -40, 0, 0])
        Panel_com('set_mode', [1 0 ]);     
        
            Panel_com('start');%% start the pattern
            ii
            pause(2);
