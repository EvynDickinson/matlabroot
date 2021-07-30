
timevec = 0.0001:0.0001:8;


trialnums{1} = [20:24];
trialnums{2} = [25:29];
trialnums{3} = [35:39];
trialnums{4} = [40:44];
trialnums{5} = [45:49];

odornames{1} = 'MD1';
odornames{2} = 'CVA1';
odornames{3} = '1Hex2';
odornames{4} = 'PO';
odornames{5} = 'phenylacetic acid 0.2%';



figure
for i = 1:length(trialnums)
    subplot(5,1,i)
    hold on
    for j = 1:length(trialnums{i})
        load(['Raw_WCwaveform_2018-06-06_E1_' num2str(trialnums{i}(j))]);
        plot(timevec,voltage1)
    end
    xlim([3 8])
    ylim([-80 -50])
    box off
    ylabel('Membrane voltage (mV)')
    xlabel('time (s)')
    plot([4 6], [-52 -52], 'k', 'LineWidth', 2)
    title(odornames{i})
end



trialnums1 = [37:40];
figure
subplot(2,1,1)
for i = 1:4
    %subplot(8,1,i)
    hold on
    load(['Raw_WCwaveform_2018-06-04_E1_' num2str(trialnums1(i))]);
    plot(timevec,voltage1)
    xlim([3 8])
    ylim([-90 -50])
    box off
    ylabel('Membrane voltage (mV)')
    xlabel('time (s)')
end
plot([4 6], [-55 -55], 'b')


trialnums1 = [48:53];
subplot(2,1,2)
for i = 1:6
    %subplot(8,1,i)
    hold on
    load(['Raw_WCwaveform_2018-06-04_E1_' num2str(trialnums1(i))]);
    plot(timevec,voltage1)
    xlim([3 8])
    ylim([-90 -50])
    box off
    ylabel('Membrane voltage (mV)')
    xlabel('time (s)')
end
plot([4 6], [-55 -55], 'b')
