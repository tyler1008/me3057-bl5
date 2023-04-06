clc
clear all
close all


load 'LowMassLowDamping.mat'
load 'LowMassHighDamping.mat'
load 'HighMassLowDamping.mat'
load 'HighMassHighDamping.mat'
load 'Time_s_lmld.mat'

%% part1 for low mass low damping
% segment data with given function
ind = segmentData(time_s,D2_m_lmld);

%create segment struct to store 7 vectors of data
i = ind(1,1):ind(1,2);      % segement 1
seg(1).t = time_s(i); 
seg(1).D1 = D1_m_lmld(i); 
seg(1).D2 = D2_m_lmld(i);
seg(1).F1 = F1_N_lmld(i);
seg(1).F2 = F2_N_lmld(i);

i = ind(2,1):ind(2,2);      % segment 2
seg(2).t = time_s(i); 
seg(2).D1 = D1_m_lmld(i); 
seg(2).D2 = D2_m_lmld(i);
seg(2).F1 = F1_N_lmld(i);
seg(2).F2 = F2_N_lmld(i);

i = ind(3,1):ind(3,2);      % segment 3
seg(3).t = time_s(i); 
seg(3).D1 = D1_m_lmld(i); 
seg(3).D2 = D2_m_lmld(i);
seg(3).F1 = F1_N_lmld(i);
seg(3).F2 = F2_N_lmld(i);

i = ind(4,1):ind(4,2);        % segment 4
seg(4).t = time_s(i); 
seg(4).D1 = D1_m_lmld(i); 
seg(4).D2 = D2_m_lmld(i);
seg(4).F1 = F1_N_lmld(i);
seg(4).F2 = F2_N_lmld(i);

i = ind(5,1):ind(5,2);          % segment 5
seg(5).t = time_s(i); 
seg(5).D1 = D1_m_lmld(i); 
seg(5).D2 = D2_m_lmld(i);
seg(5).F1 = F1_N_lmld(i);
seg(5).F2 = F2_N_lmld(i);

i = ind(6,1):ind(6,2);          % segment 6
seg(6).t = time_s(i); 
seg(6).D1 = D1_m_lmld(i); 
seg(6).D2 = D2_m_lmld(i);
seg(6).F1 = F1_N_lmld(i);
seg(6).F2 = F2_N_lmld(i);

i = ind(7,1):ind(7,2);          % segment 7
seg(7).t = time_s(i); 
seg(7).D1 = D1_m_lmld(i); 
seg(7).D2 = D2_m_lmld(i);
seg(7).F1 = F1_N_lmld(i);
seg(7).F2 = F2_N_lmld(i);

i = ind(8,1):ind(8,2);          % segment 7
seg(8).t = time_s(i); 
seg(8).D1 = D1_m_lmld(i); 
seg(8).D2 = D2_m_lmld(i);
seg(8).F1 = F1_N_lmld(i);
seg(8).F2 = F2_N_lmld(i);

i = ind(9,1):ind(9,2);          % segment 7
seg(9).t = time_s(i); 
seg(9).D1 = D1_m_lmld(i); 
seg(9).D2 = D2_m_lmld(i);
seg(9).F1 = F1_N_lmld(i);
seg(9).F2 = F2_N_lmld(i);


i = ind(10,1):ind(10,2);          % segment 7
seg(10).t = time_s(i); 
seg(10).D1 = D1_m_lmld(i); 
seg(10).D2 = D2_m_lmld(i);
seg(10).F1 = F1_N_lmld(i);
seg(10).F2 = F2_N_lmld(i);

i = ind(11,1):ind(11,2);          % segYment 7
seg(11).t = time_s(i); 
seg(11).D1 = D1_m_lmld(i); 
seg(11).D2 = D2_m_lmld(i);
seg(11).F1 = F1_N_lmld(i);
seg(11).F2 = F2_N_lmld(i);

% initialize empty vectors
gain_v = []
phaseLag_v = []
freq_v = []

for i = 1:11
    segD1_0 = seg(i).D1 - mean(seg(i).D1);       % zero segment vectors
    segD2_0 = seg(i).D2 - mean(seg(i).D2);
    
    [gain,phaseLag] = func21(seg(i).t,segD1_0,segD2_0); % calc gain and phase lag using func21
    [fit22] = fit(seg(i).t,segD1_0,'sin1') ;            
    freq = fit22.b1 / (2*pi);
    % append vectors
    gain_v = [gain_v, gain];
    phaseLag_v = [phaseLag_v,  phaseLag];
    freq_v = [freq_v,  freq];
end


box on; grid on; hold off;
%% Part2 for Lowmass high damping
ind2 = segmentData(time_s_lmhd,D2_N_lmhd);

%create segment struct to store 7 vectors of data
i = ind2(1,1):ind2(1,2);      % segement 1
seg(1).t = time_s(i); 
seg(1).D1 = D1_N_lmhd(i); 
seg(1).D2 = D2_N_lmhd(i);
seg(1).F1 = F1_N_lmhd(i);
seg(1).F2 = F2_N_lmhd(i);

i = ind2(2,1):ind2(2,2);      % segment 2
seg(2).t = time_s(i); 
seg(2).D1 = D1_N_lmhd(i); 
seg(2).D2 = D2_N_lmhd(i);
seg(2).F1 = F1_N_lmhd(i);
seg(2).F2 = F2_N_lmhd(i);

i = ind2(3,1):ind2(3,2);      % segment 3
seg(3).t = time_s(i); 
seg(3).D1 = D1_N_lmhd(i); 
seg(3).D2 = D2_N_lmhd(i);
seg(3).F1 = F1_N_lmhd(i);
seg(3).F2 = F2_N_lmhd(i);

i = ind2(4,1):ind2(4,2);        % segment 4
seg(4).t = time_s(i); 
seg(4).D1 = D1_N_lmhd(i); 
seg(4).D2 = D2_N_lmhd(i);
seg(4).F1 = F1_N_lmhd(i);
seg(4).F2 = F2_N_lmhd(i);

i = ind2(5,1):ind2(5,2);          % segment 5
seg(5).t = time_s(i); 
seg(5).D1 = D1_N_lmhd(i); 
seg(5).D2 = D2_N_lmhd(i);
seg(5).F1 = F1_N_lmhd(i);
seg(5).F2 = F2_N_lmhd(i);

i = ind2(6,1):ind2(6,2);          % segment 6
seg(6).t = time_s(i); 
seg(6).D1 = D1_N_lmhd(i); 
seg(6).D2 = D2_N_lmhd(i);
seg(6).F1 = F1_N_lmhd(i);
seg(6).F2 = F2_N_lmhd(i);

i = ind2(7,1):ind2(7,2);          % segment 7
seg(7).t = time_s(i); 
seg(7).D1 = D1_N_lmhd(i); 
seg(7).D2 = D2_N_lmhd(i);
seg(7).F1 = F1_N_lmhd(i);
seg(7).F2 = F2_N_lmhd(i);

i = ind2(8,1):ind2(8,2);          % segment 7
seg(8).t = time_s(i); 
seg(8).D1 = D1_N_lmhd(i); 
seg(8).D2 = D2_N_lmhd(i);
seg(8).F1 = F1_N_lmhd(i);
seg(8).F2 = F2_N_lmhd(i);

i = ind2(9,1):ind2(9,2);          % segment 7
seg(9).t = time_s(i); 
seg(9).D1 = D1_N_lmhd(i); 
seg(9).D2 = D2_N_lmhd(i);
seg(9).F1 = F1_N_lmhd(i);
seg(9).F2 = F2_N_lmhd(i);


i = ind2(10,1):ind2(10,2);          % segment 7
seg(10).t = time_s(i); 
seg(10).D1 = D1_N_lmhd(i); 
seg(10).D2 = D2_N_lmhd(i);
seg(10).F1 = F1_N_lmhd(i);
seg(10).F2 = F2_N_lmhd(i);

i = ind2(11,1):ind2(11,2);          % segYment 7
seg(11).t = time_s(i); 
seg(11).D1 = D1_N_lmhd(i); 
seg(11).D2 = D2_N_lmhd(i);
seg(11).F1 = F1_N_lmhd(i);
seg(11).F2 = F2_N_lmhd(i);

% initialize empty vectors
gain_v2 = []
phaseLag_v2 = []
freq_v2 = []

for i = 1:11
    segD12_0 = seg(i).D1 - mean(seg(i).D1);       % zero segment vectors
    segD22_0 = seg(i).D2 - mean(seg(i).D2);
    
    [gain2,phaseLag2] = func21(seg(i).t,segD12_0,segD22_0) % calc gain and phase lag using func21
    [fit22] = fit(seg(i).t,segD12_0,'sin1');             
    freq2 = fit22.b1 / (2*pi);
    % append vectors
    gain_v2 = [gain_v2, gain2];
    phaseLag_v2 = [phaseLag_v2,  phaseLag2];
    freq_v2 = [freq_v2,  freq2];
end


%% Part3 for High mas low damping
ind3 = segmentData(T1_s_hmld,D2_m_hmld);

%create segment struct to store 7 vectors of data
i = ind3(1,1):ind3(1,2);      % segement 1
seg(1).t = time_s(i); 
seg(1).D1 = D1_m_hmld(i); 
seg(1).D2 = D2_m_hmld(i);
seg(1).F1 = F1_N_hmld(i);
seg(1).F2 = F2_N_hmld(i);

i = ind3(2,1):ind3(2,2);      % segment 2
seg(2).t = time_s(i); 
seg(2).D1 = D1_m_hmld(i); 
seg(2).D2 = D2_m_hmld(i);
seg(2).F1 = F1_N_hmld(i);
seg(2).F2 = F2_N_hmld(i);

i = ind3(3,1):ind3(3,2);      % segment 3
seg(3).t = time_s(i); 
seg(3).D1 = D1_m_hmld(i); 
seg(3).D2 = D2_m_hmld(i);
seg(3).F1 = F1_N_hmld(i);
seg(3).F2 = F2_N_hmld(i);

i = ind3(4,1):ind3(4,2);        % segment 4
seg(4).t = time_s(i); 
seg(4).D1 = D1_m_hmld(i); 
seg(4).D2 = D2_m_hmld(i);
seg(4).F1 = F1_N_hmld(i);
seg(4).F2 = F2_N_hmld(i);

i = ind3(5,1):ind3(5,2);          % segment 5
seg(5).t = time_s(i); 
seg(5).D1 = D1_m_hmld(i); 
seg(5).D2 = D2_m_hmld(i);
seg(5).F1 = F1_N_hmld(i);
seg(5).F2 = F2_N_hmld(i);

i = ind3(6,1):ind3(6,2);          % segment 6
seg(6).t = time_s(i); 
seg(6).D1 = D1_m_hmld(i); 
seg(6).D2 = D2_m_hmld(i);
seg(6).F1 = F1_N_hmld(i);
seg(6).F2 = F2_N_hmld(i);

i = ind3(7,1):ind3(7,2);          % segment 7
seg(7).t = time_s(i); 
seg(7).D1 = D1_m_hmld(i); 
seg(7).D2 = D2_m_hmld(i);
seg(7).F1 = F1_N_hmld(i);
seg(7).F2 = F2_N_hmld(i);

i = ind3(8,1):ind3(8,2);          % segment 7
seg(8).t = time_s(i); 
seg(8).D1 = D1_m_hmld(i); 
seg(8).D2 = D2_m_hmld(i);
seg(8).F1 = F1_N_hmld(i);
seg(8).F2 = F2_N_hmld(i);

i = ind3(9,1):ind3(9,2);          % segment 7
seg(9).t = time_s(i); 
seg(9).D1 = D1_m_hmld(i); 
seg(9).D2 = D2_m_hmld(i);
seg(9).F1 = F1_N_hmld(i);
seg(9).F2 = F2_N_hmld(i);


i = ind3(10,1):ind3(10,2);          % segment 7
seg(10).t = time_s(i); 
seg(10).D1 = D1_m_hmld(i); 
seg(10).D2 = D2_m_hmld(i);
seg(10).F1 = F1_N_hmld(i);
seg(10).F2 = F2_N_hmld(i);

i = ind3(11,1):ind3(11,2);          % segYment 7
seg(11).t = time_s(i); 
seg(11).D1 = D1_m_hmld(i); 
seg(11).D2 = D2_m_hmld(i);
seg(11).F1 = F1_N_hmld(i);
seg(11).F2 = F2_N_hmld(i);

% initialize empty vectors
gain_v3 = []
phaseLag_v3 = []
freq_v3 = []

for i = 1:11
    segD13_0 = seg(i).D1 - mean(seg(i).D1);       % zero segment vectors
    segD23_0 = seg(i).D2 - mean(seg(i).D2);
    
    [gain3,phaseLag3] = func21(seg(i).t,segD13_0,segD23_0) % calc gain and phase lag using func21
    [fit23] = fit(seg(i).t,segD13_0,'sin1');             
    freq3 = fit23.b1 / (2*pi);
    % append vectors
    gain_v3 = [gain_v3, gain3];
    phaseLag_v3 = [phaseLag_v3,  phaseLag3];
    freq_v3 = [freq_v3,  freq3];
end
%% Part4 for Highmass high damping
ind4 = segmentData(time_s_hmhd,D2_N_hmhd);

%create segment struct to store 7 vectors of data
i = ind4(1,1):ind4(1,2);      % segement 1
seg(1).t = time_s_hmhd(i); 
seg(1).D1 = D1_N_hmhd(i); 
seg(1).D2 = D2_N_hmhd(i);
seg(1).F1 = F1_N_hmhd(i);
seg(1).F2 = F2_N_hmhd(i);

i = ind4(2,1):ind4(2,2);      % segment 2
seg(2).t = time_s_hmhd(i); 
seg(2).D1 = D1_N_hmhd(i); 
seg(2).D2 = D2_N_hmhd(i);
seg(2).F1 = F1_N_hmhd(i);
seg(2).F2 = F2_N_hmhd(i);

i = ind4(3,1):ind4(3,2);      % segment 3
seg(3).t = time_s_hmhd(i); 
seg(3).D1 = D1_N_hmhd(i); 
seg(3).D2 = D2_N_hmhd(i);
seg(3).F1 = F1_N_hmhd(i);
seg(3).F2 = F2_N_hmhd(i);

i = ind4(4,1):ind4(4,2);        % segment 4
seg(4).t = time_s_hmhd(i); 
seg(4).D1 = D1_N_hmhd(i); 
seg(4).D2 = D2_N_hmhd(i);
seg(4).F1 = F1_N_hmhd(i);
seg(4).F2 = F2_N_hmhd(i);

i = ind4(5,1):ind4(5,2);          % segment 5
seg(5).t = time_s_hmhd(i); 
seg(5).D1 = D1_N_hmhd(i); 
seg(5).D2 = D2_N_hmhd(i);
seg(5).F1 = F1_N_hmhd(i);
seg(5).F2 = F2_N_hmhd(i);

i = ind4(6,1):ind4(6,2);          % segment 6
seg(6).t = time_s_hmhd(i); 
seg(6).D1 = D1_N_hmhd(i); 
seg(6).D2 = D2_N_hmhd(i);
seg(6).F1 = F1_N_hmhd(i);
seg(6).F2 = F2_N_hmhd(i);

i = ind4(7,1):ind4(7,2);          % segment 7
seg(7).t = time_s_hmhd(i); 
seg(7).D1 = D1_N_hmhd(i); 
seg(7).D2 = D2_N_hmhd(i);
seg(7).F1 = F1_N_hmhd(i);
seg(7).F2 = F2_N_hmhd(i);

i = ind4(8,1):ind4(8,2);          % segment 7
seg(8).t = time_s_hmhd(i); 
seg(8).D1 = D1_N_hmhd(i); 
seg(8).D2 = D2_N_hmhd(i);
seg(8).F1 = F1_N_hmhd(i);
seg(8).F2 = F2_N_hmhd(i);

i = ind4(9,1):ind4(9,2);          % segment 7
seg(9).t = time_s_hmhd(i); 
seg(9).D1 = D1_N_hmhd(i); 
seg(9).D2 = D2_N_hmhd(i);
seg(9).F1 = F1_N_hmhd(i);
seg(9).F2 = F2_N_hmhd(i);


i = ind4(10,1):ind4(10,2);          % segment 7
seg(10).t = time_s_hmhd(i); 
seg(10).D1 = D1_N_hmhd(i); 
seg(10).D2 = D2_N_hmhd(i);
seg(10).F1 = F1_N_hmhd(i);
seg(10).F2 = F2_N_hmhd(i);

i = ind4(11,1):ind4(11,2);          % segYment 7
seg(11).t = time_s_hmhd(i); 
seg(11).D1 = D1_N_hmhd(i); 
seg(11).D2 = D2_N_hmhd(i);
seg(11).F1 = F1_N_hmhd(i);
seg(11).F2 = F2_N_hmhd(i);

% initialize empty vectors
gain_v4 = []
phaseLag_v4 = []
freq_v4 = []

for i = 1:11
    segD14_0 = seg(i).D1 - mean(seg(i).D1);       % zero segment vectors
    segD24_0 = seg(i).D2 - mean(seg(i).D2);
    
    [gain4,phaseLag4] = func21(seg(i).t,segD14_0,segD24_0) % calc gain and phase lag using func21
    [fit24] = fit(seg(i).t,segD14_0,'sin1');             
    freq4 = fit24.b1 / (2*pi);
    % append vectors
    gain_v4 = [gain_v4, gain4];
    phaseLag_v4 = [phaseLag_v4,  phaseLag4];
    freq_v4 = [freq_v4,  freq4];
end



%% plot part

% create bode plot
subplot(2,1,1);
hold on
plot(freq_v, 20*log10(gain_v))
plot(freq_v2, 20*log10(gain_v2))
plot(freq_v3, 20*log10(gain_v3))
plot(freq_v4, 20*log10(gain_v4))
box on; grid on
title('Gain  vs. Frequency','FontSize',20)
xlabel('Frequency (Hz)','FontSize',15)
xlim([1,6])
ylabel('Gain (dB)','FontSize',15)

legend('Low Mass Low Damper',' Low Mass High Damper','High Mass Low Damper',' High Mass High Damper')
hold off

subplot(2,1,2);
hold on
plot(freq_v, phaseLag_v)
plot(freq_v2, phaseLag_v2)
plot(freq_v3, phaseLag_v3)
plot(freq_v4, phaseLag_v4)
box on; grid on
title('Phase Lag vs. Frequency','FontSize',20)
xlabel('Frequency (Hz)','FontSize',15)
xlim([1,6])
ylabel('Phase Angle (deg)','FontSize',15)
legend('Low Mass Low Damper',' Low Mass High Damper','High Mass Low Damper',' High Mass High Damper')
hold off




%%
function [amp,freq, phase, offset, sinFit, stat] = sinFunc(time_v,data_v)
    dataM = mean(data_v);                                        % average data
    dataZERO = data_v - dataM;                                   % subtract average from each datapoint
    offset = time_v(end);
    timeZERO = time_v - offset;                                  % zero time data by subtracting offset
    [sinFit, stat] = fit(timeZERO, dataZERO, 'sin1');            
    amp = sinFit.a1;
    freq = sinFit.b1 / (2*pi);
    phase = rad2deg(sinFit.c1);
end

function [G, pLag] = func21(time_v, motor_v, cart_v)
    % use sinFunc defined above to obtain amplitude and phase
    [ampM, ~, phaseM, ~, ~, ~] = sinFunc(time_v,motor_v);    
    [ampC, ~, phaseC, ~, ~, ~] = sinFunc(time_v,cart_v);
    G = ampC/ampM;
    pLag = mod((phaseC - phaseM), -360);
end