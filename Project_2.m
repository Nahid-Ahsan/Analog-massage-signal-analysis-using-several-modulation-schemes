clc; close all; clear all;
% ------------------------------------------------------------
%% Sampling
% Define Message Signal
x = @(t) 5*cos(20*pi*t) + 4*sin(10*pi*t);
Fs = 100; %Sampling Frequency
t = 0:0.0001:0.5; % Time Grid
ts = 0:1/Fs:0.1; % Sampling Time Grid
xs = x(ts); % Sampled Signal
figure(1)
subplot(2,1,1), plot(t, x(t), 'linewidth', 1);
title('Sampling Process');
ylabel('Message signal, x(t)');
subplot(2,1,2), stem(ts, xs, '.', 'linewidth', 1.5);
ylabel('Sampled Message Signal, x_s(n)');
xlabel('Time (sec)')
saveas(gcf, 'project_2fig1.png')
% ------------------------------------------------------------
compressed = compand(xs,255,max(xs),'mu/compressor')
% ------------------------------------------------------------
%% Quantization (Rounding)
L = 256; b = log2(L); %Quantization Levels and bits
delta = (max(xs)-min(xs))/(L-1); % Quantization Step
xq = round((xs/delta)*delta); % Quantized Signal
SQNR = 10*log10(mean(xs.^2)/mean((xs-xq).^2));
disp(['SQNR = ', num2str(SQNR), ' dB']);
figure(2)
subplot(2,1,1), stem(ts, xs, '.', 'linewidth', 1.5);
title('Quantization Process');
ylabel('Sampled Signal, x_s(n)');
subplot(2,1,2), stem(ts, xq, '.', 'linewidth', 1.5);
ylabel('Quantized Signal, x_q(n)');
saveas(gcf, 'project_2fig2.png')
%% PCM Encoding
xb = de2bi(round((xq/delta)+ 41)); %Assigning bits for each Q-level
xbs = reshape(xb, [1, length(xq)*b]); %PCM Bitstream
figure(3)
stairs(xbs, 'linewidth', 1);
title('PCM Encoding');
ylabel('PCM Bitstream');
xlim([0, length(xbs)+1]);
ylim([-0.05, 1.05]);
saveas(gcf, 'project_2fig3.png')
expanded = compand(xs,255,max(xs),'mu/expander')
%%
x = xbs;
bp=.000001;                                                    % bit period
disp(' Binary information at Transmitter :');
disp(x);

bit=[]; 
for n=1:1:length(x)
    if x(n)==1;
       se=ones(1,100);
    else x(n)==0;
        se=zeros(1,100);
    end
     bit=[bit se];
end
t1=bp/100:bp/100:100*length(x)*(bp/100);
figure(4)
plot(t1,bit,'lineWidth',2.5);grid on;
axis([ 0 bp*length(x) -.5 1.5]);
ylabel('amplitude(volt)');
xlabel(' time(sec)');
title('transmitting information as digital signal');
saveas(gcf, 'project_2fig4.png')

A=5;                                          % Amplitude of carrier signal
br=1/bp;                                                         % bit rate
f1=br*8;                           % carrier frequency for information as 1
f2=br*2;                           % carrier frequency for information as 0
t2=bp/99:bp/99:bp;                 
ss=length(t2);
fsk=[];
for (i=1:1:length(xbs))
    if (xbs(i)==1)
        y=A*cos(2*pi*f1*t2);
    else
        y=A*cos(2*pi*f2*t2);
    end
    fsk=[fsk y];
end
t3=bp/99:bp/99:bp*length(xbs);
figure(5)
plot(t3,fsk);
xlabel('time(sec)');
ylabel('amplitude(volt)');
title('waveform for binary FSK modulation coresponding binary information');
saveas(gcf, 'project_2fig5.png')