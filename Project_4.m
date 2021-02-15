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
saveas(gcf, 'project_4fig1.png')
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
saveas(gcf, 'project_4fig2.png')
%% PCM Encoding
xb = de2bi(round((xq/delta)+ 41)); %Assigning bits for each Q-level
xbs = reshape(xb, [1, length(xq)*b]); %PCM Bitstream
figure(3)
stairs(xbs, 'linewidth', 1);
title('PCM Encoding');
ylabel('PCM Bitstream');
xlim([0, length(xbs)+1]);
ylim([-0.05, 1.05]);
saveas(gcf, 'project_4fig3.png')
expanded = compand(xs,255,max(xs),'mu/expander')
%%
x = xbs;
bp= 1;                                                    % bit period
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
saveas(gcf, 'project_4fig4.png')

A=5;                                          % Amplitude of carrier signal
br=1/bp;                                                         % bit rate
f1=br*8;                           % carrier frequency for information as 1
f2=br*2;                           % carrier frequency for information as 0
t2=bp/88:bp/88:bp;                 
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
t3=bp/88:bp/88:bp*length(x);
figure(5)
plot(t3,fsk);
xlabel('time(sec)');
ylabel('amplitude(volt)');
title('waveform for binary FSK modulation coresponding binary information');
saveas(gcf, 'project_4fig5.png')

%Noise generator SNR=Eb/No=20log(Signalrms/Noiserms)
dB = 10;
vn = sqrt(10^(dB/10)); % set noise amplitude = sqrt(Pn) = sqrt(10^(dB/10))
noise = vn*(randn(size(t3))); % noise generator
figure(6)
plot(noise, 'linewidth', 1.25);
grid on;
title('Noise Level');
saveas(gcf, 'project_4fig6.png')

% Noisy Signal
fskn=(fsk+noise); %modulated carrier plus noise
figure(7)
plot(t3,fskn, 'linewidth', 1.25);
title('Modulated Carrier Waveform Plus Noise');
grid on;
saveas(gcf, 'project_4fig7.png')
%%
mn=[];
for n=ss:ss:length(fsk)
  t=bp/88:bp/88:bp;
  y1=cos(2*pi*f1*t);                    % carrier siignal for information 1
  y2=cos(2*pi*f2*t);                    % carrier siignal for information 0
  mm=y1.*fsk((n-(ss-1)):n);
  mmm=y2.*fsk((n-(ss-1)):n);
  t4=bp/88:bp/88:bp;
  z1=trapz(t4,mm)                                             % intregation 
  z2=trapz(t4,mmm)                                            % intregation 
  zz1=round(2*z1/bp)
  zz2= round(2*z2/bp)
  if(zz1>A/2)      % logic lavel= (0+A)/2 or (A+0)/2 or 2.5 ( in this case)
    a=1;
  else(zz2>A/2)
    a=0;
  end
  mn=[mn a];
end
disp(' Binary information at Reciver :');
disp(mn);
%XXXXX Representation of binary information as digital signal which achived 
%after demodulation XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
bit=[];
for n=1:length(mn);
    if mn(n)==1;
       se=ones(1,100);
    else mn(n)==0;
        se=zeros(1,100);
    end
     bit=[bit se];
end
t4=bp/100:bp/100:100*length(mn)*(bp/100);
figure(8)
plot(t4,bit,'LineWidth',2.5);grid on
axis([ 0 bp*length(mn) -.5 1.5])
ylabel('amplitude(volt)')
xlabel(' time(sec)')
title('recived information as digital signal after binary FSK demodulation')
saveas(gcf, 'project_4fig8.png')
%%
figure(9)
subplot(2,1,1)
stairs(xbs, 'linewidth', 1);
title('PCM Encoding');
ylabel('PCM Bitstream');
xlim([0, length(xbs)+1]);
ylim([-0.05, 1.05]);

subplot(2,1,2)
plot(t4,bit,'LineWidth',2.5);grid on
axis([ 0 bp*length(mn) -.5 1.5])
ylabel('amplitude(volt)')
xlabel(' time(sec)')
title('recived information as digital signal after binary FSK demodulation')
saveas(gcf, 'project_4fig9.png')