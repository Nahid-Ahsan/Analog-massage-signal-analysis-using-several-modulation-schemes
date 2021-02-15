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
saveas(gcf, 'project_1fig1.png')
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
saveas(gcf, 'project_1fig2.png')
%% PCM Encoding
xb = de2bi(round((xq/delta)+ 41)); %Assigning bits for each Q-level
xbs = reshape(xb, [1, length(xq)*b]); %PCM Bitstream
figure(3)
stairs(xbs, 'linewidth', 1);
title('PCM Encoding');
ylabel('PCM Bitstream');
xlim([0, length(xbs)+1]);
ylim([-0.05, 1.05]);
saveas(gcf, 'project_1fig3.png')
expanded = compand(xs,255,max(xs),'mu/expander')