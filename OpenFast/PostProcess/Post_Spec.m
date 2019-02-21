% Spectral analysis

fastout = Post_LoadFastOut('/Users/nabbas/Documents/TurbineModels/DTU_10MW/DTU10MWRWT_NAUTILUS_GoM_FAST_v1.00/Baseline/DTU_10MW_NAUTILUS_GoM.out')
% 

% data = simout;
signal = fastout;
channel = 'TwrBsFyt';
% x = data.Wave1Elev;
x = signal.(channel)(end/4:end);
% % x = data.TTDspFA;
% x = data.PtfmPitch; %(1001:end);
% x = data.GenSpeed;
% x = p2.data;

fs = 1/(signal.Time(2) - signal.Time(1));
% fs = 80;
N = length(x);
xdft = fft(x(2:end));
xdft = xdft(1:N/2+1);
psd = (1/(fs*N)) * abs(xdft).^2;
psd(2:end-1) = 2*psd(2:end-1);
freq = 0:fs/N:fs/2;


psdN = (1/(2*pi*N)) * abs(xdft).^2;
psdN(2:end-1) = 2*psdN(2:end-1);
freqN = 0:(2*pi)/N:pi;


figure(1000);
% plot(freq,10*log10(psd))
semilogx(freq,10*log10(psd))
% plot(freqN/pi,10*log10(psdN))
hold on
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')

figure(1001);
myplot(signal.Time, signal.(channel))
hold on
%%

% Fs = 1000;
% t = 0:1/Fs:1-1/Fs;
% x = cos(2*pi*100*t); %+ randn(size(t));
% N = length(x);
% xdft = fft(x);
% xdft = xdft(1:N/2+1);
% psdx = (1/(Fs*N)) * abs(xdft).^2;
% psdx(2:end-1) = 2*psdx(2:end-1);
% freq = 0:Fs/N:Fs/2;
% 
% plot(freq,10*log10(psdx))
% grid on
% title('Periodogram Using FFT')
% xlabel('Frequency (Hz)')
% ylabel('Power/Frequency (dB/Hz)')