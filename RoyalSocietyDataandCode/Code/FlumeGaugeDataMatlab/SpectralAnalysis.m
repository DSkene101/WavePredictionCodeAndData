close all
clear all
clc

set(0,'defaultTextInterpreter','latex');

for loadnum=[113,123,133]


load(['00',num2str(loadnum),'.mat'])

switch loadnum
            case 111
                Tp=1.1329;
                Hs=0.0113;
            case 112
                Tp=1.1329;
                Hs=0.0226;
            case 113
                Tp=1.1329;
                Hs=0.0339;
            case 121
                Tp=2.7089;
                Hs=0.01695;
            case 122
                Tp=2.7089;
                Hs=0.0339;
            case 123
                Tp=2.7089;
                Hs=0.05085;
            case 131
                Tp=1.652;
                Hs=0.0226;
            case 132
                Tp=1.652;
                Hs=0.0452;
            case 133
                Tp=1.652;
                Hs=0.0678;
            case 114
                Tp=1.1329;
                Hs=0.0452;
            case 115
                Tp=1.1329;
                Hs=0.0565;
            case 116
                Tp=1.1329;
                Hs=0.00904;
            case 117
                Tp=1.1329;
                Hs=0.0678;
            case 118
                Tp=1.1329;
                Hs=0.0791;
end


safeportion=3;
aa=round(length(WaveGauge)/safeportion);
bb=aa*(safeportion-1);
WaveTime=WaveTime(aa:bb);
WaveGauge=WaveGauge(aa:bb,:);

figure()
plot(WaveTime,WaveGauge(:,1))

fs=128; %Hz
L=length(WaveTime);
w=2*pi*fs*(0:floor(L/2))/L;

Hs_OLD=Hs
Hs=4*(std(WaveGauge(:,1)))

x=WaveGauge(:,1:8)
N = length(x);
xdft = fft(x);
xdft = xdft(1:N/2+1,:);
psdx = (1/(fs*N)) * abs(xdft).^2;
psdx(2:end-1,:) = 2*psdx(2:end-1,:)*2*pi;
psdx=mean(psdx.').'

Yavg=psdx

% Y = fft(WaveGauge(:,1));
% Yavg=Y
% 
% Yavg=mean(Y.')

P=Yavg

avgpoints=13

P=movmean(P,avgpoints);

figure()
plot(w,P(1:L/2+1)) 

[ S, Amp, Phase ] = JONSWAP_TUCKER_NONRAND( w, Hs, Tp)
S=S/(w(2)-w(1));

hold on
plot(w,S)

xlim([0,2*pi/Tp*3])
ylim([0,max(max(P),max(S))*1.2])

title(['Case: ',num2str(loadnum)])
xlabel('$\omega$ (rad\,s$^{-1})$ ')
ylabel('$S(\omega)/\Delta \omega$ (m$^{2}$\,rad\,s$^{-1}$)')

end
