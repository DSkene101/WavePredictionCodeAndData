function [] = spacialtemporalplot(c,twindowwidth,dist,forwardtime,fignumber)
%SPACIALTEMPORALPLOT Summary of this function goes here
%   Detailed explanation goes here

%min and max phase speeds
cmin=min(c);
cmax=max(c);

%point in time and space where slowest and fastest waves meet
tmeet=cmin*twindowwidth/(cmax-cmin);
xmeet=cmax*tmeet;

tfastest=0:tmeet/100:tmeet;
tslowest=-twindowwidth:tmeet/100:tmeet;

yfastest=cmax*tfastest;
yslowest=cmin*(tslowest+twindowwidth);

figure(fignumber)

subplot(2,1,1);
xpoints=[-twindowwidth,0,tmeet];
ypoints=[0,0,xmeet];
h=patch(xpoints,ypoints,[0,0,0]);
h.FaceAlpha=0.15;

hold on
plot(forwardtime,dist,'x')
xlabel('t (s)')
ylabel('x (m)')
title('full prediction zone')

subplot(2,1,2)
xpoints=[-twindowwidth,0,tmeet];
ypoints=[0,0,xmeet];
h=patch(xpoints,ypoints,[0,0,0]);
h.FaceAlpha=0.15;
hold on
plot(forwardtime,dist,'x')
axis([0,tmeet,0,xmeet])
xlabel('t (s)')
ylabel('x (m)')
title('prediction zone with t>0')
grid on

plot(tfastest,yfastest,'k')
plot(tslowest,yslowest,'k')


end

