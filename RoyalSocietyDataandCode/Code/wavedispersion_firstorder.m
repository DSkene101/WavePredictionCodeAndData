function [Wavenum] = wavedispersion_firstorder(Omega,depth,g)
%WAVEDISPERSION Summary of this function goes here
%   Detailed explanation goes here
Wavenum=Omega*0;
for i=1:length(Omega)
    dispersionfunc=@(k) k*tanh(k*depth)-Omega(i)^2/g;
    Wavenum(i)=fzero(dispersionfunc,Omega(i)^2/g);
    if Wavenum(i)<0
        Wavenum(i)=-Wavenum(i);
    end
end
end

