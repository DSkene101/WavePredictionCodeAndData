function [wavenums,c,dist,timebacktrack,timebacktrackindex] = disttimespeed_firstorder_v2(omegas,depth,g,xs,x1,t,forwardtime)
%DISTTIMESPEED Summary of this function goes here
%   Detailed explanation goes here
%get wavenumbers for each frequency
    wavenums=wavedispersion_firstorder(omegas,depth,g)

    %get phase speed of these frequencies
    c=omegas./wavenums %sqrt(g./wavenums.*tanh(wavenums*depth))
    
%     cg=g*(tanh(depth*wavenums)+depth*wavenums.*(sech(depth*wavenums)).^2)./(2*sqrt(g*wavenums.*tanh(depth*wavenums)))

%     c=cg
    
    
    
    %determine distances and therefore phase speed backtracking time between probes
    dist=xs-x1;

    timebacktrack=zeros(length(dist),length(c))
    
    
    
    for i=1:length(dist)
        timebacktrack(i,:)=dist(i)./c-forwardtime(i);
    end

    %determine index to go back in time series for each phase
    timebacktrackindex=timebacktrack*0;
    [a,b]=size(timebacktrack);
    for i=1:a*b
        [~, timebacktrackindex(i)] = min(abs(t-timebacktrack(i)));
    end
    
    
end

