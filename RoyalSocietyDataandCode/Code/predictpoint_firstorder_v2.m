function [tpredict,etapredict] = predictpoint_firstorder_v2(wavenums,depth,omegas,forwardtime,focalprobe,pagenow,timebacktrackindex,fftstorage,fftcounter,tnow)
%PREDICTPOINT Summary of this function goes here
%   Detailed explanation goes here

        tpredictindex=pagenow-timebacktrackindex(focalprobe,:);

        
        
        for j=1:length(tpredictindex)
            if tpredictindex(j)<1
                tpredictindex(j)=length(fftcounter)+tpredictindex(j);
            end
        end
        
        
        
        complexamplitudes=zeros(length(omegas),1);
        

        
        
        for j=1:length(complexamplitudes)
            complexamplitudes(j)=fftstorage(j,tpredictindex(j))./fftcounter(tpredictindex(j));
        end
        
        
        
%         etapredict=real(exp(1i*(omegas*(forwardtime(focalprobe))))*complexamplitudes);
        etapredict=real(sum(complexamplitudes));
        tpredict=tnow+forwardtime(focalprobe);
        

end

