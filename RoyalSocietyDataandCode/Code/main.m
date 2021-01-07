clear
close all
clc

for testload=[113]%[111,112,113,114,115,116,117,118,121,122,123,131,132,133] %test to load
for samplemultipler=15:1:60 %length of observation (times peak period)
for paddingintegerloop=1  %outdated
for minomegamultiple=0.25 %bandpass lower
for maxomegamultiple=4 %bandpass upper
for avg=0 %outdated
    
    
        switch testload
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
    

        GaugeXs=[0,2.04,4.063,6.127,8.193,10.237,10.811,11.743] %gauge locations

        depth=1.1;

        load(['./FlumeGaugeDataMatlab/00',num2str(testload),'.mat'])
        WaveTime=1:length(WaveGauge);
        [AA,BB]=size(WaveGauge)
        for i=1:BB
            WaveGauge(:,i)=WaveGauge(:,i)-mean(WaveGauge(1:1000,i));
        end
        WaveTime=WaveTime-WaveTime(1);
        fs=128;
        WaveTime=WaveTime/fs;
        
        

        %Set Conditions
        
        depth=1.1;
        g=9.81;
        
        %derive key values
        omegakey=2*pi/Tp;
        wavenumkey=wavedispersion_firstorder(omegakey,depth,g);
        wavelengthkey=2*pi/wavenumkey;

        %Remove entries to speed up algoritm
        divider=2;
        cutout=1:divider:length(WaveTime);
        WaveTime=WaveTime(cutout);
        WaveGauge=WaveGauge(cutout,:);
        fs=128/divider  %1/(WaveTime(10)-WaveTime(9));

        %determine time to predict each place in the future
        x1=0;
        forwardtime=ones(length(GaugeXs(2:end)),1).*[0:0.1:5];
        [AA,BB]=size(forwardtime)
        
        xs=ones(length(GaugeXs(2:end)),BB).*GaugeXs(2:end).'
        forwardtime=forwardtime.'
        forwardtime=forwardtime(:)
        xs=xs.'
        xs=xs(:).'        

        %Vector to store probes
%         probe=WaveGauge(:,[Gauge1s,Gauge2s,Gauge3s]);
        probe=WaveGauge
        
        %vector for time
        % t=WaveTime;
        t=0:length(WaveGauge)
        t=t/fs

        %set sampling window
        samplemultiples=samplemultipler;
        twindowwidth=Tp*samplemultiples;

        %get time indicies for sampling
        tstartindex=1;
        tstart=t(tstartindex);
        tendindex=min(find(t>=(tstart+twindowwidth)));
        trimtime=t(tstartindex:tendindex);
        twidth=tendindex-tstartindex;

        %set padding for fft (1 is no padding)
        paddinginteger=paddingintegerloop;

        %get frequency for single sided fft
        N = length(trimtime)*paddinginteger
        N_2 = ceil(N/2);
        bin_vals = [0 : N-1];
        fax_rads = bin_vals*fs/N*2*pi;
        fax_rads=fax_rads(1:N_2)
        omegas=fax_rads

        %Precondition for band-pass filter
        %removes frequencies above and below multiple of Omega_peak
        omegaminTpmultiple=minomegamultiple%0.25;
        OmegamaxTpmultiple=maxomegamultiple%4;

        omegamin=1/Tp*2*pi*omegaminTpmultiple;
        omegatrimstart=min(find(omegas>=omegamin));
        if isempty(omegatrimstart)
            omegatrimstart=2;
        end

        omegamax=1/Tp*2*pi*OmegamaxTpmultiple;
        omegatrimend=min(find(omegas>=omegamax));
        if isempty(omegatrimend)
            omegatrimend=length(omegas);
        end
        omegas=omegas(omegatrimstart:omegatrimend);
        [~,omegaPindex]=min(abs(omegas-omegakey));

        
        
        %get spacial-temporal relation and relevant variables
        [wavenums,c,dist,timebacktrack,timebacktrackindex] = disttimespeed_firstorder_v2(omegas,depth,g,xs,x1,t,forwardtime)

        
        dummy=(timebacktrack-t(timebacktrackindex))*fs;
        
        bottomportion=timebacktrackindex*0;
        bottomindex=timebacktrackindex*0;
        topportion=timebacktrackindex*0;
        topindex=timebacktrackindex*0;
        
        [a,b]=size(dummy);
        for j=1:a
            for i=1:b
                if dummy(j,i)<0
                    bottomindex(j,i)=timebacktrackindex(j,i)-1;
                    bottomportion(j,i)=-dummy(j,i);

                    topindex(j,i)=timebacktrackindex(j,i);
                    topportion(j,i)=1+dummy(j,i);
                elseif dummy(j,i)>0
                    bottomindex(j,i)=timebacktrackindex(j,i);
                    bottomportion(j,i)=1-dummy(j,i);

                    topindex(j,i)=timebacktrackindex(j,i)+1;
                    topportion(j,i)=dummy(j,i);

                end
            end
        end
        
        
        
        
        %show whether points are in a causal prediction zone
        spacialtemporalplot(c,twindowwidth,dist,forwardtime,1)

        toplot=0;  
        
        if toplot==1
            %create objects for animation plot
            figure(2)
            numrows=length(xs)+1
            subplot(numrows,1,1)
            measureline=plot(0,0)
            for i=1:length(xs)
                subplot(numrows,1,i+1);
                pred(i)=plot(0,0);
                reality(i)=plot(0,0);
            end
        end
        

        %determine maximum time that needs to be known backwards
        maxtimehistory=max(max(timebacktrack))

        %get maximum amount of fft samples for this time
        storagelength=ceil(maxtimehistory*fs*paddinginteger*100)%ceil(maxtimehistory*fs)+1

        %variable to count how many fft samples are made of each phase
        fftcounter=ones(1,storagelength)

        %variable to store time of sample
        timestorage=zeros(1,storagelength);

        %variable to store fft samples
        fftstorage=zeros(length(omegas),storagelength);

        %column to write sample to
        pagenow=0

        %algoritm uses sophsticated memory rotation where it overrides values in
        %fftstroage when they are no longer needed. This preconditions this memory
        %'rotation'.
        focals=storagelength-length(trimtime)+1:storagelength

        TimeReversal=zeros(length(omegas),length(trimtime));%%%
        for i=1:length(omegas)
            TimeReversal(i,:)=exp(1i*(omegas(i)*trimtime));
        end

        
        
        predictprobe=zeros(length(probe),length(xs));
        predictt=zeros(length(probe),length(xs));
        predictindex=zeros(1,length(xs));

        %get start point (because fft looks back in time)
        startingindex=twidth+1

        startplottingindex=startingindex+ceil(twidth/2)+1;

        %buffer probes and time with zeros to make buildup easier
        backwardsbufferidex=floor(samplemultipler*Tp*fs);

        [a,b]=size(probe)

        probe=[zeros(backwardsbufferidex,b);probe];

        backtimething=-backwardsbufferidex:0;
        backtimething=backtimething(1:end-1)/fs;
        t=[backtimething,t];
%         t=t-t(1)



        
        
        %set controls to video
        tovideo=0;
        if tovideo==1
            videoname='video'
            v=VideoWriter(videoname,'MPEG-4')
            fps=20
            v.FrameRate=fps
            open(v)
        end

        toplot=0;
        plotskips=5;
        plotincrements=0;
        runstatus=0;

        if toplot==0
             close all
        else
                    animationfigure=figure(2)
            % speedfigure=figure(3)
            set(animationfigure,'units','normalized','outerposition',[1 0 1 1])
%             set(animationfigure,'units','normalized','outerposition',[0.1 0.1 0.7 0.7])
            drawnow;
        end

        
        tic
        for i=startingindex:length(t)-twidth-1
            %FFT Process

        %     if t(i)>92.77
        %         keyboard
        %     end

            %move to next analysis point
            pagenow=pagenow+1;
            if pagenow>storagelength
                pagenow=1;
            end

            %work out which part of storages need to be updated
            focals=focals+1;
            for Q=1:length(focals)
                if focals(Q)>storagelength
                    focals(Q)=1;
                end
            end

            %reset current time point
            fftstorage(:,pagenow)=fftstorage(:,pagenow)*0;

            %new observation time interval
            trimtimeindex=(i-twidth:i);
            trimtime=t(trimtimeindex);

            %data observed at this time step
            observeddata=(probe(trimtimeindex,1));
        %     subplot(4,1,1)
        %     hold on
        %     if exist('Q')
        %         delete(Q)
        %     end
        %     Q=plot(t(trimtimeindex),observeddata,'g');
        %     hold off
        %     pause(0.001)


            %fft your observation      
            [fftcomplex] = trimmed_fft_zeropadding(observeddata,N,N_2,tstartindex,tendindex,omegatrimstart,omegatrimend,paddinginteger);

            fftcounter(pagenow)=0;
            if avg==1
            fftcounter(focals)=fftcounter(focals)+1;

            fftstorage(:,focals)=fftstorage(:,focals)+TimeReversal.*fftcomplex;
            end
            
            if avg==0
            fftstorage(:,focals)=TimeReversal.*fftcomplex;
            fftcounter(focals)=1;%fftcounter(focals)+1;
            end
            
            
            
%             if toplot==1
%             t(i)
%             i-startplottingindex
            if runstatus>500
                'avg:'
                avg
                dummy1=length(t)-i;
                dummy2=(i-startingindex)/(length(t)-startingindex);
                dummy3=t(i)-t(startingindex);
                runstatus=0;
                dummy4=toc;
                paddinginteger
                fprintf(' Countdown: %i\n Portion Completed: %f\n Time Simulated: %f\n Time Run: %f\n\n',dummy1,dummy2,dummy3,dummy4)
                testload
                samplemultipler
            end
            runstatus=runstatus+1;

            %start predicting and plotting once buffer has been built
            if i>startplottingindex

                %Prediction Process
                for focalprobe=1:length(xs)
                    predictindex(focalprobe)=predictindex(focalprobe)+1;

                    [a,b] = predictpoint_firstorder_v2(wavenums,depth,omegas,forwardtime,focalprobe,pagenow,timebacktrackindex,fftstorage,fftcounter,t(i));

                    predictt(predictindex(focalprobe),focalprobe)=a;
                    predictprobe(predictindex(focalprobe),focalprobe)=b;
                end


                %Visual Process
                plotincrements=plotincrements+1;

                if plotincrements>plotskips && toplot==1
                    set(0, 'currentfigure', animationfigure)
                    viewbacktpmultiples=10;
                    viewbackindicies=find(t<t(i)-Tp*viewbacktpmultiples,1,'last');
                    if isempty(viewbackindicies)
                        viewbackindicies=1;
                    end
                    xlimits=[t(i)-viewbacktpmultiples*Tp,t(i)+max(max(max(forwardtime)))];
                    ylimits=[-Hs*0.9,Hs*0.9];
                    ylabel('\eta (m)')

                    subplot(numrows,1,1)
                    measureline=plot(t(viewbackindicies:i),probe(viewbackindicies:i,1),'k');
                    xlim(xlimits)
                    ylim(ylimits)

                    for focalprobe=1:length(xs)
                        plotanimation(predictprobe,predictt,predictindex,focalprobe,probe,t,i,xlimits,ylimits,pred,reality,numrows,viewbackindicies,dist,Hs);
                        pause(0.001)

                    end
                    plotincrements=0;

                    if tovideo==1
                        frame=getframe(animationfigure);
                        writeVideo(v,frame);
                    end

                end
            end


        end

        % if tovideo==1
        %     close(v)
        % end

clearvars fftcounter
clearvars fftstorage
clearvars fftcomplex
clearvars bin_vals
clearvars bottomindex
clearvars bottomportion
clearvars c
clearvars cutout
clearvars dummy
clearvars fax_rads
clearvars focals
clearvars observeddata
clearvars Omega
clearvars S
clearvars timebacktrack
clearvars timebacktrackindex
clearvars TimeReversal
clearvars timestorage
clearvars topindex
clearvars topportion
clearvars trimtime
clearvars trimtimeindex
clearvars wavenum
        
        runtime=toc
        savename=['./AVG_',num2str(avg),'_samplewidth_',num2str(samplemultiples),'_pad_',num2str(paddingintegerloop), ...
            '_divider_',num2str(divider), ...
            '_minO_',num2str(minomegamultiple),'_maxO_',num2str(maxomegamultiple), ...
            '_',num2str(testload),'.mat']
        close all
        save(savename)

end
end
end
end
end
end





