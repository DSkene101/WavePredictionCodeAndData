function [fftcomplex] = trimmed_fft_zeropadding(signals,N,N_2,trimstart,trimend,omegaminindex,omegamaxindex,paddinginteger)

    
    signals=signals(trimstart:trimend,:);
    
    [rowsignals,numsignals]=size(signals);
    rowsignals=rowsignals*paddinginteger;
    fftsignals=zeros(rowsignals,numsignals);
    for i=1:length(numsignals)
        fftsignals(:,i) = fft(signals(:,i),length(signals(:,i))*paddinginteger);
    end
    fftsignals=fftsignals(1:N_2,:)/N*2;
    fftcomplex=fftsignals(omegaminindex:omegamaxindex,:);
    
end

