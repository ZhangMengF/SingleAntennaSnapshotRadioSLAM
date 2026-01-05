function sigma2 = CalcuFreqDomNoiseSigma2(LoFs_UEs,PathGains_UEs,SubcarrierNum,SubcarrierSpacing,SNR_dB)
    
    % Function: Calculate the variance of the Gaussian noise in the frequency-domain channel based on "SNR_dB"

    LightSpeed = 3e8;
    f_vec = (0:SubcarrierNum-1).'*SubcarrierSpacing;
    UENum = length(LoFs_UEs);
    APNum = length(LoFs_UEs{1});
    FreqDomainChannel = zeros(UENum,APNum,SubcarrierNum);
    FirstPathLoAs = zeros(UENum,APNum);
    for UEidx = 1:UENum
        for APidx = 1:APNum
            FirstPathLoA = min(LoFs_UEs{UEidx}{APidx});
            FirstPathLoAs(UEidx,APidx) = FirstPathLoA;
            Delays_currentUE_AP = (LoFs_UEs{UEidx}{APidx}-FirstPathLoA)/LightSpeed;
            PathGains_currentUE_AP = PathGains_UEs{UEidx}{APidx};
            PathNum = length(Delays_currentUE_AP);
            FreqResp = zeros(SubcarrierNum,1);
            for Pathidx = 1:PathNum               
                FreqResp = FreqResp+PathGains_currentUE_AP(Pathidx)*exp(-sqrt(-1)*2*pi*f_vec*Delays_currentUE_AP(Pathidx));
            end    
            FreqDomainChannel(UEidx,APidx,:) = FreqResp;
        end
    end
    sigma2 = mean(abs(FreqDomainChannel).^2,'all')/10^(SNR_dB/10);
    
end