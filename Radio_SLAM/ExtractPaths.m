function [ExtractedLoFs_UEs, ExtractedPathPowers_UEs] ...
            = ExtractPaths(LoFs_UEs,PathGains_UEs,SampledUEidxs_inDataset,SubcarrierNum,SubcarrierSpacing,fdnoise_sigma2,MaxDelaySpread,OverSamplingFactor,delta_OMP,MaxIterNum_OMP)     
    % Input: 
        % LoFs_UEs - Raytraced dataset of multipath delays(converted to distance)
        % PathGains_UEs - Raytraced dataset of multipath gains
        % SampledUEidxs_inDataset - Indexs of selected UE samples
        % SubcarrierNum - Number of subcarrier
        % SubcarrierSpacing - Subcarrier spacing
        % fdnoise_sigma2 - The variance of Gaussian noise in the frequency domain channel
        % MaxDelaySpread - Maximum value for the set delay spread
        % OverSamplingFactor - OverSampling factor of OMP sensing matrix
        % delta_OMP - Threshold of OMP residual energy ratio 
        % MaxIterNum_OMP - Maximum number of iterations for OMP
    % Output: 
        % ExtractedLoFs_UEs - Extracted multipath delays (converted to distance)
        % ExtractedPathPowers_UEs - Extracted multipath powers 
    % Function: Select the UEs with the specified indexs from the Raytraced CIR dataset (these UEs' pilots share one OFDM symbol and occupy uniformly orthogonal subcarriers) and extract the multipath delay and power from their frequency-domain channels.

    LightSpeed = 3e8;
    TsOfOFDM = 1/(SubcarrierSpacing*SubcarrierNum);
    UENum = length(SampledUEidxs_inDataset);
    APNum = length(LoFs_UEs{1});

    f_vec_UEs = zeros(SubcarrierNum/UENum,UENum);
    for UEidx = 1:UENum
        f_vec_UEs(:,UEidx) = (UEidx-1:UENum:SubcarrierNum-1).'*SubcarrierSpacing;
    end
    
    % Generate frequency-domain channels using Raytraced multipath
    FreqDomainChannel = zeros(UENum,APNum,SubcarrierNum/UENum);
    FirstPathLoAs = zeros(UENum,APNum);
    for UEidx = 1:UENum
        UEidx_inDataset = SampledUEidxs_inDataset(UEidx);
        for APidx = 1:APNum
            FirstPathLoA = min(LoFs_UEs{UEidx_inDataset}{APidx});
            FirstPathLoAs(UEidx,APidx) = FirstPathLoA;
            Delays_currentUE_AP = (LoFs_UEs{UEidx_inDataset}{APidx}-FirstPathLoA)/LightSpeed;
            PathGains_currentUE_AP = PathGains_UEs{UEidx_inDataset}{APidx};
            PathNum = length(Delays_currentUE_AP);
            FreqResp = zeros(SubcarrierNum/UENum,1);
            for Pathidx = 1:PathNum               
                FreqResp = FreqResp+PathGains_currentUE_AP(Pathidx)*exp(-1i*2*pi*f_vec_UEs(:,UEidx)*Delays_currentUE_AP(Pathidx));
            end    
            FreqDomainChannel(UEidx,APidx,:) = FreqResp;
        end
    end

    % Add noise to the frequency-domain channels    
    FreqDomainChannel_noised = FreqDomainChannel + sqrt(fdnoise_sigma2/2)*randn(UENum,APNum,SubcarrierNum/UENum) ...
                                +sqrt(-1)*sqrt(fdnoise_sigma2/2)*randn(UENum,APNum,SubcarrierNum/UENum);

    % Extract CIR using OMP
    DelaysResolution = TsOfOFDM/OverSamplingFactor;
    ColumnNum = floor(MaxDelaySpread/DelaysResolution);
    Taus_tick = (0:ColumnNum-1).'*DelaysResolution;
    MeasureMatrix = zeros(SubcarrierNum/UENum,ColumnNum);

    ExtractedLoFs_UEs = cell(UENum,1);
    ExtractedPathPowers_UEs = cell(UENum,1);
    for UEidx = 1:UENum
        ExtractedLoFs = cell(APNum,1);
        ExtractedPathPowers = cell(APNum,1);
        for APidx = 1:APNum
            % Construct sensing matrix
            j2pi_f = sqrt(-1)*2*pi*f_vec_UEs(:,UEidx);
            for cidx = 1:ColumnNum
                MeasureMatrix(:,cidx) = exp(-j2pi_f*Taus_tick(cidx));
            end
            [Taus_grid,PathGains_OMP] = OMP(MeasureMatrix,squeeze(FreqDomainChannel_noised(UEidx,APidx,:)),delta_OMP,MaxIterNum_OMP);
            Est_Taus = Taus_tick(Taus_grid);
            Est_Amps = PathGains_OMP;      
            [~,SortedIdxs] = sort(abs(Est_Amps),'descend');
            Est_Taus = Est_Taus(SortedIdxs);
            Est_Amps = Est_Amps(SortedIdxs);
            ExtractedLoFs{APidx} = FirstPathLoAs(UEidx,APidx)+Est_Taus*LightSpeed;
            ExtractedPathPowers{APidx} = abs(Est_Amps).^2;
        end
        ExtractedLoFs_UEs{UEidx} = ExtractedLoFs;
        ExtractedPathPowers_UEs{UEidx} = ExtractedPathPowers;
    end
    
end