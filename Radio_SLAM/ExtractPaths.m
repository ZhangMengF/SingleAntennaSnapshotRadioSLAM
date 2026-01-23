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
    DelaysResolution = TsOfOFDM/OverSamplingFactor;
    DistRes_m = DelaysResolution * LightSpeed; 
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
        % The pilot transmission time of the UE is uniformly and randomly distributed within one delay resolution bin.
        TransmitTimeOffSet = rand*DelaysResolution;
        for APidx = 1:APNum
            FirstPathLoA = min(LoFs_UEs{UEidx_inDataset}{APidx})+TransmitTimeOffSet;
            FirstPathLoAs(UEidx,APidx) = FirstPathLoA;
            % Adjust the delays used for frequency-domain channel calculation to account for the first arriving path falling off-grid (within a delay grid bin)
            Delays_currentUE_AP = (LoFs_UEs{UEidx_inDataset}{APidx}-FirstPathLoA)/LightSpeed+mod(FirstPathLoA/LightSpeed,DelaysResolution);            
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
            % Calculate the absolute time-of-arrival (ToA) for each path based on the multipath delays extracted by OMP.
            ExtractedLoFs{APidx} = FirstPathLoAs(UEidx,APidx)-mod(FirstPathLoAs(UEidx,APidx),DistRes_m)+Est_Taus*LightSpeed;
            ExtractedPathPowers{APidx} = abs(Est_Amps).^2;

        end
        ExtractedLoFs_UEs{UEidx} = ExtractedLoFs;
        ExtractedPathPowers_UEs{UEidx} = ExtractedPathPowers;
                                   
    end
    
    % Randomly select 6 (UE, AP) pairs to visualize Ground Truth vs Extracted Paths
    % figure('Name', 'Random 6 CIR Check: GT vs OMP', 'Color', 'w', 'Position', [100, 100, 1200, 600]);
    % 
    % NumPlots = 6;
    % for i = 1:NumPlots
    %     % 1. Randomly select a UE from the current batch and an AP
    %     rand_UE_local_idx = randi(UENum);
    %     rand_AP_idx = randi(APNum);
    % 
    %     % Map local UE index back to dataset index to find Ground Truth
    %     rand_UE_dataset_idx = SampledUEidxs_inDataset(rand_UE_local_idx);
    % 
    %     % 2. Get Ground Truth Data
    %     GT_Distances = LoFs_UEs{rand_UE_dataset_idx}{rand_AP_idx};
    %     GT_Powers_dB = 10*log10(abs(PathGains_UEs{rand_UE_dataset_idx}{rand_AP_idx}).^2);
    % 
    %     % 3. Get Extracted Data (OMP results)
    %     Est_Distances = ExtractedLoFs_UEs{rand_UE_local_idx}{rand_AP_idx};
    %     Est_Powers_dB = 10*log10(ExtractedPathPowers_UEs{rand_UE_local_idx}{rand_AP_idx});
    % 
    %     % 4. Plotting
    %     subplot(2, 3, i); % 2 rows, 3 cols layout
    %     hold on; box on; grid on;
    % 
    %     % Plot GT (Red solid stems)
    %     if ~isempty(GT_Distances)
    %         stem(GT_Distances, GT_Powers_dB, 'r', 'LineWidth', 1.5, ...
    %             'MarkerFaceColor', 'r', 'DisplayName', 'GT');
    %     end
    % 
    %     % Plot Extracted (Blue dashed stems)
    %     if ~isempty(Est_Distances)
    %         stem(Est_Distances, Est_Powers_dB, 'b--', 'LineWidth', 1.2, ...
    %             'Marker', 'x', 'DisplayName', 'OMP');
    %     end
    % 
    %     % Labels and Title
    %     xlabel('Distance (m)');
    %     ylabel('Power (dB)');
    %     title(sprintf('UE (DatasetID %d) @ AP %d', rand_UE_dataset_idx, rand_AP_idx));
    % 
    %     if i == 1
    %         legend('Location', 'best');
    %     end
    % end
    % sgtitle('Multipath Extraction Check: Ground Truth vs OMP');

end