clc
clear

% Load parameters
CommonParas

% Independent variable---
Paras = [120e3,150e3,180e3,210e3,240e3]; % Subcarrier spacing
Para_length = length(Paras);

% Num of Monte Carlo trials
MCNum = 5;

% Dependent variables---
% Perfect association
PerfAssociLocaRate = zeros(Para_length,1);
PerfAssociLocaErr = zeros(Para_length,1);
PerfAssociMapRate = zeros(Para_length,1);
PerfAssociMapErr = zeros(Para_length,1);

OptimizedPerfAssociLocaRate = zeros(Para_length,1);
OptimizedPerfAssociLocaErr = zeros(Para_length,1);
OptimizedPerfAssociMapRate = zeros(Para_length,1);
OptimizedPerfAssociMapErr = zeros(Para_length,1);

% Proposed SLAM
LocaRate_SNR1 = zeros(Para_length,1);
LocaErr_SNR1 = zeros(Para_length,1);
MapRate_SNR1 = zeros(Para_length,1);
MapErr_SNR1 = zeros(Para_length,1);

LocaRate_SNR2 = zeros(Para_length,1);
LocaErr_SNR2 = zeros(Para_length,1);
MapRate_SNR2 = zeros(Para_length,1);
MapErr_SNR2 = zeros(Para_length,1);

% Load Raytraced dataset---
[SpaceBound,pAs,UEPositions_DataSet,LoAs_UEs_DataSet,PathGains_UEs_DataSet,TrueTargetPathIndexs_UEs_DataSet] ...
    = LoadRaytracedPaths();
APNum = size(pAs,2);

% Calculate the variance of the Gaussian noise in the frequency-domain channel based on SNR
fdnoise1_sigma2 = CalcuFreqDomNoiseSigma2(LoAs_UEs_DataSet,PathGains_UEs_DataSet,SubcarrierNum,SubcarrierSpacing,SNR1_dB);
fdnoise2_sigma2 = 10^((SNR1_dB-SNR2_dB)/10)*fdnoise1_sigma2;

SampledUEidxs_MC_inDataset = zeros(UENum,MCNum);

% Intermediate dependent variables---
% Perfect association
PerfAssociLocatedUEs_MC_Paras = zeros(3,UENum,MCNum,Para_length);
PerfAssociLocateSuccessFlags = false(UENum,MCNum,Para_length);
PerfAssociMap_MC_Paras = zeros(4,MCNum,Para_length);

OptimizedPerfAssociLocatedUEs_MC_Paras = zeros(3,UENum,MCNum,Para_length);
OptimizedPerfAssociLocateSuccessFlags = false(UENum,MCNum,Para_length);
OptimizedPerfAssociMap_MC_Paras = zeros(4,MCNum,Para_length);

% Proposed SLAM
LocatedUEs_MC_Paras_SNR1 = zeros(3,UENum,MCNum,Para_length);
LocateSuccessFlags_SNR1 = false(UENum,MCNum,Para_length);
Map_MC_Paras_SNR1 = zeros(4,MCNum,Para_length);

LocatedUEs_MC_Paras_SNR2 = zeros(3,UENum,MCNum,Para_length);
LocateSuccessFlags_SNR2 = false(UENum,MCNum,Para_length);
Map_MC_Paras_SNR2 = zeros(4,MCNum,Para_length);

tic
fprintf('Num of Monte Carlo trials%i, length of independent variable%i \n', MCNum, Para_length)
parfor (MCidx = 1:MCNum,6) % Randomly sample UEs
% The use of "parfor" needs "Parallel Computing Toolbox"
% for MCidx = 1:MCNum
    % Randomly select UENum positions from all the raytraced UE positions.
    SampledUEidxs_inDataset = randperm(length(LoAs_UEs_DataSet),UENum);    
    SampledUEidxs_MC_inDataset(:,MCidx) = SampledUEidxs_inDataset;

    for Paraidx = 1:Para_length

        SubcarrierSpacing = Paras(Paraidx);
        TsOfOFDM_L = LightSpeed/(SubcarrierSpacing*SubcarrierNum);        

        % Extract paths
        [ExtractedLoAs_SNR1_UEs, ExtractedPathPowers_SNR1_UEs] ...
            = ExtractPaths(LoAs_UEs_DataSet,PathGains_UEs_DataSet,SampledUEidxs_inDataset,SubcarrierNum,SubcarrierSpacing,fdnoise1_sigma2,MaxDelaySpread,OverSamplingFactor,delta_OMP,MaxIterNum_OMP);    

        [ExtractedLoAs_SNR2_UEs, ExtractedPathPowers_SNR2_UEs] ...
            = ExtractPaths(LoAs_UEs_DataSet,PathGains_UEs_DataSet,SampledUEidxs_inDataset,SubcarrierNum,SubcarrierSpacing,fdnoise2_sigma2,MaxDelaySpread,OverSamplingFactor,delta_OMP,MaxIterNum_OMP);

        % Prepare multipath delay observations for perfect association
        LoAs4PerfAssoci_UEs = LoAs_UEs_DataSet(SampledUEidxs_inDataset);
        for i = 1:UENum
            for j = 1:APNum
                CurrentLoAs = LoAs4PerfAssoci_UEs{i}{j};
                % Quantify delays
                LoAModRemains = mod(CurrentLoAs,TsOfOFDM_L);
                LoAs4PerfAssoci_UEs{i}{j} = CurrentLoAs+(LoAModRemains>TsOfOFDM_L/2).*(TsOfOFDM_L-LoAModRemains)+(LoAModRemains<=TsOfOFDM_L/2).*(-LoAModRemains);   
            end
        end        

        PerfAssociSolvedWallxy_UEs = nan(4,UENum);
        OptimizedPerfAssociSolvedWallxy_UEs = nan(4,UENum);
        SolvedWallxy_SNR1_UEs = nan(4,UENum);
        SolvedWallxy_SNR2_UEs = nan(4,UENum);

        for UEidx = 1:UENum

            % Perfect association
            LoAs4PerfAssoci = LoAs4PerfAssoci_UEs{UEidx};
            TrueTargetPathIndexs = TrueTargetPathIndexs_UEs_DataSet{SampledUEidxs_inDataset(UEidx)};
            [PerfAssociLocatedUE,PerfAssociLocatedUE_Optimized,PerfAssociWallxy,PerfAssociWallxy_Optimized,~,~] ...
                = PerfectAssociationSLAM(LoAs4PerfAssoci,TrueTargetPathIndexs,MinANum,pAs,true);
            
            if ~isempty(PerfAssociLocatedUE)
                PerfAssociLocatedUEs_MC_Paras(:,UEidx,MCidx,Paraidx) = PerfAssociLocatedUE;
                PerfAssociLocateSuccessFlags(UEidx,MCidx,Paraidx) = true;
            end
            PerfAssociSolvedWallxy_UEs(:,UEidx) = PerfAssociWallxy;

            if ~isempty(PerfAssociLocatedUE_Optimized)
                OptimizedPerfAssociLocatedUEs_MC_Paras(:,UEidx,MCidx,Paraidx) = PerfAssociLocatedUE_Optimized;
                OptimizedPerfAssociLocateSuccessFlags(UEidx,MCidx,Paraidx) = true;
            end
            OptimizedPerfAssociSolvedWallxy_UEs(:,UEidx) = PerfAssociWallxy_Optimized;

            % Proposed SLAM
            % SNR1
            LoAObservations_SNR1 = ExtractedLoAs_SNR1_UEs{UEidx};
            PathPowerObservations_SNR1 = ExtractedPathPowers_SNR1_UEs{UEidx};
            [LocatedUE_SNR1,~,~,~,Wallxy_SNR1] = ...
                SLAM(LoAObservations_SNR1,PathPowerObservations_SNR1,pAs,MinANum,PowerCheckThres,RandNumPerSeed,SolveWallTols);            
            if ~isempty(LocatedUE_SNR1)
                LocatedUEs_MC_Paras_SNR1(:,UEidx,MCidx,Paraidx) = LocatedUE_SNR1;
                LocateSuccessFlags_SNR1(UEidx,MCidx,Paraidx) = true;
            end                       
            SolvedWallxy_SNR1_UEs(:,UEidx) = Wallxy_SNR1;
            
            % SNR2
            LoAObservations_SNR2 = ExtractedLoAs_SNR2_UEs{UEidx};
            PathPowerObservations_SNR2 = ExtractedPathPowers_SNR2_UEs{UEidx};
            [LocatedUE_SNR2,~,~,~,Wallxy_SNR2] = ...
                SLAM(LoAObservations_SNR2,PathPowerObservations_SNR2,pAs,MinANum,PowerCheckThres,RandNumPerSeed,SolveWallTols);
            if ~isempty(LocatedUE_SNR2)
                LocatedUEs_MC_Paras_SNR2(:,UEidx,MCidx,Paraidx) = LocatedUE_SNR2;
                LocateSuccessFlags_SNR2(UEidx,MCidx,Paraidx) = true;
            end                       
            SolvedWallxy_SNR2_UEs(:,UEidx) = Wallxy_SNR2;

        end

        PerfAssociMap_MC_Paras(:,MCidx,Paraidx) = mean(PerfAssociSolvedWallxy_UEs,2,"omitnan");
        OptimizedPerfAssociMap_MC_Paras(:,MCidx,Paraidx) = mean(OptimizedPerfAssociSolvedWallxy_UEs,2,"omitnan");
        Map_MC_Paras_SNR1(:,MCidx,Paraidx) = mean(SolvedWallxy_SNR1_UEs,2,"omitnan");
        Map_MC_Paras_SNR2(:,MCidx,Paraidx) = mean(SolvedWallxy_SNR2_UEs,2,"omitnan");

        fprintf('Monte Carlo try%i-parameter%i finished~ take %.1f s\n', MCidx, Paraidx,toc)
    end

end

UEPositions_GroundTruth = zeros(3,UENum,MCNum);
for MCidx = 1:MCNum
    UEPositions_GroundTruth(:,:,MCidx) = UEPositions_DataSet(:,SampledUEidxs_MC_inDataset(:,MCidx));
end

MapGroudTruth = [SpaceBound(1)/2;SpaceBound(2)/2;-SpaceBound(1)/2;-SpaceBound(2)/2];

% Calculate dependent variables according to Intermediate dependent variables
for Paraidx = 1:Para_length
    
    % Localization errors------------------------------------
    derrs_SNR1 = zeros(UENum,MCNum);
    derrs_SNR2 = zeros(UENum,MCNum);
    PerfAssoci_derrs = zeros(UENum,MCNum);
    OptimizedPerfAssoci_derrs = zeros(UENum,MCNum);    
    for UEidx = 1:UENum 
        for MCidx = 1:MCNum
            CurrentGroundTruth = UEPositions_GroundTruth(:,UEidx,MCidx);
         
            if LocateSuccessFlags_SNR1(UEidx,MCidx,Paraidx)
                derrs_SNR1(UEidx,MCidx) = norm(LocatedUEs_MC_Paras_SNR1(:,UEidx,MCidx,Paraidx)-CurrentGroundTruth);
            end
            
            if LocateSuccessFlags_SNR2(UEidx,MCidx,Paraidx)
                derrs_SNR2(UEidx,MCidx) = norm(LocatedUEs_MC_Paras_SNR2(:,UEidx,MCidx,Paraidx)-CurrentGroundTruth);
            end

            if PerfAssociLocateSuccessFlags(UEidx,MCidx,Paraidx)
                PerfAssoci_derrs(UEidx,MCidx) = norm(PerfAssociLocatedUEs_MC_Paras(:,UEidx,MCidx,Paraidx)-CurrentGroundTruth);
            end

            if OptimizedPerfAssociLocateSuccessFlags(UEidx,MCidx,Paraidx)
                OptimizedPerfAssoci_derrs(UEidx,MCidx) = norm(OptimizedPerfAssociLocatedUEs_MC_Paras(:,UEidx,MCidx,Paraidx)-CurrentGroundTruth);
            end
        end
    end
    
    % Perfect association    
    PerfAssociLocaRate(Paraidx) = nnz(PerfAssoci_derrs)/(UENum*MCNum);
    PerfAssociLocaErr(Paraidx) = sum(PerfAssoci_derrs(:))/nnz(PerfAssoci_derrs);
    OptimizedPerfAssociLocaRate(Paraidx) = nnz(OptimizedPerfAssoci_derrs)/(UENum*MCNum);
    OptimizedPerfAssociLocaErr(Paraidx) = sum(OptimizedPerfAssoci_derrs(:))/nnz(OptimizedPerfAssoci_derrs);
    % Proposed SLAM
    LocaRate_SNR1(Paraidx) = nnz(derrs_SNR1)/(UENum*MCNum);
    LocaErr_SNR1(Paraidx) = sum(derrs_SNR1(:))/nnz(derrs_SNR1);
    LocaRate_SNR2(Paraidx) = nnz(derrs_SNR2)/(UENum*MCNum);
    LocaErr_SNR2(Paraidx) = sum(derrs_SNR2(:))/nnz(derrs_SNR2);

    % Mapping errors------------------------------------
    % Perfect association
    PerfAssociMapRate(Paraidx) = nnz(~isnan(PerfAssociMap_MC_Paras(:,:,Paraidx)))/(4*MCNum);
    PerfAssociMapErr(Paraidx) = mean(abs(PerfAssociMap_MC_Paras(:,:,Paraidx)-MapGroudTruth),"all",'omitnan');
    OptimizedPerfAssociMapRate(Paraidx) = nnz(~isnan(OptimizedPerfAssociMap_MC_Paras(:,:,Paraidx)))/(4*MCNum);
    OptimizedPerfAssociMapErr(Paraidx) = mean(abs(OptimizedPerfAssociMap_MC_Paras(:,:,Paraidx)-MapGroudTruth),"all",'omitnan');
    % Proposed SLAM
    MapRate_SNR1(Paraidx) = nnz(~isnan(Map_MC_Paras_SNR1(:,:,Paraidx)))/(4*MCNum);
    MapErr_SNR1(Paraidx) = mean(abs(Map_MC_Paras_SNR1(:,:,Paraidx)-MapGroudTruth),"all",'omitnan');
    MapRate_SNR2(Paraidx) = nnz(~isnan(Map_MC_Paras_SNR2(:,:,Paraidx)))/(4*MCNum);
    MapErr_SNR2(Paraidx) = mean(abs(Map_MC_Paras_SNR2(:,:,Paraidx)-MapGroudTruth),"all",'omitnan');
end
fprintf('Data processing finished~\n')

%%
% Localization performance
figure
hold on
yyaxis left
ylim([0,3.5])
ylabel('Mean Distance Error (m)')
hl1 = plot(Paras*SubcarrierNum,LocaErr_SNR1,'-','DisplayName',"Proposed, SNR = "+num2str(SNR1_dB)+" dB",'Color','blue','LineWidth',1,'Marker','.','MarkerSize',10);
hl2 = plot(Paras*SubcarrierNum,LocaErr_SNR2,'-','DisplayName',"Proposed, SNR = "+num2str(SNR2_dB)+" dB",'Color','blue','LineWidth',1,'Marker','o','MarkerSize',5);
hl3 = plot(Paras*SubcarrierNum,PerfAssociLocaErr,'-','DisplayName','Perfect Association','Color','red','LineWidth',1,'Marker','+','MarkerSize',6);
hl4 = plot(Paras*SubcarrierNum,OptimizedPerfAssociLocaErr,'-','DisplayName','Perfect Association (Optimized AP Selection)','Color','black','LineWidth',1,'Marker','*','MarkerSize',6);

yyaxis right
ylim([0,1])
ylabel('Success Rate')
plot(Paras*SubcarrierNum,LocaRate_SNR1,'-','Color','blue','LineWidth',1,'Marker','.','MarkerSize',10)
plot(Paras*SubcarrierNum,LocaRate_SNR2,'-','Color','blue','LineWidth',1,'Marker','o','MarkerSize',5)
plot(Paras*SubcarrierNum,PerfAssociLocaRate,'-','Color','red','LineWidth',1,'Marker','+','MarkerSize',6)
plot(Paras*SubcarrierNum,OptimizedPerfAssociLocaRate,'-','Color','black','LineWidth',1,'Marker','*','MarkerSize',6)
xlabel('Bandwidth (Hz)')
legend([hl1,hl2,hl3,hl4])
grid on

% Mapping performance
figure
hold on
plot(Paras*SubcarrierNum,MapRate_SNR1,'DisplayName',"Proposed, SNR = "+num2str(SNR1_dB)+" dB",'Color','blue','LineWidth',1,'Marker','.','MarkerSize',10)
plot(Paras*SubcarrierNum,MapRate_SNR2,'DisplayName',"Proposed, SNR = "+num2str(SNR2_dB)+" dB",'Color','blue','LineWidth',1,'Marker','o','MarkerSize',5)
plot(Paras*SubcarrierNum,PerfAssociMapRate,'DisplayName','Perfect Association','Marker','+','MarkerSize',6)
plot(Paras*SubcarrierNum,OptimizedPerfAssociMapRate,'DisplayName','Perfect Association (Optimized AP Selection)','Marker','*','MarkerSize',6)
xlabel('Bandwidth (Hz)')
ylabel('Finishing rate')
legend
grid on
title('Mapping Performance - Rate')

figure
hold on
plot(Paras*SubcarrierNum,MapErr_SNR1,'DisplayName',"Proposed, SNR = "+num2str(SNR1_dB)+" dB",'Color','blue','LineWidth',1,'Marker','.','MarkerSize',10)
plot(Paras*SubcarrierNum,MapErr_SNR2,'DisplayName',"Proposed, SNR = "+num2str(SNR2_dB)+" dB",'Color','blue','LineWidth',1,'Marker','o','MarkerSize',5)
plot(Paras*SubcarrierNum,PerfAssociMapErr,'DisplayName','Perfect Association','Color','red','LineWidth',1,'Marker','+','MarkerSize',6)
plot(Paras*SubcarrierNum,OptimizedPerfAssociMapErr,'DisplayName','Perfect Association (Optimized AP Selection)','Color','black','LineWidth',1,'Marker','*','MarkerSize',6)
xlabel('Bandwidth (Hz)')
ylabel('Mean Absolute Error (m)')
legend
grid on