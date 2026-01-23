clc
clear

rng(2026);

CommonParas

Paras = [120e3,150e3,180e3,210e3,240e3];% Subcarrier spacing
Para_length = length(Paras);

MCNum = 100;

% The average distance error of TDoA positioning of path combinations
AverDerrOfPathCombs = zeros(4,UENum,MCNum,Para_length);
% 1-Proposed, 2-Ground-Truth association, 3-Random association, 4-Constrained association

[SpaceBound,pAs,UEPositions_DataSet,LoAs_UEs_DataSet,PathGains_UEs_DataSet,TrueTargetPathIndexs_UEs_DataSet] ...
    = LoadRaytracedPaths();
APNum = size(pAs,2);

fdnoise_sigma2 = CalcuFreqDomNoiseSigma2(LoAs_UEs_DataSet,PathGains_UEs_DataSet,SubcarrierNum,SubcarrierSpacing,SNR1_dB);

tic
fprintf('Num of Monte Carlo trials-%i, length of independent variable-%i \n', MCNum, Para_length)

for MCidx = 1:MCNum

    % Randomly select UENum positions from all the raytraced UE positions.
    SampledUEidxs_inDataset = randperm(length(LoAs_UEs_DataSet),UENum);
    
    for Paraidx = 1:Para_length

        SubcarrierSpacing = Paras(Paraidx);
        TsOfOFDM_L = LightSpeed/(SubcarrierSpacing*SubcarrierNum);

        [ExtractedLoAs_UEs, ExtractedPathPowers_UEs] ...
            = ExtractPaths(LoAs_UEs_DataSet,PathGains_UEs_DataSet,SampledUEidxs_inDataset,SubcarrierNum,SubcarrierSpacing,fdnoise_sigma2,MaxDelaySpread,OverSamplingFactor,delta_OMP,MaxIterNum_OMP);    
                     
        % Prepare multipath delay observations for Ground-Truth association
        LoAs4PerfAssoci_UEs = LoAs_UEs_DataSet(SampledUEidxs_inDataset);
        for i = 1:UENum
            for j = 1:APNum
                CurrentLoAs = LoAs4PerfAssoci_UEs{i}{j};
                % Quantify delays
                LoAModRemains = mod(CurrentLoAs,TsOfOFDM_L);
                LoAs4PerfAssoci_UEs{i}{j} = CurrentLoAs+(LoAModRemains>TsOfOFDM_L/2).*(TsOfOFDM_L-LoAModRemains)+(LoAModRemains<=TsOfOFDM_L/2).*(-LoAModRemains);
            end
        end         
    
        for UEidx = 1:UENum
                
            UEMirrorUEGroundTruthPositions = zeros(3,6);
            pU = UEPositions_DataSet(:,SampledUEidxs_inDataset(UEidx));
            UEMirrorUEGroundTruthPositions(:,1) = pU;                
            UEMirrorUEGroundTruthPositions(:,2) = [SpaceBound(1)-pU(1);pU(2);pU(3)];    
            UEMirrorUEGroundTruthPositions(:,3) = [pU(1);SpaceBound(2)-pU(2);pU(3)];    
            UEMirrorUEGroundTruthPositions(:,4) = [-SpaceBound(1)-pU(1);pU(2);pU(3)];
            UEMirrorUEGroundTruthPositions(:,5) = [pU(1);-SpaceBound(2)-pU(2);pU(3)];
            UEMirrorUEGroundTruthPositions(:,6) = [pU(1);pU(2);-pU(3)];  
            
            % Ground-Truth association
            LoAs4PerfAssoci = LoAs4PerfAssoci_UEs{UEidx};
            TrueTargetPathIndexs = TrueTargetPathIndexs_UEs_DataSet{SampledUEidxs_inDataset(UEidx)};
            [~,~,~,~,PerfAssociLocatedUEandMirror,~] ...
                = PerfectAssociationSLAM(LoAs4PerfAssoci,TrueTargetPathIndexs,MinANum,pAs,false);
            AverDerrOfPathCombs(2,UEidx,MCidx,Paraidx) = CalcuAverDerr(PerfAssociLocatedUEandMirror,UEMirrorUEGroundTruthPositions);

            % Proposed
            ExtractedLoAs = ExtractedLoAs_UEs{UEidx};
            ExtractedPathPowers = ExtractedPathPowers_UEs{UEidx};
            [~,LocatedPoints_AdjacentMatrix,HighBeliefLocatedPoints,APSelectedProbabilities,~] = ...
                SLAM(ExtractedLoAs,ExtractedPathPowers,pAs,MinANum,PowerCheckThres,RandNumPerSeed,SolveWallTols);
            AverDerrOfPathCombs(1,UEidx,MCidx,Paraidx) = CalcuAverDerr(HighBeliefLocatedPoints,UEMirrorUEGroundTruthPositions);          
            
            % Constrained association
            AverDerrOfPathCombs(4,UEidx,MCidx,Paraidx) = CalcuAverDerr(LocatedPoints_AdjacentMatrix,UEMirrorUEGroundTruthPositions);            

            % Random association
            ExtractedLoAs = cellfun(@(x) x(1:6), ExtractedLoAs, 'UniformOutput', false);
            LocatedPoints_RandomSearch = LocateRandomSearchPathCombs(pAs,MinANum,ExtractedLoAs,APSelectedProbabilities,100);
            AverDerrOfPathCombs(3,UEidx,MCidx,Paraidx) = CalcuAverDerr(LocatedPoints_RandomSearch,UEMirrorUEGroundTruthPositions);
    
        end
        fprintf('Monte Carlo try%i-parameter%i finished~ take %.1f s\n', MCidx, Paraidx,toc)    

    end

end

function AverDerr = CalcuAverDerr(EstimatedPoints,GroundTruthPoints)
    
    PointNum = size(EstimatedPoints,2);
    MatchedDerrs = zeros(PointNum,1);
    for i = 1:PointNum
        MatchedDerrs(i) = min(vecnorm(EstimatedPoints(:,i)-GroundTruthPoints));
    end
    AverDerr = mean(MatchedDerrs);

end

function ValidPoints = LocateRandomSearchPathCombs(pAs,MinANum,LoAObservations,APSelectedProbabilities,TotalNum)    
    
    APNum = size(pAs,2);
    ValidPointNum = 0;
    ValidPoints = zeros(3,TotalNum);
    while ValidPointNum < TotalNum 
        % Generate AP selection
        if isempty(APSelectedProbabilities)
            ValidPoints = [];
            return
        else
            APSelectedFlags = rand(APNum,1)<APSelectedProbabilities;
        end
        if nnz(APSelectedFlags)>=MinANum
            % Generate path combinations
            LoAObservations_selected = LoAObservations(APSelectedFlags);
            LComb = cellfun(@(x) x(randi([1,6])), LoAObservations_selected, 'UniformOutput', true);
            [p_esti,~] = LDoALoca(pAs(:,APSelectedFlags),LComb,[]);
            if abs(p_esti(1))>0
                ValidPointNum = ValidPointNum+1;
                ValidPoints(:,ValidPointNum) = p_esti;
            end
        end
    end

end

%%
MeanDerrs = squeeze(mean(AverDerrOfPathCombs,[2 3],"omitnan")).';
figure
hold on
plot(Paras*SubcarrierNum,MeanDerrs(:,3),'DisplayName','Random Association','Color','m','LineWidth',1,'Marker','^','MarkerSize',6)
plot(Paras*SubcarrierNum,MeanDerrs(:,4),'DisplayName','Constrained Association','Color','c','LineWidth',1,'Marker','square','MarkerSize',6)
plot(Paras*SubcarrierNum,MeanDerrs(:,1),'DisplayName','Proposed','Color','blue','LineWidth',1,'Marker','.','MarkerSize',10)
plot(Paras*SubcarrierNum,MeanDerrs(:,2),'DisplayName','Ground-Truth Association','Color','red','LineWidth',1,'Marker','+','MarkerSize',6)
legend
grid on
xlabel('Bandwidth (Hz)')
ylabel('Mean Distance Error (m)')