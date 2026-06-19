%MAIN_PERFORMANCE_DA Evaluate path-association quality versus confidence retention.
% The experiment compares proposed consistency-based association with
% ground-truth and random association using mean XY-plane source-position error.
clc
clear

rng(0);

ScriptDir = fileparts(mfilename('fullpath'));
addpath(ScriptDir)

FixedSubcarrierSpacing = 120e3;
% Fractions of height-valid candidates retained after confidence ranking.
GammaGrid = [1,0.75,0.5,0.25,0.05];
GammaPercentGrid = 100*GammaGrid;
Para_length = length(GammaGrid);

[BaseParas,funcs] = Base('SubcarrierSpacing',FixedSubcarrierSpacing);
UENum = BaseParas.UENum;

MCNum = 100;

% Rows: 1) Proposed consistency-based association, 2) ground-truth association,
% 3) random association.
AverDerrOfPathCombs = nan(3,UENum,MCNum,Para_length);

[SpaceBound,pAs,UEPositions_DataSet,LoAs_UEs_DataSet,PathGains_UEs_DataSet,TrueTargetPathIndexs_UEs_DataSet] = ...
    funcs.LoadRaytracedPaths();
APNum = size(pAs,2);

% Fix the frequency-domain noise variance for all Monte Carlo trials.
fdnoise_sigma2 = funcs.CalcuFreqDomNoiseSigma2( ...
    LoAs_UEs_DataSet, ...
    PathGains_UEs_DataSet, ...
    BaseParas.SubcarrierNum, ...
    BaseParas.SubcarrierSpacing, ...
    BaseParas.SNR1_dB);

fprintf('Scheme 2 DA performance: MC=%i, SCS=%.0f kHz, gamma count=%i \n', ...
    MCNum, FixedSubcarrierSpacing/1e3, Para_length)

parfor MCidx = 1:MCNum
% PARFOR requires Parallel Computing Toolbox.
    worker_tic = tic;

    SampledUEidxs_inDataset = randperm(length(LoAs_UEs_DataSet),UENum);
    AverDerrOfPathCombs_MC = nan(3,UENum,1,Para_length);

    CurrentParas = BaseParas;
    TsOfOFDM_L = CurrentParas.LightSpeed/(CurrentParas.SubcarrierSpacing*CurrentParas.SubcarrierNum);

    [ExtractedLoAs_UEs, ExtractedPathPowers_UEs] = funcs.ExtractPaths( ...
        LoAs_UEs_DataSet, ...
        PathGains_UEs_DataSet, ...
        SampledUEidxs_inDataset, ...
        CurrentParas.SubcarrierNum, ...
        CurrentParas.SubcarrierSpacing, ...
        fdnoise_sigma2, ...
        CurrentParas.MaxDelaySpread, ...
        CurrentParas.OverSamplingFactor, ...
        CurrentParas.delta_OMP, ...
        CurrentParas.MaxIterNum_OMP);

    % Quantize ray-traced arrival distances to the bandwidth-limited delay grid.
    LoAs4PerfAssoci_UEs = LoAs_UEs_DataSet(SampledUEidxs_inDataset);
    for i = 1:UENum
        for j = 1:APNum
            CurrentLoAs = LoAs4PerfAssoci_UEs{i}{j};
            if isempty(CurrentLoAs)
                continue
            end
            LoAModRemains = mod(CurrentLoAs,TsOfOFDM_L);
            LoAs4PerfAssoci_UEs{i}{j} = CurrentLoAs ...
                + (LoAModRemains>TsOfOFDM_L/2).*(TsOfOFDM_L-LoAModRemains) ...
                + (LoAModRemains<=TsOfOFDM_L/2).*(-LoAModRemains);
        end
    end

    for UEidx = 1:UENum
        % True source set: UE, four wall-induced mirrors, and the floor mirror.
        UEMirrorUEGroundTruthPositions = zeros(3,6);
        pU = UEPositions_DataSet(:,SampledUEidxs_inDataset(UEidx));
        UEMirrorUEGroundTruthPositions(:,1) = pU;
        UEMirrorUEGroundTruthPositions(:,2) = [SpaceBound(1)-pU(1);pU(2);pU(3)];
        UEMirrorUEGroundTruthPositions(:,3) = [pU(1);SpaceBound(2)-pU(2);pU(3)];
        UEMirrorUEGroundTruthPositions(:,4) = [-SpaceBound(1)-pU(1);pU(2);pU(3)];
        UEMirrorUEGroundTruthPositions(:,5) = [pU(1);-SpaceBound(2)-pU(2);pU(3)];
        UEMirrorUEGroundTruthPositions(:,6) = [pU(1);pU(2);-pU(3)];

        % Ground-truth association uses the known path labels.
        LoAs4PerfAssoci = LoAs4PerfAssoci_UEs{UEidx};
        TrueTargetPathIndexs = TrueTargetPathIndexs_UEs_DataSet{SampledUEidxs_inDataset(UEidx)};
        PerfAssociLocatedUEandMirror = LocateGroundTruthAssociatedPoints( ...
            LoAs4PerfAssoci, ...
            TrueTargetPathIndexs, ...
            CurrentParas.MinANum, ...
            pAs, ...
            CurrentParas.z0_UE, ...
            funcs.LDoALoca);
        GroundTruthDerr = CalcuAverDerr(PerfAssociLocatedUEandMirror,UEMirrorUEGroundTruthPositions);
        AverDerrOfPathCombs_MC(2,UEidx,1,:) = GroundTruthDerr;

        % Generate all proposed candidates once, then retain each top-gamma subset.
        ExtractedLoAs = ExtractedLoAs_UEs{UEidx};
        ExtractedPathPowers = ExtractedPathPowers_UEs{UEidx};
        [~,LocatedPoints,~,APSelectedProbabilities,~,CandidateInfo] = funcs.SLAM( ...
            ExtractedLoAs, ...
            ExtractedPathPowers, ...
            pAs, ...
            CurrentParas.MinANum, ...
            CurrentParas.PowerCheckThres, ...
            CurrentParas.RandNumPerSeed, ...
            CurrentParas.SolveWallTols, ...
            CurrentParas.sigma_t, ...
            CurrentParas.sigma_P, ...
            1, ...
            2, ...
            CurrentParas.z0_UE, ...
            false, ...
            CurrentParas.beta, ...
            "UEBased", ...
            "CandidatesOnly");

        for Paraidx = 1:Para_length
            HighConfidenceLocatedPoints = SelectHighConfidencePoints( ...
                LocatedPoints, ...
                CandidateInfo.Beliefs, ...
                CandidateInfo.ValidIndexs, ...
                GammaGrid(Paraidx));
            AverDerrOfPathCombs_MC(1,UEidx,1,Paraidx) = CalcuAverDerr( ...
                HighConfidenceLocatedPoints, ...
                UEMirrorUEGroundTruthPositions);
        end

        % Randomly associate paths while matching the proposed AP-use frequencies.
        ExtractedLoAs_Random = cellfun(@(x) x(1:min(length(x),6)), ExtractedLoAs, 'UniformOutput', false);
        LocatedPoints_RandomSearch = LocateRandomSearchPathCombs( ...
            pAs, ...
            CurrentParas.MinANum, ...
            ExtractedLoAs_Random, ...
            APSelectedProbabilities, ...
            100, ...
            funcs.LDoALoca);
        RandomDerr = CalcuAverDerr(LocatedPoints_RandomSearch,UEMirrorUEGroundTruthPositions);
        AverDerrOfPathCombs_MC(3,UEidx,1,:) = RandomDerr;
    end

    AverDerrOfPathCombs(:,:,MCidx,:) = AverDerrOfPathCombs_MC;
    fprintf('Monte Carlo try %i finished, take %.1f s\n', MCidx, toc(worker_tic))
end

MeanDerrs = squeeze(mean(AverDerrOfPathCombs,[2 3],"omitnan")).';
% Plot the average nearest-source XY error for the three association methods.
Fig = figure('Color','w');
hold on
plot(GammaPercentGrid,MeanDerrs(:,3),'DisplayName','Random Association','Color','m','LineWidth',2,'Marker','^','MarkerSize',6)
plot(GammaPercentGrid,MeanDerrs(:,1),'DisplayName',sprintf('Consistency-based Association\n(Top confidence)'),'Color','blue','LineWidth',2,'Marker','.','MarkerSize',12)
plot(GammaPercentGrid,MeanDerrs(:,2),'DisplayName','Ground-Truth Association','Color','red','LineWidth',2,'Marker','+','MarkerSize',6)
legend('Location','best')
grid on
set(gca,'XDir','reverse','XTick',sort(GammaPercentGrid))
ylim([0,6])
xlabel('Top-confidence selection ratio (%)')
ylabel('Mean XY-Plane Error (m)')
Resolution_m = BaseParas.LightSpeed/(BaseParas.SubcarrierSpacing*BaseParas.SubcarrierNum);
Bandwidth = BaseParas.SubcarrierSpacing*BaseParas.SubcarrierNum;
title({ ...
    'Data Association Performance - Scheme 2', ...
    sprintf('MC=%i, SCS=%.0f kHz, BW=%.2f MHz, sigma_P=%.1f dB, sigma_t=%.2f res', ...
        MCNum, BaseParas.SubcarrierSpacing/1e3, Bandwidth/1e6, BaseParas.sigma_P, BaseParas.sigma_t/Resolution_m)}, ...
    'Interpreter','none', 'FontSize', 18)

FigureDir = fullfile(ScriptDir,'Figures');
if ~exist(FigureDir,'dir')
    mkdir(FigureDir)
end
SavePathNoExt = fullfile(FigureDir,sprintf('Performance_DA_Scheme2_GammaScan_%s',char(datetime('now', 'Format', 'MMddHHmmss'))));
savefig(Fig,[SavePathNoExt,'.fig'])
exportgraphics(Fig,[SavePathNoExt,'.png'],'Resolution',300)

fprintf('Scheme 2 DA gamma-scan figure saved to %s\n', FigureDir)

function AverDerr = CalcuAverDerr(EstimatedPoints,GroundTruthPoints)
%CALCUAVERDERR Compute mean nearest-source error in the XY plane.
% Input:
%   EstimatedPoints - 3-by-N estimated source positions.
%   GroundTruthPoints - 3-by-K true UE and mirror-source positions.
% Output:
%   AverDerr - Mean distance from each estimate to its nearest true source.

if isempty(EstimatedPoints)
    AverDerr = NaN;
    return
end

EstimatedPoints = EstimatedPoints(:,all(isfinite(EstimatedPoints),1));
GroundTruthPoints = GroundTruthPoints(:,all(isfinite(GroundTruthPoints),1));
PointNum = size(EstimatedPoints,2);
if PointNum == 0 || isempty(GroundTruthPoints)
    AverDerr = NaN;
    return
end

MatchedDerrs = zeros(PointNum,1);
for i = 1:PointNum
    MatchedDerrs(i) = min(vecnorm(EstimatedPoints(1:2,i)-GroundTruthPoints(1:2,:),2,1));
end
AverDerr = mean(MatchedDerrs);

end

function HighConfidencePoints = SelectHighConfidencePoints(LocatedPoints,Beliefs,ValidIndexs,gamma)
%SELECTHIGHCONFIDENCEPOINTS Retain the top confidence-ranked valid candidates.
% Input:
%   LocatedPoints - 3-by-N candidate positions.
%   Beliefs - N positioning-confidence values.
%   ValidIndexs - Candidate indices passing the height and finite-value checks.
%   gamma - Fraction of valid candidates to retain.
% Output:
%   HighConfidencePoints - Positions of the retained candidates.

if isempty(LocatedPoints) || isempty(Beliefs) || isempty(ValidIndexs)
    HighConfidencePoints = [];
    return
end

gamma = max(0,min(1,gamma));
[~,SortedLocalIdxs] = sort(Beliefs(ValidIndexs),'descend');
SortedValidIndexs = ValidIndexs(SortedLocalIdxs);
TopNum = max(1,ceil(gamma*length(SortedValidIndexs)));
HighConfidencePoints = LocatedPoints(:,SortedValidIndexs(1:TopNum));

end

function LocatedPoints = LocateGroundTruthAssociatedPoints(LoFObservations,TargetPathIndexs,MinANum,pAs,z0,LDoALocaFunc)
%LOCATEGROUNDTRUTHASSOCIATEDPOINTS Position sources using known path labels.
% Input:
%   LoFObservations - Per-AP cells of quantized arrival distances in meters.
%   TargetPathIndexs - AP-by-target matrix of ground-truth path indices.
%   MinANum - Minimum associated AP count required for positioning.
%   pAs - 3-by-M AP coordinate matrix in meters.
%   z0 - Optional known source height.
%   LDoALocaFunc - Handle to the arrival-difference positioning function.
% Output:
%   LocatedPoints - 3-by-target position matrix; failed targets remain NaN.

TargetNum = size(TargetPathIndexs,2);
LocatedPoints = nan(3,TargetNum);

for TargetIdx = 1:TargetNum
    ActiveAPIdxs = find(TargetPathIndexs(:,TargetIdx));
    if numel(ActiveAPIdxs) < MinANum
        continue
    end

    LoFs = nan(numel(ActiveAPIdxs),1);
    IsValid = true;
    for k = 1:numel(ActiveAPIdxs)
        APIdx = ActiveAPIdxs(k);
        PathIdx = TargetPathIndexs(APIdx,TargetIdx);
        if PathIdx < 1 || PathIdx > numel(LoFObservations{APIdx})
            IsValid = false;
            break
        end
        LoFs(k) = LoFObservations{APIdx}(PathIdx);
    end
    if ~IsValid
        continue
    end

    [LocatedPoints(:,TargetIdx),~] = LDoALocaFunc( ...
        pAs(:,ActiveAPIdxs), ...
        LoFs, ...
        [], ...
        1, ...
        1, ...
        1, ...
        1, ...
        2, ...
        z0, ...
        false);
end

end

function ValidPoints = LocateRandomSearchPathCombs(pAs,MinANum,LoAObservations,APSelectedProbabilities,TotalNum,LDoALocaFunc)
%LOCATERANDOMSEARCHPATHCOMBS Generate positions from random path associations.
% Input:
%   pAs - 3-by-M AP coordinate matrix in meters.
%   MinANum - Minimum number of selected APs.
%   LoAObservations - Per-AP cells of candidate arrival distances.
%   APSelectedProbabilities - Independent AP-selection probabilities.
%   TotalNum - Requested number of valid random positions.
%   LDoALocaFunc - Handle to the arrival-difference positioning function.
% Output:
%   ValidPoints - Randomly associated finite source positions.
% Function: Sample AP subsets and one path per selected AP until TotalNum
% valid positions are obtained or the attempt limit is reached.

APNum = size(pAs,2);
ValidPointNum = 0;
ValidPoints = zeros(3,TotalNum);

if isempty(APSelectedProbabilities) || any(~isfinite(APSelectedProbabilities))
    ValidPoints = [];
    return
end

MaxAttempts = max(1000,TotalNum*1000);
AttemptCounter = 0;
while ValidPointNum < TotalNum && AttemptCounter < MaxAttempts
    AttemptCounter = AttemptCounter + 1;
    APSelectedFlags = rand(APNum,1)<APSelectedProbabilities;
    APSelectedFlags = APSelectedFlags & cellfun(@(x) ~isempty(x), LoAObservations);

    if nnz(APSelectedFlags)>=MinANum
        LoAObservations_selected = LoAObservations(APSelectedFlags);
        LComb = cellfun(@(x) x(randi(length(x))), LoAObservations_selected, 'UniformOutput', true);
        [p_esti,~] = LDoALocaFunc(pAs(:,APSelectedFlags),LComb,[],1,1,1,1,2,[],false);
        if all(isfinite(p_esti)) && abs(p_esti(1))>0
            ValidPointNum = ValidPointNum+1;
            ValidPoints(:,ValidPointNum) = p_esti;
        end
    end
end

ValidPoints = ValidPoints(:,1:ValidPointNum);

end
