%MAIN_PERFORMANCE_LOSDETECTION_LOCA Evaluate UE localization versus LoS availability.
% The experiment compares proposed Stage 1 localization with ground-truth
% association over randomly masked LoS-path availability probabilities.
clc
clear

MCNum = 400;
OnlyObserXYLocaErr = true;

tic;

% Localization uses one parameter configuration; wall-scoring settings are inactive.
Scheme2ParamSets = repmat(struct('SolveWallTols',[],'gamma',[],'beta',[]),1,1);

Scheme2ParamSets(1).SolveWallTols = [1,1,1];Scheme2ParamSets(1).gamma = 0.1;Scheme2ParamSets(1).beta = 0.3;
LoSDetectionProbScan_Loca(Scheme2ParamSets,MCNum,OnlyObserXYLocaErr)

toc;

function LoSDetectionProbScan_Loca(Scheme2ParamSets,MCNum,OnlyObserXYLocaErr)
%LOSDETECTIONPROBSCAN_LOCA Run the LoS-availability localization experiment.
% Input:
%   Scheme2ParamSets - Proposed-method parameter configurations.
%   MCNum - Number of Monte Carlo trials per availability probability.
%   OnlyObserXYLocaErr - Plot XY error when true, otherwise 3-D error.
% Output:
%   None. The function saves a dual-axis performance figure.
% Function: Mask LoS paths, extract multipath observations, compare localization
% success/error with ground-truth association, and aggregate Monte Carlo results.

tic_start = tic;

ScriptDir = fileparts(mfilename('fullpath'));
addpath(ScriptDir)

rng(10);

% Load shared system and algorithm parameters.
[Paras,funcs] = Base('SubcarrierSpacing',120e3);
UENum = Paras.UENum;

Bandwidth = Paras.SubcarrierSpacing*Paras.SubcarrierNum;
TargetWorkerNum = 10;

NominalLoSDetectionProbability = 0.87;
LoSDetectionProbabilities = [0.35,0.50,0.65,0.75,0.87];

ConfigNum = length(Scheme2ParamSets);
Para_length = length(LoSDetectionProbabilities);

LocaRate_GroundTruth = zeros(Para_length,1);
LocaErr_GroundTruth = nan(Para_length,1);

LocaRate_Scheme2 = zeros(Para_length,ConfigNum);
LocaErr_Scheme2 = nan(Para_length,ConfigNum);

LocaErrXY_GroundTruth = nan(Para_length,1);
LocaErrXY_Scheme2 = nan(Para_length,ConfigNum);

% Recalculate noise variance for each masked dataset to preserve the target SNR.
fdnoise_sigma2_LoSPs = zeros(Para_length,1);
for Paraidx = 1:Para_length
    [~,pAs,~,LoAs_NoMask,PathGains_NoMask] = funcs.LoadRaytracedPaths(LoSDetectionProbabilities(Paraidx),NominalLoSDetectionProbability);
    fdnoise_sigma2_LoSPs(Paraidx) = funcs.CalcuFreqDomNoiseSigma2(LoAs_NoMask,PathGains_NoMask,Paras.SubcarrierNum,Paras.SubcarrierSpacing,Paras.SNR1_dB);
end
APNum = size(pAs,2);

fprintf('LoS detection probability scan (Loca only): MCNum=%i, fixed SCS=%.0f kHz\n', MCNum, Paras.SubcarrierSpacing/1e3)

CurrentPool = gcp('nocreate');
if isempty(CurrentPool)
    parpool('local',TargetWorkerNum);
elseif CurrentPool.NumWorkers ~= TargetWorkerNum
    warning('Current parallel pool has %i workers; this script will use the existing pool.', CurrentPool.NumWorkers)
end

for Paraidx = 1:Para_length
    LoSDetectionProbability = LoSDetectionProbabilities(Paraidx);

    % Per-trial counters avoid shared writes inside PARFOR.
    LocaSuccessCount_GT_MC = zeros(MCNum,1);
    LocaErrSum_GT_MC = zeros(MCNum,1);

    LocaSuccessCount_S2_MC = zeros(MCNum,ConfigNum);
    LocaErrSum_S2_MC = zeros(MCNum,ConfigNum);

    LocaErrXYSum_GT_MC = zeros(MCNum,1);
    LocaErrXYSum_S2_MC = zeros(MCNum,ConfigNum);

    parfor MCidx = 1:MCNum
        rng(2026+1000*Paraidx + MCidx)

        tic_worker = tic;

        LocalLocaSuccessCount_GT = 0;
        LocalLocaErrSum_GT = 0;

        LocalLocaSuccessCount_S2 = zeros(1,ConfigNum);
        LocalLocaErrSum_S2 = zeros(1,ConfigNum);

        LocalLocaErrXYSum_GT = 0;
        LocalLocaErrXYSum_S2 = zeros(1,ConfigNum);

        [~,~,UEPositions_DataSet,LoAs_UEs_DataSet,PathGains_UEs_DataSet,TrueTargetPathIndexs_UEs_DataSet] = ...
            funcs.LoadRaytracedPaths(LoSDetectionProbability,NominalLoSDetectionProbability);
        
        SampledUEidxs_inDataset = randperm(length(LoAs_UEs_DataSet),UENum);

        [ExtractedLoAs_UEs, ExtractedPathPowers_UEs] = funcs.ExtractPaths( ...
            LoAs_UEs_DataSet, ...
            PathGains_UEs_DataSet, ...
            SampledUEidxs_inDataset, ...
            Paras.SubcarrierNum, ...
            Paras.SubcarrierSpacing, ...
            fdnoise_sigma2_LoSPs(Paraidx), ...
            Paras.MaxDelaySpread, ...
            Paras.OverSamplingFactor, ...
            Paras.delta_OMP, ...
            Paras.MaxIterNum_OMP);

        % Quantize ray-traced arrival distances for the ground-truth baseline.
        LoAs4PerfAssoci_UEs = LoAs_UEs_DataSet(SampledUEidxs_inDataset);
        TsOfOFDM_L = Paras.LightSpeed/(Paras.SubcarrierSpacing*Paras.SubcarrierNum);
        for UEidx = 1:UENum
            for APidx = 1:APNum
                CurrentLoAs = LoAs4PerfAssoci_UEs{UEidx}{APidx};
                if isempty(CurrentLoAs)
                    continue
                end
                LoAModRemains = mod(CurrentLoAs,TsOfOFDM_L);
                LoAs4PerfAssoci_UEs{UEidx}{APidx} = CurrentLoAs+(LoAModRemains>TsOfOFDM_L/2).*(TsOfOFDM_L-LoAModRemains)+(LoAModRemains<=TsOfOFDM_L/2).*(-LoAModRemains);
            end
        end

        for UEidx = 1:UENum
            pU = UEPositions_DataSet(:,SampledUEidxs_inDataset(UEidx));

            % Ground-truth association benchmark.
            LoAs4PerfAssoci = LoAs4PerfAssoci_UEs{UEidx};
            TrueTargetPathIndexs = TrueTargetPathIndexs_UEs_DataSet{SampledUEidxs_inDataset(UEidx)};
            [LocatedUE_GT,~,~] = funcs.PerfectAssociationSLAMWLS(LoAs4PerfAssoci,TrueTargetPathIndexs,Paras.MinANum,pAs,Paras.z0_UE);
            if ~isempty(LocatedUE_GT)
                LocalLocaSuccessCount_GT = LocalLocaSuccessCount_GT + 1;
                LocalLocaErrSum_GT = LocalLocaErrSum_GT + norm(LocatedUE_GT-pU);
                LocalLocaErrXYSum_GT = LocalLocaErrXYSum_GT + norm(LocatedUE_GT(1:2)-pU(1:2));
            end

            % Proposed Stage 1 UE localization.
            for ConfigIdx = 1:ConfigNum
                CurrentConfig = Scheme2ParamSets(ConfigIdx);
                [LocatedUE_S2,~,~,~,~,~] = funcs.SLAM( ...
                    ExtractedLoAs_UEs{UEidx}, ...
                    ExtractedPathPowers_UEs{UEidx}, ...
                    pAs, ...
                    Paras.MinANum, ...
                    Paras.PowerCheckThres, ...
                    Paras.RandNumPerSeed, ...
                    CurrentConfig.SolveWallTols, ...
                    Paras.sigma_t, ...
                    Paras.sigma_P, ...
                    CurrentConfig.gamma, ...
                    2, Paras.z0_UE, Paras.UseStage1, CurrentConfig.beta, ...
                    "UEBased", "Stage1Only");

                if ~isempty(LocatedUE_S2)
                    LocalLocaSuccessCount_S2(ConfigIdx) = LocalLocaSuccessCount_S2(ConfigIdx) + 1;
                    LocalLocaErrSum_S2(ConfigIdx) = LocalLocaErrSum_S2(ConfigIdx) + norm(LocatedUE_S2-pU);
                    LocalLocaErrXYSum_S2(ConfigIdx) = LocalLocaErrXYSum_S2(ConfigIdx) + norm(LocatedUE_S2(1:2)-pU(1:2));
                end
            end
        end

        LocaSuccessCount_GT_MC(MCidx) = LocalLocaSuccessCount_GT;
        LocaErrSum_GT_MC(MCidx) = LocalLocaErrSum_GT;

        LocaSuccessCount_S2_MC(MCidx,:) = LocalLocaSuccessCount_S2;
        LocaErrSum_S2_MC(MCidx,:) = LocalLocaErrSum_S2;

        LocaErrXYSum_GT_MC(MCidx) = LocalLocaErrXYSum_GT;
        LocaErrXYSum_S2_MC(MCidx,:) = LocalLocaErrXYSum_S2;

        fprintf('LoS p=%.2f, MC %i/%i finished, take %.1f s\n', LoSDetectionProbability, MCidx, MCNum, toc(tic_worker))
    end

    LocaSuccessCount_GT = sum(LocaSuccessCount_GT_MC);
    LocaErrSum_GT = sum(LocaErrSum_GT_MC);

    LocaSuccessCount_S2 = sum(LocaSuccessCount_S2_MC,1).';
    LocaErrSum_S2 = sum(LocaErrSum_S2_MC,1).';

    LocaErrXYSum_GT = sum(LocaErrXYSum_GT_MC);
    LocaErrXYSum_S2 = sum(LocaErrXYSum_S2_MC,1).';

    % Success rate uses all UE trials; error averages only successful estimates.
    TotalLocaTrials = UENum*MCNum;

    LocaRate_GroundTruth(Paraidx) = LocaSuccessCount_GT/TotalLocaTrials;
    LocaErr_GroundTruth(Paraidx) = funcs.SafeDivide(LocaErrSum_GT,LocaSuccessCount_GT);
    LocaErrXY_GroundTruth(Paraidx) = funcs.SafeDivide(LocaErrXYSum_GT,LocaSuccessCount_GT);

    for ConfigIdx = 1:ConfigNum
        LocaRate_Scheme2(Paraidx,ConfigIdx) = LocaSuccessCount_S2(ConfigIdx)/TotalLocaTrials;
        LocaErr_Scheme2(Paraidx,ConfigIdx) = funcs.SafeDivide(LocaErrSum_S2(ConfigIdx),LocaSuccessCount_S2(ConfigIdx));
        LocaErrXY_Scheme2(Paraidx,ConfigIdx) = funcs.SafeDivide(LocaErrXYSum_S2(ConfigIdx),LocaSuccessCount_S2(ConfigIdx));
    end
end

FigureDir = fullfile(ScriptDir,'Figures');
if ~exist(FigureDir,'dir')
    mkdir(FigureDir)
end

Resolution_m = Paras.LightSpeed/(Paras.SubcarrierSpacing*Paras.SubcarrierNum);
BaseTitle = sprintf('MC=%i, sigma_P=%.1f dB, sigma_t=%.2f res, BW=%.2f MHz, SNR=%i dB, nominal LoS p=%.2f,z0=%.1f', ...
    MCNum, Paras.sigma_P, Paras.sigma_t/Resolution_m, Bandwidth/1e6, Paras.SNR1_dB, NominalLoSDetectionProbability, Paras.z0_UE);

if OnlyObserXYLocaErr
    LocaErr_GT_Plot = LocaErrXY_GroundTruth;
    LocaErr_S2_Plot = LocaErrXY_Scheme2;
    LocaErrLabel = 'Mean XY-Plane Error (m)';
else
    LocaErr_GT_Plot = LocaErr_GroundTruth;
    LocaErr_S2_Plot = LocaErr_Scheme2;
    LocaErrLabel = 'Average Localization Error (XYZ) (m)';
end

% Plot localization error and success rate on separate vertical axes.
funcs.PlotDualAxisPerformance( ...
    LoSDetectionProbabilities, ...
    LocaRate_GroundTruth, ...
    LocaErr_GT_Plot, ...
    LocaRate_Scheme2, ...
    LocaErr_S2_Plot, ...
    Scheme2ParamSets, ...
    {'Localization vs. LoS detection probability (beta scan)', BaseTitle}, ...
    'Success Rate', ...
    LocaErrLabel, ...
    fullfile(FigureDir,sprintf('LoSDetection_Localization_%s', char(datetime('now', 'Format', 'MMddHHmmss')))))

fprintf('\n LoS detection probability scan (Loca) figures saved to %s, total time: %.1f s\n', FigureDir, toc(tic_start))
end
