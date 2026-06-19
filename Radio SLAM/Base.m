function [Paras,funcs] = Base(varargin)
%BASE Initialize shared simulation parameters and public function handles.
% Input:
%   varargin - Optional name-value pairs that override fields in Paras.
% Output:
%   Paras - Structure containing system, extraction, and SLAM parameters.
%   funcs - Structure of function handles used by the experiment scripts.
% Function: Apply parameter overrides and update sigma_t from the configured
% bandwidth unless sigma_t is explicitly provided.

% System Parameters---------------------
Paras.UENum = 8;
Paras.LightSpeed = 3e8;
Paras.SNR1_dB = 15;
Paras.SNR2_dB = 10;
Paras.SubcarrierNum = 2048;
Paras.SubcarrierSpacing = 120e3;

% Method Parameters---------------------
% Extract path---
Paras.delta_OMP = 1e-6;
Paras.MaxIterNum_OMP = 10;
Paras.MaxDelaySpread = 200e-9;
Paras.OverSamplingFactor = 1;
% SLAM
Paras.MinANum = 5;
Paras.SolveWallTols = [1,1,1];
Paras.PowerCheckThres = 10; % Linear power ratio equivalent to a 10 dB difference
Paras.RandNumPerSeed = 300;

Paras.sigma_P = 4; % dB
Paras.sigma_t = 2*Paras.LightSpeed/(Paras.SubcarrierSpacing*Paras.SubcarrierNum); % Distance-domain (m), consistent with LoA-based LDoFOffsets and Vt
Paras.gamma = 0.15; % Proportion parameter for high-confidence candidate selection
Paras.beta = 1.0; % Stage 1 chi-square gate probability: e'*R_t^{-1}*e < chi2_{d-3}(beta)

Paras.z0_UE = [];% Empty means that LDoALoca estimates height without a z-coordinate prior
Paras.UseStage1 = true;

SigmaTOverridden = false;
if mod(nargin,2) ~= 0
    error('CommonParas:InvalidInput','Parameters must be provided as name-value pairs.')
end

for InputIdx = 1:2:nargin
    FieldName = char(varargin{InputIdx});
    if ~isfield(Paras,FieldName)
        error('CommonParas:UnknownParameter','Unknown parameter "%s".',FieldName)
    end
    Paras.(FieldName) = varargin{InputIdx+1};
    SigmaTOverridden = SigmaTOverridden || strcmp(FieldName,'sigma_t');
end

if ~SigmaTOverridden
    Paras.sigma_t = 2*Paras.LightSpeed/(Paras.SubcarrierSpacing*Paras.SubcarrierNum);
end

funcs.PlotDualAxisPerformance = @PlotDualAxisPerformance;
funcs.SafeDivide = @SafeDivide;
funcs.PerfectAssociationSLAMWLS = @PerfectAssociationSLAMWLS;
funcs.CalcuWallWeights = @CalcuWallWeights;
funcs.FuseWallEstimatesAcrossUEs = @FuseWallEstimatesAcrossUEs;
funcs.LoadRaytracedPaths = @LoadRaytracedPaths;
funcs.CalcuFreqDomNoiseSigma2 = @CalcuFreqDomNoiseSigma2;
funcs.ExtractPaths = @ExtractPaths;
funcs.SLAM = @SLAM;
funcs.LDoALoca = @LDoALoca;
end


function [LocatedUE,Wallxy,WallWeights] = PerfectAssociationSLAMWLS(LoFObservations,TargetPathIndexs,MinANum,pAs,z0)
%PERFECTASSOCIATIONSLAMWLS Run localization and mapping with known path associations.
% Input:
%   LoFObservations - Per-AP cells of path arrival distances in meters.
%   TargetPathIndexs - AP-by-6 path indices for the UE and its five mirrors.
%   MinANum - Minimum number of associated AP paths required for positioning.
%   pAs - 3-by-M AP coordinate matrix in meters.
%   z0 - Optional known source height; use [] to estimate height.
% Output:
%   LocatedUE - Estimated 3-D UE position, or [] if UE positioning fails.
%   Wallxy - Wall coordinates ordered as [x+; y+; x-; y-].
%   WallWeights - Inverse-residual weights for the four wall estimates.
% Function: Position ground-truth-associated sources and estimate walls from
% UE/mirror-UE midpoints.

if nargin < 5 || isempty(z0)
    z0 = [];
end

LocatedUEandMirror = zeros(3,6);
GeoMSEs = nan(6,1);
SuccessFlags = false(6,1);
for TargetIdx = 1:6
    ActiveAPIdxs = find(TargetPathIndexs(:,TargetIdx));
    if length(ActiveAPIdxs) >= MinANum
        SuccessFlags(TargetIdx) = true;
        LoFs = arrayfun(@(i) LoFObservations{i}(TargetPathIndexs(i,TargetIdx)), ActiveAPIdxs);
        [LocatedUEandMirror(:,TargetIdx),~,Metrics] = LDoALoca(pAs(:,ActiveAPIdxs),LoFs,[],1,1,1,1,2,z0);
        GeoMSEs(TargetIdx) = Metrics.GeoMSE;
    end
end

[LocatedUE,Wallxy,WallWeights] = SolveGroundTruthWallxyWLS(LocatedUEandMirror,GeoMSEs,SuccessFlags);

end

function [LocatedUE,Wallxy,WallWeights] = SolveGroundTruthWallxyWLS(LocatedUEandMirror,GeoMSEs,SuccessFlags)
%SOLVEGROUNDTRUTHWALLXYWLS Convert associated UE/mirror positions to wall estimates.
% Input:
%   LocatedUEandMirror - 3-by-6 positions of the UE and five mirror sources.
%   GeoMSEs - Geometric residual MSE for each positioned source.
%   SuccessFlags - Six-element flags indicating successful source positioning.
% Output:
%   LocatedUE - Estimated UE position, or [] when the UE estimate is missing.
%   Wallxy - UE/mirror midpoint wall coordinates [x+; y+; x-; y-].
%   WallWeights - Inverse-MSE weights for the four wall estimates.
% Function: Form wall estimates from successful wall-induced mirror sources.

Wallxy = nan(4,1);
WallWeights = zeros(4,1);
if ~SuccessFlags(1)
    LocatedUE = [];
    return
end

LocatedUE = LocatedUEandMirror(:,1);
for TargetIdx = 2:5
    if ~SuccessFlags(TargetIdx)
        continue
    end

    DirectionIdx = TargetIdx - 1;
    if TargetIdx == 2 || TargetIdx == 4
        Wallxy(DirectionIdx) = (LocatedUEandMirror(1,1)+LocatedUEandMirror(1,TargetIdx))/2;
    else
        Wallxy(DirectionIdx) = (LocatedUEandMirror(2,1)+LocatedUEandMirror(2,TargetIdx))/2;
    end
    WallWeights(DirectionIdx) = 1/max(GeoMSEs(TargetIdx),eps);
end

end

function WallWeights = CalcuWallWeights(CandidateInfo)
%CALCUWALLWEIGHTS Extract accumulated wall-estimate weights from candidate data.
% Input:
%   CandidateInfo - SLAM diagnostic structure containing mirror sets,
%                   candidate indices, and residual-based weight statistics.
% Output:
%   WallWeights - Four accumulated weights ordered as [x+; y+; x-; y-].
% Function: Use precomputed wall weight sums when available; otherwise sum
% inverse candidate residual statistics over each selected mirror set.

WallWeights = zeros(4,1);
if isfield(CandidateInfo,'WallWeightSums') && ~isempty(CandidateInfo.WallWeightSums)
    CandidateWallWeights = CandidateInfo.WallWeightSums(:);
    CandidateWallWeights(~isfinite(CandidateWallWeights) | CandidateWallWeights < 0) = 0;
    if any(CandidateWallWeights > 0)
        WallWeights = CandidateWallWeights;
        return
    end
end

if isempty(CandidateInfo) || isempty(CandidateInfo.WallPoolIndexs)
    return
end

if isfield(CandidateInfo,'WallWeightStats') && ~isempty(CandidateInfo.WallWeightStats)
    WeightStats = CandidateInfo.WallWeightStats;
elseif isfield(CandidateInfo,'GeoMSEs') && ~isempty(CandidateInfo.GeoMSEs)
    WeightStats = CandidateInfo.GeoMSEs;
else
    return
end

for DirectionIdx = 1:4
    MirrorSet = CandidateInfo.MirrorSets{DirectionIdx};
    if isempty(MirrorSet)
        continue
    end

    CandidateIndexs = CandidateInfo.WallPoolIndexs(MirrorSet);
    CandidateIndexs = CandidateIndexs(isfinite(CandidateIndexs) & CandidateIndexs >= 1 & CandidateIndexs <= numel(WeightStats));
    if isempty(CandidateIndexs)
        continue
    end
    Weights = 1 ./ max(WeightStats(CandidateIndexs),eps);
    Weights(~isfinite(Weights)) = 0;
    WallWeights(DirectionIdx) = sum(Weights);
end

end

function Stat = GetWallWeightStat(Metrics)
%GETWALLWEIGHTSTAT Select a positive residual statistic for wall weighting.
% Input:
%   Metrics - Positioning metrics returned by LDoALoca.
% Output:
%   Stat - Residual statistic, lower-bounded by eps; inf if unavailable.
% Function: Prefer reduced chi-square, then normalized geometric MSE, and
% finally geometric MSE as the inverse-precision weighting statistic.

if isfield(Metrics,'ReducedGeoChiSquare') && isfinite(Metrics.ReducedGeoChiSquare)
    Stat = Metrics.ReducedGeoChiSquare;
elseif isfield(Metrics,'WeightedGeoMSE') && isfinite(Metrics.WeightedGeoMSE)
    Stat = Metrics.WeightedGeoMSE;
elseif isfield(Metrics,'GeoMSE') && isfinite(Metrics.GeoMSE)
    Stat = Metrics.GeoMSE;
else
    Stat = inf;
end

Stat = max(Stat,eps);
end

function FusedWallxy = FuseWallEstimatesAcrossUEs(Wallxy_UEs,WallWeights_UEs)
%FUSEWALLESTIMATESACROSSUES Fuse per-UE wall estimates by weighted averaging.
% Input:
%   Wallxy_UEs - 4-by-K wall estimates, ordered as [x+; y+; x-; y-].
%   WallWeights_UEs - 4-by-K nonnegative weights for those estimates.
% Output:
%   FusedWallxy - Four fused wall coordinates; NaN where no valid estimate exists.
% Function: Independently fuse each wall using finite, strictly positive weights.

FusedWallxy = nan(4,1);
for DirectionIdx = 1:4
    WallEstimates = Wallxy_UEs(DirectionIdx,:);
    WallWeights = WallWeights_UEs(DirectionIdx,:);
    ValidFlags = isfinite(WallEstimates) & isfinite(WallWeights) & WallWeights > 0;
    if any(ValidFlags)
        FusedWallxy(DirectionIdx) = sum(WallWeights(ValidFlags).*WallEstimates(ValidFlags))/sum(WallWeights(ValidFlags));
    end
end

end

function sigma2 = CalcuFreqDomNoiseSigma2(LoFs_UEs,PathGains_UEs,SubcarrierNum,SubcarrierSpacing,SNR_dB)
%CALCUFREQDOMNOISESIGMA2 Determine frequency-domain noise variance for a target SNR.
% Input:
%   LoFs_UEs - UE/AP cell array of ray-traced path arrival distances in meters.
%   PathGains_UEs - Matching UE/AP cell array of complex path gains.
%   SubcarrierNum - Total number of OFDM subcarriers.
%   SubcarrierSpacing - OFDM subcarrier spacing in hertz.
%   SNR_dB - Target average per-subcarrier SNR in decibels.
% Output:
%   sigma2 - Complex frequency-domain noise variance.
% Function: Reconstruct noiseless channels and divide their average power by
% the linear target SNR.

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

function [ExtractedLoFs_UEs, ExtractedPathPowers_UEs] ...
            = ExtractPaths(LoFs_UEs,PathGains_UEs,SampledUEidxs_inDataset,SubcarrierNum,SubcarrierSpacing,fdnoise_sigma2,MaxDelaySpread,OverSamplingFactor,delta_OMP,MaxIterNum_OMP)     
%EXTRACTPATHS Simulate pilot observations and extract multipath components by OMP.
% Input:
%   LoFs_UEs - Dataset of UE/AP path arrival distances in meters.
%   PathGains_UEs - Dataset of matching complex path gains.
%   SampledUEidxs_inDataset - Dataset indices of the UEs used in this trial.
%   SubcarrierNum - Total number of OFDM subcarriers shared by the UEs.
%   SubcarrierSpacing - OFDM subcarrier spacing in hertz.
%   fdnoise_sigma2 - Complex frequency-domain noise variance.
%   MaxDelaySpread - Maximum extracted delay in seconds.
%   OverSamplingFactor - Delay-grid oversampling factor.
%   delta_OMP - Relative residual-energy stopping threshold for OMP.
%   MaxIterNum_OMP - Maximum number of OMP iterations.
% Output:
%   ExtractedLoFs_UEs - Per-UE/per-AP extracted arrival distances in meters.
%   ExtractedPathPowers_UEs - Matching extracted powers in linear scale,
%                             sorted from strongest to weakest.
% Function: Generate orthogonal interleaved-pilot channels, add noise, and
% recover grid-based multipath delays and powers; empty AP data are supported.

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
    
    % Generate frequency-domain channels using raytraced multipath.
    FreqDomainChannel = zeros(UENum,APNum,SubcarrierNum/UENum);
    FirstPathLoAs = nan(UENum,APNum);
    for UEidx = 1:UENum
        UEidx_inDataset = SampledUEidxs_inDataset(UEidx);
        TransmitTimeOffSet = rand*DelaysResolution;
        for APidx = 1:APNum
            if isempty(LoFs_UEs{UEidx_inDataset}{APidx})
                continue
            end

            FirstPathLoA = min(LoFs_UEs{UEidx_inDataset}{APidx})+TransmitTimeOffSet;
            FirstPathLoAs(UEidx,APidx) = FirstPathLoA;
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

    % Add noise to the frequency-domain channels.
    FreqDomainChannel_noised = FreqDomainChannel + sqrt(fdnoise_sigma2/2)*randn(UENum,APNum,SubcarrierNum/UENum) ...
                                +sqrt(-1)*sqrt(fdnoise_sigma2/2)*randn(UENum,APNum,SubcarrierNum/UENum);

    % Extract CIR using OMP.
    ColumnNum = floor(MaxDelaySpread/DelaysResolution);
    Taus_tick = (0:ColumnNum-1).'*DelaysResolution;
    MeasureMatrix = zeros(SubcarrierNum/UENum,ColumnNum);

    ExtractedLoFs_UEs = cell(UENum,1);
    ExtractedPathPowers_UEs = cell(UENum,1);
    for UEidx = 1:UENum
        ExtractedLoFs = cell(APNum,1);
        ExtractedPathPowers = cell(APNum,1);
        for APidx = 1:APNum
            if isnan(FirstPathLoAs(UEidx,APidx))
                ExtractedLoFs{APidx} = [];
                ExtractedPathPowers{APidx} = [];
                continue
            end

            % Construct sensing matrix.
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
            ExtractedLoFs{APidx} = FirstPathLoAs(UEidx,APidx)-mod(FirstPathLoAs(UEidx,APidx),DistRes_m)+Est_Taus*LightSpeed;
            ExtractedPathPowers{APidx} = abs(Est_Amps).^2;

        end
        ExtractedLoFs_UEs{UEidx} = ExtractedLoFs;
        ExtractedPathPowers_UEs{UEidx} = ExtractedPathPowers;
    end
    
end

function [p_esti,Belief,Metrics] = LDoALoca(pAs,LoAs,PathPowers_dB,sigma_t,sigma_P,V_t,V_P,ConfidenceScheme,z0,UseGeoCond)
%LDOALOCA Position one source from arrival-distance differences and score it.
% Input:
%   pAs - 3-by-M AP coordinate matrix in meters.
%   LoAs - M path arrival distances in meters.
%   PathPowers_dB - M received path powers in decibels; [] disables power scoring.
%   sigma_t - Distance-domain standard deviation of LDoF residuals in meters.
%   sigma_P - Standard deviation of compensated powers in decibels.
%   V_t - Empirical arrival-distance support length in meters.
%   V_P - Empirical received-power support length in decibels.
%   ConfidenceScheme - 1 for the baseline score or 2 for the LLRT score.
%   z0 - Optional known source height; use [] to estimate height.
%   UseGeoCond - Optional flag adding the geometric conditioning term to
%                Scheme 2; defaults to true.
% Output:
%   p_esti - Estimated 3-D source position.
%   Belief - Total positioning confidence.
%   Metrics - Structure of residuals, confidence components, and diagnostics.
% Function: Apply a coplanar-AP modification of Chan positioning and evaluate
% geometric and, when supplied, power consistency.

if nargin < 10 || isempty(UseGeoCond)
    UseGeoCond = true;
end

M = size(pAs,2);
pAsMassCenter = sum(pAs, 2) / M;
[~, RefA_id] = min(sum((pAs - pAsMassCenter).^2, 1));

NotRefA_ids = 1:M;
NotRefA_ids(RefA_id) = [];

LDoFs = LoAs(NotRefA_ids) - LoAs(RefA_id);

pA_ref = pAs(:,RefA_id);
pAs_notRef = pAs;
pAs_notRef(:,RefA_id) = [];

DeltapA = pAs_notRef-pA_ref;
DeltapA = DeltapA(1:2,:);
DeltaK = (vecnorm(pAs_notRef).^2).'-norm(pA_ref)^2;

h = -1/2*(DeltaK-LDoFs.^2);
G = -[DeltapA.' LDoFs];

z_estimation = pinv(G.'*G)*(G.'*h);

% Calculate the coordinates of the positioned point.
if ~isempty(z0)
    % Two-step xy decoupling (Case of Known Source Height, Paper eq:xyEst):
    % Step (i) already done: z_estimation(1:2) = [x0, y0] initial estimate.
    % Step (ii): compute z-constrained r_ref at [x0, y0], then solve reduced 2D LS.
    r_ref_z0 = sqrt((z_estimation(1)-pA_ref(1))^2 + (z_estimation(2)-pA_ref(2))^2 + (z0-pA_ref(3))^2);
    G_xy = G(:,1:2);
    h_tilde = h - G(:,3)*r_ref_z0;
    xy_est = pinv(G_xy.'*G_xy)*(G_xy.'*h_tilde);
    p_esti = [xy_est(1:2); z0];
else
    % Original 3D estimation.
    if z_estimation(3)<norm(pA_ref(1:2)-z_estimation(1:2))
        p_esti = [z_estimation(1:2);1];
    else
        zaxis_esti = pA_ref(3)-sqrt(max(0,z_estimation(3)^2-norm(pA_ref(1:2)-z_estimation(1:2))^2));
        p_esti = [z_estimation(1:2);zaxis_esti];
    end
end

% Geometric residuals in the same distance domain as LoAs.
LDoFs_esti = (vecnorm(p_esti-pAs_notRef)-norm(p_esti-pA_ref)).';
LDoFOffsets = LDoFs_esti-LDoFs;
GeoAbsMean = sum(abs(LDoFOffsets))/(M-1);
GeoMSE = sum(LDoFOffsets.^2)/(M-1);

HasPower = ~isempty(PathPowers_dB);
if HasPower
    d_est_vec = max(vecnorm(p_esti - pAs), 1);
    CompensatedPowers = PathPowers_dB(:).' + 20*log10(d_est_vec);
    CompensatedPowers_mean = sum(CompensatedPowers)/M;
    PowerResiduals = CompensatedPowers - CompensatedPowers_mean;
    VarCompensatedPower = sum((CompensatedPowers - CompensatedPowers_mean).^2) / M;
else
    CompensatedPowers = [];
    PowerResiduals = [];
    VarCompensatedPower = NaN;
end

sigma_t = max(sigma_t, eps);
sigma_P = max(sigma_P, eps);
V_t = max(V_t, eps);
V_P = max(V_P, eps);

if ConfidenceScheme == 1
    Belief_Geo = -GeoAbsMean;
    if HasPower
        Belief_Power = -VarCompensatedPower;
        Belief = Belief_Geo*VarCompensatedPower;
    else
        Belief_Power = 0;
        Belief = Belief_Geo;
    end
else
    % ConfidenceScheme == 2: Full-covariance LLRT.
    % R_t = sigma_t^2*(I_d + 11'),  R_t^{-1} = (1/sigma_t^2)*(I_d - 11'/M)
    d = M - 1;
    e_hat   = LDoFOffsets(:);
    e_wsum  = sum(e_hat)^2 / M;
    e_quad  = sum(e_hat.^2) - e_wsum;   % e'*(I - 11'/M)*e

    % Residual chi-square statistics with R_t^{-1} = (1/sigma_t^2)*(I - 11'/M).
    GeoChiSquare = e_quad / sigma_t^2;
    WeightedGeoMSE = GeoChiSquare / d;
    ResidualDoF = max(d - 3, 1);
    ReducedGeoChiSquare = GeoChiSquare / ResidualDoF;

    % Profile log-LR (full covariance form):
    Belief_Geo_prof = d*log(V_t/(sqrt(2*pi)*sigma_t)) ...
        - 0.5*log(M) ...
        - e_quad / (2*sigma_t^2);

    % Geometric condition term B_loc (only when UseGeoCond=true, i.e. Stage 1).
    if UseGeoCond
        d_ref      = p_esti - pA_ref;
        d_ref_norm = max(norm(d_ref), 1e-6);
        J_xy = zeros(d, 2);
        for mi = 1:d
            d_i      = p_esti - pAs_notRef(:,mi);
            d_i_norm = max(norm(d_i), 1e-6);
            J_xy(mi,:) = d_i(1:2)'/d_i_norm - d_ref(1:2)'/d_ref_norm;
        end
        P_center = eye(d) - ones(d,d)/M;
        A = J_xy' * P_center * J_xy;   % 2x2, dimensionless
        det_A = det(A);
        if det_A > 0
            B_loc = 0.5 * log(det_A);
        else
            B_loc = -inf;
        end
    else
        B_loc = 0;  % disabled for Stage 2
    end

    Belief_Geo = Belief_Geo_prof + B_loc;

    if HasPower
        Belief_Power = M*log(V_P/(sqrt(2*pi)*sigma_P)) - sum(PowerResiduals.^2)/(2*sigma_P^2);
    else
        Belief_Power = 0;
    end
    Belief = Belief_Geo + Belief_Power;
end

Metrics = struct();
Metrics.SchemeNo = ConfidenceScheme;
Metrics.RefA_id = RefA_id;
Metrics.LDoFOffsets = LDoFOffsets;
Metrics.GeoAbsMean = GeoAbsMean;
Metrics.GeoMSE = GeoMSE;
Metrics.CompensatedPowers = CompensatedPowers;
Metrics.PowerResiduals = PowerResiduals;
Metrics.VarCompensatedPower = VarCompensatedPower;
Metrics.GeoConfidence = Belief_Geo;
Metrics.PowerConfidence = Belief_Power;
Metrics.TotalConfidence = Belief;
if ConfidenceScheme == 2
    Metrics.GeoChiSquare = GeoChiSquare;             % e'*R_t^{-1}*e
    Metrics.WeightedGeoMSE = WeightedGeoMSE;         % (1/d)*e'*R_t^{-1}*e
    Metrics.ResidualDoF = ResidualDoF;               % d minus the 3D UE estimator dimension
    Metrics.ReducedGeoChiSquare = ReducedGeoChiSquare;
    Metrics.ProfileGeoConfidence = Belief_Geo_prof;
    Metrics.GeoLocScore = B_loc;
end

end

function [SpaceBound,pAs,UEPositions,LoFObservations_UEs,PathGainObservations_UEs,TrueTargetPathIndexs_UEs,MaskInfo] = LoadRaytracedPaths(LoSDetectionProbability,NominalLoSDetectionProbability)
%LOADRAYTRACEDPATHS Load Sionna multipath data and optionally mask LoS paths.
% Input:
%   LoSDetectionProbability - Target average LoS availability probability.
%   NominalLoSDetectionProbability - LoS availability in the unmasked dataset;
%                                    defaults to 0.87.
% Output:
%   SpaceBound - Room dimensions [xLength, yLength, height] in meters.
%   pAs - 3-by-M AP coordinate matrix in meters.
%   UEPositions - 3-by-U dataset of UE coordinates.
%   LoFObservations_UEs - Per-UE/per-AP path arrival distances in meters.
%   PathGainObservations_UEs - Matching complex path gains.
%   TrueTargetPathIndexs_UEs - Per-UE AP-by-6 indices for the UE and five
%                              first-order mirror sources.
%   MaskInfo - Structure describing requested, nominal, and realized LoS masking.
% Function: Randomly remove additional LoS paths to reach probabilities no
% greater than the nominal dataset availability, while preserving other paths.

    if nargin < 2 || isempty(NominalLoSDetectionProbability)
        NominalLoSDetectionProbability = 0.87;
    end
    if nargin < 1 || isempty(LoSDetectionProbability)
        LoSDetectionProbability = NominalLoSDetectionProbability;
    end

    if LoSDetectionProbability >= NominalLoSDetectionProbability
        ExtraLoSKeepProbability = 1;
    else
        ExtraLoSKeepProbability = max(0,LoSDetectionProbability/NominalLoSDetectionProbability);
    end

    ScriptDir = fileparts(mfilename('fullpath'));
    RaytracedPaths = load(fullfile(ScriptDir,'RaytracedPaths.mat')).RaytracedPaths; 

    SpaceBound = (RaytracedPaths{1}.SpaceBound).';
    pAs = RaytracedPaths{1}.APPositions;
    UENum = length(RaytracedPaths)-1;
    UEPositions = zeros(3,UENum);        
    LoFObservations_UEs = cell(UENum,1);
    PathGainObservations_UEs = cell(UENum,1);
    
    APNum = length(RaytracedPaths{2}.CIR_by_AP);
    TrueTargetPathIndexs_UEs = cell(UENum,1);
    OriginalLoSCount = 0;
    KeptLoSCount = 0;

    for u = 1:UENum
        ThisUE = RaytracedPaths{u+1};
        UEPositions(:,u) = ThisUE.UE_Pos;
        
        LoFs_by_AP = cell(APNum,1);
        PathGains_by_AP = cell(APNum,1);
        TrueTargetPathIndexs = zeros(APNum,6);
        for m = 1:APNum
            % Each row: [Delay, Gain_real, Gain_Imag, Type1, Type2, Object1, Object2, Point1, Point2]
            PathData = ThisUE.CIR_by_AP{m};

            if ~isempty(PathData)
                IsLoSPath = all(int8(PathData(:,4:5)) == [0,0],2);
                OriginalLoSCount = OriginalLoSCount + nnz(IsLoSPath);
                KeepMask = true(size(PathData,1),1);
                if ExtraLoSKeepProbability < 1
                    KeepMask(IsLoSPath) = rand(nnz(IsLoSPath),1) < ExtraLoSKeepProbability;
                end
                KeptLoSCount = KeptLoSCount + nnz(IsLoSPath & KeepMask);
                PathData = PathData(KeepMask,:);
            end

            PathNum = size(PathData,1);            
            if PathNum>0
                LoFs_by_AP{m} = PathData(:,1);
                PathGains = PathData(:,2)+sqrt(-1)*PathData(:,3);
                PathGains_by_AP{m} = PathGains;
                for PathIdx = 1:PathNum
                    PathType = int8(PathData(PathIdx,4:5));
                    PathInteractObject = PathData(PathIdx,6);
                    if all(PathType==[0,0])
                        TrueTargetPathIndexs(m,1) = PathIdx;
                    elseif all(PathType==[1,0]) 
                        % The object IDs for the four walls and the floor in Sionna are 4, 5, 6, 7, and 8.
                        TrueTargetPathIndexs(m,PathInteractObject-2) = PathIdx;
                    end
                end
            else
                LoFs_by_AP{m} = [];
                PathGains_by_AP{m} = [];
            end
        end
        LoFObservations_UEs{u} = LoFs_by_AP;
        PathGainObservations_UEs{u} = PathGains_by_AP;
        TrueTargetPathIndexs_UEs{u} = TrueTargetPathIndexs;
    end

    MaskInfo = struct();
    MaskInfo.TargetLoSDetectionProbability = LoSDetectionProbability;
    MaskInfo.NominalLoSDetectionProbability = NominalLoSDetectionProbability;
    MaskInfo.ExtraLoSKeepProbability = ExtraLoSKeepProbability;
    MaskInfo.OriginalLoSCount = OriginalLoSCount;
    MaskInfo.KeptLoSCount = KeptLoSCount;
    MaskInfo.TotalUEAPPairs = UENum*APNum;
    MaskInfo.RealizedLoSDetectionProbability = KeptLoSCount/(UENum*APNum);

end

function [omega_now,x_now] = OMP(A,y,delta,MaxIterNum)
%OMP Recover a sparse representation using orthogonal matching pursuit.
% Input:
%   A - Sensing matrix whose columns are candidate atoms.
%   y - Observation vector.
%   delta - Relative residual-energy stopping threshold.
%   MaxIterNum - Maximum number of selected atoms.
% Output:
%   omega_now - Sorted indices of the selected atoms.
%   x_now - Least-squares coefficients on the selected support.
% Function: Iteratively select the most correlated atom and refit all active
% coefficients until the residual threshold or iteration limit is reached.

    IterCounter = 0;
    omega_now = [];
    r_now = y;
    RelativeError = zeros(MaxIterNum,1);
    while true
        % Match、Persuit
        [~,index] = max(abs(A'*r_now)./vecnorm(A));
        omega_now = union(omega_now,index,'sorted');
        % Update
        A_omega_now = A(:,omega_now);
        x_now = A_omega_now \ y;
        % Judge
        r_now = y - A_omega_now*x_now;    
        IterCounter = IterCounter+1;
        RelativeError(IterCounter) = norm(r_now)^2/norm(y)^2;
        if (RelativeError(IterCounter) < delta) || (IterCounter >= MaxIterNum)
            break;
        end
    end
    
end

function [FoundLCombs,FoundPCombs,FoundLCombsIndexs] = ...
    PathCombSearch(LoAObservations, PathPowerObservations, pAs, PowerCheckThres, MinANum, RandNumPerSeed)
%PATHCOMBSEARCH Search a consistency graph for candidate LoS path combinations.
% Input:
%   LoAObservations - Per-AP cells of path arrival distances in meters.
%   PathPowerObservations - Matching per-AP powers in linear scale.
%   pAs - 3-by-M AP coordinate matrix in meters.
%   PowerCheckThres - Maximum pairwise linear power ratio.
%   MinANum - Minimum number of distinct AP paths in a valid combination.
%   RandNumPerSeed - Number of randomized clique expansions per seed path.
% Output:
%   FoundLCombs - M-by-K arrival-distance combinations; zero marks an unused AP.
%   FoundPCombs - M-by-K matching linear-power combinations.
%   FoundLCombsIndexs - M-by-K original per-AP path indices.
% Function: Connect paths satisfying pairwise geometric and power constraints,
% grow randomized cliques from each seed, and remove duplicate combinations.

    M = length(LoAObservations);

    % Data Flattening
        % Treat each path as a node; each NodeTable row is [AP_idx, ArrivalDistance, Power, Original_Path_Idx].
    TotalNodeNum = sum(cellfun(@length, LoAObservations));
    NodeTable = zeros(TotalNodeNum, 4);   
    GlobalIdx = 1;
    for m = 1:M
        NumPaths = length(LoAObservations{m});
        if NumPaths > 0
            % Vectorized assignment to accelerate processing
            indices = GlobalIdx : (GlobalIdx + NumPaths - 1);
            NodeTable(indices, 1) = m;                                  % AP_idx
            NodeTable(indices, 2) = LoAObservations{m}(:);              % Delay
            NodeTable(indices, 3) = PathPowerObservations{m}(:);        % Power
            NodeTable(indices, 4) = (1:NumPaths).';                     % Path_idx
            GlobalIdx = GlobalIdx + NumPaths;
        end
    end    
      
    if TotalNodeNum < MinANum
        FoundLCombs = []; FoundPCombs = []; FoundLCombsIndexs = [];
        return;
    end

    % Build Adjacency Matrix
        % AdjMat(i,j) = true indicates that nodes i and j are "compatible"
        % Compatibility condition: Different APs + Close power values + Geometric verification passed
    
    AdjMat = false(TotalNodeNum, TotalNodeNum);
    
    for i = 1:TotalNodeNum
       
        for j = (i+1):TotalNodeNum
            
            % Different APs
            if NodeTable(i, 1) == NodeTable(j, 1)
                continue; 
            end
            
            % Close power values        
            if max([NodeTable(i, 3),NodeTable(j, 3)])/min([NodeTable(i, 3),NodeTable(j, 3)])>PowerCheckThres
                continue; 
            end
            
            % Geometric verification
            if abs(NodeTable(i, 2)-NodeTable(j, 2))<norm(pAs(:, NodeTable(i, 1))-pAs(:, NodeTable(j, 1)))
                AdjMat(i, j) = true;
                AdjMat(j, i) = true;
            end
        end
    end

    % Search path combinations according to adjacency matrix
    FoundCliquesCell = cell(RandNumPerSeed*TotalNodeNum,1);
    FoundCount = 0;
    
    for SeedNode = 1:TotalNodeNum
        
        Candidates = find(AdjMat(SeedNode, :));
        NeighborNum = length(Candidates);

        if NeighborNum < (MinANum - 1)
            continue;
        end
        
        for RandIdx = 1:RandNumPerSeed

            % Preallocate the clique and track its active length.
            CurrentClique = zeros(1, NeighborNum + 1);
            CurrentClique(1) = SeedNode;
            clq_len = 1; % Number of active entries in CurrentClique
            
            % Track the intersection of neighbors shared by the current clique.
            ValidMask = AdjMat(SeedNode, :); 

            if RandIdx == 1 
                RandSeqCandidates = Candidates;
            else
                % Randomize the expansion order after the deterministic first trial.
                RandSeqCandidates = Candidates(randperm(NeighborNum));
            end

            for cand = RandSeqCandidates
                % A valid candidate must be adjacent to every current clique node.
                if ValidMask(cand)
                    
                    % Append through the preallocated active-length index.
                    clq_len = clq_len + 1;
                    CurrentClique(clq_len) = cand;
                    
                    % Update the shared-neighbor mask after accepting a node.
                    ValidMask = ValidMask & AdjMat(cand, :);
                    
                    if clq_len >= MinANum
                        FoundCount = FoundCount + 1;
                        % 仅在最后存储时截取有效部分
                        FoundCliquesCell{FoundCount} = CurrentClique(1:clq_len);            
                    end
                end
            end

        end    
    end
    
    % Remove any possible duplicate path combinations
    CliqueMatrix = zeros(FoundCount,M);
    for FoundIdx = 1:FoundCount 
        CliqueMatrix(FoundIdx,:) = [sort(FoundCliquesCell{FoundIdx}),zeros(1,M-length(FoundCliquesCell{FoundIdx}))];
    end
    [~, UniqueIdxs, ~] = unique(CliqueMatrix, 'rows');
    FoundCount = length(UniqueIdxs);
    FoundCliquesCell = FoundCliquesCell(UniqueIdxs);    

    % Formatting Output
    if FoundCount == 0
        FoundLCombs = []; FoundPCombs = []; FoundLCombsIndexs = [];
        return;
    end
    
    FoundLCombs = zeros(M, FoundCount);
    FoundPCombs = zeros(M, FoundCount);
    FoundLCombsIndexs = zeros(M, FoundCount);

    for k = 1:FoundCount
        CurrentNodesIdx = FoundCliquesCell{k};
        
        % Map the node index back to the path
        for node_id = CurrentNodesIdx
            ap_idx = NodeTable(node_id, 1);
            FoundLCombs(ap_idx,k) = NodeTable(node_id, 2);
            FoundPCombs(ap_idx,k) = NodeTable(node_id, 3);
            FoundLCombsIndexs(ap_idx,k) = NodeTable(node_id, 4);
        end
    end

end

function [LocatedUE,LocatedPoints,HighConfidenceLocatedPoints,APSelectedProbabilities,Wallsxy,CandidateInfo] = ...
    SLAM(LoFObservations,PathPowerObservations,pAs,MinANum,PowerCheckThres,RandNumPerSeed,SolveWallTols,sigma_t,sigma_P,gamma,ConfidenceScheme,z0,UseStage1,beta,GeomScoringMode,RunMode)
%SLAM Run candidate association, source positioning, UE localization, and mapping.
% Input:
%   LoFObservations - Per-AP cells of arrival distances, sorted by path power.
%   PathPowerObservations - Matching per-AP powers in linear scale.
%   pAs - 3-by-M AP coordinate matrix in meters.
%   MinANum - Minimum AP count for a valid path combination.
%   PowerCheckThres - Maximum pairwise linear power ratio.
%   RandNumPerSeed - Number of randomized clique expansions per seed path.
%   SolveWallTols - Symmetry tolerances [x, y, z] in meters.
%   sigma_t - Distance-domain standard deviation of LDoF residuals in meters.
%   sigma_P - Standard deviation of compensated powers in decibels.
%   gamma - Fraction of height-valid candidates retained by confidence.
%   ConfidenceScheme - 1 for the baseline score or 2 for the LLRT score.
%   z0 - Optional known source height; use [] to estimate height.
%   UseStage1 - Flag enabling UE localization before wall positioning.
%   beta - Lower-tail probability used by the Stage 1 chi-square gate.
%   GeomScoringMode - 1/"UEBased" or 2/"Joint".
%   RunMode - "Full", "Stage1Only", or "CandidatesOnly".
% Output:
%   LocatedUE - Final 3-D UE estimate, or [] if localization fails.
%   LocatedPoints - 3-by-K positions of all searched path combinations.
%   HighConfidenceLocatedPoints - Top-gamma height-valid candidate positions.
%   APSelectedProbabilities - Fraction of searched combinations using each AP.
%   Wallsxy - Four wall coordinates ordered as [x+; y+; x-; y-].
%   CandidateInfo - Structure containing candidates, scores, residuals, and
%                   wall-selection diagnostics.
% Function: Stage 1 exhaustively searches strongest/shortest-delay candidates
% for the UE and applies a chi-square gate. Stage 2 performs randomized
% multi-path association and maps walls with a fixed or jointly selected UE.

if nargin < 17 || isempty(RunMode)
    RunMode = "Full";
end
if nargin < 16 || isempty(GeomScoringMode)
    GeomScoringMode = "UEBased";
end
if nargin < 15 || isempty(beta)
    beta = 1.0;
end
if nargin < 14 || isempty(UseStage1)
    UseStage1 = false;
end
if nargin < 12 || isempty(ConfidenceScheme)
    ConfidenceScheme = 2;
end
if ~isnumeric(ConfidenceScheme) || ~isscalar(ConfidenceScheme) || ~ismember(ConfidenceScheme,[1,2])
    error('SLAM:InvalidConfidenceScheme','ConfidenceScheme must be the numeric flag 1 or 2.')
end
if nargin < 10 || isempty(gamma)
    gamma = 0.15;
end
if nargin < 9 || isempty(sigma_P)
    sigma_P = 4;
end
if nargin < 8 || isempty(sigma_t)
    sigma_t = 1;
end
RunMode = NormalizeSLAMRunMode(RunMode);
GeomScoringMode = NormalizeGeomScoringMode(GeomScoringMode);

LocatedUE = [];
LocatedPoints = [];
HighConfidenceLocatedPoints = [];
APSelectedProbabilities = [];
Wallsxy = nan(4,1);
CandidateInfo = InitCandidateInfo();
CandidateInfo.GeomScoringMode = GeomScoringMode;

% Use the 6 paths with the highest power.
LoFObservations = cellfun(@(x) x(1:min(length(x),6)), LoFObservations, 'UniformOutput', false);
PathPowerObservations = cellfun(@(x) x(1:min(length(x),6)), PathPowerObservations, 'UniformOutput', false);

% ----- Stage 1: UE localization using strongest/shortest-delay paths -----
UE_Stage1 = [];
Belief_Stage1 = -inf;
if UseStage1 && RunMode ~= "CandidatesOnly"
    [UE_Stage1, Belief_Stage1, ~] = RunStage1( ...
        LoFObservations, PathPowerObservations, pAs, PowerCheckThres, MinANum, ...
        sigma_t, sigma_P, ConfidenceScheme, z0, beta);
end
CandidateInfo.UE_Stage1 = UE_Stage1;

if RunMode == "Stage1Only"
    LocatedUE = UE_Stage1;
    if isempty(UE_Stage1)
        CandidateInfo.SelectionMode = "Stage1Only_Failed";
    else
        CandidateInfo.SelectionMode = "Stage1Only";
    end
    return
end

% ----- Stage 2: Search path combinations (full multi-path search) -----
[FoundLCombs,FoundPCombs,FoundLCombsIndexs] = PathCombSearch(LoFObservations, PathPowerObservations, pAs, PowerCheckThres, MinANum, RandNumPerSeed);
Num = size(FoundLCombs,2);

CandidateInfo.FoundLCombs = FoundLCombs;
CandidateInfo.FoundPCombs = FoundPCombs;
CandidateInfo.FoundLCombsIndexs = FoundLCombsIndexs;

if Num == 0
    return
end

APSelectedProbabilities = sum(logical(FoundLCombsIndexs),2)/Num;
[V_t,V_P] = CalcuFeasibleRanges(LoFObservations,PathPowerObservations);

% Run TDoA positioning for each path combination obtained through the search.
LocatedPoints = zeros(3,Num);
Beliefs = zeros(Num,1);
GeoConfidences = zeros(Num,1);
PowerConfidences = zeros(Num,1);
GeoMSEs = zeros(Num,1);
WallWeightStats = inf(Num,1);

for i = 1:Num
    CurrentLComb = FoundLCombs(:,i);
    CurrentPComb = FoundPCombs(:,i);
    ActiveAPIdxs = find(CurrentLComb);

    % Stage 2: UseGeoCond=false — do NOT apply B_loc here.
    % Mirror UEs need fair confidence for wall estimation.
    [LocatedPoints(:,i),Beliefs(i),Metrics] = LDoALoca( ...
        pAs(:,ActiveAPIdxs), ...
        CurrentLComb(ActiveAPIdxs), ...
        10*log10(max(CurrentPComb(ActiveAPIdxs),realmin)), ...
        sigma_t,sigma_P,V_t,V_P,ConfidenceScheme,z0,false);

    GeoConfidences(i) = Metrics.GeoConfidence;
    PowerConfidences(i) = Metrics.PowerConfidence;
    GeoMSEs(i) = Metrics.GeoMSE;
    WallWeightStats(i) = GetWallWeightStat(Metrics);
end

CandidateInfo.Beliefs = Beliefs;
CandidateInfo.GeoConfidences = GeoConfidences;
CandidateInfo.PowerConfidences = PowerConfidences;
CandidateInfo.GeoMSEs = GeoMSEs;
CandidateInfo.WallWeightStats = WallWeightStats;
CandidateInfo.V_t = V_t;
CandidateInfo.V_P = V_P;

% Filter the points located between the ceiling plane and the floor plane.
CeilingHeight = pAs(3,1);
ValidFlags = LocatedPoints(3,:)>0 & LocatedPoints(3,:)<CeilingHeight & all(isfinite(LocatedPoints),1);
ValidIndexs = find(ValidFlags);
CandidateInfo.ValidIndexs = ValidIndexs;

if isempty(ValidIndexs)
    return
end

% Sort in descending order of confidence and keep the top gamma proportion.
[~,SortedLocalIdxs] = sort(Beliefs(ValidIndexs),'descend');
SortedValidIndexs = ValidIndexs(SortedLocalIdxs);
TopNum = max(1,ceil(max(0,min(1,gamma))*length(SortedValidIndexs)));
HighConfidenceIndexs = SortedValidIndexs(1:TopNum);
HighConfidenceLocatedPoints = LocatedPoints(:,HighConfidenceIndexs);
CandidateInfo.HighConfidenceIndexs = HighConfidenceIndexs;

if RunMode == "CandidatesOnly"
    CandidateInfo.SelectionMode = "CandidatesOnly";
    return
end

if ~isempty(UE_Stage1)
    LocatedUE = UE_Stage1;
end

if GeomScoringMode == "UEBased"
    if isempty(UE_Stage1)
        CandidateInfo.SelectionMode = "UEBased_Stage1Failed";
        return
    end

    % Mode 1: UEbased GeomScoring. The UE is fixed to the Stage-1 estimate.
    PoolPoints  = [UE_Stage1, LocatedPoints(:, HighConfidenceIndexs)];
    PoolBeliefs = [Belief_Stage1;  Beliefs(HighConfidenceIndexs)];
    PoolWallWeightStats = [inf; WallWeightStats(HighConfidenceIndexs)];
    WallPoolIndexs = [nan; HighConfidenceIndexs(:)];
    [Wallsxy,MirrorSets,DirectionScores,N_walls,WallWeightSums] = LocateWalls( ...
        PoolPoints, 1, pAs, SolveWallTols, PoolBeliefs, PoolWallWeightStats, ~isempty(z0));

    CandidateInfo.SelectionMode = "UEBased";
    CandidateInfo.WallPoolIndexs = WallPoolIndexs;
    CandidateInfo.MirrorSets = MirrorSets;
    CandidateInfo.DirectionScores = DirectionScores;
    CandidateInfo.N_walls = N_walls;
    CandidateInfo.WallWeightSums = WallWeightSums;
    CandidateInfo.FixedUEWallxy = Wallsxy;
    CandidateInfo.FixedUEMirrorSets = MirrorSets;
    CandidateInfo.FixedUEDirectionScores = DirectionScores;
    CandidateInfo.FixedUEWallWeightSums = WallWeightSums;
else
    % Mode 2: Joint UE-and-mirror GeomScoring. Ignore the Stage-1 UE for
    % wall mapping and search the UE together with mirror candidates.
    [SelectedIndex,GeomScores,WallCounts] = SelectUEByGeomScoring( ...
        LocatedPoints(:,HighConfidenceIndexs), ...
        Beliefs(HighConfidenceIndexs), ...
        WallWeightStats(HighConfidenceIndexs), ...
        pAs, ...
        SolveWallTols,~isempty(z0));
    SelectedIndex = HighConfidenceIndexs(SelectedIndex);
    CandidateInfo.SelectionMode = "JointGeomScoring";
    CandidateInfo.GeomFallbackScores = GeomScores;
    CandidateInfo.GeomFallbackWallCounts = WallCounts;

    if isempty(SelectedIndex)
        return
    end

    WallPoolIndexs = unique([SelectedIndex, HighConfidenceIndexs], 'stable');
    UELocalIndex = find(WallPoolIndexs == SelectedIndex,1);
    [Wallsxy,MirrorSets,DirectionScores,N_walls,WallWeightSums] = LocateWalls( ...
        LocatedPoints(:,WallPoolIndexs), ...
        UELocalIndex, ...
        pAs, ...
        SolveWallTols, ...
        Beliefs(WallPoolIndexs), ...
        WallWeightStats(WallPoolIndexs),~isempty(z0));

    CandidateInfo.SelectedIndex = SelectedIndex;
    CandidateInfo.WallPoolIndexs = WallPoolIndexs;
    CandidateInfo.MirrorSets = MirrorSets;
    CandidateInfo.DirectionScores = DirectionScores;
    CandidateInfo.N_walls = N_walls;
    CandidateInfo.WallWeightSums = WallWeightSums;

    if isempty(UE_Stage1) && UseStage1
        LocatedUE = [];
    elseif ~UseStage1
        LocatedUE = LocatedPoints(:,SelectedIndex);
    end

end

end

function RunMode = NormalizeSLAMRunMode(RunMode)
%NORMALIZESLAMRUNMODE Convert supported run-mode aliases to canonical strings.
% Input:
%   RunMode - Scalar string or character vector naming a SLAM run mode.
% Output:
%   RunMode - "Full", "Stage1Only", or "CandidatesOnly".
% Function: Remove spaces, underscores, and hyphens, resolve aliases, and
% reject unsupported modes.

RunMode = string(RunMode);
if ~isscalar(RunMode)
    error('SLAM:InvalidRunMode','RunMode must be a scalar string or character vector.')
end

RunModeKey = lower(regexprep(char(RunMode),'[\s_\-]',''));
switch RunModeKey
    case 'full'
        RunMode = "Full";
    case {'stage1only','stage1'}
        RunMode = "Stage1Only";
    case {'candidatesonly','candidateonly','candidates','candidate'}
        RunMode = "CandidatesOnly";
    otherwise
        error('SLAM:InvalidRunMode','RunMode must be "Full", "Stage1Only", or "CandidatesOnly".')
end

end

function GeomScoringMode = NormalizeGeomScoringMode(GeomScoringMode)
%NORMALIZEGEOMSCORINGMODE Normalize geometric-scoring mode representations.
% Input:
%   GeomScoringMode - Logical, numeric, string, or character mode selector.
% Output:
%   GeomScoringMode - Canonical string "UEBased" or "Joint".
% Function: Map false/1 and supported UE-based aliases to "UEBased", and
% true/2 and supported joint-search aliases to "Joint".

if islogical(GeomScoringMode)
    if GeomScoringMode
        GeomScoringMode = "Joint";
    else
        GeomScoringMode = "UEBased";
    end
    return
end

if isnumeric(GeomScoringMode)
    if isscalar(GeomScoringMode) && GeomScoringMode == 1
        GeomScoringMode = "UEBased";
    elseif isscalar(GeomScoringMode) && GeomScoringMode == 2
        GeomScoringMode = "Joint";
    else
        error('SLAM:InvalidGeomScoringMode','GeomScoringMode must be 1/"UEBased" or 2/"Joint".')
    end
    return
end

GeomScoringMode = string(GeomScoringMode);
if ~isscalar(GeomScoringMode)
    error('SLAM:InvalidGeomScoringMode','GeomScoringMode must be a scalar string or number.')
end

ModeKey = lower(regexprep(char(GeomScoringMode),'[\s_\-]',''));
switch ModeKey
    case {'1','mode1','uebased','uebasedgeomscoring','fixedue','fixeduebased'}
        GeomScoringMode = "UEBased";
    case {'2','mode2','joint','jointgeomscoring','jointueandmirror','jointueandmirrorgeomscoring'}
        GeomScoringMode = "Joint";
    otherwise
        error('SLAM:InvalidGeomScoringMode','GeomScoringMode must be 1/"UEBased" or 2/"Joint".')
end

end

function CandidateInfo = InitCandidateInfo()
%INITCANDIDATEINFO Create an empty SLAM candidate-diagnostics structure.
% Input:
%   None.
% Output:
%   CandidateInfo - Structure with initialized candidate, confidence, wall,
%                   and selection fields.
% Function: Provide consistent default outputs for all early-return paths in SLAM.
CandidateInfo = struct();
CandidateInfo.UE_Stage1 = [];
CandidateInfo.FoundLCombs = [];
CandidateInfo.FoundPCombs = [];
CandidateInfo.FoundLCombsIndexs = [];
CandidateInfo.Beliefs = [];
CandidateInfo.GeoConfidences = [];
CandidateInfo.PowerConfidences = [];
CandidateInfo.GeoMSEs = [];
CandidateInfo.WallWeightStats = [];
CandidateInfo.V_t = [];
CandidateInfo.V_P = [];
CandidateInfo.ValidIndexs = [];
CandidateInfo.HighConfidenceIndexs = [];
CandidateInfo.SelectedIndex = [];
CandidateInfo.WallPoolIndexs = [];
CandidateInfo.MirrorSets = cell(4,1);
CandidateInfo.DirectionScores = zeros(4,1);
CandidateInfo.N_walls = 0;
CandidateInfo.WallWeightSums = zeros(4,1);
CandidateInfo.FixedUEWallxy = nan(4,1);
CandidateInfo.FixedUEMirrorSets = cell(4,1);
CandidateInfo.FixedUEDirectionScores = zeros(4,1);
CandidateInfo.FixedUEWallWeightSums = zeros(4,1);
CandidateInfo.GeomScoringMode = "";
CandidateInfo.SelectionMode = "";
CandidateInfo.GeomFallbackScores = [];
CandidateInfo.GeomFallbackWallCounts = [];
end

function [UE_Stage1, Belief_Stage1, GeoMSE_Stage1] = RunStage1(LoFObservations, PathPowerObservations, pAs, PowerCheckThres, MinANum, sigma_t, sigma_P, ConfidenceScheme, z0, beta)
%RUNSTAGE1 Localize the UE from strongest-power and shortest-delay candidates.
% Input:
%   LoFObservations - Per-AP cells of retained arrival distances in meters.
%   PathPowerObservations - Matching powers in linear scale, sorted descending.
%   pAs - 3-by-M AP coordinate matrix in meters.
%   PowerCheckThres - Maximum pairwise linear power ratio.
%   MinANum - Minimum AP count for a valid UE path combination.
%   sigma_t - Distance-domain standard deviation of LDoF residuals in meters.
%   sigma_P - Standard deviation of compensated powers in decibels.
%   ConfidenceScheme - 1 for the baseline score or 2 for the LLRT score.
%   z0 - Optional known UE height; use [] to estimate height.
%   beta - Lower-tail probability used by the residual chi-square gate.
% Output:
%   UE_Stage1 - Selected 3-D UE estimate, or [] if no candidate passes.
%   Belief_Stage1 - Confidence of the selected candidate, or -inf on failure.
%   GeoMSE_Stage1 - Geometric residual MSE of the selected candidate.
% Function: Exhaustively associate the strongest and minimum-delay candidate
% from each AP, reject invalid-height and chi-square-failing estimates, and
% select the remaining candidate with the highest geometry-corrected confidence.

if nargin < 10 || isempty(beta)
    beta = 1.0;
end

UE_Stage1    = [];
Belief_Stage1  = -inf;
GeoMSE_Stage1  = inf;

% Extract two Stage-1 candidates from each AP:
%   1) the strongest-power path (first entry after the top-power sorting);
%   2) the minimum-delay path among the retained paths.
% If they are the same path, keep it only once.
M_AP = numel(LoFObservations);
LoFObs_S1 = cell(M_AP, 1);
PwrObs_S1 = cell(M_AP, 1);
for ap = 1:M_AP
    loas = LoFObservations{ap};
    pwrs = PathPowerObservations{ap};
    if isempty(loas)
        LoFObs_S1{ap} = [];
        PwrObs_S1{ap} = [];
    else
        [~, minIdx] = min(loas);
        candidateIdxs = unique([1, minIdx], 'stable');
        LoFObs_S1{ap} = loas(candidateIdxs);
        PwrObs_S1{ap} = pwrs(candidateIdxs);
    end
end

% Exhaustive LPC search: enumerate all valid AP subsets.
[FoundLCombs_S1, FoundPCombs_S1] = Stage1PathCombSearch(LoFObs_S1, PwrObs_S1, pAs, PowerCheckThres, MinANum);
Num_S1 = size(FoundLCombs_S1, 2);
if Num_S1 == 0
    return
end

[V_t,V_P] = CalcuFeasibleRanges(LoFObservations,PathPowerObservations);
CeilingHeight = pAs(3,1);

LocatedPoints_S1    = zeros(3, Num_S1);
Beliefs_S1          = -inf(Num_S1, 1);
GeoMSEs_S1          = inf(Num_S1, 1);
GeoChiSquares_S1    = inf(Num_S1, 1);
ResidualDoFs_S1     = ones(Num_S1, 1);

for i = 1:Num_S1
    CurrentLComb = FoundLCombs_S1(:,i);
    CurrentPComb = FoundPCombs_S1(:,i);
    ActiveAPIdxs = find(CurrentLComb);
    ResidualDoFs_S1(i) = max(numel(ActiveAPIdxs) - 4, 1);
    % Stage 1: UseGeoCond=true so B_loc is included in confidence.
    [LocatedPoints_S1(:,i), Beliefs_S1(i), Metrics_S1] = LDoALoca( ...
        pAs(:,ActiveAPIdxs), ...
        CurrentLComb(ActiveAPIdxs), ...
        10*log10(max(CurrentPComb(ActiveAPIdxs), realmin)), ...
        sigma_t, sigma_P, V_t, V_P, ConfidenceScheme, z0, true);
    GeoMSEs_S1(i) = Metrics_S1.GeoMSE;
    if ConfidenceScheme == 2
        GeoChiSquares_S1(i) = Metrics_S1.GeoChiSquare;
        ResidualDoFs_S1(i) = Metrics_S1.ResidualDoF;
    else
        GeoChiSquares_S1(i) = Metrics_S1.GeoMSE * max(numel(ActiveAPIdxs) - 1, 1) / max(sigma_t^2, eps);
    end
end

% Keep only finite points between the floor and ceiling planes.
ValidFlags_S1 = LocatedPoints_S1(3,:) > 0 & LocatedPoints_S1(3,:) < CeilingHeight & all(isfinite(LocatedPoints_S1), 1);
ValidIdxs_S1  = find(ValidFlags_S1);
if isempty(ValidIdxs_S1)
    return
end

% Apply the residual chi-square gate with residual DoF d-3.
GateThresholds_S1 = arrayfun(@(nu) ChiSquareInv(beta, nu), ResidualDoFs_S1(ValidIdxs_S1));
PassGate_S1 = GeoChiSquares_S1(ValidIdxs_S1) < GateThresholds_S1;
FilteredIdxs_S1 = ValidIdxs_S1(PassGate_S1);
if isempty(FilteredIdxs_S1)
    return
end
% FilteredIdxs_S1 = ValidIdxs_S1;

% Among candidates that pass the gate, output the one with the
% highest total confidence (includes B_loc).
[~, BestLocalIdx] = max(Beliefs_S1(FilteredIdxs_S1));
BestIdx_S1    = FilteredIdxs_S1(BestLocalIdx);
UE_Stage1     = LocatedPoints_S1(:, BestIdx_S1);
Belief_Stage1  = Beliefs_S1(BestIdx_S1);
GeoMSE_Stage1  = GeoMSEs_S1(BestIdx_S1);

end

function x = ChiSquareInv(probability, dof)
%CHISQUAREINV Compute a lower-tail chi-square quantile without a toolbox call.
% Input:
%   probability - Requested cumulative probability.
%   dof - Chi-square degrees of freedom.
% Output:
%   x - Lower-tail chi-square quantile.
% Function: Clip invalid endpoint values and evaluate the quantile through
% the inverse regularized lower incomplete gamma function.

probability = min(max(probability, eps), 1 - eps);
dof = max(dof, eps);
x = 2 * gammaincinv(probability, dof/2, 'lower');

end

function [FoundLCombs, FoundPCombs] = Stage1PathCombSearch(LoFObs_S1, PwrObs_S1, pAs, PowerCheckThres, MinANum)
%STAGE1PATHCOMBSEARCH Exhaustively enumerate valid UE path combinations.
% Input:
%   LoFObs_S1 - Per-AP cells containing strongest/minimum-delay distances.
%   PwrObs_S1 - Matching per-AP powers in linear scale.
%   pAs - 3-by-M AP coordinate matrix in meters.
%   PowerCheckThres - Maximum pairwise linear power ratio.
%   MinANum - Minimum number of APs in a valid combination.
% Output:
%   FoundLCombs - M-by-K arrival-distance combinations; zero marks an unused AP.
%   FoundPCombs - M-by-K matching linear-power combinations.
% Function: Enumerate every eligible AP subset and within-AP path choice,
% retain combinations satisfying all pairwise power and geometric constraints,
% and remove duplicates.

M = length(LoFObs_S1);
FoundLCombs = [];
FoundPCombs = [];

% Identify APs that have a path measurement.
AvailMask = ~cellfun(@isempty, LoFObs_S1);
AvailAPIdxs = find(AvailMask);
AvailNum = numel(AvailAPIdxs);

if AvailNum < MinANum
    return
end

LoA_avail = LoFObs_S1(AvailAPIdxs);
Pwr_avail = PwrObs_S1(AvailAPIdxs);

% Enumerate all AP subsets and all within-AP candidate choices.
for SubsetSize = MinANum:AvailNum
    SubsetIdxsList = nchoosek(1:AvailNum, SubsetSize);  % C(AvailNum,SubsetSize) x SubsetSize
    for s = 1:size(SubsetIdxsList, 1)
        SubLocal = SubsetIdxsList(s, :);
        CandidateCounts = cellfun(@numel, LoA_avail(SubLocal));
        ChoiceCells = arrayfun(@(n) 1:n, CandidateCounts, 'UniformOutput', false);
        ChoiceGrids = cell(1, SubsetSize);
        [ChoiceGrids{:}] = ndgrid(ChoiceCells{:});
        ChoiceMatrix = zeros(numel(ChoiceGrids{1}), SubsetSize);
        for ChoiceDim = 1:SubsetSize
            ChoiceMatrix(:, ChoiceDim) = ChoiceGrids{ChoiceDim}(:);
        end

        for ChoiceIdx = 1:size(ChoiceMatrix, 1)
            CurrentLoAs = zeros(SubsetSize, 1);
            CurrentPwrs = zeros(SubsetSize, 1);
            for ii = 1:SubsetSize
                CurrentLoAs(ii) = LoA_avail{SubLocal(ii)}(ChoiceMatrix(ChoiceIdx, ii));
                CurrentPwrs(ii) = Pwr_avail{SubLocal(ii)}(ChoiceMatrix(ChoiceIdx, ii));
            end

            IsValid = true;
            for ii = 1:SubsetSize
                for jj = (ii+1):SubsetSize
                    PwrRatio = max(CurrentPwrs(ii), CurrentPwrs(jj)) / max(min(CurrentPwrs(ii), CurrentPwrs(jj)), realmin);
                    if PwrRatio > PowerCheckThres
                        IsValid = false;
                        break
                    end

                    ai = AvailAPIdxs(SubLocal(ii));
                    aj = AvailAPIdxs(SubLocal(jj));
                    if abs(CurrentLoAs(ii) - CurrentLoAs(jj)) >= norm(pAs(:,ai) - pAs(:,aj))
                        IsValid = false;
                        break
                    end
                end
                if ~IsValid, break, end
            end
            if ~IsValid, continue, end

            % Build output columns (sparse M-vector format, matching PathCombSearch output).
            LComb = zeros(M, 1);
            PComb = zeros(M, 1);
            GlobalAPIdxs = AvailAPIdxs(SubLocal);
            LComb(GlobalAPIdxs) = CurrentLoAs;
            PComb(GlobalAPIdxs) = CurrentPwrs;
            FoundLCombs = [FoundLCombs, LComb]; %#ok<AGROW>
            FoundPCombs = [FoundPCombs, PComb]; %#ok<AGROW>
        end
    end
end

if ~isempty(FoundLCombs)
    [~, UniqueIdxs] = unique([FoundLCombs; FoundPCombs].', 'rows', 'stable');
    FoundLCombs = FoundLCombs(:, UniqueIdxs);
    FoundPCombs = FoundPCombs(:, UniqueIdxs);
end
end

function [V_t,V_P] = CalcuFeasibleRanges(LoFObservations,PathPowerObservations)
%CALCUFEASIBLERANGES Estimate empirical supports for the clutter likelihood.
% Input:
%   LoFObservations - Per-AP cells of path arrival distances in meters.
%   PathPowerObservations - Matching path powers in linear scale.
% Output:
%   V_t - Max-minus-min arrival-distance support length in meters.
%   V_P - Max-minus-min received-power support length in decibels.
% Function: Pool all retained paths across APs and calculate positive support
% lengths, using 1 when the corresponding observation set is empty.
AllLoFs = [];
AllPowers_dB = [];
for APidx = 1:length(LoFObservations)
    if ~isempty(LoFObservations{APidx})
        AllLoFs = [AllLoFs; LoFObservations{APidx}(:)];
    end
    if ~isempty(PathPowerObservations{APidx})
        AllPowers_dB = [AllPowers_dB; 10*log10(max(PathPowerObservations{APidx}(:),realmin))];
    end
end

if isempty(AllLoFs)
    V_t = 1;
else
    V_t = max(AllLoFs)-min(AllLoFs);
end

if isempty(AllPowers_dB)
    V_P = 1;
else
    V_P = max(AllPowers_dB)-min(AllPowers_dB);
end

V_t = max(V_t,eps);
V_P = max(V_P,eps);
end

function [SelectedLocalIndex,GeomScores,WallCounts] = SelectUEByGeomScoring(Points,Confidences,WeightStats,pAs,Tol,OnlyXYEst,DirectionMask)
%SELECTUEBYGEOMSCORING Select a UE hypothesis using mirror-geometry consensus.
% Input:
%   Points - 3-by-N candidate source positions.
%   Confidences - N candidate confidence scores.
%   WeightStats - N positive residual statistics for wall-coordinate weighting.
%   pAs - 3-by-M AP coordinate matrix in meters.
%   Tol - Symmetry tolerances [x, y, z] in meters.
%   OnlyXYEst - If true, omit the mirror-height consistency check.
%   DirectionMask - Optional enabled-wall mask or list of direction indices.
% Output:
%   SelectedLocalIndex - Local index of the highest-scoring UE hypothesis.
%   GeomScores - UE confidence plus best mirror confidence in each direction.
%   WallCounts - Number of supported wall directions for each UE hypothesis.
% Function: Treat each candidate as the UE, identify geometrically consistent
% mirror candidates, and select the hypothesis with maximum consensus score.
CandidateNum = size(Points,2);
GeomScores = -inf(CandidateNum,1);
WallCounts = zeros(CandidateNum,1);
SelectedLocalIndex = [];

if nargin < 7 || isempty(DirectionMask)
    DirectionMask = true(4,1);
else
    DirectionMask = NormalizeWallDirectionMask(DirectionMask);
end
if CandidateNum == 0 || ~any(DirectionMask)
    return
end

for CandidateIdx = 1:CandidateNum
    [~,~,DirectionScores,N_walls] = LocateWalls(Points,CandidateIdx,pAs,Tol,Confidences,WeightStats,OnlyXYEst,DirectionMask);
    GeomScores(CandidateIdx) = Confidences(CandidateIdx) + sum(DirectionScores(DirectionMask));
    WallCounts(CandidateIdx) = N_walls;
end

if any(isfinite(GeomScores))
    [~,SelectedLocalIndex] = max(GeomScores);
end
end

function DirectionMask = NormalizeWallDirectionMask(DirectionMask)
%NORMALIZEWALLDIRECTIONMASK Convert wall-direction selections to a logical mask.
% Input:
%   DirectionMask - Four-element binary mask or vector of direction indices.
% Output:
%   DirectionMask - 4-by-1 logical mask ordered as [x+; y+; x-; y-].
% Function: Validate direction indices and normalize supported mask formats.

if islogical(DirectionMask) || (isnumeric(DirectionMask) && numel(DirectionMask) == 4 && all(ismember(DirectionMask(:),[0,1])))
    DirectionMask = logical(DirectionMask(:));
else
    DirectionIndexs = DirectionMask(:);
    DirectionMask = false(4,1);
    DirectionIndexs = DirectionIndexs(isfinite(DirectionIndexs) & DirectionIndexs == round(DirectionIndexs) & DirectionIndexs >= 1 & DirectionIndexs <= 4);
    DirectionMask(DirectionIndexs) = true;
end

end

function [Wallxy, MirrorSets, DirectionScores, N_walls, DirectionWeightSums] = LocateWalls(Points, UEIndex, pAs, Tol, Confidences, WeightStats, OnlyXYEst, DirectionMask)
%LOCATEWALLS Identify mirror candidates and estimate wall coordinates.
% Input:
%   Points - 3-by-N candidate source positions.
%   UEIndex - Index of the candidate treated as the UE.
%   pAs - 3-by-M AP coordinate matrix in meters.
%   Tol - Symmetry tolerances [x, y, z] in meters.
%   Confidences - N candidate confidence scores; defaults to zeros.
%   WeightStats - N positive residual statistics; defaults to ones.
%   OnlyXYEst - If true, omit the mirror-height consistency check.
%   DirectionMask - Optional enabled-wall mask or direction-index list.
% Output:
%   Wallxy - Estimated wall coordinates ordered as [x+; y+; x-; y-].
%   MirrorSets - Local candidate indices assigned to each wall direction.
%   DirectionScores - Maximum mirror confidence in each direction.
%   N_walls - Number of directions with positive mirror support.
%   DirectionWeightSums - Sum of midpoint-fusion weights per direction.
% Function: Apply tangential-coordinate, height, direction, and outward-midpoint
% symmetry gates, then fuse accepted UE/mirror midpoints for each wall.

N = size(Points, 2);
Wallxy = nan(4,1);
MirrorSets = cell(4,1);
DirectionScores = zeros(4,1);
DirectionWeightSums = zeros(4,1);
N_walls = 0;

if N == 0 || isempty(UEIndex) || isnan(UEIndex) || UEIndex < 1 || UEIndex > N
    return
end
if nargin < 5 || isempty(Confidences)
    Confidences = zeros(1,N);
end
if nargin < 6 || isempty(WeightStats)
    WeightStats = ones(1,N);
end
if nargin < 8 || isempty(DirectionMask)
    DirectionMask = true(4,1);
else
    DirectionMask = NormalizeWallDirectionMask(DirectionMask);
end

Confidences = Confidences(:).';
WeightStats = WeightStats(:).';

UE = Points(:,UEIndex);
xMax = max(pAs(1,:));
xMin = min(pAs(1,:));
yMax = max(pAs(2,:));
yMin = min(pAs(2,:));

for i = 1:N
    if i == UEIndex
        continue;
    end
    if ~isfinite(Confidences(i)) || Confidences(i) <= 0
        continue;
    end

    P = Points(:, i);
    Mid = (P + UE) / 2;
    if OnlyXYEst
        SameZ = 1;
    else
        SameZ = abs(P(3) - UE(3)) < Tol(3);
    end

    if SameZ && abs(P(2) - UE(2)) < Tol(2)
        if DirectionMask(1) && P(1) > UE(1) && Mid(1) > xMax - Tol(1)
            MirrorSets{1}(end+1) = i; % x+
        elseif DirectionMask(3) && P(1) < UE(1) && Mid(1) < xMin + Tol(1)
            MirrorSets{3}(end+1) = i; % x-
        end
    end

    if SameZ && abs(P(1) - UE(1)) < Tol(1)
        if DirectionMask(2) && P(2) > UE(2) && Mid(2) > yMax - Tol(2)
            MirrorSets{2}(end+1) = i; % y+
        elseif DirectionMask(4) && P(2) < UE(2) && Mid(2) < yMin + Tol(2)
            MirrorSets{4}(end+1) = i; % y-
        end
    end
end

for DirectionIdx = 1:4
    if ~DirectionMask(DirectionIdx)
        continue;
    end
    CurrentSet = MirrorSets{DirectionIdx};
    if isempty(CurrentSet)
        DirectionScores(DirectionIdx) = 0;
        continue;
    end

    DirectionScores(DirectionIdx) = max(Confidences(CurrentSet));
    [Wallxy(DirectionIdx),DirectionWeightSums(DirectionIdx)] = FuseWallCoordinate(Points, UEIndex, CurrentSet, DirectionIdx, WeightStats);
end

N_walls = nnz(DirectionScores > 0);

end

function [WallCoordinate,WeightSum] = FuseWallCoordinate(Points, UEIndex, MirrorSet, DirectionIdx, WeightStats)
%FUSEWALLCOORDINATE Fuse midpoint observations for one wall direction.
% Input:
%   Points - 3-by-N candidate source positions.
%   UEIndex - Index of the UE hypothesis.
%   MirrorSet - Candidate indices assigned to the wall direction.
%   DirectionIdx - Direction index: 1=x+, 2=y+, 3=x-, or 4=y-.
%   WeightStats - Positive residual statistics for inverse weighting.
% Output:
%   WallCoordinate - Weighted UE/mirror midpoint coordinate.
%   WeightSum - Sum of the weights used in the fusion.
% Function: Weight midpoint observations by inverse residual statistics and
% fall back to uniform weights if all inverse weights are invalid.

UE = Points(:,UEIndex);
MidPoints = (Points(:,MirrorSet) + UE) / 2;
if DirectionIdx == 1 || DirectionIdx == 3
    Observations = MidPoints(1,:);
else
    Observations = MidPoints(2,:);
end

Weights = 1 ./ max(WeightStats(MirrorSet), eps);
Weights(~isfinite(Weights)) = 0;
if sum(Weights) <= 0
    Weights = ones(size(Observations));
end

WeightSum = sum(Weights);
WallCoordinate = sum(Weights .* Observations) / WeightSum;

end

function PlotDualAxisPerformance(X,GroundTruthRate,GroundTruthErr,Scheme2Rate,Scheme2Err,Scheme2ParamSets,TitleLines,RateLabel,ErrLabel,SavePathNoExt)
%PLOTDUALAXISPERFORMANCE Plot and save success-rate and error curves.
% Input:
%   X - Shared horizontal-axis sample values.
%   GroundTruthRate - Ground-truth-association success rates.
%   GroundTruthErr - Ground-truth-association errors.
%   Scheme2Rate - Proposed-method success-rate columns for each configuration.
%   Scheme2Err - Proposed-method error columns for each configuration.
%   Scheme2ParamSets - Configuration structures used to build legend labels.
%   TitleLines - Figure title text.
%   RateLabel - Right-axis label.
%   ErrLabel - Left-axis label.
%   SavePathNoExt - Output path without a file extension.
% Output:
%   None.
% Function: Create a dual-y-axis comparison plot and save FIG and PNG files.
Fig = figure('Color','w');
hold on
Colors = lines(size(Scheme2Rate,2)+1);

yyaxis right
ylim([0,1])
ylabel(RateLabel)
hGT = plot(X,GroundTruthRate,'--o','DisplayName','GT-Association SLAM, success rate','Color',Colors(1,:),'LineWidth',1.2,'MarkerSize',5);
Handles = hGT;
for ConfigIdx = 1:size(Scheme2Rate,2)
    Handles(end+1) = plot(X,Scheme2Rate(:,ConfigIdx),'--o', ...
        'DisplayName',[ConfigLabel(Scheme2ParamSets(ConfigIdx)),', success rate'], ...
        'Color',Colors(ConfigIdx+1,:), ...
        'LineWidth',1.2, ...
        'MarkerSize',5);
end

yyaxis left
ylim([0, inf]);
ylabel(ErrLabel)
Handles(end+1) = plot(X,GroundTruthErr,'-s','DisplayName','GT-Association SLAM, error','Color',Colors(1,:),'LineWidth',1.2,'MarkerSize',5);
for ConfigIdx = 1:size(Scheme2Err,2)
    Handles(end+1) = plot(X,Scheme2Err(:,ConfigIdx),'-s', ...
        'DisplayName',[ConfigLabel(Scheme2ParamSets(ConfigIdx)),', error'], ...
        'Color',Colors(ConfigIdx+1,:), ...
        'LineWidth',1.2, ...
        'MarkerSize',5);
end

xlabel('LoS-path Availability Probability')
grid on
title(TitleLines,'Interpreter','none')
legend(Handles,'Location','best')
savefig(Fig,[SavePathNoExt,'.fig'])
exportgraphics(Fig,[SavePathNoExt,'.png'],'Resolution',300)
end

function Text = ConfigLabel(Config)
%CONFIGLABEL Format a parameter structure as a plot legend label.
% Input:
%   Config - Structure containing beta and optional tolerance, gamma, and
%            geometric-scoring mode fields.
% Output:
%   Text - Formatted configuration label.
% Function: Append available configuration values in a compact display string.
Text = sprintf('\\beta=%g', Config.beta);
if isfield(Config,'SolveWallTols') && ~isempty(Config.SolveWallTols)
    Text = sprintf('%s, tol=[%g %g %g]', Text, Config.SolveWallTols(1), Config.SolveWallTols(2), Config.SolveWallTols(3));
end
if isfield(Config,'gamma') && ~isempty(Config.gamma)
    Text = sprintf('%s, \\gamma=%g', Text, Config.gamma);
end
if isfield(Config,'GeomScoringMode') && ~isempty(Config.GeomScoringMode)
    Text = sprintf('%s, GS=%s', Text, char(NormalizeGeomScoringMode(Config.GeomScoringMode)));
elseif isfield(Config,'UseGeomScoringFallback') && ~isempty(Config.UseGeomScoringFallback)
    Text = sprintf('%s, GS=%s', Text, char(NormalizeGeomScoringMode(Config.UseGeomScoringFallback)));
end
end

function Value = SafeDivide(Numerator,Denominator)
%SAFEDIVIDE Divide two scalars while handling a zero denominator.
% Input:
%   Numerator - Scalar numerator.
%   Denominator - Scalar denominator.
% Output:
%   Value - Numerator/Denominator, or NaN when Denominator is zero.
% Function: Avoid divide-by-zero warnings in aggregated performance metrics.
if Denominator == 0
    Value = NaN;
else
    Value = Numerator/Denominator;
end
end
