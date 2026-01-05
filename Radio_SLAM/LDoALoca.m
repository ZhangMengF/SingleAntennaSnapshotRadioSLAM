function [p_esti,Belief] = LDoALoca(pAs,LoAs,PathPowers_dB)
    
    % Input:
        % pAs - AP xyz-coordinates (each column represents one AP coordinate)
        % LoAs - The arrival time (converted into distance) of each AP's path
        % PathPowers_dB - Power (dB) of each AP's path
    % Output:
        % p_esti - xyz-coordinates of the positioned point
        % Belief - positioning belief
    % Function: Run TDoA positioning

    M = size(pAs,2);
    pAsMassCenter = mean(pAs,2);
    [~,RefA_id] = min(arrayfun(@(i) norm(pAs(:,i)-pAsMassCenter),1:M));% 找在中间位置的AP做为参考AP

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

    % Calculate the coordinates of the positioned point
    if z_estimation(3)<norm(pA_ref(1:2)-z_estimation(1:2))
        p_esti = [z_estimation(1:2);1];
    else
        zaxis_esti = pA_ref(3)-sqrt(z_estimation(3)^2-norm(pA_ref(1:2)-z_estimation(1:2))^2);
        p_esti = [z_estimation(1:2);zaxis_esti];
    end
            
    % Calculate positioning belief
    LDoFs_esti = (vecnorm(p_esti-pAs_notRef)-norm(p_esti-pA_ref)).';
    LDoFOffsets = LDoFs_esti-LDoFs;
    Belief_Geo = -mean(abs(LDoFOffsets));

    if isempty(PathPowers_dB)
        Belief = Belief_Geo;
    else
        d_est_vec = vecnorm(p_esti - pAs);
        CompensatedPowers = PathPowers_dB(:)' + 20*log10(d_est_vec);
        VarCompensatedPower = var(CompensatedPowers);           
        Belief = Belief_Geo*VarCompensatedPower;
    end
    
end