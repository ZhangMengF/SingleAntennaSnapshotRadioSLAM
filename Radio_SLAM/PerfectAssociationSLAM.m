function [LocatedUE,LocatedUE_Optimized,Wallxy,Wallxy_Optimized,LocatedUEandMirror_valid,LocatedUEandMirror_Optimized_valid] = ...
            PerfectAssociationSLAM(LoFObservations,TargetPathIndexs,MinANum,pAs,DoOptimizationFlag)
    % Input:
        % LoFObservations - The set of all arrival times (converted into distances) of the paths on the APs
        % TargetPathIndexs - The arrival time indices of the paths corresponding to the UE and its images across the APs
        % MinANum - The minimum number of arrival times required for TDoA positioning to operate
        % pAs - AP xyz-coordinates (each column represents one AP coordinate)
        % DoOptimizationFlag - The flag to denote whether to optimize AP selection
    % Output:
        % LocatedUE - Estimation of the xyz-coordinates of the UE
        % LocatedUE_Optimized - Estimation of the xyz-coordinates of the UE, optimized AP selection
        % Wallxy - Estimation of the position of the walls 
        % Wallxy_Optimized - Estimation of the position of the walls, optimized AP selection
        % LocatedUEandMirror_valid - Estimation of the xyz-coordinates of the UE and UE mirrors
        % LocatedUEandMirror_Optimized_valid - Estimation of the xyz-coordinates of the UE and UE mirrors, optimized AP selection
    % Function: The indices of the paths corresponding to the UE and its images are known at each AP, namely perfect association. Directly execute TDoA positioning to estimate the positions of the UE and walls.

    % Use the arrival times of the paths extracted from all APs to run TDoA positioning.
    LocatedUEandMirror = zeros(3,6);
    SuccessFlags = false(6,1);
    for TargetIdx = 1:6
        ActiveAPIdxs = find(TargetPathIndexs(:,TargetIdx));
        if length(ActiveAPIdxs) >= MinANum
            SuccessFlags(TargetIdx) = true;
            LoFs = arrayfun(@(i) LoFObservations{i}(TargetPathIndexs(i,TargetIdx)), ActiveAPIdxs);
            [LocatedUEandMirror(:,TargetIdx),~] ...
                = LDoALoca(pAs(:,ActiveAPIdxs),LoFs,[]);       
        end

    end
    
    [LocatedUE,Wallxy] = SolveWallxy(LocatedUEandMirror,SuccessFlags);
    LocatedUEandMirror_valid = LocatedUEandMirror(:,SuccessFlags);

    LocatedUEandMirror_Optimized = zeros(3,6);
    if DoOptimizationFlag
        % Select APs (≥MinANum), run TDoA positioning on their paths' arrival times, and output the highest-belief result
        for TargetIdx = 1:6
    
            ActiveAPIdxs = find(TargetPathIndexs(:,TargetIdx));
            ActiveAPNum = length(ActiveAPIdxs);
    
            if ActiveAPNum <= MinANum
                LocatedUEandMirror_Optimized(:,TargetIdx) = LocatedUEandMirror(:,TargetIdx);
            else
                CombNum = sum(arrayfun(@(i) nchoosek(ActiveAPNum,i), MinANum:ActiveAPNum));
                LocatedPs = zeros(3,CombNum);
                Beliefs = zeros(CombNum,1);
                PCounter = 1;
                for UsedAPNum = MinANum:ActiveAPNum
        
                    LocalAPIdxCombs = nchoosek(1:ActiveAPNum,UsedAPNum);
        
                    for LocalAPIdxs = LocalAPIdxCombs.'
                        SelectedActiveAPIdxs = ActiveAPIdxs(LocalAPIdxs);
                        LoFs = arrayfun(@(i) LoFObservations{i}(TargetPathIndexs(i,TargetIdx)), SelectedActiveAPIdxs);
                        [LocatedPs(:,PCounter),Beliefs(PCounter)] ...
                            = LDoALoca(pAs(:,SelectedActiveAPIdxs),LoFs,[]);
                        PCounter = PCounter+1;
                    end
        
                end
                [~,MaxIdx] = max(Beliefs);
                LocatedUEandMirror_Optimized(:,TargetIdx) = LocatedPs(:,MaxIdx);
            end
    
        end    
    end
    [LocatedUE_Optimized,Wallxy_Optimized] = SolveWallxy(LocatedUEandMirror_Optimized,SuccessFlags);
    LocatedUEandMirror_Optimized_valid = LocatedUEandMirror_Optimized(:,SuccessFlags);

end

function [LocatedUE,Wallxy] = SolveWallxy(LocatedUEandMirror,SuccessFlags)
    
    SolvedWallxy = nan(2,2);
    if SuccessFlags(1)
        LocatedUE = LocatedUEandMirror(:,1);
        for TargetIdx = 2:5
            if SuccessFlags(TargetIdx)
                if TargetIdx==2
                    SolvedWallxy(1,1) = (LocatedUEandMirror(1,1)+LocatedUEandMirror(1,TargetIdx))/2;
                elseif TargetIdx==3
                    SolvedWallxy(2,1) = (LocatedUEandMirror(2,1)+LocatedUEandMirror(2,TargetIdx))/2;
                elseif TargetIdx==4
                    SolvedWallxy(1,2) = (LocatedUEandMirror(1,1)+LocatedUEandMirror(1,TargetIdx))/2;
                elseif TargetIdx==5
                    SolvedWallxy(2,2) = (LocatedUEandMirror(2,1)+LocatedUEandMirror(2,TargetIdx))/2;
                end
            end
        end
    else
        LocatedUE = [];
    end    
    Wallxy = SolvedWallxy(:);

end