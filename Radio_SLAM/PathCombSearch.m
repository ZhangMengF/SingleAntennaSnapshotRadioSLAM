function [FoundLCombs,FoundPCombs,FoundLCombsIndexs] = ...
    PathCombSearch(LoAObservations, PathPowerObservations, pAs, PowerCheckThres, MinANum, RandNumPerSeed)
    
    % Input:
        % LoAObservations - The set of all arrival times (converted into distances) of the paths on the APs
        % PathPowerObservations - The set of all powers of the paths on the APs
        % pAs - AP xyz-coordinates (each column represents one AP coordinate)
        % PowerCheckThres - Power check threshold
        % MinANum - The minimum number of paths included in a valid path combination
        % RandNumPerSeed - Control the total number of search attempts
    % Output:
        % FoundLCombs - The arrival times (converted into distances) of the searched path combinations
        % FoundPCombs - The powers of the searched path combinations
        % FoundLCombsIndexs - The order (from largest to smallest power) of the paths in the searched path combination on the AP
    % Function: Based on the arrival times and power values of the paths on the APs, search for path combinations that satisfy the power and geometric consistency constraints.

    M = length(LoAObservations);

    % Data Flattening
        % Treat each path on each AP as a node and establish a unified index, each row in NodeTable: [AP_idx, Delay, Power, Original_Path_Idx]    
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

            CurrentClique = SeedNode;
            if RandIdx == 1 
                RandSeqCandidates = Candidates;
            else
                RandSeqCandidates = Candidates(randperm(NeighborNum));
            end

            for cand = RandSeqCandidates
                % Check whether 'cand' is adjacent to all the nodes in 'CurrentClique'.
                if all(AdjMat(cand, CurrentClique))           
                    CurrentClique = [CurrentClique, cand];
                    if length(CurrentClique) >= MinANum
                        FoundCount = FoundCount + 1;
                        FoundCliquesCell{FoundCount} = CurrentClique;            
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