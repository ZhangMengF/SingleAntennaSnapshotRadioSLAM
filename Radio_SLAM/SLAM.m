function [LocatedUE,LocatedPoints,HighBeliefLocatedPoints,APSelectedProbabilities,Wallsxy] = ...
    SLAM(LoFObservations,PathPowerObservations,pAs,MinANum,PowerCheckThres,RandNumPerSeed,SolveWallTols)
    % Input:
        % LoAObservations - The set of all arrival times (converted into distances) of the paths on the APs
        % PathPowerObservations - The set of all powers of the paths on the APs
        % pAs - AP xyz-coordinates (each column represents one AP coordinate)
        % MinANum - The minimum number of paths included in a valid path combination
        % PowerCheckThres - Power check threshold       
        % RandNumPerSeed - Control the total number of search attempts
        % SolveWallTols - The allowable position error between the UE and the mirrored UE when determining the position of the walls
    % Output:
        % LocatedUE - Estimation of the xyz coordinates of UE
        % LocatedPoints - The xyz coordinates of the TDoA positioning of the searched path combinations through geometric and power consistency
        % HighBeliefLocatedPoints - The part of the points in 'LocatedPoints' with the highest belief
        % APSelectedProbabilities - Probabilities of AP selection of the searched path combinations through geometric and power consistency
        % Wallsxy - Estimation of wall positions
    % Function:
        % Locate a UE and surrounding walls, using the multipath arrival times and powers extracted from the UE's pilot on the APs

    % Use the 6 paths with the highest power.
    LoFObservations = cellfun(@(x) x(1:min(length(x),6)), LoFObservations, 'UniformOutput', false);          
    PathPowerObservations = cellfun(@(x) x(1:min(length(x),6)), PathPowerObservations, 'UniformOutput', false);
    
    % Serach path combinations
    [FoundLCombs,FoundPCombs,FoundLCombsIndexs] = PathCombSearch(LoFObservations, PathPowerObservations, pAs, PowerCheckThres, MinANum, RandNumPerSeed);    
    Num = size(FoundLCombs,2);
    APSelectedProbabilities = sum(logical(FoundLCombsIndexs),2)/Num;
    
    % Run TDoA positioning for each path combination obtained through the search.
    LocatedPoints = zeros(3,Num);
    Beliefs = zeros(Num,1);
    for i = 1:Num
        CurrentLComb = FoundLCombs(:,i);
        CurrentPComb = FoundPCombs(:,i);
        ActiveAPIdxs = find(CurrentLComb);

        [LocatedPoints(:,i),Beliefs(i)] ...
            = LDoALoca(pAs(:,ActiveAPIdxs),CurrentLComb(ActiveAPIdxs),10*log10(CurrentPComb(ActiveAPIdxs)));
    end

    % Filter the points located between the ceiling plane and the floor plane
    CeilingHeight = pAs(3,1);
    BaseSpaceFilteredIndexs = find(arrayfun(@(i) LocatedPoints(3,i)>0 && LocatedPoints(3,i)<CeilingHeight, 1:Num));
    
    % Sort in descending order of belief
    [~,SortedIdxs] = sort(Beliefs(BaseSpaceFilteredIndexs),'descend');
    
    % Output the points with the 15% highest belief
    FinedIndexs = BaseSpaceFilteredIndexs(SortedIdxs(1:floor(0.15*length(SortedIdxs))));
    HighBeliefLocatedPoints = LocatedPoints(:,FinedIndexs);

    % Obtain the UE positioning result
    LocatedUEIdx = find(arrayfun(@(i) all(FoundLCombsIndexs(:,FinedIndexs(i))<2),1:length(FinedIndexs)),1);
    LocatedUE = LocatedPoints(:, FinedIndexs(LocatedUEIdx));
    
    % Obtain the wall positioning result
    Wallsxy = LocateWalls(LocatedPoints(:, FinedIndexs),LocatedUEIdx,[max(pAs(1,:)),max(pAs(2,:))],SolveWallTols);

end