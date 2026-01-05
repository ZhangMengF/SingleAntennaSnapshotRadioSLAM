function [SpaceBound,pAs,UEPositions,LoFObservations_UEs,PathGainObservations_UEs,TrueTargetPathIndexs_UEs] = LoadRaytracedPaths()
    
    % Output: 
        % SpaceBound - The true values of the positions of the four walls
        % pAs - AP xyz-coordinates (each column represents one AP coordinate)
        % UEPositions - The dataset of xyz-coordinates of all UE position samples
        % LoFObservations_UEs - The dataset of all arrival times (converted into distances) of the paths on the APs
        % PathGainObservations_UEs - The dataset of all gains of the paths on the APs
        % TrueTargetPathIndexs_UEs - The dataset of arrival time indices of the paths corresponding to the UE and its images across the APs
    % Function: Load Sionna raytraced multipath data

    RaytracedPaths = load('RaytracedPaths.mat').RaytracedPaths; 

    SpaceBound = (RaytracedPaths{1}.SpaceBound).';
    pAs = RaytracedPaths{1}.APPositions;
    UENum = length(RaytracedPaths)-1;
    UEPositions = zeros(3,UENum);        
    LoFObservations_UEs = cell(UENum,1);
    PathGainObservations_UEs = cell(UENum,1);
    
    APNum = length(RaytracedPaths{2}.CIR_by_AP);
    TrueTargetPathIndexs_UEs = cell(UENum,1);
    for u = 1:UENum
        % Extract the current UE's data
        ThisUE = RaytracedPaths{u+1};
        UEPositions(:,u) = ThisUE.UE_Pos;
        
        LoFs_by_AP = cell(APNum,1);
        PathGains_by_AP = cell(APNum,1);
        TrueTargetPathIndexs = zeros(APNum,6);
        for m = 1:APNum
            % Obtain the data matrix of the m-th AP
            % Each row: [Delay, Gain_real, Gain_Imag, Type1, Type2, Object1, Object2, Point1, Point2]
            PathData = ThisUE.CIR_by_AP{m};            
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
                        % The Object_ids for the four walls and the floor in Sionna are 4, 5, 6, 7, and 8, respectively.
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

end