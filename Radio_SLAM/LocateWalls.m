function Wallxy = LocateWalls(Points, UEIndex, xyBounds, Tol)

    N = size(Points, 2);
    
    if isempty(UEIndex)
        Wallxy = nan(4,1);
        return
    else
        UE = Points(:,UEIndex);
    end

    x_cuts = []; y_cuts = [];
    
    for i = 1:N
        if i == UEIndex
            continue;
        end
        P = Points(:, i);        
        Mid = (P + UE) / 2;
            
        if abs(P(2) - UE(2)) < Tol(2) && abs(P(3) - UE(3)) < Tol(3) && abs(Mid(1))>xyBounds(1)
            x_cuts(end+1) = Mid(1);
        elseif abs(P(1) - UE(1)) < Tol(1) && abs(P(3) - UE(3)) < Tol(3) && abs(Mid(2))>xyBounds(2)
            y_cuts(end+1) = Mid(2);
        end
    end
  

    xywalls = nan(2,2);
    xywalls(1,1) = mean(x_cuts(x_cuts>0));
    xywalls(1,2) = mean(x_cuts(x_cuts<0));
    xywalls(2,1) = mean(y_cuts(y_cuts>0));
    xywalls(2,2) = mean(y_cuts(y_cuts<0));
    
    Wallxy = xywalls(:);
end