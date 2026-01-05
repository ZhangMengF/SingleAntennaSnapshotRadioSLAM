function [omega_now,x_now] = OMP(A,y,delta,MaxIterNum)

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