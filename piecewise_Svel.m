function speed = piecewise_Svel(i, j, t, s, Svel)
    % Convert time to time index
    t_index = max(ceil(t / 10), 1); % Assuming each time step represents 10 minutes
    
    % Extract the corresponding Svel value
%     speed = Svel(i, j, t_index, s);
    speed = max(Svel(i, j, t_index, s),0);
    
    if ~isreal(speed)
        speed = 0;
    end
        
end
