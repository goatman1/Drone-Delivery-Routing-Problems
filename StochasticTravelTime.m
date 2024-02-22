function travelTime = StochasticTravelTime(i, j, t0, s, Svel, Distance, timeBreakpoints)
    if i ~= j
        t = t0;
        d = Distance(i,j);
        tPrime = t + (d / piecewise_Svel(i, j, t, s, Svel))/60;
    
        % Find the time breakpoint after which t falls
        k = find(timeBreakpoints > t, 1);
    
        % If k is empty, t is after the last breakpoint
        if isempty(k)
            k = length(timeBreakpoints);
        end

        tk = timeBreakpoints(k);
    
        % Go directly to the time breakpoint before tPrime
        if k > 0 && tPrime ~= Inf
            while tPrime > tk
                dCovered = piecewise_Svel(i, j, t, s, Svel) * 60 * (tk - t);
                d = d - dCovered;
                t = tk;
                tPrime = t + (d / piecewise_Svel(i, j, t, s, Svel))/60;
                k = k+1;
                if k <= length(timeBreakpoints)
                    tk = timeBreakpoints(k);
                else 
                    fprintf('pair (%s) at time %s will be traveling out of opeartion period in scenario %s \n;', num2str([i,j]), num2str(t0),num2str(s))
                    break
                end
            end
        end
    
        % Calculate the remaining time from the last breakpoint to the destination
        travelTime = tPrime - t0;
    else
        travelTime = 0;
    end
end
