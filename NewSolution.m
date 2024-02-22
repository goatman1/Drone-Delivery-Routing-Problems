function S2 = NewSolution(S1, rule)
    Length = length(S1); % Get the number of elements in S1
    S2 = S1; % Initialize S2 as a copy of S1

    % Generate two distinct random indices between 2 and Length-1
    i = randi([2, Length-1]);
    j = randi([2, Length-1]);
    % Ensure i is not equal to j
    while i == j
        j = randi([2, Length-1]);
    end
    % Ensure i is less than j
    if i > j
        temp = i;
        i = j;
        j = temp;
    end
    
    % Apply the chosen exchange rule
    switch rule
        case 'simple' % Simple exchange (swap two elements)
            S2([i, j]) = S2([j, i]);
            
        case '2-opt' % 2-opt swap (reverse a subsequence)
            S2(i:j) = S2(j:-1:i);
            
        case 'cross' % Cross exchange (swap two subsequences)
            k = min(j+randi([1, Length-j]), Length); % Index for the second part of cross
            if k - j > i - 1 % Ensure subsequences are non-empty and non-overlapping
                S2([i:j, (j+1):k]) = S2([(j+1):k, i:j]); % Swap subsequences
            end
            
        otherwise
            error('Unknown rule specified');
    end
end
