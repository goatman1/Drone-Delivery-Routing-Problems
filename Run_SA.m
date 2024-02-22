%% Multiple run SA
% Initialize variables
num_runs = 10; % Number of runs
results = zeros(num_runs, 1); % Initialize array to store results

% Run the algorithm multiple times
for i = 1:num_runs
    % Call your algorithm or function here and store the result
    result = Main_SA_SVRP(); % Replace 'your_algorithm()' with your actual algorithm or function
    
    % Store the result in the array
    results(i) = result;
end

% Compute mean and standard deviation
mean_result = mean(results);
std_result = std(results);

% Display the mean and standard deviation
fprintf('Mean result: %.2f\n', mean_result);
fprintf('Standard deviation: %.2f\n', std_result);
