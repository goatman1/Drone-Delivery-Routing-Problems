clear;clc;

% load('City2.mat')	      %需求点经纬度，用于画实际路径的XY坐标
load('Nodes.mat')
load('Distance.mat')	  %距离矩阵
%% Generate random travel time matrix
City = nodes;
CityNum=size(City,1)-1;
NodeNum=size(City,1);
S = 100; % # of scenarios
p = 1.0/S; % Prob of a scenario 
u = 20; %drone air speed (m/s)
tau = 180; %launh/land time (s)
m_w = 7; % Speed Mean
sd_w = 2; % Speed STD
m_theta = pi/3; % Direction Mean
sd_theta = pi/3; % Direction STD
w = m_w + sd_w*randn(S, 1); %wind speed 
theta = m_theta + sd_theta*randn(S, 1); %wind direction


% Generate random direction and distance matrix
% Distance_core = randi([1e3,1e4],[NodeNum, NodeNum]); % generate random integers in the interval (1000, 10000)
Distance_core = dists;
phi_core = rand(NodeNum, NodeNum) * 2*pi; % generate random values in the interval (0, 2π)
% create a symmetric matrix with the random values and set the diagonal elements to zero
Distance = tril(Distance_core,-1) + tril(Distance_core,-1)';
phi = tril(phi_core,-1) + tril(phi_core,-1)';

% Initialize t to be an empty array
STime = zeros(NodeNum, NodeNum, S);

% Loop over each index in the 3-dimensional array t
for i = 1:NodeNum
    for j = 1:NodeNum
        for s = 1:S
            % Calculate the value of t at index (i,j,s) using the provided equation
            if i ~= j
               numerator = Distance(i,j);
               denominator = sqrt(u^2 - (w(s)*sin(phi(i,j)-theta(s)))^2) + w(s)*cos(phi(i,j)-theta(s));
               STime(i,j,s) = numerator / denominator + tau;
            else 
               STime(i,j,s) = 0;
            end
        end
    end
end

% Initialize t to be an empty array
STime2py = cell(S,1);

% Loop over each index in the 3-dimensional array t
for i = 1:NodeNum
    for j = 1:NodeNum
        for s = 1:S
            % Calculate the value of t at index (i,j,s) using the provided equation
            if i ~= j
               numerator = Distance(i,j);
               denominator = sqrt(u^2 - (w(s)*sin(phi(i,j)-theta(s)))^2) + w(s)*cos(phi(i,j)-theta(s));
               STime2py{s}(i,j) = numerator / denominator + tau;
            else 
               STime2py{s}(i,j) = 0;
            end
        end
    end
end
save('STime2.mat', 'STime','STime2py');
