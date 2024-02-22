function consump=Energy(Demand)
% input: load (lb) from i to j 

% output: energy consumption from i to j

% Linear Approx Coefficients
alpha = 0.217; % (kg/kW)
beta = 0.185; % (kW)
q = 1.5; % self-weight (kg)

consump = alpha*(Demand/2.205 + q) + beta; %kW

% Non Linear Real Consump