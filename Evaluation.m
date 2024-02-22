function ttlE=Evaluation(route,Svel,Distance,Demand,Travelcon,Capacity,S,Pr,timeBreakpoints)
%% 计算各个体的能量消耗 适应度函数  
% 输入：
% route         种群矩阵
% Svel      random ground speed
% Distance  Distance matrix
% Demand        各点需求量
% Travelcon     行程约束
% Capacity      容量约束
% S             scenarios number
% Pr          scenarios probability
% timeBreakpoints   opearation time period

% 输出：
% ttlE	此个体路径的能量

len=size(route,2);
%相关数据初始化
EnTraveled=zeros(S,1);  % 汽车已经行驶的能量
delivery=0; % 汽车已经送货量，即已经到达点的需求量之和置零
ttlE=0; %此方案所有车辆的总能量
unit = 0.2; % Electricity Price
C = 500; % Drone Cost
t = zeros(S,length(Demand)); % initial time row: scenario, col: #vehicle
num_v = 1; % initial drone number
s = 1:S; % scenario vector
t_s = ones(S,1); % fixed time cost for delivery action 1 min

for j=2:len
    traveltime = arrayfun(@(x,y)StochasticTravelTime(route(j-1),route(j),x,y,Svel,Distance,timeBreakpoints),t(:,num_v),s'); % Compute travel time in different scenarios
    EnTraveled = EnTraveled + Energy(Capacity - delivery)*traveltime; % Compute energy in different scenarios
    delivery = delivery+Demand(route(j)); % cumulative delivered payloads
    t(:,num_v) = t(:,num_v) + traveltime + t_s; % cumulative travel time + service time
    if min(EnTraveled) > Travelcon || delivery > Capacity
        ttlE = Inf;  % penalize infeasible solution
        break
    end
    if route(j)==1 && Pr*traveltime ~= 0 %if it is depot
        num_v = num_v + 1; % add new drone
        ttlE=ttlE + Pr*EnTraveled*(60^-1)*unit + C; % total cost = Previous Energy Cost + Drone Cost
        EnTraveled=zeros(S,1); %已行驶距离置零
        delivery=0; %已配送置零
    elseif route(j)==1
        ttlE=ttlE + Pr*EnTraveled*(60^-1)*unit; % total cost = Previous Energy Cost + Drone Cost
        EnTraveled=zeros(S,1); %已行驶距离置零
        delivery=0; %已配送置零
    end
end
