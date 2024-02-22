% function mindisever=Main_SA_SVRP()
clc 
clear 
close all
tic % current time
rng(0); % Sets the seed to 0
%% SA Algorithm for TDSVRP
%Input：
%City           City Coordinates
%Demand        
%Travelcon      Energy Limit
%Capacity       Capacity Limit
%T0             Initial Temp
%Tend           End Temp
%L              Chain Length 
%q              Rate of Temp Loss
%S              Scenarios Number


%Output：
%bestroute      minimal energy routes
%mindisever     minimal energy
%% Load Data
% load('City.mat')	      %(x,y) coordinates from Google Map
% load('Demand.mat')       % Demand
% load('STime.mat')       % Stochastic travel time
% load('Capacity.mat')     % Capacity Limit
load('STD_Stag_Rand_H80_C50_S100_Z1_para_9am_1pm.mat') % Parameters load


%% Initialize Model Parameters
City = nodes';
CityNum=size(City,1)-1;    % Customer Num
S = size(Svel,4); % # of scenarios
s = 1:S; %scenarios vector
Pr = 1.0/S*ones(1,S); % Prob of each scenario 
Travelcon = 0.222*60; % Battery Capacity (kW*min) = 22.2V*10000mAh/1e6*60
Capacity = 5; % Load Capacity (lb)
timeBreakpoints = 0:10:60; % time range
%% Initialize Algo Hyper Parameters
T0=1000;        % Initial Temp
Tend=1e-3;      % End Temp
L=200;          % Chain Length
q=0.95;          % Rate of Temp Loss
count=0;        % 
rule = 'cross'; % random exchange rule: simple; 2-opt; cross

%% Initial Solution (Must feasible)
TSProute=[0,randperm(CityNum)]+1; % Initialize Random TSP Route
S1=ones(1,CityNum*2+1); %
    
EnergyTraveled=zeros(S,1);  % used energy 
delivery=0; % delivered freights
t = zeros(S,CityNum); % initial time row: scenario, col: #vehicle
num_v = 1; % initial vehicle number
k=1;
for j=2:CityNum+1
	k=k+1;   % 
    EnergyNow = Energy(Capacity - delivery)*arrayfun(@(x,y)StochasticTravelTime(S1(k-1),TSProute(j),x,y,Svel,Distance,timeBreakpoints),t(:,num_v),s');
    EnergyNext = Energy(Capacity - delivery - Demand(TSProute(j)))*arrayfun(@(x,y)StochasticTravelTime(TSProute(j),1,x,y,Svel,Distance,timeBreakpoints),t(:,num_v),s');
    
	if min(EnergyTraveled+EnergyNow+EnergyNext) > Travelcon || delivery+Demand(TSProute(j)) > Capacity
        num_v = num_v + 1;
        S1(k)=1; % Back to depot  
        % Send a new vehicle
        EnergyTraveled=Energy(Capacity - delivery)*arrayfun(@(x,y)StochasticTravelTime(1,TSProute(j),x,y,Svel,Distance,timeBreakpoints),t(:,num_v),s'); % Initial Energy = Depot to Current Node
        delivery=Demand(TSProute(j)); % Initialize Capacity
        k=k+1;
        t(:,num_v) = t(:,num_v) + arrayfun(@(x,y)StochasticTravelTime(1,TSProute(j),x,y,Svel,Distance,timeBreakpoints),t(:,num_v),s');
        S1(k)=TSProute(j);  % Add Node to Route
	else %
        EnergyTraveled=EnergyTraveled+EnergyNow;
        delivery=delivery+Demand(TSProute(j)); 
        t(:,num_v) = t(:,num_v) + arrayfun(@(x,y)StochasticTravelTime(S1(k-1),TSProute(j),x,y,Svel,Distance,timeBreakpoints),t(:,num_v),s');
        S1(k)=TSProute(j); 
	end
end

%% Compute #iterations (Annealing Time), i.e., solve T0 * q^x = Tend
Time=ceil(double(log(Tend/T0)/log(q)));

bestind=zeros(1,CityNum*2+1);      %每代的最优路线矩阵初始化
BestObjByIter=zeros(Time,1);       %每代目标值矩阵初始化

%% Iteration
while T0 > Tend
    count = count+1;     %更新迭代次数
    Population = zeros(L,CityNum*2+1); %为此温度下迭代个体矩阵分配内存
    ObjByIter = zeros(L,1); %为此温度下迭代个体的目标函数值矩阵分配内存
    for k = 1:L
        %% 产生新解 (optimizable）
        S2 = NewSolution(S1, rule);
        %% Metropolis法则判断是否接受新解
        [S1,ttlE] = Metropolis(S1,S2,Svel,Distance,Demand,Travelcon,Capacity,T0,S,Pr,timeBreakpoints);  % Metropolis 抽样算法
        ObjByIter(k) = ttlE;    %此温度下每迭代一次就存储一次目标函数值
        Population(k,:) = S1;          %此温度下每迭代一次就存储一次此代最优个体
    end
    
    %% 记录每次迭代过程的最优路线
    [d0,index] = min(ObjByIter); %取出此温度下所有迭代中最优的一次
    if count == 1 || d0 < BestObjByIter(count-1) %若为第一次迭代或上次迭代比这次更满意
        BestObjByIter(count) = d0;            %如果当前温度下最优路程小于上一路程则记录当前路程
        bestind = Population(index,:);  %记录当前温度的最优路线
    else
        BestObjByIter(count) = BestObjByIter(count-1);  %如果当前温度下最优路程大于上一路程则记录上一路程的目标函数值
    end
    
    T0 = q * T0; %降温
     %% 显示此代信息
    fprintf('Iteration = %d, Min Cost = %.2f $  \n',count,BestObjByIter(count)) %输出当前迭代信息
end

%% 找出历史最短距离和对应路径
mindisever = BestObjByIter(count); %找出历史最优目标函数值
bestroute=bestind; % 取最优个体

%删去路径中多余1
for i=1:length(bestroute)-1
    if bestroute(i)==bestroute(i+1)
        bestroute(i)=0;  %相邻位都为1时前一个置零
    end
end
bestroute(bestroute==0)=[];  %删去多余零元素

bestroute=bestroute-1;  % 编码各减1，与文中的编码一致

%% 计算结果数据输出到命令行
disp('-------------------------------------------------------------')
toc %显示运行时间
fprintf('Total Energy = %s kW*s \n',num2str(mindisever))
fprintf('Total Cost = %s $ \n',num2str(mindisever))
TextOutput(Svel,Distance,Demand,bestroute,Capacity,S,Pr,timeBreakpoints)  %显示最优路径
disp('-------------------------------------------------------------')
% Save results
save('SVRP_routes.mat','bestroute','mindisever')

%% 优化过程迭代图
figure
plot(BestObjByIter,'LineWidth',2)
xlim([1 count]) %设置 x 坐标轴范围
set(gca, 'LineWidth',1)
xlabel('Iterations')
ylabel('Min Energy (kW*s)')
title('SA Optimization Process')

%% 绘制最优解的实际路线
% Add wind visual
load('Wind_Field.mat')
%%
zeta = 1; % speed up/down factor
zoom_factor = 10; % factor to increase dimensions
zoom_matrix = ones(zoom_factor); % zoom matrix
time_range = 79 - 55; % 9 am to 1 pm
zoomed_magnitude = zeros(size(magnitude, 1) * zoom_factor, size(magnitude, 2) * zoom_factor, time_range);
zoomed_direction = zeros(size(magnitude, 1) * zoom_factor, size(magnitude, 2) * zoom_factor, time_range);
zoomed_u_data = zeros(size(magnitude, 1) * zoom_factor, size(magnitude, 2) * zoom_factor, time_range);
zoomed_v_data = zeros(size(magnitude, 1) * zoom_factor, size(magnitude, 2) * zoom_factor, time_range);
for t = 1:time_range
    % Scale up the magnitude matrix using kron
    zoomed_magnitude(:, :, t) = zeta * kron(zoom_matrix, magnitude(:, :, t));
    zoomed_direction(:, :, t) = kron(zoom_matrix, direction(:, :, t));
    zoomed_u_data(:, :, t) = zeta * kron(zoom_matrix, u_data(:, :, t));
    zoomed_v_data(:, :, t) = zeta * kron(zoom_matrix, v_data(:, :, t));   
end

%%
% for t = 1:size(magnitude, 3)
%     % Scale up the magnitude matrix using kron
%     zoomed_u_data(:, :, t) = zeta * kron(zoom_matrix, u_data(:, :, t));
%     zoomed_v_data(:, :, t) = zeta * kron(zoom_matrix, v_data(:, :, t));    
% end
l_sw = size(zoomed_magnitude,1);
l_cs = size(zoomed_magnitude,2);
% Plot the points
[x, y] = meshgrid(1:l_cs, 1:l_sw);
figure(2)
% Wind Field
imagesc('CData',zoomed_magnitude(:,:,1))
colorbar
hold on 
quiver(x(1:150:end,1:150:end), y(1:150:end,1:150:end), zoomed_v_data(1:150:end,1:150:end,1), zoomed_u_data(1:150:end,1:150:end,1),'k','AutoScaleFactor', 0.5);
hold on
%%
figure(2)
DrawPath(bestroute,City)
