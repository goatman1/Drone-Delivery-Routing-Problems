function [S,ttlE]=Metropolis(S1,S2,Svel,Distance,Demand,Travelcon,Capacity,T,K,Pr,timeBreakpoints)
%% 输入
% S1        当前解
% S2        新解
% Svel      random ground speed (between two nodes)
% Distance  Distance matrix
% Demand        各点需求量
% Travelcon     行程约束
% Capacity      容量约束
% T         当前温度
% K          scenarios number
% Pr         scenarios probability
% timeBreakpoints   opearation time period
%% 输出
% S         下一个当前解
% ttlDis	下一个当前解的路线距离

%%
ttlE1 = Evaluation(S1,Svel,Distance,Demand,Travelcon,Capacity,K,Pr,timeBreakpoints);  %计算路线expected energy consumption
ttlE2 = Evaluation(S2,Svel,Distance,Demand,Travelcon,Capacity,K,Pr,timeBreakpoints);  %计算路线energy consumption
dC = ttlE2 - ttlE1;   %计算路线enery difference

if dC < 0       %如果能力降低 接受新路线
    S = S2;
    ttlE = ttlE2;
elseif exp(-dC/T) >= rand   %以exp(-dC/T)概率接受新路线
    S = S2;
    ttlE = ttlE2;
else  %不接受新路线
    S = S1;
    ttlE = ttlE1;
end