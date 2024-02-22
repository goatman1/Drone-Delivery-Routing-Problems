function TextOutput(Svel,Distance,Demand,route,Capacity,S,Pr,timeBreakpoints)
%% 输出路径函数
%输入：route 路径
%输出：p 路径文本形式

%% 总路径
len=length(route); %路径长度
disp('Best Route:')

p=num2str(route(1)); %配送中心位先进入路径首位
for i=2:len
    p=[p,' -> ',num2str(route(i))]; %路径依次加入下一个经过的点
end
disp(p)




%% 子路径

route=route+1; %路径值全体+1，为方便下面用向量索引

Vnum=1; %
EnergyTraveled=zeros(S,1);  % 汽车已经行驶的距离
delivery=0;       % 汽车已经送货量，即已经到达点的需求量之和
subpath='0'; %子路径路线
unit = 0.2; % electricity rate 0.2$/KJ
ttlE=0; %此方案所有车辆的总能量
t = zeros(S,length(Demand)); % initial time row: scenario, col: #vehicle
s = 1:S; % scenario vector
t_s = ones(S,1); % fixed time cost for delivery action 1 min

for j=2:len
    traveltime = arrayfun(@(x,y)StochasticTravelTime(route(j-1),route(j),x,y,Svel,Distance,timeBreakpoints),t(:,Vnum),s'); 
    EnergyTraveled = EnergyTraveled+Energy(Capacity - delivery)*traveltime; %每两点间距离累加
    delivery = delivery+Demand(route(j)); %累加可配送量
    subpath=[subpath,' -> ',num2str(route(j)-1)]; %子路径路线输出
    t(:,Vnum) = t(:,Vnum) + traveltime + t_s; % cumulative travel time + service time
	if route(j)==1 %若此位是配送中心
        disp('-------------------------------------------------------------')
        fprintf('Route of Vehichle No.%d: %s  \n',Vnum,subpath)%输出：每辆车 路径 
        fprintf('Total Expected Travel Time: %.2f mins, Energy traveled: %.2f kW*s, Expected Energy Cost: %.2f $, load rate: %.2f%%;  \n',mean(t(:,Vnum)), Pr*EnergyTraveled,Pr*EnergyTraveled*(60^-1)*unit,delivery/Capacity*100)%输出：expected energy cost 满载率
        ttlE=ttlE + Pr*EnergyTraveled*(60^-1); 
        Vnum=Vnum+1; %车辆数累加
        EnergyTraveled=zeros(S,1); %已行驶距离置零
        delivery=0; %已配送置零
        subpath='0'; %子路径重置
	end
end
fprintf('Total Energy Consumption: %.2f kW*h',ttlE)
