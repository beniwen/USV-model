clc;
clear all;
close all;


% 定义初始参数
initial_params = [-20000, -350000, -350000, -350000]'; % 初始参数值
% 定义模型结构
FileName = 'vehicle_model'; % 模型文件
Order = [8 4 8]; % 输出、输入和状态向量的维数 [ny nu nx]
InitialStates = zeros(8, 1);                % 系统状态的初始值
Ts = 0; % 连续时间系统

% 生成一个系统辨识对象
nlgr = idnlgrey(FileName, Order, initial_params,InitialStates, Ts, ...
    'Name', 'Extended Vehicle Model', ...
    'InputName', {'tau_x', 'tau_y', 'tau_p', 'tau_r'}', ...
    'InputUnit', {'N', 'N', 'Nm', 'Nm'}', ...
    'OutputName', {'u', 'v', 'p', 'r', 'phi', 'psi', 'x', 'y'}', ...
    'OutputUnit', {'m/s', 'm/s', 'rad/s', 'rad/s', 'rad', 'rad', 'm', 'm'}', ...
    'TimeUnit', 's');

% 设置系统状态向量的名称和单位            
nlgr = setinit(nlgr, 'Name', {'u', 'v', 'p', 'r', 'phi', 'psi', 'x', 'y'}');
nlgr = setinit(nlgr, 'Unit', {'m/s', 'm/s', 'rad/s', 'rad/s', 'rad', 'rad', 'm', 'm'}');

% 设置待辨识参数的名称
param_names = {'Xudot', 'Yvdot', 'Kpdot', 'Nrdot'};
nlgr = setpar(nlgr, 'Name', param_names);
for i = 2:4
    nlgr.Parameters(i).Fixed = false; % 将参数设置为可变参数
    nlgr.Parameters(i).Minimum = -8e5; % 设置参数的下界
    nlgr.Parameters(i).Maximum = -3.5e5;  % 设置参数的上界
end
% % Xudot = -17400.0; Yvdot = -393000; Kpdot = -774000; Nrdot = -38.7e4; % 附加质量
    nlgr.Parameters(1).Fixed = false;
    nlgr.Parameters(1).Minimum = -2e4; % 设置参数的下界
    nlgr.Parameters(1).Maximum = -1e4;  % 设置参数的上界

% 打印出nlgr的信息
disp(nlgr);

% 导入数据
dataset = load('vehicle_example.mat');
z = iddata(dataset.Output, dataset.Input, 0.01, 'Name', 'Extended Vehicle Data'); 
z.InputName = {'tau_x', 'tau_y', 'tau_p', 'tau_r'}; 
z.InputUnit = {'N', 'N', 'Nm', 'Nm'};
z.OutputName = {'u', 'v', 'p', 'r', 'phi', 'psi', 'x', 'y'}; 
z.OutputUnit = {'m/s', 'm/s', 'rad/s', 'rad', 'rad', 'rad', 'm', 'm'};
z.Tstart = 0; 
z.TimeUnit = 's'; 

% 打印出数据集的信息
disp(z);
% 处理缺失值
z = misdata(z);
% 设置初始状态
X0init = zeros(8,1);
% 使用遗传算法进行参数辨识
options = optimoptions('ga', 'Display', 'iter', 'MaxGenerations', 100, 'PopulationSize', 50);

% 定义适应度函数
fitnessFunction = @(params) evaluate_model(params, nlgr, dataset.Input, dataset.Output);

% 运行遗传算法
optimal_params = ga(fitnessFunction, length(initial_params), [], [], [], [], ...
    [-2e4, -8e5, -8e5, -8e5], [-1e4, -3.5e5, -3.5e5, -3.5e5], [], options);

% 输出优化结果
disp('优化后的参数:');
disp(optimal_params);

function fitness = evaluate_model(params, nlgr, Input, Output)
    % 更新模型参数
    nlgr.Parameters(1).Value = params(1); % Xudot
    nlgr.Parameters(2).Value = params(2); % Yvdot
    nlgr.Parameters(3).Value = params(3); % Kpdot
    nlgr.Parameters(4).Value = params(4); % Nrdot
    
    % 创建 iddata 对象
    data = iddata(Output, Input);
    
    % 使用 idnlgrey 模型进行模拟
    sim_out = sim(nlgr, data); % 使用 iddata 对象
    
    % 计算适应度（例如，均方误差）
    fitness = norm(sim_out.OutputData - Output) ^ 2;
end
