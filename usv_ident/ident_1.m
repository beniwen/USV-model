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
% nlgr = setinit(nlgr, 'Value', num2cell(X0init)); 
% % 对比系统真实输出和模型输出（参数辨识前）
% figure(1)
% compare(z, nlgr, [], compareOptions('InitialCondition', X0init));

% 定义PSO参数
nParticles = 30;
nIterations = 100;
paramDim = length(initial_params);  % 参数维度

% 粒子群初始化
particles = rand(nParticles, paramDim) .* (maxBound - minBound) + minBound;  % 在参数范围内初始化
velocities = zeros(nParticles, paramDim);
pBest = particles;  % 个体最优
gBest = particles(1, :);  % 全局最优
gBestScore = inf;

% 优化循环
for iter = 1:nIterations
    for i = 1:nParticles
        % 计算当前粒子的目标函数值
        score = objective_function(particles(i, :), nlgr, z);
        
        % 更新个体最优
        if score < objective_function(pBest(i, :), nlgr, z)
            pBest(i, :) = particles(i, :);
        end
        
        % 更新全局最优
        if score < gBestScore
            gBest = particles(i, :);
            gBestScore = score;
        end
    end
    
    % 更新粒子速度和位置
    for i = 1:nParticles
        velocities(i, :) = 0.5 * velocities(i, :) + rand() * (pBest(i, :) - particles(i, :)) ...
                           + rand() * (gBest - particles(i, :));
        particles(i, :) = particles(i, :) + velocities(i, :);
    end
    
    % 打印当前迭代的最优值
    disp(['Iteration ' num2str(iter) ', Best Error: ' num2str(gBestScore)]);
end

% 将最佳参数更新到模型
nlgr = setpar(nlgr, 'Value', num2cell(gBest));

% 对比系统真实输出和模型输出
figure;
compare(z, nlgr);
function error = objective_function(params, nlgr, z)
    % 设定模型的当前参数
    nlgr = setpar(nlgr, 'Value', num2cell(params));
    
    % 使用模型生成输出
    [~, y_sim] = sim(nlgr, z.SamplingInstants);
    
    % 计算模型输出与真实数据之间的误差
    y_real = z.OutputData;
    
    % 计算误差 (可以使用 MSE 或其他衡量指标)
    error = sum((y_sim(:) - y_real(:)).^2);  % 二次误差 (MSE)
end
