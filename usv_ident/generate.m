clc;
clear all;
close all;

% 仿真时间和采样时间
period = 0.01; % 采样时间
T = 100; % 仿真总时间

% 无人车扩展模型的初始状态
u = 0; v = 0; p = 0; r = 0; phi = 0; psi = 0; x = 0; y = 0;
state = [u; v; p; r; phi; psi; x; y];

 % 假设实际系统的参数（用于生成数据）
    m = 350*1025; % 实际质量
    Xudot = -17400.0; Yvdot = -393000; 
    Kpdot = -774000; Nrdot = -38.7e4; % 附加质量

% 数据保存
Input = [];
Output = [];

% 控制输入（例如正弦输入）
for t = 0:period:T
    tau = [500000; 0; 0; 50000000]; % 控制输入力和力矩
    
    % 系统动力学模型迭代
    dstate = vehicle_model(t, state, tau,  Xudot, Yvdot, Kpdot, Nrdot);
    state = state + dstate * period;
    
    % 保存数据
    Input = [Input; tau'];
    Output = [Output; state'];
end

% 保存数据
save('vehicle_example.mat', 'Input', 'Output');


