%% 微电网两阶段鲁棒优化经济调度方法(成功)
%  start date: 2024.4.7
%  completion date: 2024.4.10
clc;
clear;
yalmip('clear')

% 迭代次数
itr = 10;

%% 参数
% 微型燃气轮机
    P_Gmin = 80;
    P_Gmax = 800;
    a = 0.67;

% 储能
    P_Smax = 500;
    E_Smax = 1800;
    E_Smin = 400;
    E_S0 = 1000;
    K_S = 0.38;
    yita = 0.95;
% 需求响应负荷
    K_DR = 0.32;
    D_DR = 2940;
    DR_max = 200;
    DR_min = 50;
    P_DRstar = ...
    [   8.2e+001
        7.1e+001
        6.1e+001
        5.1e+001
        7.1e+001
        7.3e+001
        9.1e+001
        1.02e+002
        1.22e+002
        1.54e+002
        1.70e+002
        2.00e+002
        1.39e+002
        1.03e+002
        1.01e+002
        1.22e+002
        1.40e+002
        1.50e+002
        1.90e+002
        2.00e+002
        2.00e+002
        1.91e+002
        1.01e+002
        8.2e+001];
% 配电网交互功率
    P_Mmax = 1500;

%% 不确定集 [u_pv, u_L]
u_pv_est = [
    0
    0
    0
    0
    0
    4.17074232453333e+000
    5.94883875758584e+001
    2.48146877320564e+002
    5.36796660947943e+002
    8.81018103602796e+002
    1.09189870176427e+003
    1.18055770363352e+003
    1.13031213991243e+003
    9.30100120349269e+002
    8.79862750620951e+002
    8.24069853788441e+002
    5.90516477607354e+002
    3.95835403170051e+002
    1.28924282385476e+002
    0
    0
    0
    0
    0];

u_L_est = [
    3.47699739155637e+002
    3.24362863491450e+002
    3.08909284127137e+002
    2.89508260071491e+002
    3.37179016520143e+002
    4.67769297652401e+002
    5.39194280745822e+002
    5.55273886581007e+002
    6.06834122307023e+002
    7.41418220461791e+002
    8.44399574920297e+002
    9.11819147908415e+002
    7.38490967056323e+002
    6.95405274852671e+002
    6.79951695488359e+002
    7.27634045019805e+002
    7.71392136025505e+002
    7.91442372717612e+002
    8.43048980774804e+002
    8.90800888803014e+002
    8.35838083276978e+002
    5.75735677712298e+002
    4.61561201816250e+002
    4.02639358516085e+002];

%% 设置价值向量
lambda_t = [
    0.45
    0.45
    0.45
    0.45
    0.45
    0.45
    0.45
    0.9
    1.35
    1.35
    1.35
    0.9
    0.9
    0.9
    0.9
    0.9
    0.9
    0.9
    1.35
    1.35
    1.35
    1.35
    1.35
    0.45];  % 日前交易电价

%% 约束矩阵

E24 = eye(24,24);  % 24*24 单位矩阵
Z24 = zeros(24,24);    % 24*24 零矩阵
L24 = tril(ones(24,24));    % 24*24 下三角形为1阵

% 价值矩阵 c  240*1
    c = [
               a*ones(24,1);
        K_S*yita*ones(24,1);
        K_S/yita*ones(24,1);
                zeros(24,1);
            K_DR*ones(48,1);
                   lambda_t;
                  -lambda_t;
                zeros(48,1);
        ];
% D y >= d    144*240  144*1
    D = [
        E24,        Z24,           Z24,       Z24;
       -E24,        Z24,           Z24,       Z24;
        Z24,  yita.*L24,  -1/yita.*L24,       Z24;
        Z24, -yita.*L24,   1/yita.*L24,       Z24;
        Z24,        Z24,           Z24,       E24;
        Z24,        Z24,           Z24,      -E24;
        ];

    D = [D, zeros(24*6, 24*6)];

    d = [
                 P_Gmin*ones(24,1);
                -P_Gmax*ones(24,1);
        (E_Smin - E_S0)*ones(24,1);
        (E_S0 - E_Smax)*ones(24,1);
                 DR_min*ones(24,1);
                -DR_max*ones(24,1);
        ];

% K y == k  50*240  50*1
    K = [
         zeros(1,24),   yita*ones(1,24),   -1/yita*ones(1,24), zeros(1,24*7);
         zeros(1,24*3), ones(1,24),     zeros(1,24*6);
         Z24,               Z24,                  Z24,          E24,   E24, -E24,  Z24,  Z24,  Z24,  Z24;
        -E24,               E24,                 -E24,          E24,   Z24,  Z24, -E24,  E24, -E24,  E24;
        ];

    O = [
        0;
        D_DR;
        P_DRstar;
        zeros(24,1)
        ];
 
% Fx + Gy >= h  96*48   96*240   96*1
    F = [
         P_Smax*E24,            Z24;
        -P_Smax*E24,            Z24;
                Z24,     P_Mmax*E24;
                Z24,    -P_Mmax*E24;
        ];

    G = [
        Z24,  Z24, -E24, Z24, Z24, Z24,  Z24,  Z24, Z24, Z24;
        Z24, -E24,  Z24, Z24, Z24, Z24,  Z24,  Z24, Z24, Z24;
        Z24,  Z24,  Z24, Z24, Z24, Z24, -E24,  Z24, Z24, Z24;
        Z24,  Z24,  Z24, Z24, Z24, Z24,  Z24, -E24, Z24, Z24;
        ];

    h = [
               zeros(24,1);
        -P_Smax*ones(24,1);
               zeros(24,1);
        -P_Mmax*ones(24,1);
        ];


% I_u y == u_hat 
    I_u = [
        zeros(24,192), E24, Z24;
        zeros(24,192), Z24, E24;
        ];


%% 初始主问题设置

% 迭代次数
k = 1;

u_pv_k(:,k) = u_pv_est; %- u_pv_est*0.15;
u_L_k(:,k) = u_L_est; %+ u_L_est*0.1;

% [U_S, U_M]
x = binvar(48,1);

% [P_G, P_Sch, P_Sdis, P_DR, P_DR1, P_DR2, P_Mbuy, P_Msell, P_pv, P_L]
y_k = sdpvar(240,itr,'full');

% 第一阶段目标函数
alpha = sdpvar(1);
MP_Obj = alpha;

% 第一阶段约束
% MP_cons = [
%     alpha >= c'*y_k(:,1);
%     D*y_k(:,1) >= d;
%     K*y_k(:,1) == O;
%     F*x + G*y_k(:,1) >= h
%     I_u*y_k(:,1) == u_hat + du;
%     y_k(:,1) >= 0
%     ];

MP_cons = [
            alpha >= c'*y_k(:,k);
            
            D*y_k(:,k) >= d;
            K*y_k(:,k) == O;
            F*x + G*y_k(:,k) >= h
            I_u*y_k(:,k) == [u_pv_est; u_L_est];
            y_k(:,k) >= 0
            ];


%% 循环

% 上下界
UB = inf;
LB = -inf;
lb = []; ub = [];

% 大 M 法用到的 M
M = 1e5;

% 不确定性调节参数
Gamma_pv = 6;
Gamma_L = 12;

% 求解器设置
opt = sdpsettings('verbose', 0,'solver','gurobi');

% u_k存储
u_pv_k = zeros(24,itr);
u_L_k = zeros(24,itr);

x_kk = zeros(48,itr);

tic;
hold on;
while UB - LB >= 1e-5

    % 求解第一阶段问题
    % disp('                          ');
    % disp('**************************');
    % disp(['第',num2str(k),'次迭代 MP']);
    % disp('--------------------------');
    % disp('                          ');

    MP = optimize(MP_cons, MP_Obj, opt);

    % 更新下界
    LB = max([LB,value(MP_Obj)]);
    lb = [lb, LB];

    % 取得最优的 x
    x_k = value(x);
    x_kk(:,k) = x_k; 
    
    % 第二阶段约束矩阵
    A = [D;K;-K;G;I_u;-I_u];

    % 第二阶段问题的决策变量
    y = sdpvar(240,1);
    
    % 第二阶段问题对偶变量定义
    pi = sdpvar(size(A,1),1);

    % 不确定集的 u
    u_pv = sdpvar(24,1);
    u_L = sdpvar(24,1);

    % 第二阶段右端向量
    b = [d;O;-O;h-F*x_k;u_pv;u_L;-u_pv;-u_L];

    % 为解决双线性项引入的0-1变量与连续变量
    B_pv = binvar(24,1);
    B_L = binvar(24,1);
    

    % 为线性化互补松紧性引入的 0-1 变量
    v = binvar(size(b,1),1);
    w = binvar(240,1);

    % 第二阶段问题目标函数
    SP_Obj = b'*pi;

    % 第二阶段问题约束
    SP_cons = [

        A*y >= b;
        A'*pi <= c;

        u_pv == u_pv_est - B_pv.*u_pv_est*0.15;
        u_L == u_L_est + B_L.*u_L_est*0.1;

        sum(B_pv) <= Gamma_pv;
        sum(B_L) <= Gamma_L;


        % 互补松紧性
        A*y - b <= M*(1 - v);
        pi <= M*v;

        c - A'*pi <= M*(1 - w);
        y <= M*w;

        y >= 0;
        pi >= 0;
        ];
    
    % 求解第二阶段问题
    % disp('                          ');
    % disp('--------------------------');
    % disp(['第',num2str(k),'次迭代 SP']);
    % disp('--------------------------');
    % disp('                          ');
    SP = optimize(SP_cons, -SP_Obj, opt);

    if SP.problem == 0
        % feasible

        % 更新上界
        UB = min(UB, value(SP_Obj));
        ub = [ub, UB];

        % 获取 u_k
        u_pv_k(:,k) = u_pv_est - value(B_pv).*u_pv_est*0.15;
        u_L_k(:,k) = u_L_est + value(B_L).*u_L_est*0.1;

        % 主问题增加约束
        MP_cons = [
            MP_cons;

            alpha >= c'*y_k(:,k+1);

            D*y_k(:,k+1) >= d;
            K*y_k(:,k+1) == O;
            F*x + G*y_k(:,k+1) >= h
            I_u*y_k(:,k+1) == [u_pv_k(:,k); u_L_k(:,k)];
            
            y_k(:,k+1) >= 0

            ];
    else
        % infeasible
        disp("问题无解");
        break;
    end
    

    disp("------------------------");
    disp(['第', num2str(k), '次迭代']);
    disp(['UB = ', num2str(UB)]);
    disp(['LB = ', num2str(LB)]);
    disp("------------------------");

    k = k + 1;

    % 收敛情况
    plot(ub,'Color','r','Marker','*');
    plot(lb,'Color','b','Marker','+');
    legend("UB","LB");
    xlabel("迭代次数");
    ylabel("最优解上下界");

    if k >= itr
        break;
    end

end
hold off;
toc;

%% 可视化
k = k-1;

% 生成索引序列
index = zeros(10,1);
for i = 1 : 10
    index(i) = 1 + 24*(i-1);
end

% 取出最优解
y_opt = value(y);

% 计算各时刻功率
P_G = y_opt(1:24);
P_buysell = -y_opt(index(7):index(8)-1) + y_opt(index(8):index(9)-1);
P_es = -y_opt(index(2):index(3)-1) + y_opt(index(3):index(4)-1);
P_dr = y_opt(index(4):index(5)-1);


figure;
bar([P_G,P_buysell]);   % 微燃机输出功率及微电网购售电功率
legend('微型燃气轮机','微电网购售电');
xlabel("时间/h");
ylabel("功率/kW");
ylim([-1500,1500]);

figure;
bar(P_es);  % 储能充放电功率
legend("储能充放电");
xlabel("时间/h");
ylabel("功率/kW");
ylim([-600,600]);

figure;
bar([P_dr,P_DRstar]);
legend("实际调度功率","期望用电计划");
xlabel("时间/h");
ylabel("功率/kW");
ylim([0,300]);