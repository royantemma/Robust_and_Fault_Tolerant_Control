%% Robust Control Assignment

%% Controller Design
clc; close all; clear;
load('ECP_values.mat');
w = logspace(-6,4,100);

%% Get model

s = tf("s");

% Physical system parameters
J_1 = ECP_values(1);            % Disk 1 inertia kgm^
%J_1 = 0.0325;
J_2 = ECP_values(2);            % Disk 2 inertia kgm^2
J_3 = ECP_values(3);            % Disk 3 inertia kgm^2
k_1 = ECP_values(4);            % Shaft 1-2 stiffness Nm/rad
k_2 = ECP_values(5);            % Shaft 2-3 stiffness Nm/rad
b_1 = mean(ECP_values([6 7]));  % Disk 1 damping and friction Nms/rad
b_2 = mean(ECP_values([8 9]));  % Disk 2 damping and friction Nms/rad
b_3 = mean(ECP_values([10 11]));% Disk 3 damping and friction Nms/rad
T_Cp = ECP_values(12);          % Disk 1 Coulomb friction in positive direction
T_Cm = ECP_values(13);          % Disk 1 Coulomb friction in negative direction
atan_scale = 100;               % Sign approximation factor
w_th = 0.75;                    % Threshold angular velocity rad/s

A = [ 0 1 0 0 0 0
      -k_1/J_1 -b_1/J_1 k_1/J_1 0 0 0
      0 0 0 1 0 0
      k_1/J_2 0 -(k_1+k_2)/J_2 -b_2/J_2 k_2/J_2 0
      0 0 0 0 0 1
      0 0 k_2/J_3 0 -k_2/J_3 -b_3/J_3];

B = [ 0 0
      1/J_1 0
      0 0
      0 1/J_2
      0 0
      0 0];

C = [  1 0 0 0 0 0
       0 0 1 0 0 0
       0 0 0 0 1 0 ];

D = zeros(3,2);

x_0 = zeros(1,6);
C_td = [0 0 0 0 1 0];
B_td = B(:,1);
D_td = zeros(1,1);

[num_td, den_td] = ss2tf(A,B_td,C_td,D_td);

G_td = tf(num_td, den_td);

%% Question 9 Mixed sensitivity setup
close all;
% It is going to be difficult and take some time - start simple and
% increase the order of the weight function 2nd or maybe higher

% Weights
omega_B = 1;
A = 10^-9;
M = 2;

W_P = 2.5 * (s/M + omega_B)*(s/M + 0.5)/((s+omega_B*A)*(s+0.001)) ; %cutband filter
W_K = 0.5 * (s/2 + 1)^2 / ((s/20 + 1)^2);
W_T = [];

% Optimal Controller design
[K, CL, Gam, INFO] = mixsyn(G_td, W_P, W_K, W_T);
Gam

% Get sensitivity and complementary sensitivity function
S = (1 + G_td*K)^-1;
T = G_td*K*(1 + G_td*K)^-1; 
L = G_td*K;
u = K*S;

% Plots
figure(Name="Sensitivity vs Weight");
bodemag(S, 'b', 1/W_P, 'r', w); 
grid on;
title('Sensitivity vs Weight');
legend('S','1/W_P');

figure(Name="Controller Bode Plot");
bodemag(K*S, 'b', 1/W_K, 'r', w); 
grid on;
title('Control Effort');
legend('KS', '1/W_K');

figure(Name="Loop Transfer L Bode Plot");
bodemag(L, 'g', w); 
grid on;
title('Loop Transfer L');
legend('L');

figure(Name="Controller Bode Plot");
bodemag(K, 'g', w); 
grid on;
title('K');
legend('K (Controller)');

figure(Name="Step Response");
step(T); 
grid on;
title('Closed-loop Step Response');

figure;
step(u);
grid on;
title('Control Input');

figure(Name="Nyquist");
nyquist(S)

%% Relaxed design

alpha = 1.05;   % 5% relaxation

W_P2 = W_P / alpha;
W_T2 = W_T / alpha;
W_K2 = W_K / alpha;

[K2, CL2, gam2] = mixsyn(G_td, W_P2, W_K2, W_T2);

% Get sensitivity and complementary sensitivity function
S2 = (1 + G_td*K2)^-1;
T2 = G_td*K2*(1 + G_td*K2)^-1;

% Plot both S & T in one bode plot

figure(Name="S and T Bode Plot");
bodemag(S, 'r', T, 'b', S2, 'y--', T2, 'g--'); 
grid on;
title('S and T');
legend('S (Sensitivity)', 'T (Complementary Sensitivity)', 'S2 (Sensitivity)', 'T2 (Complementary Sensitivity)');

% Plot K bode plot
figure(Name="Controller Bode Plot");
bodemag(K, 'g', K2, 'b--'); 
grid on;
title('K');
legend('K (Controller)','K2 (Controller)');

figure(Name="Step Response");

step(T, 'r', T2, 'b--'); grid on;

%% Question 10 Reduction of controller order
close all;

n_red = 6;
K_red = balred(K, n_red);

T_red  = feedback(G_td*K_red,1);

% Check stability
isstable(T)
isstable(T_red)

% Poles and zeros comparison
pole(K_red)
zero(K_red)

figure;
pzmap(K); 
hold on;
pzmap(K_red);
grid on;
legend('Full K','Reduced K');
title('Pole-Zero Map Comparison');

% Controllers comparison
figure;
bodemag(K, K_red, w);
grid on;
legend('Full K','Reduced K');
title('Controller Magnitude Comparison');

% Closed-Loop comparison
figure;
bodemag(T, T_red, w);
grid on;
legend('Full CL','Reduced CL');
title('Closed-loop Frequency Response');

figure;
step(T, T_red);
grid on;
legend('Full','Reduced');
title('Step Response Comparison');

%% Question 11 H2 controller design
close all;
% H2 design (approximate using mixsyn with large gamma)
gamRange = [50,50];

[K_H2,~,gamma_H2] = mixsyn(G_td, W_P, W_K, W_T,gamRange);
gamma_H2

T_H2 = feedback(G_td*K_H2,1);

% Comparisons
% Controllers
figure; bodemag(K, K_H2); grid on;
legend('Hinf','H2');
title('Controller Comparison');

% Closed-loop transfer functions
figure; bodemag(T, T_H2); grid on;
legend('Hinf','H2');
title('Closed-loop Comparison');

% Step responses
figure; step(T, T_H2); grid on;
legend('Hinf','H2');
title('Step Response Comparison');

%% Question 13 Robustness analysis
close all;

L_red = G_td * K_red;
S_red = feedback(1,L_red);

% Condition : |WI * T| < 1  → max WI = 1/|T|

WI_max = 1/T;

figure;
bodemag(WI_max,w)
grid on;
title('Maximum multiplicative uncertainty |W_I|');


%% Question 14 Model uncertainty
close all;

J1_new = 0.0325;

A_new = [ 0 1 0 0 0 0
         -k_1/J1_new -b_1/J1_new k_1/J1_new 0 0 0
          0 0 0 1 0 0
          k_1/J_2 0 -(k_1+k_2)/J_2 -b_2/J_2 k_2/J_2 0
          0 0 0 0 0 1
          0 0 k_2/J_3 0 -k_2/J_3 -b_3/J_3];

B_new = [ 0
      1/J1_new
      0
      0
      0
      0];

[num_new,den_new] = ss2tf(A_new,B_new,C_td,D_td);
G_new = tf(num_new,den_new);

T_new = feedback(G_new*K_red,1);

figure;
bodemag(w, 1/T_new);
grid on;
title('New |W_I|');

%% Question 15 Hinf design with uncertainty
close all;

% Uncertainty weight (given)
WI = 0.833*s / (s + 0.089);

figure;
bodemag(WI, 'g', WI_max, 'r', 1/T_new, 'b', w);
grid on;
title('Multiplicative uncertainty |W_I|');

% Weight functions
omega_B = 1;
A = 10^-9;
M = 2;

W_P = 1.1 * (s/M + omega_B)*(s/M + 0.5)/((s+omega_B*A)*(s+0.001));
W_K = 0.5 * (s/2 + 1)^2 / ((s/20 + 1)^2);
W_T = WI;

% Optimal Controller design
[K_unc,CL_unc,gamma_unc] = mixsyn(G_td, W_P, W_K, W_T);
gamma_unc

% Get sensitivity and complementary sensitivity function
L_unc = G_td * K_unc;
S_unc = feedback(1,L_unc);
T_unc = feedback(L_unc,1);

% Plots
figure;
bodemag(S_unc, 1/W_P);
grid on;
title('Sensitivity vs weight');

figure;
bodemag(K_unc*S_unc, 1/W_K);
grid on;
title('Control effort');

figure;
bodemag(T_unc, 1/W_T);
grid on;
title('Complementary sensitivity');

figure;
bodemag(K_unc);
grid on;
title('Controller');

figure;
step(T_unc);
grid on;
title('Step response');