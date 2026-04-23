%% Robust Control Assignment

%% Controller Design
clc; close all; clear;
load('ECP_values.mat');

%% Get model

s = tf("s");

% Physical system parameters
J_1 = ECP_values(1);            % Disk 1 inertia kgm^2
J_2 = ECP_values(2);            % Disk 2 inertia kgm^2
J_3 = ECP_values(3);            % Disk 3 inertia kgm^2
k_1 = ECP_values(4);            % Shaft 1-2 stiffness Nm/rad
k_2 = ECP_values(5);            % Shaft 2-3 stiffness Nm/rad
b_1 = mean(ECP_values([6 7]));  % Disk 1 damping and friction Nms/rad
b_2 = mean(ECP_values([8 9]));  % Disk 2 damping and friction Nms/rad
b_3 = mean(ECP_values([10 11]));% Disk 3 damping and friction Nms/rad
T_Cp = ECP_values(12);          % Disk 1 Coulomb friction in positive direction
T_Cm = ECP_values(13);          % Disk 1 Coulomb friction in negative direction


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

C_td = [0 0 0 0 1 0];
B_td = B(:,1);
D_td = zeros(1,1);

[num_td, den_td] = ss2tf(A,B_td,C_td,D_td);

G_td = tf(num_td, den_td);

%% Mixed sensitivity setup
close all;
% It is going to be difficult and take some time - start simple and
% increase the order of the weight function 2nd or maybe higher

% Weights
omega_B = 1;
A = 10^-6;
M = 2;

W_P = 2.5 * (s/M + omega_B)*(s/M + 0.5)/((s+omega_B*A)*(s+0.5));
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
bodemag(S, 'b', 1/W_P, 'r'); 
grid on;
title('Sensitivity vs Weight');
legend('S','1/W_P');

figure(Name="Controller Bode Plot");
bodemag(K*S, 'b', 1/W_K, 'r'); 
grid on;
title('Control Effort');
legend('KS', '1/W_K');

figure(Name="Loop Transfer L Bode Plot");
bodemag(L, 'g'); 
grid on;
title('Loop Transfer L');
legend('L');

figure(Name="Controller Bode Plot");
bodemag(K, 'g', G_td, 'b'); 
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

% Relaxed design

alpha = 1.05;   % 5% relaxation

W_P2 = W_P / alpha;
W_T2 = W_T / alpha;
W_K2 = W_K;

[K2, CL2, gam2] = mixsyn(G_td, W_P2, W_K2, W_T2);

% Get sensitivity and complementary sensitivity function
S2 = (1 + G_td*K2)^-1;
T2 = G_td*K2*(1 + G_td*K2)^-1;

% Plot both S & T in one bode plot

figureName="S and T Bode Plot";
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

%% Reduction of controller order
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
bodemag(K, K_red);
grid on;
legend('Full K','Reduced K');
title('Controller Magnitude Comparison');

% Closed-Loop comparison
figure;
bodemag(T, T_red);
grid on;
legend('Full CL','Reduced CL');
title('Closed-loop Frequency Response');

figure;
step(T, T_red);
grid on;
legend('Full','Reduced');
title('Step Response Comparison');

%% H2 controller design
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