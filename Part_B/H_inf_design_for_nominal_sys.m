%% Robust Control Assignment
%% Controller Design
clc; close all; clear;
load('ECP_values.mat');

% H_inf controller design for the nominal system

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

omega_B = 3;
A = 10^-9;
M = 1;
W_P = (s/M + omega_B)/(s+omega_B*A);
W_K = [];
W_T = [];

figure(Name="W_P Bode Plot");
bodemag(1/W_P, 'b'); 
grid on;
title('W_P');
legend('1/W_p');

%%

% Controller design
[K, CL, Gam, INFO] = mixsyn(G_td, W_P, W_K, W_T);
% Get sensitivity and complementary sensitivity function
S = (1 + G_td*K)^-1;
T = G_td*K*(1 + G_td*K)^-1;

%[mag_S, ~, w_S] = bode(S);
%[mag_T, ~, w_T] = bode(CL);
%[mag_K, ~, w_K] = bode(K);

% Plot both S & T in one bode plot

figureName="S and T Bode Plot";
bodemag(S, 'r', T, 'b'); 
grid on;
title('S and T');
legend('S (Sensitivity)', 'T (Complementary Sensitivity)');

% Plot K bode plot
figure(Name="Controller Bode Plot");
bodemag(K, 'g'); 
grid on;
title('K');
legend('K (Controller)');

figure(Name="Step Response");

step(T); grid on;

figure(Name="Nyquist");

nyquist(S)

