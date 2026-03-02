clear all;
close all;
clc;
load('ECP_values.mat');
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
atan_scale = 100;               % Sign approximation factor
w_th = 0.75;                    % Threshold angular velocity rad/s

% The system states are [theta_1;omega_1;theta_2;omega_2;theta_3;omega_3]
x_0 = [0;0;0;0;0;0];            % Initial conditions
T_s = 0.004;
simTime = 200;
% Sampling period
sigma_meas = 0.0093*eye(3);     % Measurements covariance matrix


% Utilities 
s = tf("s");


%% State space representation
close all;
A = [ 0 1 0 0 0 0
      -k_1/J_1 -b_1/J_1 k_1/J_1 0 0 0
      0 0 0 1 0 0
      k_1/J_2 0 -(k_1+k_2)/J_2 -b_2/J_2 k_2/J_2 0
      0 0 0 0 0 1
      0 0 k_2/J_3 0 -k_2/J_3 -b_3/J_3];

B = [ 0 0
      1 0
      0 0
      0 1
      0 0
      0 0];

C = [  1 0 0 0 0 0
       0 0 1 0 0 0
       0 0 0 0 1 0 ];

D = zeros(3,2);

E_x = [ 0; 1; 0; 0; 0; 0 ];

E_y = [ 0; 0; 0 ];

F_x = [ 0 0 0 0 0 0 
        1/J_1 0 0 0 0 0 
        0 0 0 0 0 0 
        0 1/J_2 0 0 0 0 
        0 0 0 0 0 0 
        0 0 0 0 0 0 ];

F_y = [ 0 0 1 0 0
        0 0 0 1 0
        0 0 0 0 1];

sys = ss(A,B,C,D);
% step(sys)
sys_d = c2d(sys, T_s, 'tustin');

%%%%%%%% SKIP THAT WHEN SIMULATING IN OPEN LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discrete time
F_d = sys_d.A;
G_d = sys_d.B;

% State-feedback LQR design
Q_c = diag([2 0 2 0 2.5 0.0024]);
R_c = diag([10 10]);
K_c = [];

% Scaling of reference
C_ref = [];

% Kalman filter with friction estimation - DO NOT MODIFY
F_aug = [F_d G_d(:,1);zeros(1,6) 1];
G_aug = [G_d;0 0];
C_aug = [C zeros(3,1)];
% Kalman gain
L_aug = dlqe(F_aug,eye(7),C_aug,1e-3*eye(7),sigma_meas(1,1).^2*eye(3));
L_o = L_aug(1:6,:);
L_d = L_aug(7,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Residual filter design
% Design of characteristic polynomial
f_c = 1; % Hz - Desired cut-off frequency
w_n = 2*pi*f_c;
w_n2 = w_n^2;
zeta = 1 / sqrt(2);

G = w_n^2 / ((s^2 + 2*zeta*w_n*s + w_n^2)*(s+zeta*w_n)^2);
% figure;
% bode(G)
% figure;
% step(G)

% Residuals with Filter in laplace domain
G_r1y2 = (s^2 + b_2/J_2*s + (k_1+k_2)/J_2)*G;
G_r1y1 = -(k_1/J_2)*G;
G_r1y3 = -(k_2/J_2)*G;
G_r1u2 = -(1/J_2)*G;

G_r2y3 = (s^2 + b_3/J_3*s + (k_2)/J_3)*G;
G_r2y2 = -(k_2/J_3)*G;

% Discretize 
G_r1y2_d = c2d(G_r1y2,T_s, 'tustin');
G_r1y1_d = c2d(G_r1y1,T_s, 'tustin');
G_r1y3_d = c2d(G_r1y3,T_s, 'tustin');
G_r1u2_d = c2d(G_r1u2,T_s, 'tustin');
G_r2y3_d = c2d(G_r2y3,T_s, 'tustin');
G_r2y2_d = c2d(G_r2y2,T_s, 'tustin');

% Linearcombinations of residuals
a = -(1/J_2);
b = -(k_1/J_2);
c = (s^2 + b_2/J_2*s + (k_1+k_2)/J_2);
d = -(k_2/J_2);
e = -(k_2/J_3);
f = (s^2 + b_3/J_3*s + (k_2)/J_3);

G_r3y1 = b*G;
G_r3y3 = (d-c*f/e)*G;
G_r3u2 = a*G;

G_r4y1 = -f/d*b*G;
G_r4y2 = (e-f*c/d)*G;
G_r4u2 = -f/d*a*G;

G_r3y1_d = c2d(G_r3y1,T_s, 'tustin');
G_r3y3_d = c2d(G_r3y3,T_s, 'tustin');
G_r3u2_d = c2d(G_r3u2,T_s, 'tustin');
G_r4y1_d = c2d(G_r4y1,T_s, 'tustin');
G_r4y2_d = c2d(G_r4y2,T_s, 'tustin');
G_r4u2_d = c2d(G_r4u2,T_s, 'tustin');

%% Question 3: Experimental Data
close all;
% Load in data
data = load("ECP502Data.mat");
t = data.t;
u_1 = data.u_1;
u_2 = data.u_2;
y_meas = data.y_meas;

% Plotting
close all;
figure;
subplot(4,1,1)
plot(t,y_meas)
legend({'y_1','y_2','y_3'})
title('Measured y-values')
ylabel('Sensor Measurements')
xlabel('Time (s)')
grid on;
xlim([0 T_s*length(u_1)])
subplot(4,1,2)
hold on;
plot(t,u_1)
plot(t,u_2)
hold off;
legend({'u_1','u_2'})
title('Given Inputs')
xlim([0 T_s*length(u_1)])
ylabel('Input Torque (Nm)')
xlabel('Time (s)')
grid on;

subplot(4,1,3)
r_1 = lsim(G_r1y2_d,y_meas(:,2)) + lsim(G_r1y1_d,y_meas(:,1),t) + lsim(G_r1y3_d,y_meas(:,3),t) + lsim(G_r1u2_d,u_2,t); 
r_2 = lsim(G_r2y3_d,y_meas(:,3),t) + lsim(G_r2y2_d,y_meas(:,2),t);
hold on;
grid on;
plot(t,r_1)
plot(t,r_2)
hold off;
legend({'r_1','r_2'})
title('Residuals')
ylabel('Residual Value')
xlabel('Time (s)')
xlim([0 T_s*length(u_1)])

subplot(4,1,4)
r_3 = lsim(G_r3y1_d,y_meas(:,1)) + lsim(G_r3y3_d,y_meas(:,3),t) + lsim(G_r3u2_d,u_2,t); 
r_4 = lsim(G_r4y1_d,y_meas(:,1)) + lsim(G_r4y2_d,y_meas(:,2),t) + lsim(G_r4u2_d,u_2,t); 
hold on;
grid on;
plot(t,r_3)
plot(t,r_4)
hold off;
legend({'r_3','r_4'})
title('New Residuals')
ylabel('Residual Value')
xlabel('Time (s)')
xlim([0 T_s*length(u_1)])
%% Strong and weak detectability
H_rf = tf(0);

%% GLR
f_m = [0;-0.025;0];     % Sensor fault vector (added to [y1;y2;y3])
h = 0;                  % Put the threshold from GLR here

%% Virtual actuator
% Failure in actuator 2
% Do the desing first in continuous time
va_eig_d = [];  % Discrete time eigenvalues
va_eig = log(va_eig_d)/T_s;     % Continuous time eigenvalues
% Then discretise your VA

B_change = [1 0;0 0];

%% Simulation for sensor fault (f_u = 0)
simTime = 45;                   % Simulation duration in seconds
f_u_time = 25;                  % Actuator fault occurence time
detect_time = f_u_time + 3.75;
f_u = [0;0];                    % Actuator fault vector (added to [u1;u2])
u_fault = 0;                    % Disable VA meachanism
f_m_time = 8.5;                 % Sensor fault occurence time
sim('threeDiskOscillatorRig');

%% Simulation for actuator fault (f_m = 0)
f_u = [0;-0.1];                 % Actuator fault vector (added to [u1;u2])
u_fault = 1;                    % Enable VA meachanism
f_m = [0;0;0];                  % Sensor fault vector (added to [y1;y2;y3])
sim('threeDiskOscillatorRig');

%% Plot settings
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultAxesFontSize',20);
set(0,'DefaultLineLineWidth', 2);

