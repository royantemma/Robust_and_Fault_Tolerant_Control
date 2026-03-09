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

F_x = [ 0 0 0 0 0;
        0 0 0 1/J_1 0; 
        0 0 0 0 0; 
        0 0 0 0 1/J_2; 
        0 0 0 0 0; 
        0 0 0 0 0];

F_y = [1 0 0 0 0
       0 1 0 0 0
       0 0 1 0 0];

sys = ss(A,B,C,D);
% step(sys)
sys_d = c2d(sys, T_s,'tustin');
sys_d_cont = c2d(sys, T_s);
%%%%%%%% SKIP THAT WHEN SIMULATING IN OPEN LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discrete time
F_d = sys_d_cont.A;
G_d = sys_d_cont.B;

% State-feedback LQR design
Q_c = diag([2 0 2 0 2.5 0.0024]);
R_c = diag([10 10]);
K_c = dlqr(F_d,G_d,Q_c,R_c);
% Check closed loop stability
disp("LQR - CL poles")
eig(F_d-G_d*K_c)

cl_sys = ss(F_d-G_d*K_c,G_d,C,0,T_s);
pzmap(cl_sys)

% Scaling of reference
C_ref = pinv(C(3,:)/(eye(6) - F_d + G_d*K_c)*G_d*K_c); % reference scaling (V in the slides)


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
f_c = 10; % Hz - Desired cut-off frequency
w_n = 2*pi*f_c;
w_n2 = w_n^2;
zeta = 1 / sqrt(2);

G = w_n^2 / (s^2 + 2*zeta*w_n*s + w_n^2);
G_lc = G /((s+zeta*w_n)^2);
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

G_r3y1 = b*G_lc;
G_r3y3 = (d-c*f/e)*G_lc;
G_r3u2 = a*G_lc;

G_r4y1 = -f/d*b*G_lc;
G_r4y2 = (e-f*c/d)*G_lc;
G_r4u2 = -f/d*a*G_lc;

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

u1_ts = timeseries(u_1, t);
%u1_ts = setinterpmethod(u1_ts,'zoh');
u2_ts = timeseries(u_2, t);
%u2_ts = setinterpmethod(u2_ts,'zoh');
y_meas_ts = timeseries(y_meas, t);
y_meas_ts = setinterpmethod(y_meas_ts,'zoh');

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
r_1 = lsim(G_r1y2_d,y_meas(:,2),t) + lsim(G_r1y1_d,y_meas(:,1),t) + lsim(G_r1y3_d,y_meas(:,3),t) + lsim(G_r1u2_d,u_2,t); 
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
r_3 = lsim(G_r3y1_d,y_meas(:,1),t) + lsim(G_r3y3_d,y_meas(:,3),t) + lsim(G_r3u2_d,u_2,t); 
r_4 = lsim(G_r4y1_d,y_meas(:,1),t) + lsim(G_r4y2_d,y_meas(:,2),t) + lsim(G_r4u2_d,u_2,t); 
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
clc;
syms s
% Constuct symbolic tf's
[n_yu1,d_yu1] = ss2tf(A,B,C,D,1);
[n_yu2,d_yu2] = ss2tf(A,B,C,D,2);
[n_yd,d_yd] = ss2tf(A,E_x,C,E_y);
[n_yf1,d_yf1] = ss2tf(A,F_x,C,F_y,1);
[n_yf2,d_yf2] = ss2tf(A,F_x,C,F_y,2);
[n_yf3,d_yf3] = ss2tf(A,F_x,C,F_y,3);
[n_yf4,d_yf4] = ss2tf(A,F_x,C,F_y,4);
[n_yf5,d_yf5] = ss2tf(A,F_x,C,F_y,5);



function H = tf_from_numden(n,d,s)
    
    [m, ~] = size(n);
    H = sym(zeros(m,1));
    for i = 1:m
        num = poly2sym(n(i,:), s);
        
        if size(d, 1) > 1
            den = poly2sym(d(i,:), s);
        else
            den = poly2sym(d, s);
        end
        
        H(i) = num / den;
    end
end

H_yu1 = tf_from_numden(n_yu1,d_yu1,s);
H_yu2 = tf_from_numden(n_yu2,d_yu2,s);
H_yu = [H_yu1 H_yu2];

H_yd = tf_from_numden(n_yd,d_yd,s);


H_yf1 = tf_from_numden(n_yf1,d_yf1,s);
H_yf2 = tf_from_numden(n_yf2,d_yf2,s);
H_yf3 = tf_from_numden(n_yf3,d_yf3,s);
H_yf4 = tf_from_numden(n_yf4,d_yf4,s);
H_yf5 = tf_from_numden(n_yf5,d_yf5,s);

H_yf = [H_yf1 H_yf2 H_yf3 H_yf4 H_yf5];

H = [H_yu H_yd;
     eye(size(H_yu,2)) zeros(size(H_yu,2),size(H_yd,2))];

F = null(H')';

% V_ry = F(:,1:size(H_yu,1));
% H_yf = V_ry * 

weak_detectability = zeros(5,1);
strong_detectability = zeros(5,1);
for i = 1:5
    weak_detectability(i) = rank([H_yd H_yf(:,i)]) > rank(H_yd);
    strong_detectability(i) = any(limit(F*[H_yf(:,i); 0; 0],s,inf) ~= 0);
end
disp("Case Disturbance: Unknown")
weak_detectability
strong_detectability

disp("Case Disturbance: Known (Modelled as an input")

H_yu_dasinput = [H_yu1 H_yu2 H_yd];
H_yd_dasinput = zeros(3,1);
H_dasinput = [H_yu_dasinput H_yd_dasinput;
     eye(size(H_yu_dasinput,2)) zeros(size(H_yu_dasinput,2),size(H_yd_dasinput,2))];

F_dasinput = null(H_dasinput')';

weak_detectability = zeros(5,1);
strong_detectability = zeros(5,1);
for i = 1:5
    weak_detectability(i) = rank([H_yd H_yf(:,i)]) > rank(H_yd);
    strong_detectability(i) = any(limit(F_dasinput*[H_yf(:,i); 0; 0;0;],s,inf) ~= 0);
end

weak_detectability
strong_detectability

%{
% Finding H_rf
% H_rf = V_ry * H_yf
H_yf = C * (s*eye(6)-A)^-1 * F_x + F_y;

G_r2y1 = 0;
G_r3y2 = 0;
G_r4y3 = 0;
G_r1u1 = 0;
G_r2u1 = 0;
G_r2u2 = 0;
G_r3u1 = 0;
G_r4u1 = 0;


V_ry = [G_r1y1 G_r1y2 G_r1y3;
        G_r2y1 G_r2y2 G_r2y3;
        G_r3y1 G_r3y2 G_r3y3;
        G_r4y1 G_r4y2 G_r4y3;];

H_rf = V_ry*H_yf;

H_rf = minreal(H_rf);

% Investigate Weak Detectability
% None column of Hyf should be able to be written as a linear
% combination of the columns in H_yd
H_yd = C * (s*eye(6)-A)^-1 * E_x + E_y;

% F(s) = [V_ry  V_ru]
V_ru = [G_r1u1 G_r1u2;
        G_r2u1 G_r2u2;
        G_r3u1 G_r3u2;
        G_r4u1 G_r4u2;];
F = [V_ry V_ru];
weak_detectability = zeros(5,1);
strong_detectability = zeros(5,1);


for i = 1:5
    %H_yf_i = H_yf(:,i);
    %weak_detectability(i) = rank([H_yd H_yf(:,i)]) > rank(H_yd);
    %disp(dcgain([H_yf(:,i);0;0]));
    strong_detectability(i) = any(dcgain(F*[H_yf(:,i); 0; 0]) ~= 0);
end
strong_detectability
%}

%% Variance of r2
H = [G_r2y2, G_r2y3];

sys_r = ss(G_r2y2);

A_r = sys_r.A;
B_r = sys_r.B;
Q_w = sigma_meas(2,2)

Q = lyap(A_r,B_r*Q_w*B_r');

C_r = sys_r.C
D_r = sys_r.D

sigma_r21 = C_r*Q*C_r' + D_r*Q_w*D_r'

sys_r = ss(G_r2y3);

A_r = sys_r.A;
B_r = sys_r.B;
Q_w = sigma_meas(2,2)

Q = lyap(A_r,B_r*Q_w*B_r');

C_r = sys_r.C
D_r = sys_r.D

sigma_r22 = C_r*Q*C_r' + D_r*Q_w*D_r'

sigma_r2 = sigma_r21 + sigma_r22


%% GLR
f_m = [0;-0.025;0];     % Sensor fault vector (added to [y1;y2;y3])
%r2_nofault = out.r_nofault.Data(:,2);

% Threshold h
mu_0 = 0.0016; %mean(r2_nofault)
sigma_r2 = 0.0087; %var(r2_nofault)

P_F = 1e-4;

syms zz gg;     % zz is the integrant variable, gg is h in symbolic
% Gamma distribution parameters
k = 1/2;
theta = 2;
% Density function expression
pd_zz = 1/(theta^k*gamma(k))*zz^(k-1)*exp(-zz/theta);
p_zz = int(pd_zz,zz,2*gg,Inf);  % Integrate over the probability space
eq_1 = P_F - p_zz == 0;  % Equation to be solved
h = double(vpasolve(eq_1,gg))  % Convert to doublef_m
% check
p_calc = double(int(pd_zz,zz,2*h,Inf));
%disp(P_F - p_calc);
%disp(h - chi2inv(1 - P_F,1)/2); % Compare with the other methodw
% Window size M
f2 = f_m(2);

h = 7.5684

P_M = 0.01;

for M = 1:500
    lambda = M*f2^2/sigma_r2;
    P_D = ncx2cdf(h,1,lambda);
    if P_D <= 1-P_M
        break
    end
end

M

%% Virtual actuator
% Failure in actuator 2
% Do the desing first in continuous time
va_eig_d = [];  % Discrete time eigenvalues
va_eig = log(va_eig_d)/T_s;     % Continuous time eigenvalues
% Then discretise your VA

B_change = [1 0;0 0];

B_f = [B(:,1), zeros(size(B,1),1)];

% Discretization
fsys = ss(A,B_f,C,D);
fsys_d = c2d(sys,T_s);

G_f = fsys_d.B;
if rank(B_f) == rank([B B_f])
    disp('Perfect static matching for actuator fault');
else
    disp('Imperfect static matching for actuator fault');
end

fS_c = ctrb(A,B_f);
fS_o = obsv(A,C);
if rank(fS_c) == size(A,1)
    disp('Faulty system is controllable');
else
    disp('Faulty system is not controllable');
end

% Continuous time
va_eig = log(eig(F_d - G_d*K_c))/T_s;
M_va = place(A,B_f,va_eig);
A_D = A - B_f*M_va;
N_D = pinv(B_f)*B;
B_D = B - B_f*N_D;
C_D = C;

% Discrete time
va_eig_d = exp(va_eig*T_s);
M_d_va = place(F_d,G_f,va_eig_d);
F_D = F_d - G_f*M_d_va;
N_D_d = pinv(G_f)*G_f;
G_D = G_d - G_f*N_D_d;
C_D = C;



%% Simulation for sensor fault (f_u = 0)
simTime = 45;                   % Simulation duration in seconds
f_u_time = 25;                  % Actuator fault occurence time
detect_time = f_u_time + 3.75;
f_u = [0;0];                    % Actuator fault vector (added to [u1;u2])
u_fault = 0;                    % Disable VA meachanism
f_m_time = 8.5;                 % Sensor fault occurence time
%sim('threeDiskOscillatorRig');

%% Simulation for actuator fault (f_m = 0)
f_u = [0;-0.1];                 % Actuator fault vector (added to [u1;u2])
u_fault = 1;                    % Enable VA meachanism
f_m = [0;0;0];                  % Sensor fault vector (added to [y1;y2;y3])
sim('threeDiskOscillatorRigNew');
% threeDiskOscillatorRigNew.slx
%% Plot settings
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultAxesFontSize',20);
set(0,'DefaultLineLineWidth', 2);



