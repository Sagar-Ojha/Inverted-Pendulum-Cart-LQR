clear; clc; close all;
%% Settings----------------------------------------------------------------
nonlinear = 1;  % 1 enables nonlinear dynamics
up = 1;         % 1 enables up stabilizing controller
animation = 1;   % 1 enables animation
anim_pause = 0.0005; % higher values make animation slower

%% ------------------------------------------------------------------------
% Initial Conditions
rng('default');
q0 = 0;%rand(1);

% Equilibrium States
q_eq = q0;      % Keep it within the xlims when running animation
x_eq_down = [q_eq; 0; 0; 0]; x_eq_up = [q_eq; 0; pi; 0];

qdot0 = 0;%rand(1);
theta0 = 140;
thetadot0 = 0;%rand(1);
error_q0 = -q0;   % Error between the positions of reference and cart

theta0 = deg2rad(theta0);
thetadot0 = deg2rad(thetadot0);

initial_velocity = [qdot0; thetadot0];
initial_position = [q0; theta0];
initial_condition = [q0; qdot0; theta0; thetadot0];

% Time parameters
t0 = 0;
tfinal = 80;
tstep = 0.005;

%% ------------------------------------------------------------------------
% System Parameters
cart_mass = 2; pendulum_mass = 1; pendulum_length = 0.3;
lin_spring_stiff = 0; rot_spring_stiff = 0; accel_gravity = 10;

%% ------------------------------------------------------------------------
% Linearized Dynamics
syms n m1 m2 k kappa l g q q_dot theta theta_dot f_ext

X = [q; q_dot; theta; theta_dot];

M = [m1+m2 m2*l*cos(theta);
     m2*l*cos(theta) m2*l^2];
G = [-k*q - m2*l*sin(theta)*theta_dot^2;
     m2*l*g*sin(theta) + kappa*theta];
F = [f_ext; 0];

M_inv_FG = M\(F - G);
X_dot = [X(2);
         M_inv_FG(1,:);
         X(4);
         M_inv_FG(2,:)];
A = jacobian(X_dot, X);
B = jacobian(X_dot, f_ext);

A_0 = subs(A, [m1, m2, l, k, kappa, g, q_dot, theta, theta_dot, f_ext], ...
              [cart_mass, pendulum_mass, pendulum_length, lin_spring_stiff, ...
               rot_spring_stiff, accel_gravity, 0,0,0,0]);
B_0 = subs(B, [m1, m2, l, k, kappa, g, q_dot, theta, theta_dot, f_ext], ...
              [cart_mass, pendulum_mass, pendulum_length, lin_spring_stiff, ...
               rot_spring_stiff, accel_gravity, 0,0,0,0]);

A_pi = subs(A, [m1, m2, l, k, kappa, g, q_dot, theta, theta_dot, f_ext], ...
               [cart_mass, pendulum_mass, pendulum_length, lin_spring_stiff, ...
                rot_spring_stiff, accel_gravity, 0,pi,0,0]);
B_pi = subs(B, [m1, m2, l, k, kappa, g, q_dot, theta, theta_dot, f_ext], ...
               [cart_mass, pendulum_mass, pendulum_length, lin_spring_stiff, ...
                rot_spring_stiff, accel_gravity, 0,pi,0,0]);

A_0 = double(A_0);
B_0 = double(B_0);
A_pi = double(A_pi);
B_pi = double(B_pi);

% Augmented Matrices for Integral Action
C = [1 0 0 0];
Aa_0  = [A_0; -C];   Aa_0 = cat(2, Aa_0, zeros(5,1));
Aa_pi = [A_pi; -C]; Aa_pi = cat(2, Aa_pi, zeros(5,1));
Ba_0  = [B_0; 0]; Ba_pi  = [B_pi; 0]; Br = [0;0;0;0;1];
D = 0;

%% ------------------------------------------------------------------------
R1 = 1*eye(5);
R1(5,5) = 1e-1;
R2 = 10^(1);
[Ka_pi_stabilize,~,~] = lqr(Aa_pi, Ba_pi,R1,R2);
% For faster cart position tracking, compute gains using place
% Ka_pi_stabilize = place(Aa_pi, Ba_pi, 7.5*[-0.1 -0.2 -0.3 -0.5 -0.25]);
fprintf('\nLQR-I Gain\n');
disp(Ka_pi_stabilize);
fprintf('\nClosed Loop Poles for Tracking and Upwards Stabilization\n');
disp(eig(Aa_pi - Ba_pi*Ka_pi_stabilize));

% Gain for Tracking and Destabilizing Up Position--------------------------
Ka_pi_destabilize = place(Aa_pi, Ba_pi, [-0.01 -0.02 0.05 0 -1.1]);
fprintf('\nUpward Destabilizing Gain\n');
disp(Ka_pi_destabilize);
fprintf('\nClosed Loop Poles for Tracking and Upwards Destabilization\n');
disp(eig(Aa_pi - Ba_pi*Ka_pi_destabilize));

% Gain for Tracking and Stabilizing Down Position--------------------------
[Ka_0_stabilize,~,~] = lqr(Aa_0, Ba_0,R1,R2);
% For faster cart position tracking, compute gains using place
% Ka_0_stabilize = place(Aa_0, Ba_0, [-1 -2 -3 -4 -5]);
fprintf('\nLQR-I Gain\n');
disp(Ka_0_stabilize);
fprintf('\nClosed Loop Poles for Tracking and Downwards Stabilization\n');
disp(eig(Aa_0 - Ba_0*Ka_0_stabilize));

% Gain for Tracking and Destabilizing Down Position------------------------
Ka_0_destabilize = place(Aa_0, Ba_0, [-1 -2 4 3 -2.25]);
fprintf('\nDownward Destabilizing Gain\n');
disp(Ka_0_destabilize);
fprintf('\nClosed Loop Poles for Downwards Destabilization\n');
disp(eig(Aa_0 - Ba_0*Ka_0_destabilize));

%% ------------------------------------------------------------------------
% % State space to TF
% Kx = -Ka_pi_stabilize(1:4);
% Kq = -Ka_pi_stabilize(5);
% [b,a] = ss2tf(A_pi+B_pi*Kx, B_pi*Kq,C,D);
% G = tf(b,a)
% pzmap(G)
% 
% % Root locus
% G = G/tf('s');
% rlocus(G)

%% ------------------------------------------------------------------------
% Observer Gains
V1 = 1*eye(4);
% V1(2,2) = 0.00001;
V2 = 10^(-2);
% C = [1 0 0 0];
[F,~,~] = lqr(A_0', -C', V1, V2);
F = F';
F = (place(A_0', -C', 5*[-10 -10.01 -40.02 -20.03]))';

%% ------------------------------------------------------------------------
% Tracking Cart Reference Position
Ka_pair = [Ka_pi_destabilize; Ka_0_stabilize];
if (up == 1)
    x_eq(3) = pi;
    Ka_pair(1,:) = Ka_0_destabilize;
    Ka_pair(2,:) = Ka_pi_stabilize;
end

%% ------------------------------------------------------------------------
% Simple Pendulum on a Cart with Integral Action
tspan = 0:tstep:tfinal;

if (nonlinear == 1)
[time, states] = ode45(@(time,states) ...
    pendulum_cart_tracking(time,states,Ka_pair,[x_eq_down;0],[x_eq_up;0],up, ...
    cart_mass, pendulum_mass, pendulum_length, lin_spring_stiff, ...
    rot_spring_stiff, accel_gravity), tspan, [initial_condition;error_q0]);
else
[time, states] = ode45(@(time,states) ...
    linear_pendulum_cart(time,states,Ka_pair,Aa_0,Ba_0,Aa_pi,Ba_pi,Br,up, ...
    [x_eq_up;0], [x_eq_down;0]), tspan, [initial_condition;error_q0]);
end

%% ------------------------------------------------------------------------
% Animation
if (animation == 1)
    animate_pendulum_cart(anim_pause, states, time, cart_mass, ...
                          pendulum_mass, pendulum_length);
end

% Plots
figure(2);
subplot(2,2,1);
plot(time, states(:,1));
xlabel('Time (s)');
ylabel('$q$ (m)', 'interpreter', 'latex');
title('Cart Position')

subplot(2,2,2);
plot(time, states(:,2), '--');
xlabel('Time (s)');
ylabel('$\dot{q}$ (m/s)', 'interpreter', 'latex');
title('Cart Velocity')

subplot(2,2,3);
plot(time, rad2deg(states(:,3)));
xlabel('Time (s)');
ylabel('$\theta$ (deg)', 'interpreter', 'latex');
title('Pendulum Position')

subplot(2,2,4);
plot(time, rad2deg(states(:,4)), '--');
xlabel('Time (s)');
ylabel('$\dot{\theta}$ (deg/s)', 'interpreter', 'latex');
title('Pendulum Velocity')