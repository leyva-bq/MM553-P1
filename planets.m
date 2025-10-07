clc; clear;

% define epsilon and t final 
e = 1e-2;
tf = 10;

% calculate num of steps
nstep = tf/e;

% define constants
G = 1;
M = 1;
m = 0.1;

% define V & F
V = @(x,y) - G*m*M * (x^2 + y^2)^(-1/2);
Fx = @(x,y) - G*m*M * x * (x^2 + y^2)^(-3/2);
Fy = @(x,y) - G*m*M * y * (x^2 + y^2)^(-3/2);

% define p & q
p(1,:) = [1 -1];
q(1,:) = [10 11];
t(1) = 0;

% define Hamiltonian
H(1) = dot(p,p)/(2*m) + V(q(1,1), q(1,2));

%% euler-cromer iteration
for i=1:nstep
    p(i+1,1) = p(i,1) + e * Fx(q(i,1), q(i,2)); % Px
    p(i+1,2) = p(i,2) + e * Fy(q(i,1), q(i,2)); % Py
    q(i+1,:) = q(i,:) + e * p(i+1,:)/m; % Q
    H(i+1,:) = dot(p(i,:),p(i,:))/(2*m) + V(q(i,1), q(i,2));
    t(i+1) = t(i) + e;
end

%% RK2
RK2_p(1) = p(1);
RK2_q(1) = q(1);
RK2_H(1) = H(1);
RK2_t(1) = t(1);

for i=1:nstep
    RK2_k1_q = RK2_p(i)/m;
    RK2_k1_p = F(RK2_q(i));

    RK2_q_mid = RK2_q(i) + e * 0.5 * RK2_k1_q;
    RK2_p_mid = RK2_p(i) + e * 0.5 * RK2_k1_p;

    RK2_k2_q = RK2_p_mid / m;
    RK2_k2_p = F(RK2_q_mid);

    RK2_p(i+1) = RK2_p(i) + e * RK2_k2_p;
    RK2_q(i+1) = RK2_q(i) + e * RK2_k2_q;
    
    RK2_H(i+1) = RK2_p(i+1) * RK2_p(i+1) / (2*m) + V(RK2_q(i+1));
    RK2_t(i+1) = RK2_t(i) + e;
end

%% PLOTS
plot(RK2_t, RK2_p, 'x', t, p, '.');
title('Momentum (p vs t)')
pause;
plot(RK2_t, RK2_q, 'x', t, q, '.');
title('Position (q vs t)')
pause;
plot(RK2_p, RK2_q, 'x', p, q, '.');
title('Phase diagram(p vs q)')
pause;
plot(RK2_t, RK2_H, 'x', t, H, '.');
title('Hamiltonian (H vs t)')
pause;

%% SIMULATION
for i=1:nstep
    %hold on
    scatter(q(i,1), q(i,2));
    xlim([-15 15]);
    ylim([-15 15]);
    pause(0.01);
    %hold off
end