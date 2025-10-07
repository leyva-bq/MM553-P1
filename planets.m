clc; clear;

% define epsilon and t final 
e = 1e-2;
tf = 10;

% calculate num of steps
nstep = tf/e;

% define mass and k constant
m = 0.35;
k = 0.4;

% define V & F
V = @(x) 0.5*k*x^2;
F = @(x) -k*x;

% define p & q
p(1) = 0.2;
q(1) = 0.5;
t(1) = 0;

% define Hamiltonian
H(1) = p(1) * p(1)/(2*m) + V(q(1));

%% euler-cromer iteration
for i=1:nstep
    p(i+1) = p(i) - k*q(i) * e;
    q(i+1) = q(i) + e * p(i+1)/m; % this p(i+1) change is euler vs euler-cromer
    H(i+1) = p(i+1) * p(i+1) / (2*m) + V(q(i+1));
    t(i+1) = t(i) + e;
end

RK2_p(1) = p(1);
RK2_q(1) = q(1);
RK2_H(1) = H(1);
RK2_t(1) = t(1);

%% RK2
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
pause;
plot(RK2_t, RK2_q, 'x', t, q, '.');
pause;
plot(RK2_p, RK2_q, 'x', p, q, '.');
pause;
plot(RK2_t, RK2_H, 'x', t, H, '.');
pause;