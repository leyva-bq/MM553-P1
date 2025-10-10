clear;

% define epsilon and t final 
e = 1e-3;
tf = 5;
nstep = tf/e;
T = 0:e:tf;

% define constants
G = 3;
M = 30;
m = 1.1e-1;

% define V & F
V = @(q) - G*m*M * (q(1)^2 + q(2)^2)^(-1/2);
%Fx = @(q) - G*m*M * q(1) * (q(1)^2 + q(2)^2)^(-3/2);
%Fy = @(q) - G*m*M * q(2) * (q(2)^2 + q(2)^2)^(-3/2);
F = @(q) [(- G*m*M * q(1) * (q(1)^2 + q(2)^2)^(-3/2))...
          (- G*m*M * q(2) * (q(1)^2 + q(2)^2)^(-3/2))];

% define p(t, n, p_n) & q(t, n, q_n)
P(1,:) = [0.3 -0.8];
Q(1,:) = [1 1];

% define Hamiltonian
H(1,:) = dot(P(1,:),P(1,:))/(2*m) + V(Q(1,:));

[EC_P,EC_Q,EC_H] = euler_cromer(P,Q,F,V,H,m,e,nstep);

%% PLOTS
plot3(EC_P(:,1), EC_P(:,2), T, 'x-');
pause;

plot3(EC_Q(:,1), EC_Q(:,2), T, 'x-');
pause;

plot(EC_H);
