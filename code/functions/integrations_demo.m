clear; clc;

% specific problem
k = 0.4; % spring constant
F = @(q) -k * q;
V = @(q) 0.5 * k * q^2;

% initial conditions
P(1) = 0.2;
Q(1) = 0.5;
M = 0.35;
H(1) = 1/2 * P(1).^2 / M + V(Q(1));

% time stuff
e = 1e-3;
tf = 100;
nstep = tf/e;
T = 0:e:tf;

[EC_P,EC_Q,EC_H] = euler_cromer(P,Q,F,V,H,M,e,nstep);
[LF_P,LF_Q,LF_H] = leapfrog(P,Q,F,V,H,M,e,nstep);
[RK2_P,RK2_Q,RK2_H] = RK2(P,Q,F,V,H,M,e,nstep);

%% PLOTS
plot(T, RK2_P, 'o-', ...
     T, LF_P, 'x-', ...
     T, EC_P, '.-');
legend('RK2', 'LF', 'EC');
title('Momenta over time');
xlabel('t');
ylabel('p');
pause;

plot(T, RK2_Q, 'o-', ...
     T, LF_Q, 'x-', ...
     T, EC_Q, '.-');
legend('RK2', 'LF', 'EC');
title('Position over time');
xlabel('t');
ylabel('q');
pause;

plot(RK2_P, RK2_Q, 'o-', ...
     LF_P, LF_Q, 'x-', ...
     EC_P, EC_Q, '.-');
legend('RK2', 'LF', 'EC');
title('Phase space diagram');
xlabel('t');
ylabel('q');
pause;

plot(T, RK2_H, 'o-', ...
     T, LF_H, '.-', ...
     T, EC_H, '.-');
legend('RK2', 'LF', 'EC');
title('Energy over time');
xlabel('t');
ylabel('H');
pause;