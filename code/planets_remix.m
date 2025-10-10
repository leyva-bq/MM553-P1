clear;

% TIME STUFF
e = 1e-3;
tf = 5;
nstep = tf/e;
T = 0:e:tf;

% CONSTANTS
G = 3;
M = 30;
m = 1.1e-1;

% DEFINE V & F
V = @(q) - G*m*M * (q(1)^2 + q(2)^2)^(-1/2);
F = @(q) [(- G*m*M * q(1) * (q(1)^2 + q(2)^2)^(-3/2))...
          (- G*m*M * q(2) * (q(1)^2 + q(2)^2)^(-3/2))];

% DEFINE P & Q
P(1,:) = [0.3 -0.8];
Q(1,:) = [1 1];

% DEFINE HAMILTONIAN
H(1,:) = dot(P(1,:),P(1,:))/(2*m) + V(Q(1,:));

% INTEGRATION METHODS
[EC_P,EC_Q,EC_H] = euler_cromer(P,Q,F,V,H,m,e,nstep);
[LF_P,LF_Q,LF_H] = leapfrog(P,Q,F,V,H,m,e,nstep);
[RK2_P,RK2_Q,RK2_H] = RK2(P,Q,F,V,H,m,e,nstep);

%% PLOTS
plot3(RK2_P(:,1), RK2_P(:,2), T, 'x-',...
      LF_P(:,1) , LF_P(:,2) , T, 'o-',...
      EC_P(:,1) , EC_P(:,2) , T, 's-');
legend('RK2', 'LF', 'EC');
title('Momenta over time');
pause;

plot3(RK2_Q(:,1), RK2_Q(:,2), T, 'x-',...
      LF_Q(:,1) , LF_Q(:,2) , T, 'o-',...
      EC_Q(:,1) , EC_Q(:,2) , T, 's-');
legend('RK2', 'LF', 'EC');
title('Position over time');
pause;

plot(T, RK2_H(:,1), 'x-',...
     T, LF_H(:,1) , 'o-',...
     T, EC_H(:,1) , 's-');
legend('RK2', 'LF', 'EC');
title('Hamiltonian over time');
pause;
