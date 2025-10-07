clc; clear;

% define epsilon and t final 
e = 1e-3;
tf = 10;

% calculate num of steps
nstep = tf/e;

% define constants
G = 3;
M = 10;
m = 1e-1;

% define V & F
V = @(x,y) - G*m*M * (x^2 + y^2)^(-1/2);
Fx = @(x,y) - G*m*M * x * (x^2 + y^2)^(-3/2);
Fy = @(x,y) - G*m*M * y * (x^2 + y^2)^(-3/2);

% define p & q
p(1,:) = [0 -0.5];
q(1,:) = [1 1];
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
RK2_p(1,:) = p(1,:);
RK2_q(1,:) = q(1,:);
RK2_H(1,:) = H(1,:);
RK2_t(1) = t(1);

for i=1:nstep
    % k1
    RK2_k1_q(i,:) = RK2_p(i,:)/m;
    RK2_k1_p(i,1) = Fx(RK2_q(i,1), RK2_q(i,2)); % Px
    RK2_k1_p(i,2) = Fy(RK2_q(i,1), RK2_q(i,2)); % Py

    % mid for k2
    RK2_q_mid(i,:) = RK2_q(i,:) + e * 0.5 * RK2_k1_q(i,:);
    RK2_p_mid(i,:) = RK2_p(i,:) + e * 0.5 * RK2_k1_p(i,:);

    % k2
    RK2_k2_q(i,:) = RK2_p_mid(i,:)/m;
    RK2_k2_p(i,1) = Fx(RK2_q_mid(i,1), RK2_q_mid(i,2)); % Px
    RK2_k2_p(i,2) = Fy(RK2_q_mid(i,1), RK2_q_mid(i,2)); % Py

    %
    RK2_p(i+1,:) = RK2_p(i,:) + e * RK2_k2_p(i,:);
    RK2_q(i+1,:) = RK2_q(i,:) + e * RK2_k2_q(i,:);
    
    RK2_H(i+1,:) = dot(RK2_p(i,:), RK2_p(i,:)) / (2*m) +...
                   V(RK2_q(i+1,1), RK2_q(i+1,2));
    RK2_t(i+1) = RK2_t(i) + e;
end

%% PLOTS
plot3(RK2_p(:,1), RK2_p(:,2), RK2_t, 'x-', p(:,1), p(:,2), t, '.-')
title('Momentum (p vs t)');
legend('RK2', 'EC');
xlabel('p_x');
ylabel('p_y');
zlabel('t');
pause;
plot(RK2_q(:,1), RK2_q(:,2), 'x-', q(:,1), q(:,2), '.-')
title('Position (q vs t)');
legend('RK2', 'EC');
xlabel('q_x');
ylabel('q_y');
zlabel('t');
pause;
plot3(RK2_p(:,1), RK2_q(:,1), RK2_t, '.', p(:,1), q(:,1), t, '.')
title('Phase space diagram(px vs qx)');
legend('RK2', 'EC');
xlabel('px');
ylabel('qy');
pause;
plot(RK2_t, RK2_H(:,1), 'x-', t, H(:,1), '.-')
title('Hamiltonian (H vs t)');
legend('RK2', 'EC');
xlabel('t');
ylabel('H');
pause;

%% SIMULATION
for i=1:nstep
    p = plot(RK2_q(i,1), RK2_q(i,2), 'o', 0,0,'o');
    p(2).MarkerSize = 10;
    xlim([-3 3]);
    ylim([-3 3]);
    pause(e);
end