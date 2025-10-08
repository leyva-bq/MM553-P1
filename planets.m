clc; clear;

% define epsilon and t final 
e = 1e-3;
tf = 3;

% calculate num of steps
nstep = tf/e;

% define constants
G = 3;
M = 10;
m = [1e-1 0.5e-1];

% define V & F
V = @(x,y,m) - G*m*M * (x^2 + y^2)^(-1/2);
Fx = @(x,y,m) - G*m*M * x * (x^2 + y^2)^(-3/2);
Fy = @(x,y,m) - G*m*M * y * (x^2 + y^2)^(-3/2);

% define p(t, n, p_n) & q(t, n, q_n)
p(1,1,:) = [0 -1];
p(1,2,:) = [0 1];
q(1,1,:) = [1 1];
q(1,2,:) = [-1 -1];
t(1) = 0;

% define Hamiltonian
p1 = squeeze(p(1,1,:))';
p2 = squeeze(p(1,2,:))';
q1 = squeeze(q(1,1,:))';
q2 = squeeze(q(1,2,:))';
H(1,:) = dot(p1,p1)/(2*m(1)) + V(q1(1), q1(2), m(1)) +...
         dot(p2,p2)/(2*m(2)) + V(q2(1), q2(2), m(2));

%% euler-cromer iteration
for n=1:2
    for i=1:nstep
        pn = squeeze(p(i,n,:))';
        px = pn(1); py = pn(2);
        qn = squeeze(q(i,n,:))';
        qx = qn(1); qy = qn(2);

        p(i+1,n,1) = px + e * Fx(qx, qy, m(n)); % Px
        p(i+1,n,2) = py + e * Fy(qx, qy, m(n)); % Py
        pn1 = squeeze(p(i+1,n,:))';
        q(i+1,n,:) = qn + e * pn1/m(n); % Q
        H(i+1,:) = H(i,:) + dot(pn,pn)/(2*m(n)) + V(qx, qy, m(n));
        if n==1
            t(i+1) = t(i) + e;
        end
    end
end

%% RK2
RK2_p(1,:,:) = p(1,:,:);
RK2_q(1,:,:) = q(1,:,:);
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
plot3(RK2_p(:,1), RK2_p(:,2), RK2_t, 'x-',...
      p(:,1,1), p(:,1,2), t, '.-',...
      p(:,2,1), p(:,2,2), t, '.-');
title('Momentum (p vs t)');
legend('RK2', 'EC');
xlabel('p_x');
ylabel('p_y');
zlabel('t');
pause;
plot(RK2_q(:,1), RK2_q(:,2), 'x-', q(:,1), q(:,2), '.-')
title('Position (qx vs qy)');
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