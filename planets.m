clc; clear;

% define epsilon and t final 
e = 1e-3;
tf = 3;

% calculate num of steps
nstep = tf/e;

% define constants
G = 3;
M = 30;
m = [1.3e-1 1.3e-1];

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
        if n==1
            H(i+1,:) = dot(pn,pn)/(2*m(n)) + V(qx, qy, m(n));
            t(i+1) = t(i) + e;
        else
            H(i+1,:) = H(i+1,:) + dot(pn,pn)/(2*m(n)) + V(qx, qy, m(n));
        end
    end
end

%% RK2
RK2_p(1,:,:) = p(1,:,:);
RK2_q(1,:,:) = q(1,:,:);
RK2_H(1,:) = H(1,:);
RK2_t(1) = t(1);

for n=1:2
    for i=1:nstep
        pn = squeeze(RK2_p(i,n,:))';
        px = pn(1); py = pn(2);
        qn = squeeze(RK2_q(i,n,:))';
        qx = qn(1); qy = qn(2);

        % k1
        RK2_k1_q(i,:) = pn/m(n);
        RK2_k1_p(i,1) = Fx(qx, qy, m(n)); % Px
        RK2_k1_p(i,2) = Fy(qx, qy, m(n)); % Py
    
        % mid for k2
        RK2_q_mid(i,:) = qn + e * 0.5 * RK2_k1_q(i,:);
        RK2_p_mid(i,:) = pn + e * 0.5 * RK2_k1_p(i,:);
    
        % k2
        RK2_k2_q(i,:) = RK2_p_mid(i,:)/m(n);
        RK2_k2_p(i,1) = Fx(RK2_q_mid(i,1), RK2_q_mid(i,2), m(n)); % Px
        RK2_k2_p(i,2) = Fy(RK2_q_mid(i,1), RK2_q_mid(i,2), m(n)); % Py
    
        %
        RK2_p(i+1,n,:) = pn + e * RK2_k2_p(i,:);
        RK2_q(i+1,n,:) = qn + e * RK2_k2_q(i,:);
        
        RK2_H(i+1,:) = RK2_H(i,:) + dot(pn,pn) / (2*m(n)) + V(qx, qy, m(n));
        if n==1
            RK2_t(i+1) = RK2_t(i) + e;
        end
    end
end

%% PLOTS
plot3(RK2_p(:,1,1), RK2_p(:,1,2), RK2_t, 'x-',...
      RK2_p(:,2,1), RK2_p(:,2,2), RK2_t, 'x-',...
      p(:,1,1), p(:,1,2), t, '.-',...
      p(:,2,1), p(:,2,2), t, '.-');
title('Momentum (p vs t)');
legend('P1 - RK2', 'P2 - RK2', ...
       'P1 - EC', 'P2 - EC');
xlabel('p_x');
ylabel('p_y');
zlabel('t');
pause;

plot3(RK2_q(:,1,1), RK2_q(:,1,2), RK2_t, 'x-',...
      RK2_q(:,2,1), RK2_q(:,2,2), RK2_t, 'x-',...
      p(:,1,1), p(:,1,2), t, '.-',...
      p(:,2,1), p(:,2,2), t, '.-');
title('Position (qx vs qy)');
legend('P1 - RK2', 'P2 - RK2', ...
       'P1 - EC', 'P2 - EC');
xlabel('q_x');
ylabel('q_y');
zlabel('t');
pause;

plot3(RK2_p(:,1,1), RK2_q(:,1,1), RK2_t, 'x-',...
      RK2_p(:,2,1), RK2_q(:,2,1), RK2_t, 'x-',...
      p(:,1,1), q(:,1,1), t, '.-',...
      p(:,2,1), q(:,2,1), t, '.-');
title('Phase space diagram(px vs qx)');
legend('P1 - RK2', 'P2 - RK2', ...
       'P1 - EC', 'P2 - EC');
xlabel('px');
ylabel('qy');
pause;

plot(RK2_t, RK2_H(:,1), 'x-', ...
     t, H(:,1), '.-');
title('Hamiltonian (H vs t)');
legend('RK2', 'EC');
xlabel('t');
ylabel('H');
pause;

%% SIMULATION
trail = 100;
for i=1:5:nstep
    p = plot(q(i,1,1), q(i,1,2), 'o', ...
             q(i,2,1), q(i,2,2), 'o', ...
             0,0,'o',...
             q(1:i,1,1), q(1:i,1,2), '-',...
             q(1:i,2,1), q(1:i,2,2), '-');

    % PLANET 1
    p(1).MarkerSize = 10;
    p(1).MarkerFaceColor = 'red';
    p(4).Color = [0.5 0 0];

    % PLANET 2
    p(2).MarkerSize = 10;
    p(2).MarkerFaceColor = 'blue';
    p(5).Color = [0 0 0.5];
    
    % SUN
    p(3).MarkerSize = 20;
    p(3).MarkerFaceColor = 'yellow';


    legend('Planet 1', 'Planet 2', 'Sun');
    xlim([-3 3]);
    ylim([-3 3]);
    pause(e);
end