clear;

% TIME STUFF
e = 1e-3;
tf = 5;
nstep = tf/e;
T = 0:e:tf;

% CONSTANTS
G = 3;
M = 30;
all_m = [1.1e-1 1.3e-1 1.2e-1];
nplanets = 3;

% DEFINE P & Q
all_P(1,1,:) = [0.3 -0.8];
all_P(1,2,:) = [-0.1 0.9];
all_P(1,3,:) = [-0.5 -0.9];

all_Q(1,1,:) = [1 1];
all_Q(1,2,:) = [-1 -1];
all_Q(1,3,:) = [1 0];

tic
for n=1:nplanets
    P = squeeze(all_P(:,n,:))';
    Q = squeeze(all_Q(:,n,:))';
    m = all_m(n);
    
    % DEFINE V & F
    V = @(q) - G*m*M * (q(1)^2 + q(2)^2)^(-1/2);
    F = @(q) [(- G*m*M * q(1) * (q(1)^2 + q(2)^2)^(-3/2))...
              (- G*m*M * q(2) * (q(1)^2 + q(2)^2)^(-3/2))];
    
    % DEFINE HAMILTONIAN
    H(1,:) = dot(P(1,:),P(1,:))/(2*m) + V(Q(1,:));
    
    % INTEGRATION METHODS
    [all_EC_P(:,n,:),...
     all_EC_Q(:,n,:),...
     all_EC_H(:,n)] = euler_cromer(P,Q,F,V,H,m,e,nstep);
    
    [all_LF_P(:,n,:),...
     all_LF_Q(:,n,:),...
     all_LF_H(:,n)] = leapfrog(P,Q,F,V,H,m,e,nstep);
    
    [all_RK2_P(:,n,:),...
     all_RK2_Q(:,n,:),...
     all_RK2_H(:,n)] = RK2(P,Q,F,V,H,m,e,nstep);
end
toc

%% REAL PLOTS
plot3(all_RK2_P(:,1,1), all_RK2_P(:,1,2), T, '^-',...
      all_LF_P(:,1,1) , all_LF_P(:,1,2) , T, 'x-',...
      all_EC_P(:,1,1) , all_EC_P(:,1,2) , T, 'o-',...
      all_RK2_P(:,2,1), all_RK2_P(:,2,2), T, '^-',...
      all_LF_P(:,2,1) , all_LF_P(:,2,2) , T, 'x-',...
      all_EC_P(:,2,1) , all_EC_P(:,2,2) , T, 'o-',...
      all_RK2_P(:,3,1), all_RK2_P(:,3,2), T, '^-',...
      all_LF_P(:,3,1) , all_LF_P(:,3,2) , T, 'x-',...
      all_EC_P(:,3,1) , all_EC_P(:,3,2) , T, 'o-');
title('Momentum (P_x vs P_y vs T)');
legend('P1 - RK2', 'P1 - LF', 'P1 - EC', ...
       'P2 - RK2', 'P2 - LF', 'P2 - EC', ...
       'P3 - RK2', 'P3 - LF', 'P3 - EC');
xlabel('p_x');
ylabel('p_y');
zlabel('t');
pause;

plot3(all_RK2_Q(:,1,1), all_RK2_Q(:,1,2), T, '^-',...
      all_LF_Q(:,1,1) , all_LF_Q(:,1,2) , T, 'x-',...
      all_EC_Q(:,1,1) , all_EC_Q(:,1,2) , T, 'o-',...
      all_RK2_Q(:,2,1), all_RK2_Q(:,2,2), T, '^-',...
      all_LF_Q(:,2,1) , all_LF_Q(:,2,2) , T, 'x-',...
      all_EC_Q(:,2,1) , all_EC_Q(:,2,2) , T, 'o-',...
      all_RK2_Q(:,3,1), all_RK2_Q(:,3,2), T, '^-',...
      all_LF_Q(:,3,1) , all_LF_Q(:,3,2) , T, 'x-',...
      all_EC_Q(:,3,1) , all_EC_Q(:,3,2) , T, 'o-');
title('Position (Q_x vs Q_y vs T)');
legend('P1 - RK2', 'P1 - LF', 'P1 - EC', ...
       'P2 - RK2', 'P2 - LF', 'P2 - EC', ...
       'P3 - RK2', 'P3 - LF', 'P3 - EC');
xlabel('q_x');
ylabel('q_y');
zlabel('t');
pause;

plot(T, all_RK2_H, '^-',...
     T, all_LF_H, 's-',...
     T, all_EC_H, 'o-');
title('Energy over time');
legend('P1 - RK2', 'P1 - LF', 'P1 - EC', ...
       'P2 - RK2', 'P2 - LF', 'P2 - EC', ...
       'P3 - RK2', 'P3 - LF', 'P3 - EC');
xlabel('t');
ylabel('H');
pause;

%% SIMULATION
trail = nstep;
for t=1:5:nstep
    if t <= trail
        range = 1:t;
    else
        range = t-trail:t;
    end

    X = squeeze(all_RK2_Q(:,:,1));
    Y = squeeze(all_RK2_Q(:,:,2));
    trail_X = squeeze(all_RK2_Q(1:t,:,1));
    trail_Y = squeeze(all_RK2_Q(1:t,:,2));
    
    p = plot(0, 0, 'o', ...
             X(t, 1), Y(t, 1), 'o', ...
             X(t, 2), Y(t, 2), 'o', ...
             X(t, 3), Y(t, 3), 'o', ...
             X(range, 1), Y(range, 1), '-', ...
             X(range, 2), Y(range, 2), '-', ...
             X(range, 3), Y(range, 3), '-', ...
             MarkerSize=10);

    % SUN
    p(1).MarkerSize = 20;
    p(1).MarkerFaceColor = 'yellow';
    p(1).MarkerEdgeColor = 'red';

    % PLANETS
    p(2).MarkerFaceColor = 'green';
    p(3).MarkerFaceColor = 'cyan';
    p(4).MarkerFaceColor = 'red';
    
    % MORE STUFF
    legend('Sun', 'Planet 1', 'Planet 2', 'Planet 3');
    title(['System of planets (t = ' num2str(t*e) ')']);
    xlim([-3 3]);
    ylim([-3 3]);
    pause(e);
end

%% KEPLER
kep = [];

for n=1:nplanets
    planet_q = squeeze(all_RK2_Q(:,n,:));
    planet_r = sqrt(planet_q(:,1).^2 + planet_q(:,2).^2);
    
    [max, I_max] = findpeaks(planet_r);
    [min, I_min] = findpeaks(-planet_r);
    
    perihelion = planet_q(I_min(1), :);
    aphelion = planet_q(I_max(1), :);
    V = perihelion - aphelion;
    a = sqrt(V * V') / 2;
    T_orbit = I_min(2) - I_min(1);
    
    kep = [kep; a^3/T_orbit^2];
    
end
bar(kep);
title('Kepler constant (a^3/T^2)');
text(1:length(kep), kep, num2str(kep), 'vert','bottom','horiz','center');