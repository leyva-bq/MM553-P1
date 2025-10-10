clear;

% TIME STUFF
e = 1e-3;
tf = 5;
nstep = tf/e;
T = 0:e:tf;

% CONSTANTS
G = 3;
M = 30;
all_m = [1.1e-1 1.3e-1];
nplanets = 2;

% DEFINE P & Q
all_P(1,1,:) = [0.3 -0.8];
all_P(1,2,:) = [-0.1 0.9];
all_Q(1,1,:) = [1 1];
all_Q(1,2,:) = [-1 -1];

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
     all_EC_H(:,n,:)] = euler_cromer(P,Q,F,V,H,m,e,nstep);
    
    [all_LF_P(:,n,:),...
     all_LF_Q(:,n,:),...
     all_LF_H(:,n,:)] = leapfrog(P,Q,F,V,H,m,e,nstep);
    
    [all_RK2_P(:,n,:),...
     all_RK2_Q(:,n,:),...
     all_RK2_H(:,n,:)] = RK2(P,Q,F,V,H,m,e,nstep);
end

% SUM HAMILTONIANS WIP
% all_EC_H(:,n,) =
% all_LF_H
% all_RK2_H

%% PLOTS
hold on
for n=1:nplanets
    plot3(all_RK2_P(:,n,1), all_RK2_P(:,n,2), T, 'x-',...
          all_LF_P(:,n,1) , all_LF_P(:,n,2) , T, 'o-',...
          all_EC_P(:,n,1) , all_EC_P(:,n,2) , T, 's-');
    legend(['P' num2str(n) ' - RK2'], ...
           ['P' num2str(n) ' - LF'], ...
           ['P' num2str(n) ' - EC']);
    pause;
end
hold off
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

%% SIMULATION
trail = nstep;
for i=1:5:nstep
    if i <= trail
        range = 1:i;
    else
        range = i-trail:i;
    end

    p1_x = q(i,1,1);
    p1_y = q(i,1,2);
    p1_r = sqrt(p1_x.^2 + p1_y.^2);
    
    p2_x = q(i,2,1);
    p2_y = q(i,2,2);
    p2_r = sqrt(p2_x.^2 + p2_y.^2);

    p = plot(q(i,1,1), q(i,1,2), 'o', ...
             q(i,2,1), q(i,2,2), 'o', ...
             0,0,'o',...
             q(range,1,1), q(range,1,2), '-',...
             q(range,2,1), q(range,2,2), '-');

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
    text(p1_x, p1_y, ['    r = ', num2str(p1_r)]);
    text(p2_x, p2_y, ['    r = ', num2str(p2_r)]);
    xlim([-3 3]);
    ylim([-3 3]);
    pause(e);
end

%% KEPLER
kep = [];

for n=1:2
    planet_q = squeeze(q(:,n,:));
    planet_r = sqrt(planet_q(:,1).^2 + planet_q(:,2).^2);
    
    [max, I_max] = findpeaks(planet_r);
    [min, I_min] = findpeaks(-planet_r);
    
    perihelion = planet_q(I_min(1), :);
    aphelion = planet_q(I_max(1), :);
    V = perihelion - aphelion;
    a = sqrt(V * V') / 2;
    T = I_min(2) - I_min(1);
    
    kep = [kep; a^3/T^2];
    
end
bar(kep);
title('Kepler constant (a^3/T^2)');
text(1:length(kep), kep, num2str(kep), 'vert','bottom','horiz','center');