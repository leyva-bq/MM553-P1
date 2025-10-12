clear;

% TIME STUFF
e = 1e-2;
tf = 10;
nstep = tf/e;
T = 0:e:tf;

% CONSTANTS
n_particles = 500; % number of particles
length = 500; % length of string

% DEFINE P & Q
Q(1,:) = linspace(1,length,n_particles);
P(1,:) = (rand(1, n_particles) * 0.5) - 0.5;

% DEFINE V & F
V_h = @(X) [ diff(X,1,2) , 0 ];
V = @(Q) 1/2 * (abs(V_h(Q)).^2)...
       + 1/3 * (abs(V_h(Q)).^3)...
       + 1/4 * (abs(V_h(Q)).^4);

F_h = @(X) (abs(X) + abs(X).^2 + abs(X).^3) * sign(X);
F = @(Q) arrayfun(F_h, [ 0 , flip(diff(flip(Q),1,2)) ]) - ...
         arrayfun(F_h, [ -diff(Q,1,2) , 0 ]);

% DEFINE H
H(1,:) = 1/2 * P.^2 + V(Q);
tic
[EC_P,EC_Q,EC_H] = euler_cromer(P,Q,F,V,H,1,e,nstep, ...
                                [n_particles, length]);
[LF_P,LF_Q,LF_H] = leapfrog(P,Q,F,V,H,1,e,nstep, ...
                                [n_particles, length]);
[RK2_P,RK2_Q,RK2_H] = RK2(P,Q,F,V,H,1,e,nstep, ...
                                [n_particles, length]);
toc

plot(T, sum(EC_H,2), 'o-', ...
     T, sum(LF_H,2), 'x-', ...
     T, sum(RK2_H,2), 's-');

%% SIMULATION
% for i=1:100:nstep
%     scatter(EC_Q(i,:), 0, 'o');
%     xlim([0.5 length+0.5]);
%     pause(e);    
% end

%% PLOTS
plot_i = 2;
plot(t, P(:,plot_i));
title('Momentum (t vs p)');
pause;

plot(t, Q(:,plot_i));
title('Position (t vs q)');
pause;

plot(P(:,plot_i), Q(:,plot_i));
title('Phase space diagram (p vs q)');
pause;

plot(t, H);
title('Energy (t vs H)');
pause;

%% 4. AVG VELOCITY^2
avg_vel_squared = sum(EC_P'.^2) / n_particles;
plot(avg_vel_squared, 'x-');
title('Average velocity squared over time (v^2(t) vs t)');
xlabel('t');
ylabel('<v^2>');

%% 5. HISTOGRAM
for t=1:1:nstep
    h = RK2_P(t,:).^2;
    histogram(h);
    %xlim([0 0.5]);
    %ylim([0 200]);
    title(['Histogram of velocity squared over time (t = ' num2str(t) ')'])
    pause(0.01);
end