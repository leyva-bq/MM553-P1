clear;

% TIME STUFF
e = 1e-3;
tf = 10;
nstep = tf/e;
T = 0:e:tf;

% CONSTANTS
n_particles = 500; % number of particles
length = 500; % length of string

% DEFINE P & Q
Q(1,:) = zeros(1,n_particles);
range = 1;
P(1,:) = (rand(1, n_particles) * 2 * range) - range;
P(1,1) = 0;
P(1,end) = 0;

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
[EC_P,EC_Q,EC_H] = euler_cromer(P,Q,F,V,H,1,e,nstep, n_particles);
[LF_P,LF_Q,LF_H] = leapfrog(P,Q,F,V,H,1,e,nstep, n_particles);
[RK2_P,RK2_Q,RK2_H] = RK2(P,Q,F,V,H,1,e,nstep, n_particles);
toc

%% ENERGY
plot(T, EC_H, 'o-', ...
     T, LF_H, 'x-', ...
     T, RK2_H, 's-');
title('Energy (T vs H)');
legend('EC', 'LF', 'RK2');
xlabel('T');
ylabel('H');

%% SIMULATION
for i=1:100:nstep
    string = 0:n_particles;
    scatter(RK2_Q(i,1:10) + string(1:10), i, 'o');
    xlim([-0.5 10.5]);
    pause(e);    
end

%% PLOTS
plot_i = 300;
plot(T, RK2_P(:,plot_i));
title('Momentum (t vs p)');
pause;

plot(T, RK2_Q(:,plot_i));
title('Position (t vs q)');
pause;

plot(RK2_P(:,plot_i), RK2_Q(:,plot_i));
title('Phase space diagram (p vs q)');
pause;

plot(T, RK2_H);
title('Energy (t vs H)');
pause;

%% 4. AVG VELOCITY^2
avg_vel_squared = sum(EC_P'.^2) / n_particles;
plot(T, avg_vel_squared, 'x-');
title('Average velocity squared over time (v^2(t) vs t)');
xlabel('t');
ylabel('v^2');

%% 5. HISTOGRAM
for t=1:1:nstep
    h = RK2_P(t,:).^2;
    histogram(h,15);
    % xlim([0 0.5]);
    %ylim([0 200]);
    title(['Histogram of velocity squared over time (t = ' num2str(t*e) ')'])
    pause(e);
end

%% MAXWELL
m = 1;
k = 1.38e-20;
temp = 100;
v = RK2_P(100,:);

for temp=0:1000:10e20
    h = sqrt(m / (2*pi*k*temp)) * exp(-m*v.^2 / (2*k*temp));
    scatter(v.^2,h);
    pause(e);
end