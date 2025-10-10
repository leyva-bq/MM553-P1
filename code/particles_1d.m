clc; clear;

% define epsilon and t final 
e = 1e-2;
t_final = 50;

% calculate num of steps
n_step = t_final/e;

N = 500; % number of particles
length = 10; % length of string
Q = linspace(1,length,N); % initial positions

% define V and F
V = @(i, j) 1/2 * (abs(i - j)^2)...
          + 1/3 * (abs(i - j)^3)...
          + 1/4 * (abs(i - j)^4);
F = @(q) (abs(q) + abs(q)^2 + abs(q)^3) * sign(q);

P = (rand(1, N) * 0.5) - 0.5; % define momenta
t(1) = 0;

% calculate potential
function potential = get_potentials(Q, V, N)
    potential = 0;
    for i=1:N-1
        potential = potential + V(Q(i), Q(i+1));
    end
end

% define Hamiltonian
H = 1/2 * sum(P.^2) + get_potentials(Q, V, N);

%% euler-cromer iteration
for i=1:n_step
    % get forces
    force = zeros(1,N);
    for j=2:N-1
        if j == 1 
            force(j) = 0; % right neighbor
        elseif j == N
            force(j) = 0; % left neighbor
        else
            L = Q(i,j-1) - Q(i,j);
            R = Q(i,j) - Q(i,j+1);
            force(j) = F(L) - F(R);
        end
    end

    P(i+1,:) = P(i,:) + e .* force;
    Q(i+1,:) = Q(i,:) + e .* P(i+1,:);
    Q(i+1,1) = 1;
    Q(i+1,N) = length;
    H(i+1,:) = 1/2 * sum(P(i+1,:).^2) + get_potentials(Q(i,:),V,N);
    t(i+1) = t(i) + e;
end

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

%% RK2
RK2_p(1) = p(1);
RK2_q(1) = q(1);
RK2_H(1) = H(1);
RK2_t(1) = t(1);

for i=1:n_step
    RK2_k1_q = RK2_p(i)/m;
    RK2_k1_p = F(RK2_q(i));

    RK2_q_mid = RK2_q(i) + e * 0.5 * RK2_k1_q;
    RK2_p_mid = RK2_p(i) + e * 0.5 * RK2_k1_p;

    RK2_k2_q = RK2_p_mid / m;
    RK2_k2_p = F(RK2_q_mid);

    RK2_p(i+1) = RK2_p(i) + e * RK2_k2_p;
    RK2_q(i+1) = RK2_q(i) + e * RK2_k2_q;
    
    RK2_H(i+1) = RK2_p(i+1) * RK2_p(i+1) / (2*m) + V(RK2_q(i+1));
    RK2_t(i+1) = RK2_t(i) + e;
end

%% SIMULATION
for i=1:10:n_step
    %hold on
    %line([1.5 -1.5], [0 0], 'Color','r');
    scatter(Q(i,:), 0, 'o');
    xlim([0.5 100+0.5]);
    pause;
    %hold off
end

%% 4. AVG VELOCITY^2
avg_vel_squared = sum(P'.^2) / N;
plot(avg_vel_squared, 'x-');
title('Average velocity squared over time (v^2(t) vs t)');
xlabel('t');
ylabel('<v^2>');

%% 5. HISTOGRAM
for t=1:5:n_step
    h = P(t,:).^2;
    histogram(h);
    %xlim([0 0.5]);
    %ylim([0 200]);
    title(['Histogram of velocity squared over time (t = ' num2str(t) ')'])
    pause(0.01);
end