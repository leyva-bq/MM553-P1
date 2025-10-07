clc; clear;

% define epsilon and t final 
e = 1e-3;
t_final = 10;

% calculate num of steps
n_step = t_final/e;

% define N, mass, and initial positions
N = 20;
% M = ones(1,N);
X = sort(rand(1,N));
X(1) = 0;
X(N) = 0;

% define V and F
V = @(i, j) 1/2 * (abs(i - j)^2)...
          + 1/3 * (abs(i - j)^3)...
          + 1/4 * (abs(i - j)^4);
F = @(q) (abs(q) + abs(q)^2 + abs(q)^3) * sign(q);

% define p and q
P = (rand(1, N) * 2) - 1;
Q = X;
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

% hold on
% line([5 0], [0 0]);
% scatter(X, 0);
% hold off

% euler-cromer iteration
for i=1:n_step
    % get forces
    force = Q(i,:);
    for j=1:N
        if j == 1
            force(j) = -F(Q(i,j) - Q(i,j+1)); % only right neighbor
        elseif j == N
            force(j) = F(Q(i,j-1) - Q(i,j)); % only left neighbor
        else
            L = Q(i,j-1) - Q(i,j);
            R = Q(i,j) - Q(i,j+1);
            force(j) = F(L) - F(R);
        end
    end

    P(i+1,:) = P(i,:) + e .* force;
    Q(i+1,:) = Q(i,:) + e .* P(i+1,:);
    % Q(i+1,1) = 0; % energy is better when i dont restrain these...
    % Q(i+1,N) = 0;
    H(i+1,:) = 1/2 * sum(P(i+1,:).^2) + get_potentials(Q(i,:),V,N);
    t(i+1) = t(i) + e;
end

plot_i = 2;
plot(t, P(:,plot_i ));
pause;
plot(t, Q(:,plot_i ));
pause;
plot(P(:,plot_i ), Q(:,plot_i));
pause;
plot(t, H);
pause;

for i=1:n_step
    %hold on
    %line([1.5 -1.5], [0 0], 'Color','r');
    scatter(Q(i,:), 0);
    xlim([-1.5 1.5]);
    pause(0.01);
    %hold off
end
%%
RK2_p(1) = p(1);
RK2_q(1) = q(1);
RK2_H(1) = H(1);
RK2_t(1) = t(1);

%% RK2
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

%% RK4
RK4_p(1) = p(1);
RK4_q(1) = q(1);
RK4_H(1) = H(1);
RK4_t(1) = t(1);

% RK4
for i=1:n_step-1
    % for k1
    RK4_k1_q = RK4_p(i)/m;
    RK4_k1_p = F(RK4_q(i));

    % for k2
    q_int = RK4_q(i) + e/2 * RK4_k1_q;
    p_int = RK4_p(i) + e/2 * RK4_k1_p;

    RK4_k2_q = p_int/m;
    RK4_k2_p = F(q_int);

    % for k3
    q_int = RK4_q(i) + e/2 * RK4_k1_q;
    p_int = RK4_p(i) + e/2 * RK4_k1_p;

    RK4_k3_q = p_int/m;
    RK4_k3_p = F(q_int);

    % for k4
    q_int = RK4_q(i) + e * RK4_k1_q;
    p_int = RK4_p(i) + e * RK4_k1_p;

    RK4_k4_q = p_int/m;
    RK4_k4_p = F(q_int);

    % last step
    RK4_q(i+1) = RK4_q(i) + e/6 * (RK4_k1_q + 2*RK4_k2_q + 2*RK4_k3_q + RK4_k4_q);
    RK4_p(i+1) = RK4_p(i) + e/6 * (RK4_k1_p + 2*RK4_k2_p + 2*RK4_k3_p + RK4_k4_p);
    
    RK4_t(i+1) = RK4_t(i) + e;
    RK4_H(i+1) = RK4_p(i+1) * RK4_p(i+1) / (2*m) + V(RK4_q(i+1));
end

%% PLOTS
plot(RK4_t, RK4_p, 'o', RK2_t, RK2_p, 'x', t, p, '.');
pause;
plot(RK4_t, RK4_q, 'o', RK2_t, RK2_q, 'x', t, q, '.');
pause;
plot(RK4_p, RK4_q, 'o', RK2_p, RK2_q, 'x', p, q, '.');
pause;
plot(RK4_t, RK4_H, 'o', RK2_t, RK2_H, 'x', t, H, '.');
pause;