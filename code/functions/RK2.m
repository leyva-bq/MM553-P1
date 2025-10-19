function [P,Q,H] = RK2(P,Q,F,V,H,M,e,nstep,NB)
%RK2 RK2 integration
%   A function for RK2 integration.
arguments (Input)
    P % MOMENTA matrix
    Q % POSITIONS / Degrees of freedom matrix
    F % FORCE function
    V % POTENTIAL force function
    H % HAMILTONIAN matrix
    M % MASS matrix
    e % EPSILON time step
    nstep % Number of steps
    NB = 0 % Fixed Boundary Conditions (optional)
             % n, boundary
end

arguments (Output)
    P % Full MOMENTA matrix with nstep entries
    Q % Full POSITIONS / DoF matrix with nstep entries
    H % Full HAMILTONIAN matrix with nstep entries
end

for t=1:nstep
    k1_Q = P(t,:)/M;
    k1_P = F(Q(t,:));

    Q_mid = Q(t,:) + e/2 * k1_Q;
    P_mid = P(t,:) + e/2 * k1_P;

    k2_Q = P_mid / M;
    k2_P = F(Q_mid);

    P(t+1,:) = P(t,:) + e * k2_P;
    
    % Check boundary
    if NB
        P(t+1,1) = 0;
        P(t+1,NB(1)) = 0;
    end
    
    Q(t+1,:) = Q(t,:) + e * k2_Q;
    
    % Check boundary
    if NB
        Q(t+1,1) = 0;
        Q(t+1,NB(1)) = 0;
    end
    
    H(t+1,:) = 1/2 * P(t+1,:).^2 / M + V(Q(t+1,:));
end

H = sum(H, 2);

end