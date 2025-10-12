function [P,Q,H] = euler_cromer(P,Q,F,V,H,M,e,nstep,NB)
%EULER_CROMER Euler-Cromer integration
%   A function for EC integration.
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
    P(t+1,:) = P(t,:) + e * F(Q(t,:));
    if NB
        P(t+1,1) = 0;
        P(t+1,NB(1)) = 0;
    end
    Q(t+1,:) = Q(t,:) + e * P(t+1,:)/M; % Change P(i+1) to P(i)
                                     % for EC to Euler, respectively
    if NB
        Q(t+1,1) = 1;
        Q(t+1,NB(1)) = NB(2);
    end
    
    H(t+1,:) = 1/2 * dot(P(t+1,:),P(t+1,:)) / M + V(Q(t+1,:));
end

H = sum(H, 2);

end