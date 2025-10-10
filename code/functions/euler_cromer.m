function [P,Q,H] = euler_cromer(P,Q,F,V,H,M,e,nstep)
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
end

arguments (Output)
    P % Full MOMENTA matrix with nstep entries
    Q % Full POSITIONS / DoF matrix with nstep entries
    H % Full HAMILTONIAN matrix with nstep entries
end

for t=1:nstep
    P(t+1,:) = P(t,:) + e * F(Q(t,:));
    Q(t+1,:) = Q(t,:) + e * P(t+1,:)/M; % Change P(i+1) to P(i)
                                     % for EC to Euler, respectively
    H(t+1,:) = 1/2 * dot(P(t+1,:),P(t+1,:)) / M + V(Q(t+1,:));
end

end