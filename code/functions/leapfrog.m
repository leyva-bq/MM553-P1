function [P,Q,H] = leapfrog(P,Q,F,V,H,M,e,nstep)
%LEAPFROG Leapfrog integration
%   A function for LF integration.
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

% First half step
P(2,:) = P(1,:) + e/2 * F(Q(1,:));

% Middle steps
for t=2:nstep
    Q(t,:) = Q(t-1,:) + e * P(t,:)/M;
    P(t+1,:) = P(t,:) + e * F(Q(t,:));
    H(t,:) = 1/2 * dot(P(t+1,:),P(t+1,:)) / M + V(Q(t,:));
end

% Last half step
Q(nstep+1,:) = Q(nstep,:) + e * P(nstep+1,:)/M;
%P(nstep+2) = P(nstep+1) + e/2 * F(Q(nstep+1)); % useless half step?
H(nstep+1,:) = 1/2 * dot(P(nstep+1,:), P(nstep+1,:)) / M + V(Q(nstep+1,:));

end