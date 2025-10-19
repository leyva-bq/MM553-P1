function [P,Q,H] = leapfrog(P,Q,F,V,H,M,e,nstep,NB)
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
    NB = 0 % Fixed Boundary Conditions (optional)
           % n, boundary
end

arguments (Output)
    P % Full MOMENTA matrix with nstep entries
    Q % Full POSITIONS / DoF matrix with nstep entries
    H % Full HAMILTONIAN matrix with nstep entries
end

% First half step
P(2,:) = P(1,:) + e/2 * F(Q(1,:));

% Check boundary
if NB
    P(2,1) = 0;
    P(2,NB(1)) = 0;
end

% Middle steps
for t=2:nstep
    Q(t,:) = Q(t-1,:) + e * P(t,:)/M;

    % Check boundary
    if NB
        Q(t,1) = 1;
        Q(t,NB(1)) = NB(2);
    end

    P(t+1,:) = P(t,:) + e * F(Q(t,:));

    % Check boundary
    if NB
        P(t+1,1) = 0;
        P(t+1,NB(1)) = 0;
    end

    H(t,:) = 1/2 * P(t+1,:).^2 / M + V(Q(t,:));
end

% Last half step
Q(nstep+1,:) = Q(nstep,:) + e * P(nstep+1,:)/M;

% Check boundary
if NB
    Q(nstep,1) = 1;
    Q(nstep,NB(1)) = NB(2);
end

% P(nstep+2) = P(nstep+1) + e/2 * F(Q(nstep+1)); % not needed, since
                                                 % nstep+2 exceeds range
H(nstep+1,:) = 1/2 * P(t+1,:).^2 / M + V(Q(nstep+1,:));

H = sum(H, 2);

end