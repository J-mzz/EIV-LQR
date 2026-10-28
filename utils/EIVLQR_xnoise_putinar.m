function out = EIVLQR_xnoise_putinar(sim, d, T, sysd, xk, Yk, UU)
% EIVLQR_XNOISE_PUTINAR  Input-constrained EIV LQR SDP with Putinar's Psatz
% (scalarized formulation), when only the state measurements are noisy.
%   out = EIVLQR_xnoise_putinar(sim, d, T, sysd, xk, Yk, UU) has the same
%   inputs as EIVLQR_xnoise_scherer. The matrix condition F(A,B) >= 0 is
%   replaced by z'*F(A,B)*z >= 0 on the sphere z'*z = 1, which allows a
%   scalar Putinar certificate; the relaxation order d must be at least 2.
%
%   xk    current state;  Yk  previous Y (Y >= Yk is enforced)
%   UU    input constraint ||UU*u||_2 <= 1 at xk

% extract data and size
X = sim.X_noise(:,1:T);     % state data
U = sim.U(:,1:T-1);         % input data
eps = sim.epsilon;          % noise bound epsilon
tol = sim.tolerance;        % delta in paper
n = size(X,1);              % dim of state
m = size(U,1);              % dim of input

% extract performance matrix C,D,E
r = sysd.r;
E = sysd.E;
C = sysd.C;
D1 = sysd.D1;       % D used in H2

% system parameters as variables
A = sdpvar(n,n,'full');
B = sdpvar(n,m,'full');

% auxiliary variable z for the scalarization
s = 2*n;
z = sdpvar(s,1,'full');
vars = recover(depends([A(:);B(:);z(:)]));
vars = vars(:).';

% constraints
g = eps*ones(2*n*T,1);            % eq.(15) inequality constraint
for i = 1:T-1                     % eq.(15) equality constraint
    h(1+n*(i-1):n*i,:) = X(:,i+1) - A*X(:,i) - B*U(:,i);
end

hz = 1 - z'*z;
cons_data = struct('ineq', g, 'eq', h, 'zeq', hz);

% define sdp F
Y = sdpvar(n);              % P
S = sdpvar(m,n,'full');     % M
Z = sdpvar(r);              % W
eta = sdpvar(1);
gamma = sdpvar(1);
C1 = [Y-E*E' A*Y+B*S; Y*A'+S'*B' Y];
C2 = [Z C*Y+D1*S; Y*C'+S'*D1' Y];
C3 = gamma-trace(Z);
F = C1-eta*eye(size(C1));
Cons = [(Y >= tol):'Y'; (Z >= tol):'Z'; (eta >= tol):'eta'; (C2 >= tol):'C2'; (C3 >= tol):'C3'];

% input constraint ||UU*u||_2 <= 1 at the current state, and Y >= Yk
beta = sdpvar(1);
C4 = [beta*eye(m), UU*S; (UU*S)', Y];
C5 = [beta, beta*xk'; (beta*xk')', Y];
Cons = [Cons; (C4 >= tol):'C4'; (C5 >= tol):'C5'; (beta >= tol):'beta'; (Y >= Yk):'Yk'];

out.C4 = C4;
out.C5 = C5;
out.beta = beta;

% scalarized condition
q = z' * F * z;

[p_psatz, cons_psatz, Gram, Coef] = psatz_xnoise_putinar(q, cons_data, d, vars, A, n, T, B);

Cons = [Cons; cons_psatz];

out.poly = p_psatz;
out.cons = Cons;
out.Gram = Gram;
out.Coef = Coef;
out.Y = Y;
out.S = S;
out.obj = gamma;
out.vars = vars;

end
