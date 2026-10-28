function out = EIVLQR_linf(sim, d, T, sysd, xk, Yk, UU)
% EIVLQR_LINF  EIV data-driven LQR SDP for l_inf-bounded noise (Scherer, (10)).
%   out = EIVLQR_linf(sim, d, T, sysd, xk, Yk, UU) builds the SDP for a state
%   feedback u = S*inv(Y)*x that minimizes the worst-case LQR cost bound
%   gamma over all (A, B) consistent with the first T noisy samples in sim,
%   using the relaxation order d.
%
%   sim   struct with fields X_noise (n x T measured states), U_noise
%         (m x T-1 measured inputs), epsilon = [eps_x; eps_u; eps_w]
%         (l_inf bounds of the state and input measurement errors and of the
%         process noise), and tolerance
%   sysd  performance matrices from generate_sys
%   xk    current state;  Yk  previous Y (Y >= Yk is enforced)
%   UU    input constraint ||UU*u||_2 <= 1 at xk; with UU = [] (or omitted)
%         there is no input constraint, which gives the unconstrained EIV
%         controller
%
%   Solve with optimize(out.cons, out.obj, opts); then K = value(out.S)/value(out.Y).

% extract data and size
X = sim.X_noise(:,1:T);     % state data
U = sim.U_noise(:,1:T-1);   % input data
eps = sim.epsilon;          % noise bound epsilon
tol = sim.tolerance;        % delta in paper
n = size(X,1);              % dim of state
m = size(U,1);              % dim of input

% extract performance matrix C,D,E
r = sysd.r;
E = sysd.E;
C = sysd.C;
D1 = sysd.D1;       % D1 used in H2

% system parameters as variables
A = sdpvar(n,n,'full');
B = sdpvar(n,m,'full');
vars = [A(:);B(:)];

% constraints
g = [];
if eps(1) ~= 0
    g = [g; eps(1)*ones(2*n*T,1)];
end
if eps(2) ~= 0
    g = [g; eps(2)*ones(2*m*(T-1),1)];
end
if eps(3) ~= 0
    g = [g; eps(3)*ones(2*n*(T-1),1)];
end
for i = 1:T-1
    h(1+n*(i-1):n*i,:) = X(:,i+1) - A*X(:,i) - B*U(:,i);
end
cons_data = struct('ineq', g, 'eq', h);

% define sdp F
Y = sdpvar(n);          % P
S = sdpvar(m,n,'full'); % M
Z = sdpvar(r);          % W
tau = sdpvar(1);
gamma = sdpvar(1);
C1 = [Y-E*E' A*Y+B*S; Y*A'+S'*B' Y];
C2 = [Z C*Y+D1*S; Y*C'+S'*D1' Y];
C3 = gamma-trace(Z);
F = C1-tau*eye(size(C1));
Cons = [Y >= tol; Z >= tol; tau >= tol; C2 >= tol; C3 >= tol];

% input constraint ||UU*u||_2 <= 1 at the current state, and Y >= Yk
if nargin >= 7 && ~isempty(UU)
    beta = sdpvar(1);
    C4 = [beta*eye(m), UU*S; (UU*S)', Y];
    C5 = [beta, beta*xk'; (beta*xk')', Y];
    Cons = [Cons; C4 >= tol; C5 >= tol; beta >= tol; Y >= Yk];

    out.C4 = C4;
    out.C5 = C5;
    out.beta = beta;
end

% get the alternative constraint
[p_psatz, cons_psatz, Gram, Coef] = psatz_linf(F, cons_data, d, vars, A, B, eps, n, m, T);
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
