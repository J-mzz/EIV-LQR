function out = EIVLQR_l2(sim, d, T, sysd, xk, Yk, UU)
% EIVLQR_L2  EIV data-driven LQR SDP for pointwise l_2-bounded noise (Scherer, (11)).
%   out = EIVLQR_l2(sim, d, T, sysd, xk, Yk, UU) builds the SDP for a state
%   feedback u = S*inv(Y)*x that minimizes the worst-case LQR cost bound
%   gamma over all (A, B) consistent with the first T noisy samples in sim,
%   using local SOC multipliers and the relaxation order d.
%
%   sim   struct with fields X_noise (n x T measured states), U_noise
%         (m x T-1 measured inputs), epsilon = [eps_x; eps_u; eps_w]
%         (pointwise l_2 bounds of the state and input measurement errors and
%         of the process noise), and tolerance
%   sysd  performance matrices from generate_sys
%   xk    current state;  Yk  previous Y (Y >= Yk is enforced)
%   UU    input constraint ||UU*u||_2 <= 1 at xk; with UU = [] (or omitted)
%         there is no input constraint, which gives the unconstrained EIV
%         controller
%
%   Solve with optimize(out.cons, out.obj, opts); then K = value(out.S)/value(out.Y).

    %#ok<*AGROW>

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Xall = sim.X_noise;
    Uall = sim.U_noise;

    eps = sim.epsilon(:);
    tol = sim.tolerance;

    epsx = eps(1);
    epsu = eps(2);
    epsw = eps(3);

    % only use a window of length T
    X = Xall(:,1:T);
    U = Uall(:,1:T-1);

    n = size(X,1);
    m = size(U,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % performance matrices
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    r  = sysd.r;
    E  = sysd.E;
    C  = sysd.C;
    D1 = sysd.D1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % system variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    A = sdpvar(n,n,'full');
    B = sdpvar(n,m,'full');
    vars = [A(:); B(:)];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % residuals h_t = x_{t+1} - A x_t - B u_t
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h = cell(T-1,1);
    for t = 1:T-1
        h{t} = X(:,t+1) - A*X(:,t) - B*U(:,t);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % controller/LQR variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y = sdpvar(n,n,'symmetric');     % P
    S = sdpvar(m,n,'full');          % M
    Z = sdpvar(r,r,'symmetric');     % W
    eta   = sdpvar(1,1);
    gamma = sdpvar(1,1);

    % main LQR matrices
    C1 = [Y - E*E',         A*Y + B*S;
          Y*A' + S'*B',     Y];

    C2 = [Z,                C*Y + D1*S;
          Y*C' + S'*D1',    Y];

    C3 = gamma - trace(Z);

    % robust matrix to certify on consistency set
    F = C1 - eta*eye(2*n);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % base constraints
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Cons_base = [];
    Cons_base = [Cons_base; Y >= tol*eye(n)];
    Cons_base = [Cons_base; Z >= tol*eye(r)];
    Cons_base = [Cons_base; eta >= tol];
    Cons_base = [Cons_base; C2 >= tol*eye(size(C2,1))];
    Cons_base = [Cons_base; C3 >= tol];

    % input constraint ||UU*u||_2 <= 1 at the current state, and Y >= Yk
    if nargin >= 7 && ~isempty(UU)
        beta = sdpvar(1,1);
        C4 = [beta*eye(m), UU*S;
              (UU*S)',     Y];
        C5 = [beta,        beta*xk';
              beta*xk,     Y];
        Cons_base = [Cons_base; C4 >= tol*eye(m+n)];
        Cons_base = [Cons_base; C5 >= tol*eye(n+1)];
        Cons_base = [Cons_base; beta >= tol];
        Cons_base = [Cons_base; Y >= Yk];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % P-satz / SOS certificate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [p_psatz, cons_psatz, Gram, Coef, mult] = ...
        psatz_l2( ...
            F, ...           % 2n x 2n polynomial matrix
            h, ...           % cell(T-1,1), each n x 1
            d, ...
            vars, ...
            A, B, ...
            epsx, epsu, epsw, ...
            n, m, T-1);

    % total constraints
    Cons = [Cons_base; cons_psatz];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % outputs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    out.poly = p_psatz;
    out.cons = Cons;
    out.cons_base = Cons_base;
    out.cons_psatz = cons_psatz;

    out.Gram = Gram;
    out.Coef = Coef;
    out.mult = mult;

    out.Y = Y;
    out.S = S;
    out.Z = Z;
    out.eta = eta;
    out.obj = gamma;

    out.vars = vars;
    out.Avar = A;
    out.Bvar = B;
    out.h = h;
    out.F = F;
    out.C1 = C1;
    out.C2 = C2;
    out.C3 = C3;
end
