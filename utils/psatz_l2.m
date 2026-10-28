function [p_psatz, cons_psatz, Gram, Coef, mult] = ...
    psatz_l2(Cmat, h, d, vars, A, B, epsx, epsu, epsw, n, m, N)
% PSATZ_L2  Scherer Positivstellensatz certificate for pointwise l_2-bounded
% EIV noise, with local second-order-cone (SOC) multipliers.
%
% Certificate:
%   Cmat - I_{2n} * (
%       sum_{t=1}^{N+1} rho_x{t} * epsx
%     + sum_{t=1}^{N}   rho_u{t} * epsu
%     + sum_{t=1}^{N} ( rho_w{t} * epsw - v_w{t}' * h_t )
%   ) >= 0
%
% subject to:
%   v_x{1}' + v_w{1}'*A = 0
%   v_x{t}' - v_w{t-1}' + v_w{t}'*A = 0,  t=2,...,N
%   v_x{N+1}' - v_w{N}' = 0
%   v_u{t}' + v_w{t}'*B = 0,                t=1,...,N
%
% and cone constraints:
%   [rho_x{t}, v_x{t}'; v_x{t}, rho_x{t} I_n] >= 0
%   [rho_u{t}, v_u{t}'; v_u{t}, rho_u{t} I_m] >= 0
%   [rho_w{t}, v_w{t}'; v_w{t}, rho_w{t} I_n] >= 0

    n_var = length(vars);
    s_p   = size(Cmat,1);

    if s_p ~= 2*n
        error('Cmat must be 2n x 2n.');
    end

    if ~iscell(h)
        if size(h,1) ~= n || size(h,2) ~= N
            error('If h is not a cell, it must be n x N.');
        end
        hcell = cell(N,1);
        for t = 1:N
            hcell{t} = h(:,t);
        end
    else
        if length(h) ~= N
            error('If h is a cell, it must have length N.');
        end
        hcell = h(:);
    end

    cons_psatz = [];
    Gram = struct();
    Coef = struct();
    mult = struct();

    %--------------------------------------------------------------
    % monomial bases
    %
    % rho_x, rho_u, rho_w : degree <= 2d   (SOS)
    % v_w                 : degree <= 2d-1
    % v_x, v_u            : degree <= 2d
    %--------------------------------------------------------------
    [pow_rho, ~] = momentPowers(0, n_var, d);
    z_rho = recovermonoms(pow_rho, vars);
    l_rho = length(z_rho);

    deg_vw = max(0, 2*d - 1);
    [pow_vw, ~] = momentPowers(0, n_var, deg_vw);
    z_vw = recovermonoms(pow_vw, vars);
    l_vw = length(z_vw);

    deg_vxu = 2*d;
    [pow_vxu, ~] = momentPowers(0, n_var, deg_vxu);
    z_vxu = recovermonoms(pow_vxu, vars);
    l_vxu = length(z_vxu);

    %--------------------------------------------------------------
    % multipliers
    %--------------------------------------------------------------
    rho_x = cell(N+1,1);
    rho_u = cell(N,1);
    rho_w = cell(N,1);

    v_x = cell(N+1,1);
    v_u = cell(N,1);
    v_w = cell(N,1);

    Gram.rho_x = cell(N+1,1);
    Gram.rho_u = cell(N,1);
    Gram.rho_w = cell(N,1);

    Coef.v_x = cell(N+1,1);
    Coef.v_u = cell(N,1);
    Coef.v_w = cell(N,1);

    % rho_x, v_x
    for t = 1:N+1
        Gram.rho_x{t} = sdpvar(l_rho, l_rho);
        rho_x{t} = z_rho' * Gram.rho_x{t} * z_rho;
        cons_psatz = [cons_psatz; Gram.rho_x{t} >= 0];

        vtemp = sdpvar(n,1);
        Coef.v_x{t} = sdpvar(l_vxu, n, 'full');
        for i = 1:n
            vtemp(i) = z_vxu' * Coef.v_x{t}(:,i);
        end
        v_x{t} = vtemp;
    end

    % rho_u, v_u
    for t = 1:N
        Gram.rho_u{t} = sdpvar(l_rho, l_rho);
        rho_u{t} = z_rho' * Gram.rho_u{t} * z_rho;
        cons_psatz = [cons_psatz; Gram.rho_u{t} >= 0];

        vtemp = sdpvar(m,1);
        Coef.v_u{t} = sdpvar(l_vxu, m, 'full');
        for i = 1:m
            vtemp(i) = z_vxu' * Coef.v_u{t}(:,i);
        end
        v_u{t} = vtemp;
    end

    % rho_w, v_w
    for t = 1:N
        Gram.rho_w{t} = sdpvar(l_rho, l_rho);
        rho_w{t} = z_rho' * Gram.rho_w{t} * z_rho;
        cons_psatz = [cons_psatz; Gram.rho_w{t} >= 0];

        vtemp = sdpvar(n,1);
        Coef.v_w{t} = sdpvar(l_vw, n, 'full');
        for i = 1:n
            vtemp(i) = z_vw' * Coef.v_w{t}(:,i);
        end
        v_w{t} = vtemp;
    end

    mult.rho_x = rho_x;
    mult.rho_u = rho_u;
    mult.rho_w = rho_w;
    mult.v_x   = v_x;
    mult.v_u   = v_u;
    mult.v_w   = v_w;

    %--------------------------------------------------------------
    % main certificate
    %--------------------------------------------------------------
    scalar_term = 0;

    for t = 1:N+1
        scalar_term = scalar_term + epsx * rho_x{t};
    end
    for t = 1:N
        scalar_term = scalar_term + epsu * rho_u{t};
    end
    for t = 1:N
        scalar_term = scalar_term + epsw * rho_w{t} - v_w{t}' * hcell{t};
    end

    p_psatz = Cmat - eye(2*n) * scalar_term;

    [cons_main, G_main] = local_matrix_sos(p_psatz, d, vars);
    cons_psatz = [cons_psatz; cons_main];
    Gram.main = G_main;

    %--------------------------------------------------------------
    % local cone constraints
    %--------------------------------------------------------------
    Gram.cone_x = cell(N+1,1);
    Gram.cone_u = cell(N,1);
    Gram.cone_w = cell(N,1);

    for t = 1:N+1
        Mx = [rho_x{t}, v_x{t}'; v_x{t}, rho_x{t}*eye(n)];
        [cons_x, Gx] = local_matrix_sos(Mx, d, vars);
        cons_psatz = [cons_psatz; cons_x];
        Gram.cone_x{t} = Gx;
    end

    for t = 1:N
        Mu = [rho_u{t}, v_u{t}'; v_u{t}, rho_u{t}*eye(m)];
        [cons_u, Gu] = local_matrix_sos(Mu, d, vars);
        cons_psatz = [cons_psatz; cons_u];
        Gram.cone_u{t} = Gu;
    end

    for t = 1:N
        Mw = [rho_w{t}, v_w{t}'; v_w{t}, rho_w{t}*eye(n)];
        [cons_w, Gw] = local_matrix_sos(Mw, d, vars);
        cons_psatz = [cons_psatz; cons_w];
        Gram.cone_w{t} = Gw;
    end

    %--------------------------------------------------------------
    % polynomial equality constraints
    %--------------------------------------------------------------
    expr = v_x{1}' + v_w{1}' * A;
    cons_psatz = [cons_psatz; local_poly_eq(expr, vars)];

    for t = 2:N
        expr = v_x{t}' - v_w{t-1}' + v_w{t}' * A;
        cons_psatz = [cons_psatz; local_poly_eq(expr, vars)];
    end

    expr = v_x{N+1}' - v_w{N}';
    cons_psatz = [cons_psatz; local_poly_eq(expr, vars)];

    for t = 1:N
        expr = v_u{t}' + v_w{t}' * B;
        cons_psatz = [cons_psatz; local_poly_eq(expr, vars)];
    end
end

function [cons_out, G] = local_matrix_sos(P, d, vars)
    n_var = length(vars);
    s_p   = size(P,1);

    [pow0, ~] = momentPowers(0, n_var, d);
    vect0 = recovermonoms(pow0, vars);
    l0 = length(vect0);

    temp = kron(vect0, eye(s_p));
    G = sdpvar(l0*s_p, l0*s_p);

    Q = temp' * G * temp;
    cons_out = [G >= 0; local_poly_eq(P - Q, vars)];
end

function cons_out = local_poly_eq(expr, vars)
    cons_out = [];
    expr = expr(:);
    for k = 1:length(expr)
        cons_out = [cons_out; coefficients(expr(k), vars) == 0];
    end
end