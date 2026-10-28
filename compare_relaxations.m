%% Scherer vs Putinar relaxations (Section VII of the paper, after Example 1)
%
% Runs five receding-horizon steps of the input-constrained EIV controller
% with the matrix-valued certificate (Scherer, r = 1) and with the
% scalarized certificate (Putinar, r = 2), and reports the solve time per
% step and the size of the largest Gram matrix.
%
% Both runs use the same data: the system of Example 1 with T = 20 samples
% (N = 19 transitions) in which only the state measurements are noisy
% (||dx||_inf <= 0.1, exact inputs, no process noise), and the input
% constraint ||U*u||_2 <= 1 with U = 0.5*I.

clear
yalmip('clear')

%% Settings
relaxations = {'scherer', 1, 20;    % name, relaxation order r, # of samples T
               'putinar', 2, 20};
steps = 5;                          % receding-horizon steps
umax = 1;                           % inputs of the data are uniform in [-umax, umax]
eps = 0.1;                          % bound on the state measurement errors (l_inf)
tol = 1e-6;                         % delta in paper
opts = sdpsettings('solver','mosek','verbose',0);

root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(root, 'utils')))

%% System
A = [0.6863    0.3968
    0.3456    1.0388];
B = [0.4170    0.0001
    0.7203    0.3023];

n = size(A,1);                 % dim of state
m = size(B,2);                 % dim of input
sysd = generate_sys(n,m,A,B);

UU = 0.5 * eye(m);             % input constraint ||UU*u||_2 <= 1

summary = cell(size(relaxations,1), 1);
for k = 1:size(relaxations,1)
    name = relaxations{k,1};
    d = relaxations{k,2};
    T = relaxations{k,3};
    fprintf('\n=== %s, r = %d, T = %d ===\n', name, d, T)

    yalmip('clear')
    rng(1, 'twister')

    %% Generate the data (noise on the state measurements only)
    U = (2*rand(m, T-1)-1)*umax;
    X = zeros(n,T);
    X(:, 1) = rand(n,1);
    noise_x = (2*rand(size(X))-1)*eps;
    for t = 1:(T-1)
        X(:,t+1) = A*X(:,t) + B*U(:,t);
    end
    X_noise = X + noise_x;
    sim = struct('U',U,'X_noise',X_noise,'epsilon',eps,'tolerance',tol);

    %% Receding-horizon control
    xk = 2 * ones(n,1);
    Yk = eye(n);

    X_log = [xk, zeros(n,steps)];
    U_log = zeros(m,steps);
    NOISE = (2*rand(n, steps)-1)*eps;
    solver_time = 0;

    for t = 1:steps
        fprintf('Iteration: %d\n', t)

        % calculate control
        switch name
            case 'scherer'
                out = EIVLQR_xnoise_scherer(sim, d, T, sysd, xk, Yk, UU);
            case 'putinar'
                out = EIVLQR_xnoise_putinar(sim, d, T, sysd, xk, Yk, UU);
        end
        tic
        sol = optimize(out.cons, out.obj, opts);
        solver_time = solver_time + toc;
        Yk = value(out.Y);
        Sk = value(out.S);
        K = Sk/Yk;

        % apply control
        uk = K*xk;
        xk = A*xk + B*uk;
        xk_noise = xk + NOISE(:,t);

        % online trajectory log
        X_log(:,t+1) = xk_noise;
        U_log(:,t) = uk;

        % add the new measurement to the data
        sim.X_noise(:,end+1) = xk_noise;
        sim.U(:,end+1) = uk;
    end

    gram_size = size(out.Gram{1}, 1);
    summary{k} = sprintf('%-8s r = %d, T = %2d: Gram size %3d, %.1f s per step (%.1f s for %d steps), max ||U u_k||_2 = %.3f', ...
        name, d, T, gram_size, solver_time/steps, solver_time, steps, max(vecnorm(UU*U_log, 2)));
    results.(name) = struct('X_log', X_log, 'U_log', U_log, 'K', K, ...
        'solver_time', solver_time, 'gram_size', gram_size, 'T', T, 'd', d);
end

fprintf('\n')
fprintf('%s\n', summary{:})
