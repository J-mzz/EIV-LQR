%% Example 2 of "Errors-in-Variables Data-Driven Receding-Horizon LQR" (Fig. 2)
%
% Same system as Example 1 (example1_linf.m), with pointwise l_2 noise
% bounds for each time step k:
%   ||dx_k||_2 <= 0.03, ||du_k||_2 <= 0.03, ||w_k||_2 <= 0.1
% Data: one trajectory with T = 20 noisy states (N = 19 transitions)
% Input constraint: ||U*u||_2 <= 1 with U = 0.3*I; relaxation order r = 1
%
% Set `method` below:
%   'proposed'       input-constrained EIV controller (Algorithm 1 with (11))
%   'unconstrained'  EIV controller without the input constraint (ref. [6])

clear
yalmip('clear')
rng(1, 'twister')

%% Settings
method = 'proposed';        % 'proposed' or 'unconstrained'

umax = 1;                   % inputs of the data are uniform in [-umax, umax]
eps = [3e-2; 3e-2; 0.1];    % noise bounds [eps_x; eps_u; eps_w], pointwise l_2
T = 20;                     % # of samples (states) in the data
d = 1;                      % relaxation order
tol = 1e-6;                 % delta in paper
steps = 10;                 % receding-horizon steps
opts = sdpsettings('solver','mosek','verbose', 0);

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

UU = 0.3 * eye(m);             % input constraint ||UU*u||_2 <= 1

switch method
    case 'proposed',      UU_design = UU;
    case 'unconstrained', UU_design = [];
    otherwise
        error('method must be ''proposed'' or ''unconstrained''.')
end

%% Generate the data (one noisy trajectory)
U = (2*rand(m, T-1)-1)*umax;
X = zeros(n,T);
X(:, 1) = rand(n,1);
noise_x = getl2noise(n,T,eps(1));
noise_u = getl2noise(m,T-1,eps(2));
noise_w = getl2noise(n,T-1,eps(3));
for t = 1:(T-1)
    X(:,t+1) = A*X(:,t) + B*U(:,t) + noise_w(:,t);
end
X_noise = X + noise_x;
U_noise = U + noise_u;
sim = struct('X_noise',X_noise,'U_noise',U_noise,'epsilon',eps,'tolerance',tol);

%% Receding-horizon control
% initialization
xk = 2 * ones(n,1);
Yk = eye(n);

X_log = [xk, zeros(n,steps)];
U_log = zeros(m,steps);
NOISE_x = getl2noise(n,steps,eps(1));
NOISE_u = getl2noise(m,steps,eps(2));
NOISE_w = getl2noise(n,steps,eps(3));
solver_time = 0;

for t = 1:steps
    fprintf('Iteration: %d\n', t)

    % calculate control
    out = EIVLQR_l2(sim, d, T, sysd, xk, Yk, UU_design);
    sol = optimize(out.cons, out.obj, opts);
    solver_time = solver_time + sol.solvertime;
    fprintf('Full problem status: %d | %s\n', sol.problem, sol.info);

    Yk = value(out.Y);
    Sk = value(out.S);
    K = Sk/Yk;

    % apply the control to the true system (the controller is evaluated at
    % the true state xk; X_log stores the measured states)
    uk = K*xk;
    xk = A*xk + B*uk + NOISE_w(:,t);
    xk_noise = xk + NOISE_x(:,t);
    uk_noise = uk + NOISE_u(:,t);

    % online trajectory log
    X_log(:,t+1) = xk_noise;
    U_log(:,t) = uk;

    % add the new measurement to the data
    sim.X_noise(:,end+1) = xk_noise;
    sim.U_noise(:,end+1) = uk_noise;
end

%% Save and plot
out_dir = fullfile(root, 'results');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
out_file = fullfile(out_dir, ['example2_' method '.mat']);
save(out_file, 'method', 'A', 'B', 'UU', 'eps', 'T', 'd', 'steps', 'sim', ...
    'X_log', 'U_log', 'K', 'NOISE_x', 'NOISE_u', 'NOISE_w', 'solver_time')
fprintf('Saved %s\n', out_file)

res = vecnorm(UU*U_log, 2);
fprintf('max_k ||U*u_k||_2 = %.4f\n', max(res))

figure()
plot(X_log(1,:), X_log(2,:),'*--')
title(['Trajectory (' method ')'])
figure()
plot(1:steps,res,'*--')
title(['||U u_k||_2 (' method ')'])
