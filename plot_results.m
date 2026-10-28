function plot_results(example, folder)
% PLOT_RESULTS  Draw Fig. 1 ('example1') or Fig. 2 ('example2') of the paper.
%   plot_results(example) uses the results reported in the paper (data/).
%   plot_results(example, 'results') uses new runs of example1_linf.m or
%   example2_l2.m; run the script once with each method first.
%
%   The figures show the measured trajectories and the normalized inputs
%   ||U*u_k||_2 of the proposed controller and of the unconstrained EIV
%   controller, with the input bound ||U*u||_2 <= 1.

root = fileparts(mfilename('fullpath'));
if nargin < 2
    folder = 'data';
end
if ~isfolder(folder)
    folder = fullfile(root, folder);   % allow paths relative to the repo
end
if ~any(strcmp(example, {'example1', 'example2'}))
    error('example must be ''example1'' or ''example2''.')
end

P = load(fullfile(folder, [example '_proposed.mat']), 'X_log', 'U_log', 'UU', 'steps');
N = load(fullfile(folder, [example '_unconstrained.mat']), 'X_log', 'U_log');
UU = P.UU;
steps = P.steps;

X_ic = P.X_log;
Ub_ic = getBound(P.U_log, UU, steps);
X_un = N.X_log;
Ub_un = getBound(N.U_log, UU, steps);

color1 = [0         0.4470    0.7410];
color3 = [0.9290    0.6940    0.1250];

%% Figure 1: Trajectories

figure()
hold on
plot(X_un(1,:), X_un(2,:),  '*--', 'LineWidth', 1, 'Color', color3)
plot(X_ic(1,:), X_ic(2,:),  '*--', 'LineWidth', 1, 'Color', color1)
legend('EIV unconstrained','Proposed','Location','northwest')
title('Trajectories')
xlabel('$x_1$', 'Interpreter', 'latex')
ylabel('$x_2$', 'Interpreter', 'latex')
grid on
axis square
xlim([-.5, 2.5])
ylim([-.5, 2.5])
hold off

ax = gca; % Get the current axes handle
ax.Title.FontSize = 14;
ax.XLabel.FontSize = 14;
ax.YLabel.FontSize = 14;
ax.Legend.FontSize = 14;

%% Figure 2: Control

figure()
hold on
plot(1:steps, Ub_un,    '*--', 'LineWidth', 1, 'Color', color3)
plot(1:steps, Ub_ic,    '*--', 'LineWidth', 1, 'Color', color1)

% red constraint boundary
yline(1, 'r-', 'LineWidth', 1.5);
% red shaded region ||Uu||_2 > 1
patch([1, steps, steps, 1], ...
      [1, 1, 1.6, 1.6], ...
      [1, 0, 0], ...
      'FaceAlpha', 0.08, ...
      'EdgeColor', 'none', ...
      'HandleVisibility', 'off');

legend('EIV unconstrained','Proposed','Input bound','Interpreter', 'latex')
title('Control')
xlabel('$t$', 'Interpreter', 'latex')
ylabel('$\|Uu\|_2$', 'Interpreter', 'latex')
grid on
axis square
xlim([1, 10])
ylim([0, 1.6])
hold off

ax = gca; % Get the current axes handle
ax.Title.FontSize = 14;
ax.XLabel.FontSize = 14;
ax.YLabel.FontSize = 14;
ax.Legend.FontSize = 14;

end

function res = getBound(U_log, UU, steps)
    res = zeros(1, steps);
    for t = 1:steps
        res(t) = norm(UU*U_log(:,t), 2);
    end
end
