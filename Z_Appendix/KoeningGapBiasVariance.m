clear; close all;


% Monte Carlo settings for error statistics over varying Hankel sizes
nTrials = 10;                        % number of Monte Carlo repetitions per noise level
l_list = [120,150];                   % different Hankel row sizes to test (x-axis)
sigma_var_list = 0:0.015:0.2;        % keep original noise list (from prefix)
T = 1000;
Ts = 1;

% Load System data
load("GaussianSystemData.mat", "sys_cell_cell")
sys_stack = reshape([sys_cell_cell{1:ceil(nTrials/20),1}],[],1);
sys_nom = sys_stack{1};

% Preallocate storage: dimensions (nHankelSizes x nNoise x nTrials)
nL = length(l_list);
nS = length(sigma_var_list);

gap_koen_trials = zeros(nL, nS, nTrials);
gm_val = zeros(1, 1, nTrials);

%% Input
[butter_nom, butter_denom] = butter(3, 0.7, 'low');
u_input = filter(butter_nom, butter_denom, rand(1,T));

w2_list = zeros(2*nTrials, T);
%% pre-generate data
for t = 1:nTrials
    
    %% Clean trajectories
    y1 = lsim(sys_nom, u_input, 1:Ts:T);
    y2 = lsim(sys_stack{t}, u_input, 1:Ts:T);

    %% Build w = [u; y]
    w1 = [u_input; y1'];

    w2_list(2*t-1:2*t, :) = [u_input; y2'];
end

%% For methods that are deterministic w.r.t noise realization we still compute per trial
% so that statistics are consistent. Main triple loop: over hankel sizes, noise levels, trials.
for iL = 1:nL
    l_cur = l_list(iL);
    fprintf("Testing hankel size l = %d\n", l_cur);
    for iS = 1:nS
        sigma_var = sigma_var_list(iS);
        fprintf("  Noise σ = %.3f\n", sigma_var);
        for t = 1:nTrials
            %% Block Hankels (for Koenings)
            %% Add noise
            w1n = w1 + sigma_var * [zeros(1,T); randn(1,T)];
            w2n = w2_list(2*t-1:2*t, :) + sigma_var * [zeros(1,T); randn(1,T)];

            gap_koen_trials(iL, iS, t) = KoeningsGap(w1n, w2n, l_cur, l_cur);
        end
    end
end
% The value of the gapmetric is independent of noise and length "l". We can
% replicate this matrix for only one trial
gap_matlab_trials = repmat(gm_val, [nL, nS, 1]);
% Compute error statistics: err = Koen - MATLAB, then mean and variance over trials
err_trials = gap_koen_trials - gap_matlab_trials;    % size nL x nS x nTrials
err_mean = mean(err_trials, 3);                      % nL x nS
err_var  = var(err_trials, 0, 3);                    % nL x nS

% Plotting: for each noise level, plot mean and variance vs hankel size l
colors = lines(nS);

% Mean plot
figMean = figure; clf; hold on; grid on; box on;
for iS = 1:nS
    plot(l_list, err_mean(:, iS), '-o', 'Color', colors(iS,:), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('\\sigma=%.3f', sigma_var_list(iS)));
end
xlabel('Hankel row size l')
ylabel('Mean error (Koen - MATLAB)')
title('Mean error vs Hankel size for different noise levels')
legend('Location','best')
tightfig(figMean)

% Variance plot
figVar = figure; clf; hold on; grid on; box on;
for iS = 1:nS
    plot(l_list, err_var(:, iS), '-s', 'Color', colors(iS,:), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('\\sigma=%.3f', sigma_var_list(iS)));
end
xlabel('Hankel row size l')
ylabel('Variance of error (Koen - MATLAB)')
title('Variance of error vs Hankel size for different noise levels')
legend('Location','best')
tightfig(figVar)

% Save statistics to workspace
stats.l_list = l_list;
stats.sigma_var_list = sigma_var_list;
stats.err_mean = err_mean;
stats.err_var = err_var;
assignin('base','gap_error_stats',stats);

%% Functions
function [L11, L21, L22, L31, L32, L33, Q1, Q2, Q3] = split_lq(L_decomp,Q_decom, T_ini, N)
    l_1 = 2*T_ini; l_2 = l_1 + N;
    
    %%%%%%%%%%%%%
    L11 = L_decomp(1:l_1,1:l_1);
    L21 = L_decomp(l_1+1:l_2, 1:l_1); L22 = L_decomp(l_1+1:l_2, l_1+1:l_2);
    L31 = L_decomp(l_2+1:end, 1:l_1); L32 = L_decomp(l_2+1:end, l_1+1:l_2); L33 = L_decomp(l_2+1:end, l_2+1:end);

    Q1 = Q_decom(1:l_1,:); Q2 = Q_decom(l_1+1:l_2,:); Q3 = Q_decom(l_2+1:end,:); 
end

function [L,Q] = lq_decomp(X)
% Computes the LQ decomposition of a matrix
[U, R] = qr(X', 0);
L = R'; Q = U';
end

function [U1, U2, S1, S2, V1, V2] = splitSVD(A, last_svd, varargin)
    % Split the SVD into 2 parts: a noise-free and a noise-based part.
    % Using varargin = "tol", you can split the svd for a specific
    % tolerance value

    [U, S, V] = svd(A);
    tol_u = last_svd >= 1:size(S,1);
    tol_v = last_svd >= 1:size(S,2);

    if ~isempty(varargin)
        if strcmpi(varargin{1}, "tol")
            tol_u = diag(S) > last_svd;
            tol_v = diag(S) > last_svd;
        end
    end
    
    U1 = U(:,tol_u); V1 = V(:,tol_v); S1 = S(tol_u,tol_v);
    U2 = U(:,~tol_u); V2 = V(:,~tol_v); S2 = S(~tol_u,~tol_v);
end