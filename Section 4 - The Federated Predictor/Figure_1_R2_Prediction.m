%% Compare the predicting capabilities in an R^2 of multiple controller
%{
- Purpose of this code: Compare multi-step prediction performance (R^2) of several
  predictors/controllers on mass-spring system datasets.
- Workflow:
  1) If available, load precomputed "Figure_1_R2_Values.mat" to skip sims.
  2) Otherwise run Monte Carlo (N_sims) trials: generate I/O data,
     build candidate predictors with varying regularization (eps_list),
     compute multi-step predictions and store results.
  3) Reorder predictions and compute R^2(k) = 1 - SS_res/SS_tot for each
     predictor and prediction horizon k (1..N).
  4) Plot selected predictors' R^2 across horizons and export figures.
  5) Produce example predicted trajectories and export illustrative plots.
- Inputs:
  * MassSpringSystemData.mat (sys_cell_cell, GapMatCell, etc.)
  * Optional precomputed Figure_1_R2_Values.mat
- Outputs:
  * PDF figures in 'Figures\' (R^2 plot and prediction visualizations)
  * Variables: R2, Predictions_reOrdered, predictions_cell for further use
- Notes:
  * PredictorMonteCarlo encapsulates predictor construction and per-trial outputs.
  * Set plot_on = 0 to disable per-trial plotting for faster runs.
%}
close all

%% In the paper, we utilize the data in the  file "Figure_1_R2_Values.mat"
if exist("Figure_1_R2_Values.mat", "file") == 0
    %% Parameters
    N_sims  = 100 ;
    l       = 18;
    N       = l;
    T_ini   = l;
    T       = 6*l;
    
    lambda_g  = 0;
    eps_list  = logspace(-3, -0.5, 8); 
    plot_on = 1;    
    
    load MassSpringSystemData.mat
    
    N_order = size(sys_cell_cell{1}{1}.A,2);
    
    N_sys = length(sys_cell_cell{1});
    opts = optimoptions(@quadprog,'Algorithm', 'interior-point-convex','Display', 'Off', 'OptimalityTolerance',1e-13, 'ConstraintTolerance', 1e-9, 'MaxIterations', 5e3);
    %% Preallocate storage
    y_real_list  = zeros(N_sims, N);
    pad_list = zeros(N_sims, 2);
    %% Monte Carlo simulations
    predictions_cell = cell(N_sims,1);
    for i = 1:N_sims
        disp("Running " + i + "/" + N_sims)
        i_select = 3; 
        if N_sims > 5
            plot_on = 0;
        end
        % Run predictor Monte Carlo
        [u_traj, y_traj, y_predicted_cell, matNames, gap_list, alpha_list, gap_pad] = ...
            PredictorMonteCarlo(N_order, N_sys, lambda_g, T, T_ini, N, sys_cell_cell{i_select}(1:N_sys), GapMatCell{i_select}(1:N_sys,1), eps_list, opts, plot_on);
        pad_list(i, :) = gap_pad';
        predictions_cell{i} = cell2mat({y_predicted_cell{:}}');
    
        y_real_list(i,:)  = y_traj(T_ini+1:T_ini+N);
    end
    
    %% Compute R² values
    N_predictors = size(predictions_cell{1},1); % Total amount of predictors (for different reg vals)
    R2 = zeros(N_predictors, N);
    Predictions_reOrdered = cell(N_predictors,1);
    
    for i = 1:N_predictors
        predictions_i = zeros(N_sims,N);
        for i_sim = 1:N_sims
            pred_i = predictions_cell{i_sim};
            predictions_i(i_sim,:) = pred_i(i,:);
        end
        Predictions_reOrdered{i} = predictions_i;
    
        for h = 1:N
            y_true_h = y_real_list(1:N_sims, h);
            y_pred_h = predictions_i(:, h);
          
            SS_res = sum((y_true_h - y_pred_h).^2);
            SS_tot = sum((y_true_h - mean(y_true_h)).^2);
            
            R2(i,h) = 1 - SS_res/SS_tot;
        end
    end

else
    load("Figure_1_R2_Values.mat")
end
%% Plot results
i_select = [1, 2, 19]; 
colorlist = {[0 0 0], [0  0.447    0.741],  [0, 0.5, 0]};

fig = figure();
for i = 1:length(i_select)
plot(R2(i_select(i),:)' * 100, '-o', 'LineWidth', 1.5, 'Color',colorlist{i}); grid on; hold on
end
ylim([75 104])
ylabel("$R^2$ (\%)")

% title('$R^2$ value of predictions')
xlabel("k-step ahead prediction")

xticks(1:2:N)
xtickangle(0)
xlim([1 N])
leg_names = {"$\{u^{\mathbf{0}}, y^{\mathbf{0}} \}$",...
    "$\{\tilde{u}^{\mathbf{0}}, \tilde{y}^{\mathbf{0}}\}$", "$\{\tilde{u}^{\mu}, \tilde{y}^{\mu}\}$"};
% legend(leg_names, "Location","southwest", "Fontsize", 18)
    tightfig(fig)

exportgraphics(gcf, 'Figures\Figure_1_R2_values_plot.pdf', 'ContentType','vector', 'BackgroundColor','none');
disp("Saved image!")
	    
%% Make a plot highlighting the prediction process of the multi-step predictor
predict_plot_fig = figure;
% % % U PLOT
subplot(2,1,2)
plot(1:(T_ini+N), u_traj, '-o', 'LineWidth', 1.5, 'Color','k'); grid on
xlim([0 N+T_ini])
x = xline(T_ini, 'r--', 'T_{ini}', 'LineWidth',2, 'LabelColor','red');
% change label size
x.LabelHorizontalAlignment = 'left';     %
x.LabelVerticalAlignment   = 'top';   %
x.FontSize = 24;      
ylabel('u(k)')
xlabel('k')

% % % Y PLOT
predicted_y = cell2mat(y_predicted_cell([1,2,20]));
y_traj_ini = y_traj(1:T_ini);
y_traj_predicted = [repmat(y_traj_ini(end),3,1), predicted_y];
subplot(2,1,1)
plot(1:T_ini, y_traj_ini, '-o', 'LineWidth', 1.5, 'Color','k'); hold on; grid on
xlim([0 N+T_ini])
ylabel('y(k)')
ylim([-4 2])
x = xline(T_ini, 'r--', 'T_{ini}', 'LineWidth',2, 'LabelColor','red');

% change label size
x.LabelHorizontalAlignment = 'left';     %
x.LabelVerticalAlignment   = 'top';   %
x.FontSize = 24;                         % change label font size

%Set fig size
predict_plot_fig.Position = [100, 100 , 1240, 720];
exportgraphics(gcf, 'Figures\PredictedPlots.pdf', 'ContentType','vector', 'BackgroundColor','none');

%
for i = 1:3
plot(T_ini:(T_ini+N), y_traj_predicted(i,:), '-o', 'LineWidth', 1.5, 'Color',colorlist{i});
    if i == 1    
        exportgraphics(gcf, 'Figures\PredictedPlots_FutureVals_OracleOnly.pdf', 'ContentType','vector', 'BackgroundColor','none');
    end
end

exportgraphics(gcf, 'Figures\PredictedPlots_FutureVals.pdf', 'ContentType','vector', 'BackgroundColor','none');
%%
function [u_traj, y_traj, y_predicted_cell, matNames, GapMat, alpha_list, gap_pad] = ...
         PredictorMonteCarlo(N_order, N_sys, lambda_g, T, T_ini, N, sys_cell, GapMat, eps_list, opts, plot_on)    
    %% Generate system data
    % disp('Generating systems...')
    % [sys_cell, GapMat] = GenerateRandomSys(N_order, N_sys, show_progress, ...
    % 'MaxGap', MaxGap, 'MinGap', MinGap, "TrueRandom",1,"Sys_nom", Sys_nom);

    % We use MegaSystemData, which uniformly created
    % 250 datasets, in gaprange [0 0.55]

    Sys_nom = sys_cell{1};
    SNR = 20;
    [y_data, u_data] = GenerateIOData(sys_cell, T, zeros(N_order,1)); 
    u_nom = u_data(1,:); 
    y_nom = y_data(1,:);

    H_orac = traj2Hankel([u_nom; y_nom], T_ini, N); 
    H_reg  = traj2Hankel([u_nom; addNoise(y_nom, SNR)], T_ini, N);
    y_data_noisy = addNoise(y_data, SNR);
    % variance = 0.2*rms(y_nom);
    % H_reg  = traj2Hankel([u_nom; (y_nom + variance.*randn(size(y_nom)))], T_ini, N);
    % y_data_noisy = (y_data + variance.*randn(size(y_data)));

    %% Initialize pruning and weights
    alpha_list = zeros(N_sys, 2 + numel(eps_list)); 

    %% Construct controllers
    OracleMats  = SaveControllerData(H_orac, T_ini, N, 1, 1); 
    RegularMats = SaveControllerData(H_reg,  T_ini, N, 1, 1); 
    
    if any(size(H_orac, 2) <= [T_ini, N])
        error("Matrix too small, increase T or decrease T_ini, N")
    end
    matNames = {"Oracle", "Regular"};
    HMats    = {OracleMats, RegularMats};
    
    prune_select           = GapMat < 0.45;
    % y_data_noisy_pruned    = y_data_noisy(prune_select, :);
    y_data_noisy_pruned_2  = y_data_noisy(prune_select, :);

    %% Build weighted controllers for each epsilon
    for k = 1:numel(eps_list)
        % GAP + Prune weight
        alpha_prune = 1 ./ (GapMat(prune_select, 1) + eps_list(k));
        alpha_prune = alpha_prune ./ sum(alpha_prune);
        
        y_fed_prune  = alpha_prune' * y_data_noisy_pruned_2;
        H_fed_prune  = traj2Hankel([u_nom; y_fed_prune], T_ini, N);
        FedPruneMats = SaveControllerData(H_fed_prune, T_ini, N, 1, 1, 1);

        HMats{end+1}   = FedPruneMats;
        matNames{end+1} = "FedPrune$_{\epsilon " + k+"}$";

        % Save weights
        alpha_list(prune_select, k+1) = alpha_prune;

        % CVX nominal alpha
        y_data_noisy_norm = y_data_noisy./norm(y_data_noisy(1,:));
        Q = y_data_noisy_norm*y_data_noisy_norm' + 1e2.*eps_list(k)*eye(N_sys);
        c =  y_data_noisy_norm*y_data_noisy_norm(1,:)';

        cvx_begin quiet
            variable beta_LSFit(N_sys) nonnegative
            minimize(beta_LSFit'*Q*beta_LSFit -2*c'*beta_LSFit + norm(beta_LSFit, 1))
            subject to
                sum(beta_LSFit) == 1;
        cvx_end

        % beta_LSFit = (Q\(c - ((ones(N_sys,1)'*(Q\c))-1)/(ones(1,N_sys)*(Q\ones(N_sys,1)))));
        % disp([beta_LSFit, a_NomMatch_anal])
        
        y_NominalAlpha = beta_LSFit'*y_data_noisy;
        if beta_LSFit(1) > 0.99
            disp('LS weighting gone wrong!')
        end
        H_LS_fit  = traj2Hankel([u_nom; y_NominalAlpha], T_ini, N);
        LSFitMats = SaveControllerData(H_LS_fit, T_ini, N, 1, 1);
        HMats{end+1} = LSFitMats;
        matNames{end+1} = "LSFit$_{\epsilon " + k+"}$";
        % if k ==8
        %     % stop!
            gap_pad(1) = 0; %PadoanGap([u_nom; y_NominalAlpha], [u_nom; y_nom], 20);
            gap_pad(2) = 0; %PadoanGap([u_nom; addNoise(y_nom, SNR)], [u_nom; y_nom], 20);
        %     % disp(gap_pad)
        % end
    end
        H_mean  = traj2Hankel([u_nom; mean(y_data + 0.01.*randn(size(y_data)), 1)], T_ini, N);
        MeanMats = SaveControllerData(H_mean, T_ini, N, 1, 1);
        HMats{end+1} = MeanMats;
                matNames{end+1} = "Mean";

        % We purpsefully take a model with a Wrong description of the
        % system to analyse predictor-bias
        if size(y_data_noisy_pruned_2,1) < 2 % If the pruning fails, we just grab the second available trajectory
            y_data_noisy_pruned_2 = y_data_noisy(1:2,:);
        end
        H_modelMismatch = traj2Hankel([u_nom; y_data(2,:)], T_ini, N);
        modelMismatchMats = SaveControllerData(H_modelMismatch, T_ini, N, 1, 1);
        HMats{end+1} = modelMismatchMats;
        matNames{end+1} = "Biased";


    N_ctrls = numel(HMats);
    y_predicted_cell = cell(N_ctrls, numel(lambda_g)); 

    %% System simulation with random inputs
    % randomSig = 1 - 2 .* rand(1, T_ini+N); 
    randomSig = multisine([2/T, 1/3], 1, T, 'PhaseResponse', 'NormalDistribution');
    u_predict_func = @(x,u,t) randomSig(t); 

    [x_traj, u_traj] = stepStates(rand(N_order, 1), T_ini+N, Sys_nom.A, Sys_nom.B, 1, 0, u_predict_func);
    y_traj = Sys_nom.C * x_traj;

    %% Predictions
    disp("Running predictions...")
    u_ini = u_traj(1:T_ini)'; 
    y_ini = y_traj(1:T_ini)';
    u_f   = u_traj(T_ini+1:T_ini+N)';
    y_f   = y_traj(T_ini+1:T_ini+N)';
    
    if plot_on
        figure()
        plot(0:N, [y_ini(end); y_f], 'ko--', LineWidth=1.4); hold on; grid on;
    end

    g_size = T-N-T_ini+1;

    for i_pred = 1:N_ctrls
        % Construct predictor matrix
         Z_mat = [HMats{i_pred}.H_p; HMats{i_pred}.U_f];
         z_ini = [u_ini; y_ini; u_f];
            
        if i_pred == 1 
                % If we call the oracle controller, compute via Pinv
                % Since our oracle uses no regularization
                g = pinv(Z_mat)*z_ini;  
                y_predicted_cell{i_pred, 1} = (HMats{i_pred}.Y_f * g)';    
        else
            for i_reg = 1:numel(lambda_g)
             % Else, solve via quadprog (for each reg value)
                lambda_g_i = lambda_g(i_reg);
                g = pinv(Z_mat)*z_ini;  
                % g = quadprog(lambda_g_i .* eye(g_size), zeros(g_size,1), [], [], Z_mat, z_ini, [],[],[],opts);
                % linprog()
                y_predicted_cell{i_pred, i_reg} = (HMats{i_pred}.Y_f * g)';  
            end
        end
            if plot_on
                if any(i_pred == [2,N_ctrls])
                    plot(0:N, [y_ini(end), y_predicted_cell{i_pred, 1}], 'o-', LineWidth=1.4)
                end
            end    
    end

    if plot_on
    legend(["real", "Regular", "Federated"], "Interpreter", "none")
    ylabel('y'); xlabel('k-step ahead prediction')

    figure()
    subplot(2,1,1)
        plot(y_data(1,:), 'k--', 'LineWidth',1.6); hold on; grid on;
        plot(y_data_noisy(1,:), 'LineWidth',1.6)
        plot(y_fed_prune, 'LineWidth',1.6)
        ylabel('y')
        legend("Ora", "Reg", "Fed")
    subplot(2,1,2)
        plot(u_data(1,:), 'k', 'LineWidth',1.6); grid on
        ylabel('u')
        xlabel('t')
    sgtitle('Utilized IO trajectories')
    end
end
