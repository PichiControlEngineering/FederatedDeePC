%% Analyse all the weightings as shown in FederatedPredictorResearch, we analyse the efficacy of the Federated, Gap-Based approach, as well as the LS-Fit approach
N_sims  = 100;
T       = 50;
T_sim = 60;
l       = 8;
   
r = @(t) 1.*sin((t-l)/2.5);
% r = @(t) 1.*square((t-l)/5);
% r = @(t) 0*t;
% r = @(t) heaviside(t-15);

lambda_g    = logspace(-4, 2, 16);

load("MassSpringSystemData.mat")
% load("GaussianSystemData.mat")


N_order = size(sys_cell_cell{1}{1}.A, 1);
N_sys = numel(sys_cell_cell{1});
plot_cell_cell = cell(N_sims, 1);
% Preallocate table
Results     = table('Size',[N_sims 8], ...
    'VariableTypes', repmat("cell",1,8), ...
    'VariableNames', ["RMSE","RMSU","RMSE_orac","RMSU_orac","J","J_norm", "GapList", "AlphaList"]);

for i = 1:N_sims
    i_select = i; %randi(10);

    disp("Running " + i + "/" + N_sims)

    % Run simulation
    [RMSE, RMSU, RMSE_orac, RMSU_orac, matNames, plot_cell, gap_list, alpha_list] = ...
        OptimalPrune(N_order, r, lambda_g, T, T_sim, l, GapMatCell{i_select}, sys_cell_cell{i_select});
    Q = 1; R = 1e-2;
    
    % Performance metrics
    J = RMSE.*Q + RMSU.*R;
    J_norm = J ./ J(1);

    % Store in table (each cell holds a vector/matrix)
    Results.RMSE{i}       = RMSE;
    Results.RMSU{i}       = RMSU;
    Results.RMSE_orac{i}  = RMSE_orac;
    Results.RMSU_orac{i}  = RMSU_orac;
    Results.J{i}          = J;
    Results.J_norm{i}     = J_norm;
    Results.AlphaList{i}  = alpha_list;
    Results.GapList{i}    = gap_list;

    plot_cell_cell{i} = plot_cell;
end

save("Data\MassSpringSimulationData3", "Results", "N_sims", "N_order", "N_sys","T","l","GapMatCell","lambda_g", "sys_cell_cell" , "plot_cell_cell", "matNames", "r")
function [RMSE, RMSU, RMSE_orac, RMSU_orac, matNames, plot_cell, gap_list, alpha_list] = OptimalPrune(N_order, r, lambda_g, T, T_sim, l, gap_mats, sys_cell)
rng('shuffle')    
%% PARAMETERS
    eps_list = logspace(-2, 1, 12);
    SNR = 20;
    N_sys   = numel(sys_cell);

    % Controller parameters
    N = l;
    T_ini = l;
    Q = 1;
    R = 1e-2;
    lambda_y = [] ;
    
    [y_data, u_data] = GenerateIOData(sys_cell, T); 
    u_nom = u_data(1,:); 
    y_nom = y_data(1,:);
    Sys_nom = sys_cell{1};

    H_orac = traj2Hankel([u_nom; y_nom], T_ini, N); 
    H_reg = traj2Hankel([u_nom; addNoise(y_nom, SNR)], T_ini, N);
    y_data_noisy = addNoise(y_data, SNR);

    %% LIST of pruning values & initializing other vars
    gap_list = gap_mats(:,1);
    alpha_list = zeros(N_sys, 2); 

    %% Construct controllers
    OracleMats  = SaveControllerData(H_orac, T_ini, N, Q, R); 
    RegularMats = SaveControllerData(H_reg,  T_ini, N, Q, R); 
    
    matNames = {"Oracle", "Regular"};
    HMats    = {OracleMats, RegularMats};
    
    prune_select = gap_list < 0.3;
    if ~any(prune_select)
        prune_select = 0.*gap_list;
        prune_select(1) = 1;
        disp('Federated selection failed! Choosing regular')
    end
    y_data_noisy_pruned  = y_data_noisy(prune_select, :);
    alpha_list(1,1:2) = 1; %Regular & oracle controllers
    %% Build weighted controllers for each epsilon
    for k = 1:numel(eps_list)
        % GAP + Prune weight
        alpha_prune = 1 ./ (gap_list(prune_select) + eps_list(k));
        alpha_prune = alpha_prune ./ sum(alpha_prune);
        
        y_fed_prune  = alpha_prune' * y_data_noisy_pruned;
        H_fed_prune  = traj2Hankel([u_nom; y_fed_prune], T_ini, N);
        FedPruneMats = SaveControllerData(H_fed_prune, T_ini, N, Q, R);

        HMats{end+1}   = FedPruneMats;
        matNames{end+1} = "FedPrune$_{\epsilon " + k+"}$";

        % Save weights
        alpha_list(prune_select, end+1) = alpha_prune;
        
        % LS fit weighting (using CVX)
        y_data_noisy_norm = y_data_noisy; %./norm(y_data_noisy(1,:));
        W = y_data_noisy_norm*y_data_noisy_norm' + 1e1.*eps_list(k)*eye(N_sys);
        c =  y_data_noisy_norm*y_data_noisy_norm(1,:)';

        if mod(k,2) == 0
            lambda_g1Norm = 0;
            matNames{end+1} = "LSFit$_{\epsilon " + k+"}$";
        else
            lambda_g1Norm = 1;
            matNames{end+1} = "LSFit$_{\epsilon " + k+", 1Norm}$";
        end

        cvx_begin quiet
            variable beta_LSFit(N_sys) nonnegative
            minimize(beta_LSFit'*W*beta_LSFit -2*c'*beta_LSFit + 1*lambda_g1Norm.*norm(beta_LSFit, 1))
            subject to
                sum(beta_LSFit) == 1;
                beta_LSFit >= zeros(N_sys,1);
        cvx_end

        y_NominalAlpha = beta_LSFit'*y_data_noisy;
        alpha_list(:, end+2) = beta_LSFit;

        if beta_LSFit(1) > 0.99
            disp("LS weighting gone wrong! eps = " + num2str(eps_list(k), 2))
        end
        H_LS_fit  = traj2Hankel([u_nom; y_NominalAlpha], T_ini, N);
        LSFitMats = SaveControllerData(H_LS_fit, T_ini, N, Q, R);
        HMats{end+1} = LSFitMats;
        
        if sqrt(eps_list(k)) < 1
            nom_weight = 1- sqrt(eps_list(k)); 
        else
           nom_weight = 0;
        end
        Beta_BSFit = [nom_weight; (1-nom_weight)/(N_sys-1).*ones(N_sys-1,1)];
        y_BSFit = Beta_BSFit'*y_data_noisy;
        H_BS_fit  = traj2Hankel([u_nom; y_BSFit], T_ini, N);
        BSFitMats = SaveControllerData(H_BS_fit, T_ini, N, Q, R);
        HMats{end+1} = BSFitMats;

        alpha_list(:, end+1) = Beta_BSFit;      
        matNames{end+1} = "BS Fit" + k;

    end
        H_mean  = traj2Hankel([u_nom; mean(y_data_noisy, 1)], T_ini, N);
        MeanMats = SaveControllerData(H_mean, T_ini, N, Q, R);
        HMats{end+1} = MeanMats;
        matNames{end+1} = "Mean";
        alpha_list(:, end+1) = ones(N_sys,1)./N_sys;
        
        % H_modelMismatch = traj2Hankel([u_nom; y_data_noisy(2,:)], T_ini, N);
        % modelMismatchMats = SaveControllerData(H_modelMismatch, T_ini, N, Q, R, 1);
        % HMats{end+1} = modelMismatchMats;
        % matNames{end+1} = "Biased";

    N_ctrls = numel(HMats);

    % Generate controllers recursively
    ctrl = cell(N_ctrls, 1); % for each lambda_g, and each predictor, generate 1 trajectory
    ctrl{1} = @(x,u,t) DeePCController(Sys_nom.C*x, u, t, T_ini, N, HMats{1}, r, 0, []);
    for i_ctrl = 2:N_ctrls % skip the oracle controller
        HM = HMats{i_ctrl};
        ctrl{i_ctrl} = @(x, u, t, i_reg) DeePCController(Sys_nom.C*x, u, t, T_ini, N, HM, r, lambda_g(i_reg), lambda_y);
    end
    %% Run system simulation
    % disp("Running systems...")
    active_time = T_ini+1:T_sim;
    y_traj = zeros(N_ctrls, T_sim); % Initialize place to save data
    y_traj_opt = zeros(N_ctrls, T_sim);
    x_traj_opt = zeros(N_order*N_ctrls, T_sim);
    u_traj_opt = zeros(N_ctrls, T_sim);
    lambda_g_opt = zeros(N_ctrls, 1);
    RMSE = zeros(N_ctrls, length(lambda_g));
    RMSU= zeros(N_ctrls, length(lambda_g));
    RMSE_orac = zeros(N_ctrls-1, length(lambda_g));
    RMSU_orac = zeros(N_ctrls-1, length(lambda_g));
    x_init = zeros(N_order, 1);

    % Run oracle seperately
    [x_orac, u_orac] = stepStates(x_init, T_sim, Sys_nom.A, Sys_nom.B, 1, [], ctrl{1}, T_ini);
    y_orac = Sys_nom.C*x_orac;
    RMSE(1,:) = ones(size(lambda_g)).*rms(y_orac(active_time) - r(active_time));
    RMSU(1,:) = ones(size(lambda_g)).*rms(u_orac(active_time));
    y_traj_opt(1,:) = y_orac;
    u_traj_opt(1,:) = u_orac;
    x_traj_opt(1:N_order, :) = x_orac;
    
    for i_ctrl = 2:N_ctrls
        % disp("Running controller " + i_ctrl + "/" + N_ctrls)
        for i_reg = 1:length(lambda_g)
            ctrl_i = ctrl{i_ctrl};
            [x_out, u_d] = stepStates(x_init, T_sim, Sys_nom.A, Sys_nom.B, 1, [], @(x, u, t) ctrl_i(x, u, t, i_reg), T_ini);
            y_traj(i_ctrl, :) = Sys_nom.C*x_out;
            
            % Save performance variables
            RMSE(i_ctrl, i_reg) = rms(y_traj(i_ctrl,active_time) - r(active_time));
            RMSE_orac(i_ctrl-1, i_reg) = rms(y_traj(i_ctrl,active_time) - y_orac(1,active_time));
            RMSU(i_ctrl, i_reg) = rms(u_d(active_time));
            RMSU_orac(i_ctrl-1, i_reg) = rms(u_d(active_time) - u_orac(active_time));
            
            if rms(y_traj(i_ctrl,active_time) - r(active_time))> 1e5
                disp("Halt")
            end

            % Check if the new regularization has improved results compared to
            % previous, then save with corresponding values
            if (RMSE(i_ctrl, i_reg) < rms(y_traj_opt(i_ctrl,active_time) - r(active_time))) || i_reg ==1
                y_traj_opt(i_ctrl, :) = y_traj(i_ctrl,:);
                u_traj_opt(i_ctrl, :) = u_d;
                x_traj_opt((N_order*(i_ctrl-1)+1:N_order*i_ctrl), :) = x_out;
                lambda_g_opt(i_ctrl) = lambda_g(i_reg);
            end
        end
    end
    
    plot_cell = cell(2, N_ctrls);
    for i = 1:N_ctrls
        plot_cell{1, i} = x_traj_opt((N_order*(i-1)+1):1:N_order*i, :);
        plot_cell{2, i} = u_traj_opt(i, :);
    end
end