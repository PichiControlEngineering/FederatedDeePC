clear; close all
ExpData_directory = convertCharsToStrings(findFolder("Experiment Data"));

%% Settings
N_experiments = 12;
Pato_NR = 7;
refType = "squareWave"; DataSets = {"Reg", "FedGap", "FedGap", "FedGap", "FedGap", "FedGap","FedGap"};
% ExtraInfo = ["FineGrid", "FineGrid2"]; % It is possible to add multiple extra info's
ExtraInfo = ["FineGrid","_7_11_12_17","_7_10_12_17","_7_10_15_17","_7_1_11_17","_7_10_1_14","_7_15_12_14", "_7_15_1_11_V4"];
% ExtraInfo = ["FineGrid", "_7_10_1_17_V2", "_7_10_11_14_V2", "_7_1_14_17_V2", "_7_10_1_11_V2", "_7_15_11_14_V2", "_7_15_12_14_V2", "_7_10_15_17_V2", "_7_15_11_17_V2", "_7_10_1_14_V2", "_7_12_14_17_V2"];
% ExtraInfo = ["FineGrid", "_7_1_12_14_V3", "_7_1_12_17_V3", "_7_10_1_17_V3", "_7_10_12_14_V3", "_7_10_12_17_V3", "_7_11_14_17_V3", "_7_15_1_14_V3", "_7_15_1_17_V3", "_7_15_11_12_V3", "_7_15_11_14_V3", "_7_15_12_14_V3","_7_15_12_17_V3"];
% ExtraInfo = ["FineGrid", "_7_15_1_11_V4", "_7_10_15_12_V4", "_7_10_1_11_V4", "_7_12_14_17_V4", "_7_11_14_17_V4", "_7_15_11_17_V4", "_7_11_12_14_V4"];

% refType = "SmoothBF";
% refType = "FASTStepUD"; DataSets = {"Fed_sys_15", "Fed_sys_15_4_16"};

% DataSets = {"Regsys_15__4_16", "Fedsys_15_4_16"};
legend_names = {"Target", "Gap"};

%% Preallocate results
N_datasets = numel(DataSets);
RMSE_list = zeros(N_datasets, N_experiments);
RMSU_list = zeros(N_datasets, N_experiments);
RMS_DU_list = zeros(N_datasets, N_experiments);

lambda_g_list = zeros(N_datasets, N_experiments);

k_start = 400;
k_end = 8800;

%% 
if isscalar(ExtraInfo)
    ExtraInfo = repmat(ExtraInfo, N_datasets, 1);
end
%% First pass: compute RMSE, RMSU, lambda_g for all (d,i)
for d = 1:N_datasets
    for i = 1:N_experiments
        
        load(ExpData_directory+"\ExperimentResults\"+Pato_NR+"\"+refType+DataSets{d}+ExtraInfo(d)+i+".mat");
        % Extract needed indices
        signalNames = logsout.getElementNames;
        idx_u   = find(strcmp(signalNames, 'u_total'));
        idx_y   = find(strcmp(signalNames, 'y1'));
        idx_ref = find(strcmp(signalNames, 'reference'));

        % Reference only once
        if d == 1 && i == 1
            r = resample(logsout{idx_ref}.Values,0:1e-3:24,'zoh').Data;
            r_mask = ~isnan(r) & (r ~= 0);
            offset_r = round(DataStruct.T_ini./1e-3);
            r_trim = r(k_start-offset_r:k_end-offset_r);
            t = ((k_start:1:k_end)-1)'*1e-3;
            Q = DataStruct.Q; 
            R = DataStruct.R;
        end

        % Save lambda_g
        lambda_g_list(d,i) = exp_settings.lambda_g;

        % Trim data
        y_i = logsout{idx_y}.Values.Data;
        u_i = logsout{idx_u}.Values.Data;
            
        % If signals are shorter than required window, pad with zeros at end
        req_len = length(k_start:k_end);
        if length(y_i) < req_len
            pad_len = k_end;
            y_i = [y_i; zeros(pad_len, size(y_i,2))];
            u_i = [u_i; zeros(pad_len, size(u_i,2))];
        end
        y_trim = y_i(k_start:k_end);
        u_trim = u_i(k_start:k_end);

        % Compute metrics
        RMSE_list(d,i) = rms(y_trim - r_trim);
        RMSU_list(d,i) = rms(u_trim);
        RMS_DU_list(d,i) = rms(u_trim(1:end-1) - u_trim(2:end));
    end
end
legend("DeePC", "F-DeePC")
%% Compute performance J
J_y = RMSE_list;
J_u = RMSU_list;
J_DU = RMS_DU_list;
J_perf = J_y;

%% Determine optimal lambda_g for each dataset
[~, opt_J_vals] = min(J_perf, [], 2);   % index per dataset

%% Plot trajectories of optimal runs
ref_track_plot = figure(1);
i_leg = 1;
i_select = [1:N_datasets];
for d = i_select    
    i_opt = opt_J_vals(d);

    % Reload optimal experiment
    load("Experiment Data\ExperimentResults\"+Pato_NR+"\"+refType+DataSets{d}+ExtraInfo(d)+i_opt+".mat");

    signalNames = logsout.getElementNames;
    idx_u   = find(strcmp(signalNames, 'u_total'));
    idx_y   = find(strcmp(signalNames, 'y1'));

    y_trim = logsout{idx_y}.Values.Data(k_start:k_end);
    u_trim = logsout{idx_u}.Values.Data(k_start:k_end);

    % Store legend
    leg_cell{i_leg} = DataSets{d} + " $\lambda_g$ = " + num2str(lambda_g_list(d,i_opt),'%1.1e');
    i_leg = i_leg + 1;

    % Plot
    subplot(3,1,1); stairs(t, y_trim, 'LineWidth',1.4); hold on; grid on
    subplot(3,1,2); stairs(t, u_trim, 'LineWidth',1.4); hold on; grid on
    subplot(3,1,3); stairs(t, y_trim-r_trim, 'LineWidth',1.4); hold on; grid on
end

%% Finalize tracking plot
subplot(3,1,1)
xlim([t(1), t(end)])
stairs(t, r_trim, 'k--','LineWidth',1.4)
ylabel('y (rad)')
% legend([legend_names, 'r(t)'], 'Location','southeast')
% xlim([0.4, 2.2])

subplot(3,1,2)
ylabel('u (V)')
xlim([t(1), t(end)])
% xlim([0.4, 2.2])

subplot(3,1,3)
ylabel('e (rad)')
xlabel('t (s)')
xlim([t(1), t(end)])
% xlim([0.4, 2.2])

ref_track_plot.Position = 0.75.*[100, 100, 800, 800];
tightfig(ref_track_plot)
%% Save plot
% file_location = current_folder+"\trajectoryThesis.pdf";
% exportgraphics(gcf, file_location, 'ContentType','vector', 'BackgroundColor','none');

%% Order datasets for plotting (easier to visualize
[~, i_jy] = sort(min(J_y,[],2));
J_y_o = J_y(i_jy,:);
% [~, i_ju] = sort(min(J_u,[],2));
J_u_o = J_u(i_jy,:);
% [~, i_j_Du] = sort(min(J_DU,[],2));
J_DU_o = J_DU(i_jy);

%% Plot cost function J vs lambda_g
lambda_g_vals = figure(2); clf
select_vals =  [1,N_datasets];
for d = select_vals
    subplot(2,1,1)
    semilogx(lambda_g_list(d,:), J_y_o(d,:), 'LineWidth', 1.4);
    ylim([0,0.5])
    hold on

    subplot(2,1,2)
    semilogx(lambda_g_list(d,:), J_u_o(d,:), 'LineWidth', 1.4);
    hold on

    % subplot(3,1,3)
    % semilogx(lambda_g_list(d,:), J_DU(d,:), 'LineWidth', 1.4);
    % hold on
end

Title_id = ["MSE", "MSU", "MSDU"];
ylabel_id = ["${\mathcal{J}}_y$", "${\mathcal{J}}_u$"];
for i = 1:2
    subplot(2,1,i)
    grid on
    if i == 2
    xlabel('$\lambda_g$')
    title(Title_id(i) + " vs. $\lambda_g$")
    end
    ylabel(ylabel_id(i))
    % ylim([0,0.5])
    % title('Cost function vs. \lambda_g', 'Interpreter','latex')
end
tightfig(lambda_g_vals)
legend(ExtraInfo(i_jy(select_vals)), 'Interpreter','none', 'Location','best')


%%
[RMSE_vals, id_best] = sort(min(RMSE_list, [], 2));
disp(ExtraInfo(id_best)'+ ", RMSE = "+ round(RMSE_vals,3))