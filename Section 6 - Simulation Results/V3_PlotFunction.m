function V3_PlotFunction(filename)
%% Use The results from OptPruneMonteCarlo and plot them
% Data processing
if filename{1}(6) == 'M'
  system_type = 'MassSpring';
  
elseif filename{1}(6) == 'G'
   system_type = 'Gaussian';
else
    system_type = "Unknown";
end
default_colors = {[0    0.447    0.741],...
          [0.850    0.325    0.098]...
          [0.9290    0.6940    0.1250]...
          [0.4940    0.1840    0.5560] ...
          [0.4660    0.6740    0.1880]};
load(filename) % Load one of the datasets
lambda_start = 4; lambda_end = 14; % length(lambda_g);

N_ctrl = length(matNames) - 1; % remove oracle
J_mat = zeros(N_ctrl, length(lambda_g), N_sims);
J_mat_orac = zeros(N_ctrl, length(lambda_g), N_sims);

Q = 1; R = 1e-2;
% Create data for swarmplot
for i = 1:N_sims
    J_norm_sel = Results.J_norm{i};
    J_mat(:,:, i) = J_norm_sel(2:end,:);
    J_mat_orac(:,:, i) = Q.*Results.RMSE_orac{i} + R.*Results.RMSU_orac{i};
end

% Create labels for swarmplot
% Vectorized label creation
reg_labels = compose('%.1e', lambda_g(lambda_start:lambda_end));  % cell array of strings

% % Remove leading zero in exponent
% reg_labels = regexprep(reg_labels, 'e0?(-?\d+)', 'e$1');

J_opt_choice = squeeze(min(J_mat,[],2)); % selecting the lambda_g that minimizes the J
% J_opt_choice_orac = squeeze(min(J_mat_orac,[],2)); % selecting the lambda_g that minimizes the J

%
MatNames_4Plots = {"Baseline", "Gap", "LS", "Simple"};

% We split in 3 groups, the Gap, LS, and Simple
ind_Gap = 2:3:N_ctrl;
ind_LS = 3:3:N_ctrl;
ind_Simple = 4:3:N_ctrl;

[gap_optimal_i] = Mean_Median_STD_compute(J_opt_choice(ind_Gap, :));
[LS_optimal_i] = Mean_Median_STD_compute(J_opt_choice(ind_LS, :));
[Simple_optimal_i] = Mean_Median_STD_compute(J_opt_choice(ind_Simple, :));

% Select, for each ctrl type, which ctrler minimizes the mean optimal
% performance given some lambda_g
Gap_ctrl_i = ind_Gap(gap_optimal_i(1));
LS_ctrl_i = ind_LS(LS_optimal_i(1));
Simple_ctrl_i = ind_Simple(Simple_optimal_i(1));
% % Defining the quantative performance variables & analysing

i_select = [1, [Gap_ctrl_i, LS_ctrl_i, Simple_ctrl_i]]; %The selected controllers, with best selected of each

%% Local optima
% % Defining the quantative performance variables & analysing
J_opt_med = median(J_opt_choice,2);
J_opt_mean = mean(J_opt_choice,2);
J_opt_std = std(J_opt_choice,[],2);

k_select = 15;
[~, ind_median] = mink(J_opt_med, k_select);
[~, ind_mean] = mink(J_opt_mean, k_select);
[~, ind_std] = mink(J_opt_std, k_select);

disp("Median Optimal values: ")
disp(vertcat(matNames{[1, ind_median']+1 }) + ":    " + J_opt_med([1, ind_median']))

disp("Mean Optimal values: ")
disp(vertcat(matNames{[1,ind_mean']+1}) + ":    " + J_opt_mean([1,ind_mean']))

disp("Best STD values: ")
disp(vertcat(matNames{[1,ind_std']+1}) + ":    " + J_opt_std([1,ind_std']))

%%
swarmplot_RMS = figure();
tiledlayout(1,length(i_select),'TileSpacing','compact','Padding','compact');
% 
i_count = 1;
for i = i_select
    nexttile

    Y = squeeze(J_mat(i, lambda_start:lambda_end, :))';   % size = [S x L]

    % J_i = J_mat(i, lambda_start:lambda_end, :);
    % J_plot = J_i; J_plot(J_i > 5) = NaN;
    % Y = squeeze(J_plot)';   % size = [S x L]

    [L, S] = deal(size(Y,2), size(Y,1));
    
    x = repelem(1:L, S)';    % vector of length S*L
    y = Y(:);                % vector of length S*L
   
    swarmchart(x, y, 0.8, default_colors{i_count}, 'filled');
    hold on
    m = median(Y, 1);
    for k = 1:L
        plot([k k], [m(k) m(k)], 'k_', 'LineWidth', 1); % horizontal tick at mean
    end
    hold off
    grid on;

    if i_count == 1
        ylabel('$\overline{\mathcal{J}}$')
    end

    y_ub = max(prctile(squeeze(J_mat(2, lambda_start:lambda_end, :))', 99));
    title(MatNames_4Plots{i_count})
    ylim([0 y_ub])
    xticks(1:L)
    xticklabels(reg_labels)
    xlabel('$\lambda_g$')

    i_count = i_count + 1;
end
swarmplot_RMS.Position = [150, 100, 1100, 600];
tightfig(swarmplot_RMS)
file_location = "Figures\"+system_type+"Swarmplot.pdf";
exportgraphics(gcf, file_location{1}, 'ContentType','vector', 'BackgroundColor','none');
disp("Saved image!")

%% Creating Histograms
% Sizes
[~, ~, nFeat] = size(J_mat);

% Find argmin across columns (same as before)
[~, col_idx] = min(mean(J_mat,3), [], 2);

% Preallocate
J_hist = zeros(N_ctrl, nFeat);
% Select row-wise
for i = 1:N_ctrl
    J_hist(i,:) = squeeze(J_mat(i, col_idx(i), :));
end

%% Plotting "Optimal Choice" histogram: What if for each experiment, we would take the time to independently find our optimal lambda_g?

% Find argmin across columns (same as before)
max_val = ceil(prctile(J_opt_choice,98,"all")); %determining maximum
hist_yBounds = [0 N_sims/5];
h_edge = 0:(max_val/50):max_val;

% Plot histograms
hist = figure();
mean_J_for_i_select = zeros(size(i_select));
for i = 1:length(i_select)
    mean_J_for_i_select(i) = median(J_opt_choice(i_select(i),:));
    subplot(length(i_select),1,i)
    histogram(J_opt_choice(i_select(i),:), h_edge, 'FaceColor',default_colors{i}); hold on
    xline(mean_J_for_i_select(i), 'k--', 'LineWidth',1.2)
    xlim([0 min(4, max_val)]); ylim(hist_yBounds)
    legend({MatNames_4Plots{i}, "Median"})
    ylabel('counts'); 
end
    xlabel('$\overline{\mathcal{J}}$')
% sgtitle("$\mathcal{J}$ for optimal $\lambda_g$, $N ="+N_sims+"$ Monte-Carlo runs")

hist.Position = [200, 100, 800, 600];
tightfig(hist)
file_location = "Figures\"+system_type+"Histogram.pdf";
exportgraphics(gcf, file_location{1}, 'ContentType','vector', 'BackgroundColor','none');
disp("Saved image!")

%% Plotting trajectories
[~, Opt_Fed_Ctrl_select] = min(mean_J_for_i_select(2:end));
ctrl_select_4_plot = [1,2, i_select(Opt_Fed_Ctrl_select+1)+1]; %we have to shift the i_select by 1 in order to compensate for the fact this one was in the shifted index (wihtout oracle)
color_select = [1,2, Opt_Fed_Ctrl_select+2];
% Compute the means
[h, w] = size(plot_cell_cell{1}{1});
plot_mat = zeros(h+1, 3*w, N_sims);
for i = 1:N_sims
    % Plot optimal trajectory
    plot_cell_i = plot_cell_cell{i}(:, ctrl_select_4_plot); %Extract trajectories
    plot_mat(:,:,i) = cell2mat(plot_cell_i);
end
obsv_state = find(sys_cell_cell{1}{1}.C ~= 0);
plot_mat_mean = mean(plot_mat([obsv_state,end],:,:), 3);
% mean_plot_cell = mat2cell(plot_mat_mean, [h 1],[w w w]);
mean_plot_cell = mat2cell(plot_mat_mean, [1 1],[w w w]);

% Plot the means
empty_leg = {"", "","","","","","","",""}; %we put legend in the caption
plot_fig = PlotFunction(mean_plot_cell, l, r, '', empty_leg, 1, "max", 1, color_select);
plot_fig.Position = [150, 150, 800, 800];
tightfig(plot_fig)
file_location = "Figures\"+system_type+"Trajectory.pdf";
exportgraphics(gcf, file_location{1}, 'ContentType','vector', 'BackgroundColor','none');
disp("Saved Trajectory image!")

end

function ind_out = Mean_Median_STD_compute(J_in)
%This function computes the mean, median, and std and outputs which indices
%they belong to
J_opt_mean = mean(J_in,2);
J_opt_med = median(J_in,2);
J_opt_std = std(J_in,[],2);

[~, ind_mean] = min(J_opt_mean);
[~, ind_med] = min(J_opt_med);
[~, ind_std] = min(J_opt_std);
ind_out = [ind_mean, ind_med, ind_std];
end