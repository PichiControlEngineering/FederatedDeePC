%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run DeePC Experiments: Du Cost penalization %
% Author: Gert Vankan                         %
% Created: 13/08/2025                         %
% Version: 11/03/2026                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preliminaries
clear; clc;

addpath(genpath('Data')); 
addpath(genpath('Reference Generating'))
addpath(genpath('F_DeePC_Functions'))

% Change the current folder to the folder of this m-file.
if(~isdeployed)
  cd(fileparts(which(mfilename)));
end

exp_path = convertCharsToStrings(findFolder("ExperimentResults")); % Look for the folder 'ExperimentResults' on matlab path
tg = slrealtime;

%% Reference definition
% Examples:
% t_points = [0,0.9,1,2,2.1,3]; p_points = pi.*[0,0,1,1,0,0]; current_ref = "\StepUD";
% t_points = [0,0.9,1,2]; p_points = pi.*[0,0,1,1]; current_ref = "\RepStep";
% t_points   = 0:0.5:2; p_points   = 2*pi.*[0,1,0,-1,0]; current_ref = "\SmoothBF";
t_points   = linspace(0, 3, 40); p_points   = 0.5*pi*(squarewave(2/3*pi*t_points - 1)+1); current_ref = "\squareWave";

%% Parameters
Ts_ctrl         = 4e-3; 
% % lambda_g_list   = logspace(0.5, -5, 32);
lambda_g_list   = logspace(0.5, -0.5, 12);

lambda_y        = 5e2;
l_val           = 10;
Current_Pato_Nr = 7;
% 
full_PATO_list  = [7, 10, 15, 1, 11, 12, 14, 17]; %MAKE SURE the first number in list corresponds to current PATO number


combi_list      = nchoosek(full_PATO_list(full_PATO_list~= Current_Pato_Nr),3);

% Choose a random number
[~,idx]=sort(rand(size(combi_list,1),1));
N_exp = 12;
random_ints = idx(1:N_exp);
full_PATO_list       = [Current_Pato_Nr.*ones(N_exp,1), combi_list(random_ints,:)];

for i_pato_list = 1:N_exp
PATO_list = full_PATO_list(i_pato_list,:);
Q               = 1;
R               = [1e-8, 8e-1];
IOscale_init    = [1, 1]; %["norm", "norm"];
RegOrFed        = 'FedGap';

s = num2str(PATO_list); 
s_new = "_"+strrep(strrep(s,'  ','_'),'__','_'); %Generate the correct string type 
ExtraInfo       = s_new+"_V4";
T_switch_val    = 0;

%% Generate data structure
if PATO_list(1) ~= Current_Pato_Nr
    error('You want to set the first index in PATO_list to the Pato to be controlled')
end

[DataStruct, pos_table] = InitializeDeePCFunc_V2(PATO_list, ...
    Ts_ctrl, l_val, Q, R, IOscale_init, RegOrFed, ...
    t_points, p_points);

%% Wrap as Simulink.Parameters
makeParam = @(val) setfield(Simulink.Parameter,'Value',val); %Make function to set values as parameter

lambda_g     = makeParam(lambda_g_list(1)); 
lambda_g.CoderInfo.StorageClass = 'ExportedGlobal';

p_tab        = makeParam(pos_table.signals.values);
t_tab        = makeParam(pos_table.time);
Ts_ctrl_sl   = makeParam(Ts_ctrl);
f_pre_r      = makeParam(DataStruct.f);
H_p          = makeParam(DataStruct.H_p);
U_f          = makeParam(DataStruct.U_f);
Y_f          = makeParam(DataStruct.Y_f);
Ts           = makeParam(DataStruct.Ts);
H            = makeParam(DataStruct.H);
T_ini        = makeParam(DataStruct.T_ini);
N            = makeParam(DataStruct.N);
IOscale      = makeParam(DataStruct.IOscale);
u_bounds     = makeParam(DataStruct.u_bounds);
y_bounds     = makeParam(DataStruct.y_bounds);
l            = makeParam(l_val);
T_switch     = makeParam(T_switch_val);


%% Experiment settings struct (for logging)
exp_settings = struct('Ts_ctrl',Ts_ctrl, ...
    'lambda_g',lambda_g.Value, ...
    'lambda_y',lambda_y, ...
    'l',l.Value, ...
    'Pato_List',PATO_list, ...
    'Q',Q,'R',R, ...
    'IOscale',IOscale.Value);

%% Connect & build model
connect(tg);
answer = 'Yes';
% answer = questdlg('Do you want to rebuild the simulink model? This might take a while.', 'Rebuild Model?', ...
                      % 'Yes', 'No', 'No');  % Default is 'No'
if strcmp(answer, 'Yes')
    slbuild('DeePC_SimulinkV4','ForceTopModelBuild',true);
end
%% START SIMULATION LOOP
for i = 1:length(lambda_g_list)
    disp("NOW RUNNING EXPERIMENT #" + i)

    % Update Simulink.Parameter before experiment
    lambda_g.Value = lambda_g_list(i);
    p_tab.Value = pos_table.signals.values;
    t_tab.Value = pos_table.time;
    % Update bookkeeping struct
    exp_settings.lambda_g = lambda_g.Value;

    % Reload model on target
    load(tg, 'DeePC_SimulinkV4');

    % Update parameter on Speedgoat
    setparam(tg, '','lambda_g', lambda_g.Value);

    % Start execution
    start(tg);
    pause(12.5);   % adjust for experiment duration

    % Save results
    save(exp_path +"\"+ Current_Pato_Nr + current_ref + RegOrFed + ExtraInfo + i, ...
        "logsout","exp_settings","DataStruct");
end
end