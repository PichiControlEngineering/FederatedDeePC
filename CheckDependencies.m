% Check if all required components are installed

if ~exist("DeePCController.m", 'file')
disp("Adding dependencies to matlab path...")
addpath(genpath(pwd))
disp('...Dependencies added successfully.');
end
check_dependencies(0); % Set to "1" to check dependencies for replicating experiments

function missing_total = check_dependencies(include_experimental)
%CHECK_DEPENDENCIES Check for required MATLAB toolboxes and third-party tools.
%   MISSING = CHECK_DEPENDENCIES(INCLUDE_EXPERIMENTAL) returns a cell array of
%   missing components. If INCLUDE_EXPERIMENTAL is true, experimental (optional)
%   products are also checked. Default is false.
%
%   Required:
%     - Robust Control Toolbox
%     - Control System Toolbox
%     - Symbolic Math Toolbox
%     - Statistics and Machine Learning Toolbox
%     - Optimization Toolbox
%
%   Experimental (optional):
%     - Simulink
%     - Simulink Real-Time
%     - Simulink Compiler
%     - MATLAB Compiler
%
%   Third-party:
%     - CVX (checked by existence of cvx_setup or cvx_install or cvx_begin)
%
%   Example:
%     missing = check_dependencies(true);

if nargin < 1 || isempty(include_experimental)
    include_experimental = false;
end

% List of MathWorks products to check (name, toolbox identifier from ver)
required_products = {
    'Robust Control Toolbox'             , 'Robust Control Toolbox';
    'Control System Toolbox'             , 'Control System Toolbox';
    'Statistics and Machine Learning Toolbox', 'Statistics and Machine Learning Toolbox';
    'Optimization Toolbox'               , 'Optimization Toolbox'
    };

experimental_products = {
    'Simulink'                           , 'Simulink';
    'Simulink Real-Time'                 , 'Simulink Real-Time';
    'Simulink Compiler'                  , 'Simulink Compiler';
    'MATLAB Compiler'                    , 'MATLAB Compiler'
    };

% Get installed products once
v = ver;
installedNames = {v.Name};

missing = {};
missing_experimental = {};

% Helper to check presence (case-sensitive match of Name)
isInstalled = @(id) any(strcmp(installedNames, id));

% Check required products
for k = 1:size(required_products,1)
    name = required_products{k,1};
    id   = required_products{k,2};
    if ~isInstalled(id)
        missing{end+1,1} = name; %#ok<AGROW>
    end
end

% Check experimental if requested
if include_experimental
    for k = 1:size(experimental_products,1)
        name = experimental_products{k,1};
        id   = experimental_products{k,2};
        if ~isInstalled(id)
            missing_experimental{end+1,1} = name; %#ok<AGROW>
        end
    end
end

% Check third-party: CVX
% We consider CVX present if any of common entry points exist on path.
cvxPresent = ~isempty(which('cvx_setup')) || ~isempty(which('cvx_install')) || ~isempty(which('cvx_begin'));
if ~cvxPresent
    missing{end+1,1} = 'CVX (third-party)';
end

% Display results
if isempty([missing, missing_experimental])
    fprintf('All required components are installed.\n');
else
    fprintf('Missing components (%d):\n', numel(missing));
    for i = 1:numel(missing)
        fprintf('  - %s\n', missing{i});
    end

    fprintf('Missing components for experimental interfacing (%d):\n', numel(missing_experimental));
    for i = 1:numel(missing_experimental)
        fprintf('  - %s\n', missing_experimental{i});
    end
end
missing_total = [missing, missing_experimental];

%Check if all internal functions exist
% Check for internal functions and add to missing list if not found
internalFunctions = {'DeePCController.m',...
                     'DeePCController_Du_Cost.m',...
                      'AddNoise.m',...
                      'GenerateIOData.m',...
                      'square.m',...
                      'squarewave.m',...
                      'StepStates.m',...
                      'Traj2Hankel.m',...
                      'multisine.m',...
                      'KoeningsGap.m',...
                      'PadoanGap.m',...
                      'subspacea.m',...
                      'PadoanGap.m',...
                      'lq_decomp.m',...
                      'PredictMat.m',...
                      'SaveControllerData.m',...
                      'SplitLQ.m',...
                      'SplitSVD.m',...
                      'PlotFunction.m',...
                      'PlotSettings.m',...
                      'tightfig.m',...
                      }; % Replace with actual function names
    for j = 1:numel(internalFunctions)
        if isempty(which(internalFunctions{j}))
            missing{end+1,1} = sprintf('Internal function: %s', internalFunctions{j}); %#ok<AGROW>
        end
    end
end