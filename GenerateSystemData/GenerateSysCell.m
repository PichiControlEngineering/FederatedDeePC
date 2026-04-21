%% Generate N random systems and save them in a cell, these cells are used for monte-carlo simulation
N_order = 2; % Model order
N_sys = 20;  % Amount of systems in each sys_cell
N_repeats = 1000; % Amount of total system sets, in total we generate N_repeats nominal systems, and N_sys-1 addiditional systems for each nominal system
show_progress = 0; % Progress counter for debugging
%%

tic()
sys_cell_cell = cell(N_repeats, 1); GapMatCell = cell(N_repeats, 1);
disp('Generating systems...')
for i = 1:N_repeats
        [sys_cell_i, GapMat_i] = GenerateRandomSys(N_order, N_sys, show_progress, ...
                            'MaxGap', 0.7, 'MinGap', 0, "TrueRandom", 1);
    sys_cell_cell{i} = sys_cell_i;
    GapMatCell{i} = GapMat_i;
    disp(i)
end
    t_sim = toc();

disp("time per generaterated system: " + t_sim/(N_sys*N_repeats))



function [sys_list, varargout] = GenerateRandomSys(n_order, n_systems, show_progress, varargin)
%% GenerateRandomSys
% Generates a list of stable random state-space systems with constraints 
% on the gap metric relative to a nominal system.
%
% Inputs:
%   n_order      - Order of the systems
%   n_systems    - Number of systems to generate
%   show_progress- Boolean flag for progress display
%   varargin     - Optional parameters:
%       'MaxGap' : Upper bound on allowed gap
%       'MinGap' : Lower bound on allowed gap
%       'Sys_nom': Nominal system
%       'KGap'   : Use Koenig's gap metric (slower)
%
% Outputs:
%   sys_list     - Cell array of generated systems
%   varargout{1} - Gap matrix (GapMat or KGapMat)

%% Step 1: Generate system matrices B, C, D
B = zeros(n_order, 1); 
C = zeros(1, n_order); 
B(1) = 1; 
C(end) = 1; 
D = 0;

sys_list = cell(n_systems, 1);
GapMat   = zeros(n_systems);
KGapMat  = zeros(n_systems);

%% Parse optional inputs
if isempty(show_progress)
    show_progress = 0;
end

%% Koenig's gap option
u_multSin = 0;
Use_KGap = strcmp(varargin{end}, 'KGap');

if Use_KGap
    u_multSin = MultiSineGenerator(400, 1/20, 1/2, 1)';
    varargin = varargin(1:end-1);
    GapMat  = zeros(n_systems, 1);
    KGapMat = zeros(n_systems, 1);
end

p = inputParser;
addParameter(p, 'MaxGap', []);
addParameter(p, 'MinGap', []);
addParameter(p, 'Sys_nom', []);
addParameter(p, 'TrueRandom', 0);
parse(p, varargin{:});

MaxGap  = p.Results.MaxGap;
MinGap  = p.Results.MinGap;
Sys_nom = p.Results.Sys_nom;
TrueRandom = p.Results.TrueRandom;

%% Step 2: Generate stable systems
A_stab_list = cell(1, n_systems);
for i_sys = 1:n_systems
    
    % Generate stable random system
    A_i =  randMat(n_order);
    while any(abs(eig(A_i)) > 1)
        A_i = randMat(n_order);
    end
    
    % If nominal system is provided, assign it to sys_list{1}
    if ~isempty(Sys_nom)
        sys_list{1} = Sys_nom;
    else
        Sys_nom = ss(A_i, B, C, D, 1);
    end
    
    % Enforce gap constraints if required
    if ~isempty(MaxGap) && i_sys > 1
        k = 1;
        gap_k = 1;
        
        if all(abs(eig(A_i)) < 1)
            gap_k = computeGap(Sys_nom, A_i, B, C, D, Use_KGap, u_multSin);
        end
        
        while any(~[all(abs(eig(A_i)) < 1), gap_k < MaxGap, gap_k > MinGap])
            if TrueRandom
                A_i = randMat(n_order);
            else
                A_i = Sys_nom.A + 0.5.*randn(n_order);
            end
            if show_progress
                disp("System: " + i_sys + " run: " + k + " gap = " + gap_k)
            end
            
            if all(abs(eig(A_i)) < 1)
                gap_k = computeGap(Sys_nom, A_i, B, C, D, Use_KGap, u_multSin);
            else
                gap_k = 1;
            end
            
            k = k + 1;
            
            % Restart entire function if too many iterations
            if k > 10000
                warning("Exceeded 10000 iterations. Rerunning GenerateRandomSys with new system...")
                Sys_nom = ss(randMat(n_order), B, C, D, 1);
                [sys_list, varargout{:}] = GenerateRandomSys(n_order, n_systems, show_progress, "MaxGap", MaxGap, "MinGap", MinGap, "Sys_nom", Sys_nom);
                return
            end
        end
    end
    
    % Save stable A
    A_stab_list{i_sys} = A_i;
    
    % Check controllability & observability
    if any([rank(ctrb(A_i, B)), rank(obsv(A_i, C))] < n_order)
        warning('System is unobservable/uncontrollable')
    end
    
    % Store system
    sys_list{i_sys} = ss(A_i, B, C, D, 1);
    
    % Compute gap matrices
    It_Max = Use_KGap * i_sys + (~Use_KGap) * 1;
    for j = 1:It_Max
        if ~Use_KGap
            GapMat(i_sys, j) = gapmetric(sys_list{1}, ss(A_i, B, C, D, 1));
        else
            KGapMat(i_sys, j) = ComputeKGap(sys_list{1}, ss(A_i, B, C, D, 1), @(x,u,t) u_multSin(t));
        end
    end
    
    % Output GapMat / KGapMat if requested
    if nargout > 1
        if Use_KGap
            varargout{1} = KGapMat;
        else
            varargout{1} = GapMat;
        end
    end
end
end


%% Helper Functions
function gap_k = computeGap(Sys_nom, A_i, B, C, D, Use_KGap, u_multSin)
    if ~Use_KGap
        gap_k = gapmetric(Sys_nom, ss(A_i, B, C, D, 1));
    else
        gap_k = ComputeKGap(Sys_nom, ss(A_i, B, C, D, 1), @(x,u,t) u_multSin(t));
    end
end


function gap = ComputeKGap(sys_1, sys_2, u_in)
    x0 = zeros(size(sys_1.A,1), 1);
    [x_1, u_1] = stepStates(x0, 200, sys_1.A, sys_1.B, 1, [], u_in);
    [x_2, u_2] = stepStates(x0, 200, sys_2.A, sys_2.B, 1, [], u_in);
    gap = KoeningsGap([u_1; sys_1.C*x_1], [u_2; sys_2.C*x_2], 50, 50);
    % gap = PadoanGap([u_1; sys_1.C*x_1], [u_2; sys_2.C*x_2], 4);
end

function RandomMat = randMat(n_o)
RandomMat = rand(n_o); % Generate equally distributed random sys
% RandomMat = RandomMat./(1.1*norm(RandomMat));

end