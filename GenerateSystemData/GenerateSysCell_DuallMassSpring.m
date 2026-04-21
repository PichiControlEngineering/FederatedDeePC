%% Generate N mass-spring-damper systems with random coefficients and save them in a big cell
N_sys = 20;  % Amount of systems in each sys_cell
N_repeats = 1000; % Amount of total system sets, in total we generate N_repeats nominal systems, and N_sys-1 addiditional systems for each nominal system
show_progress = 0; % Progress counter for debugging
%%

tic()
sys_cell_cell = cell(N_repeats, 1); GapMatCell = cell(N_repeats, 1);
disp('Generating systems...')
for i = 1:N_repeats
        [sys_cell_i, GapMat_i] = GenerateRandomMassSpring(N_sys, show_progress, ...
                            'MaxGap', 0.6, 'MinGap', 0);
    sys_cell_cell{i} = sys_cell_i;
    GapMatCell{i} = GapMat_i;
    disp(i)
end
    t_sim = toc();

disp("time per generaterated system: " + t_sim/(N_sys*N_repeats))





%% Function to generate the random mass spring
function [sys_list, varargout] = GenerateRandomMassSpring(n_systems, show_progress, varargin)
%% GenerateRandomSys
% Generates a list of stable random state-space systems with constraints 
% on the gap metric relative to a nominal system.
%
% Inputs:
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
sys_list = cell(n_systems, 1);
GapMat   = zeros(n_systems);
KGapMat  = zeros(n_systems);

%% Parse optional inputs
if isempty(show_progress)
    show_progress = 0;
end
m1 = 3;
m2 = 3;
k1 = 1;
k2 = 1.5;
%% Koenig's gap option
u_multSin = 0;
Use_KGap = strcmp(varargin{end}, 'KGap');

if Use_KGap
    u_multSin = MultiSineGenerator(400, 1/20, 1/3, 1)';
    varargin = varargin(1:end-1);
    GapMat  = zeros(n_systems, 1);
    KGapMat = zeros(n_systems, 1);
end

p = inputParser;
addParameter(p, 'MaxGap', []);
addParameter(p, 'MinGap', []);
addParameter(p, 'Sys_nom', []);
parse(p, varargin{:});

MaxGap  = p.Results.MaxGap;
MinGap  = p.Results.MinGap;
Sys_nom = p.Results.Sys_nom;

% If nominal system is provided, assign it to sys_list{1}
if ~isempty(Sys_nom)
    sys_list{1} = Sys_nom;
else
    sys_list{1} = randomMassSpring(m1, m2, k1, k2);
end

%% Step 2: Generate stable systems
for i_sys = 1:n_systems
    % Generate stable random system
    ss_i = randomMassSpring(m1, m2, k1, k2);
    
    % Enforce gap constraints if required
    if ~isempty(MaxGap) && i_sys > 1
        k = 1;
        gap_k = 1;
        
        if all(abs(pole(ss_i)) < 1)
            gap_k = computeGap(sys_list{1}, ss_i, Use_KGap, u_multSin);
        end
        
        while ~(all(abs(pole(ss_i)) < 1) && gap_k < MaxGap && gap_k > MinGap)
            ss_i = randomMassSpring(m1,m2, k1, k2);
            
            if show_progress
                disp("System: " + i_sys + " run: " + k + " gap = " + gap_k)
            end
            
            if all(abs(pole(ss_i)) < 1)
                gap_k = computeGap(sys_list{1}, ss_i, Use_KGap, u_multSin);
            else
                gap_k = 1;
                % disp('Systems Unstable!')
            end
            
            k = k + 1;
            
            % Restart entire function if too many iterations
            if k > 10000
                warning("Exceeded 10000 iterations. Rerunning GenerateRandomSys with new system...")
                A_i = randomMassSpring(m1,m2, k1, k2);
                return
            end
        end
    end
   
    % Check controllability & observability
    if any([rank(ctrb(ss_i)), rank(obsv(ss_i))] < 4)
        warning('System is unobservable/uncontrollable')
    end

    % Store system
    sys_list{i_sys} = ss_i;
    
    % Compute gap matrices
    It_Max = Use_KGap * i_sys + (~Use_KGap) * 1;
    for j = 1:It_Max
        if ~Use_KGap
            GapMat(i_sys, j) = gapmetric(sys_list{1}, ss_i);
        else
            KGapMat(i_sys, j) = ComputeKGap(sys_list{1}, ss_i, @(x,u,t) u_multSin(t));
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
function gap_k = computeGap(Sys_1, Sys_2, Use_KGap, u_multSin)
    if ~Use_KGap
        gap_k = gapmetric(Sys_1, Sys_2);
    else
        gap_k = ComputeKGap(Sys_1, Sys_2, @(x,u,t) u_multSin(t));
    end
end


function gap = ComputeKGap(sys_1, sys_2, u_in)
    x0 = zeros(size(sys_1.A,1), 1);
    [x_1, u_1] = stepStates(x0, 300, sys_1.A, sys_1.B, 1, [], u_in);
    [x_2, u_2] = stepStates(x0, 300, sys_2.A, sys_2.B, 1, [], u_in);
    gap = KoeningsGap([u_1; sys_1.C*x_1], [u_2; sys_2.C*x_2], 50, 50);
    % gap = PadoanGap([u_1; sys_1.C*x_1], [u_2; sys_2.C*x_2], 4);
end

function ss_disc = randomMassSpring(m1_avg, m2_avg, k1_avg, k2_avg)
    variance = 0.5;
    m1 = m1_avg*(1+variance*(2*rand(1)-1));  
    m2 = m2_avg*(1+variance*(2*rand(1)-1));  
    k1 = k1_avg*(1+variance*(2*rand(1)-1)); 
    k2 = k2_avg*(1+variance*(2*rand(1)-1)); 
    d = 0.5;
    
    A_i=[0,1,0,0;...
       -(k1+k2)/m1,-2*d/m1,k2/m1,d/m1;...
       0,0,0,1;...
       k2/m2,d/m2,-k2/m2,-d/m2];
    B_i = [0;0;0;1/m2];
    C = [1,0,0,0]; 
    D = 0;
    ss_disc = c2d(ss(A_i,B_i,C,D),1);

end