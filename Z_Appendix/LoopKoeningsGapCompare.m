clear; close all;

%% PARAMETERS
l = 80;                 % Hankel row size
T = 1000;               % Signal length
Ts = 1;

sigma_var_list = 0:0.015:0.2;

% System: Double Mass spring
m1=3; m2=1; k1=0.1; k2=0.1; d=5;
A=[0,1,0,0;...
   -(k1+k2)/m1,-2*d/m1,k2/m1,d/m1;...
   0,0,0,1;...
   k2/m2,d/m2,-k2/m2,-d/m2];
B=[0;0;0;1/m2]; 
C=[1,0,0,0]; 
D=0;

A_i = @(k_i, m_i) [0,1,0,0;...
               -(k_i+k2)/m_i,-2*d/m_i,k2/m_i,d/m_i;...
               0,0,0,1;...
               k2/m2,d/m2,-k2/m2,-d/m2];

sys1 = c2d(ss(A_i(3*k1,m1), B, C, D), Ts, 'foh');
sys2 = c2d(ss(A,B,C,D),Ts,'foh');

%% Preallocate storage
gap_matlab = zeros(size(sigma_var_list));
gap_padoan = zeros(size(sigma_var_list));
gap_angle  = zeros(size(sigma_var_list));
gap_koen   = zeros(size(sigma_var_list));

%% MAIN LOOP
for idx = 1:length(sigma_var_list)

    sigma_var = sigma_var_list(idx);
    fprintf("Running noise σ = %.3f\n", sigma_var);

    %% Input
    u_input = rand(1,T);

    %% Clean trajectories
    y1 = lsim(sys1, u_input, 1:Ts:T);
    y2 = lsim(sys2, u_input, 1:Ts:T);

    %% Add noise
    y1n = y1 + sigma_var * randn(T,1);
    y2n = y2 + sigma_var * randn(T,1);

    %% Build w = [u; y]
    w1 = [u_input; y1n'];
    w2 = [u_input; y2n'];

    %% Full Hankels (for Padoan & angle)
    H1_L = traj2Hankel(w1, l);
    H2_L = traj2Hankel(w2, l);

    %% Block Hankels (for Koenings)
    T_ini = l/2; N = l/2;

    [Up1, Yp1, Uf1, Yf1] = traj2Hankel(w1, T_ini, N);
    [Up2, Yp2, Uf2, Yf2] = traj2Hankel(w2, T_ini, N);

    H1 = [Up1; Yp1; Uf1; Yf1];
    H2 = [Up2; Yp2; Uf2; Yf2];

    %% ===== MATLAB gap =====
    gap_matlab(idx) = gapmetric(sys1, sys2);

    %% ===== Padoan Gap =====
    gap_padoan(idx) = PadoanGap(w1, w2, l);

    %% ===== Subspace angle gap =====
    gap_angle(idx) = sin(subspace(H1_L, H2_L));

    %% ===== Koenings Gap =====

    % system 1
    [L1,Q1] = lq_decomp(H1);
    [L11_1, L21_1, L22_1, L31_1, L32_1] = split_lq(L1, Q1, T_ini, N);
    L_U1 = L32_1 * pinv(L22_1);     % CORRECT multiplication form
    Idf1 = [eye(size(L_U1,2)); L_U1];
    [U1_1,~,~,~,~,~] = splitSVD(Idf1, 0.1, "tol");
    U1 = U1_1;

    % system 2
    [L2,Q2] = lq_decomp(H2);
    [L11_2, L21_2, L22_2, L31_2, L32_2] = split_lq(L2, Q2, T_ini, N);
    L_U2 = L32_2 * pinv(L22_2);
    Idf2 = [eye(size(L_U2,2)); L_U2];
    [U1_2,~,~,~,~,~] = splitSVD(Idf2, 0.1, "tol");
    U2 = U1_2;

    % Koenings symmetric projection gap
    d1 = norm(U1 - U2*U2'*U1);
    d2 = norm(U2 - U1*U1'*U2);
    gap_koen(idx) = max(d1, d2);

end

%% ============================================================
%% PLOTS
%% ============================================================

gapfig = figure; hold on; grid on; box on;
plot(sigma_var_list, gap_matlab, 'k--', 'LineWidth',2, 'DisplayName','$\delta_{g}$');
plot(sigma_var_list, gap_padoan, 'o-', 'LineWidth',2, 'DisplayName','$\delta_{pad}$');
plot(sigma_var_list, gap_angle,  'o-','LineWidth',2, 'DisplayName','$\delta_{\theta}$');
plot(sigma_var_list, gap_koen,  'o-', 'LineWidth',2, 'DisplayName','$\delta_{koen}$');

xlabel(' $(\sigma_e)^2$','Interpreter','latex')
ylabel('$\delta_{i}$','Interpreter','latex')
% title('Gap estimation vs noise','Interpreter','latex')
legend('Location','northwest')
tightfig(gapfig)


% Function
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