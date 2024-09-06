
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                               %
%   This code computes the passive-only endemic equilibrium of the Warwick HAT model                            %
%                                                                                                               %
%   Inputs:                                                                                                     %
%       N_H - number denoting human population size used in ODE                                                 %
%       Paras - structure containing location-specific parameters (fixed, fitted and intervention parameters)   %
%                                                                                                               %
%   Outputs:                                                                                                    %
%       meff - number denoting effective vector density, calculated by the relation Paras.R0^2 ~ meff           %
%       ICs - cell array containing endemic equilibrium of given set of parameters                              %
%                                                                                                               %
%   Note: hosts are (1) low-risk, random participants                                                           %
%                   (2) high-risk, random participants                                                          %
%                   (3) low-risk, non-participants                                                              %
%                   (4) high-risk, non-participants                                                             %
%                   (5) reservoir animals                                                                       %
%                   (6) non-reservoir animals, no dynamics and is ignored                                       %
%                                                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [meff, ICs] = GetEndemicEq(N_H, Paras)
    
    k4 = 1 - Paras.k1 - Paras.k2 - Paras.k3;
    Plow = Paras.k1 + Paras.k3;
    Phigh = Paras.k2 + k4;

    N_Hlow = N_H * Plow;
    N_Hhigh = N_H * Phigh;
    f_Hlow = Paras.f_H * Plow / (Plow + Paras.r * Phigh);
    f_Hhigh = Paras.f_H * Paras.r * Phigh / (Plow + Paras.r * Phigh);

    N_A = N_H * Paras.k_A;

    if N_Hhigh == 0
        N_Hhigh = 1;
        f_Hhigh = 0;
    end
    if N_A == 0
        N_A = 1;
        Paras.f_A = 0;
    end

    eta_H0 = Paras.eta_H * Paras.b_eta_H0;
    gamma_H0 = Paras.gamma_H * Paras.b_gamma_H0;


%%% Compute m_eff using NGM approach, where the order of matrix elements is
%%% (E_H, I1_H, I2_H) for low risk, (E_H, I1_H, I2_H) for high risk, E_A, I1_A, E1_V, E2_V, E3_V, I_V
    T = zeros(12,12);
    S = zeros(12,12);
    S_Vstar = Paras.mu_V * N_H / (Paras.alpha + Paras.mu_V); % equilibrium S_V
    G_Vfrom0 = Paras.alpha * N_H / (Paras.alpha + Paras.mu_V); % equilibrium (no infection) G_V

    meff = 1; % assign a value for computing the corresponding R_0

    %T is transmissions
    T(1,12) = Paras.alpha * meff * f_Hlow; % I_V infects E_Hlow
    T(4,12) = Paras.alpha * meff * f_Hhigh; % I_V infects E_Hhigh
    T(7,12) = Paras.alpha * meff * Paras.f_A; % I_V infects E_A
    T(9,2) = Paras.alpha * f_Hlow * Paras.p_V * (S_Vstar + Paras.epsilon * G_Vfrom0) / N_Hlow; % I1_Hlow infects E1_V
    T(9,3) = T(9,2); % I2_Hlow infects E1_V
    T(9,5) = Paras.alpha * f_Hhigh * Paras.p_V * (S_Vstar + Paras.epsilon * G_Vfrom0) / N_Hhigh; % I1_Hhigh infects E1_V
    T(9,6) = T(9,5); % I2_Hhigh infects E1_V
    T(9,8) = Paras.alpha * Paras.f_A * Paras.p_V * (S_Vstar + Paras.epsilon * G_Vfrom0) / N_A; % I1_A infects E1_V

    %S is transissions (including passive stage1 detection)
    S(1,1) = - Paras.sigma_H - Paras.mu_H; % leaves E_Hlow
    S(2,1) = Paras.sigma_H; % enters I1_Hlow from E_Hlow
    S(2,2) = - eta_H0 - Paras.phi_H - Paras.mu_H; % leaves I1_Hlow
    S(3,2) = Paras.phi_H; % enters I2_Hlow from I1_Hlow
    S(3,3) = - gamma_H0 - Paras.mu_H; % leaves I2_Hlow
    S(4,4) = - Paras.sigma_H - Paras.mu_H; % leaves E_Hhigh
    S(5,4) = Paras.sigma_H; % enters I1_Hhigh from E_Hhigh
    S(5,5) = - eta_H0 - Paras.phi_H - Paras.mu_H; % leaves I1_Hhigh
    S(6,5) = Paras.phi_H; % enters I2_Hhigh from I1_Hhigh
    S(6,6) = - gamma_H0 - Paras.mu_H; % leaves I2_Hhigh
    S(7,7) = - Paras.sigma_A - Paras.mu_A; % leaves E_A
    S(8,7) = Paras.sigma_A; % enters I_A from E_A
    S(8,8) = - Paras.mu_A; % leaves I1_A (no treatment and no staging) 
    S(9,9) = - 3 * Paras.sigma_V - Paras.mu_V; % leaves E1_V
    S(10,9) =  3 * Paras.sigma_V; % enters E2_V from E1_V
    S(10,10) = - 3 * Paras.sigma_V - Paras.mu_V; % leaves E2_V
    S(11,10) = 3 * Paras.sigma_V; % enters E3_V from E2_V
    S(11,11) = - 3 * Paras.sigma_V - Paras.mu_V; % leaves E3_V
    S(12,11) = 3 * Paras.sigma_V; % enters I_V from E3_V
    S(12,12) = - Paras.mu_V; % leaves I_V

    K = - T * inv(S);
    R0_current = max(abs(eig(K)));
    meff=(Paras.R0/R0_current)^2 * meff;


%%% Compute Endemic Equilibrium
    A = Paras.omega_H * (eta_H0 + (gamma_H0 * Paras.phi_H) / (gamma_H0 + Paras.mu_H)) * Paras.sigma_H / (Paras.mu_H * (Paras.omega_H + Paras.mu_H) * (eta_H0 + Paras.phi_H + Paras.mu_H)) - Paras.sigma_H / Paras.mu_H - 1;
    B = 27 * Paras.sigma_V^3 / (Paras.mu_V * (3 * Paras.sigma_V + Paras.mu_V)^2);
    C = (Paras.sigma_H / (eta_H0 + Paras.phi_H + Paras.mu_H)) * (1 + (Paras.phi_H / (gamma_H0 + Paras.mu_H)));
    D = Paras.sigma_H + Paras.mu_H;
    E = 3 * Paras.sigma_V + Paras.mu_V;
    F = Paras.alpha * Paras.p_V * Paras.epsilon * (1 + 3 * Paras.sigma_V / Paras.mu_V) * C;
    G = Paras.alpha * Paras.p_V * N_H * C * (Paras.mu_V + Paras.epsilon *Paras.alpha) / (Paras.alpha + Paras.mu_V);
    H = (1 + Paras.sigma_A / Paras.mu_A);
    J = Paras.sigma_A + Paras.mu_A;
    K = Paras.alpha * meff * B * H;
    L = Paras.alpha * Paras.p_V * Paras.epsilon * (1 + 3 * Paras.sigma_V / Paras.mu_V) * Paras.sigma_A / Paras.mu_A;
    M = Paras.alpha * Paras.p_V * N_H * Paras.sigma_A * (Paras.mu_V + Paras.epsilon * Paras.alpha) / ((Paras.alpha + Paras.mu_V) * Paras.mu_A);

    a = Paras.alpha * meff * B;
    b = Paras.alpha * meff * A * B;

    Z1 = b * f_Hlow * f_Hhigh * Paras.f_A * Paras.relprob * (-E * K * b + F * a * K * (f_Hlow + f_Hhigh) - L * a * b * Paras.f_A);

    Z2 = J * b * f_Hlow * f_Hhigh * N_A * (-E * b + F * a * (f_Hlow + f_Hhigh))...
       + Paras.f_A * Paras.relprob * (D * b * (f_Hlow * N_Hhigh + f_Hhigh * N_Hlow) * (E * K + a * L * Paras.f_A) - a * F * D * K * (f_Hlow^2 * N_Hhigh + f_Hhigh^2 * N_Hlow))...
       + a * b * f_Hlow * f_Hhigh * Paras.f_A * Paras.relprob * (b * M * Paras.f_A - G * K * (f_Hlow + f_Hhigh));

    Z3 = Paras.f_A * Paras.relprob * D^2 * (-E * K - a * L * Paras.f_A) * N_Hlow * N_Hhigh...
       + D * b * (f_Hlow * N_Hhigh * N_A + f_Hhigh * N_Hlow * N_A) * (E * J - M * a * Paras.f_A^2 * Paras.relprob / N_A)...
       + D * a * (f_Hlow^2 * N_Hhigh * N_A + f_Hhigh^2 * N_Hlow * N_A) * (-F * J + G * K * Paras.f_A * Paras.relprob / N_A)...
       - G * J * b * a * f_Hlow * f_Hhigh * (f_Hlow + f_Hhigh) * N_A;

    Z4 = -E * J * D^2 * N_Hlow * N_Hhigh * N_A + a * G * D * J * (f_Hlow^2 * N_Hhigh * N_A + f_Hhigh^2 * N_Hlow * N_A)... 
       + a * M * D^2 * Paras.f_A^2 * Paras.relprob * N_Hlow * N_Hhigh;

    E1_V = max([0; roots([Z1 Z2 Z3 Z4])]);

    E_Hlowstar = a * f_Hlow * E1_V / (D - b * f_Hlow * E1_V / N_Hlow);
    E_Hhighstar= a * f_Hhigh * E1_V / (D - b * f_Hhigh * E1_V / N_Hhigh);
    % splits low risk/high risk groups by participation
    E_H1 = Paras.k1 / Plow * E_Hlowstar;
    E_H3 = Paras.k3 / Plow * E_Hlowstar;
    E_H2 = Paras.k2 / Phigh * E_Hhighstar;
    E_H4 = k4 / Phigh * E_Hhighstar;

    E_H = [E_H1, E_H2, E_H3, E_H4];
    E_A = a * Paras.f_A * Paras.relprob * E1_V / (J + K * Paras.f_A * Paras.relprob * E1_V / N_A);

    % computes other equilibria from E's
    S_H = N_H * [Paras.k1, Paras.k2, Paras.k3, k4] + A * E_H;
    I1_H = (Paras.sigma_H / (eta_H0 + Paras.phi_H + Paras.mu_H)) * E_H;
    I2_H = (Paras.phi_H / (gamma_H0 + Paras.mu_H)) * I1_H;
    R_H = (eta_H0 * I1_H + gamma_H0 * I2_H) / (Paras.omega_H + Paras.mu_H);

    S_A = N_H * Paras.k_A - H * E_A;
    I_A = Paras.sigma_A / Paras.mu_A * E_A;

    S_V = S_Vstar;
    G_V = Paras.alpha * S_V / Paras.mu_V - ((3 * Paras.sigma_V + Paras.mu_V) / Paras.mu_V) * E1_V;
    E2_V = (3 * Paras.sigma_V / (3 * Paras.sigma_V + Paras.mu_V)) * E1_V;
    E3_V = (3 * Paras.sigma_V / (3 * Paras.sigma_V + Paras.mu_V)) * E2_V;
    I_V = (3 * Paras.sigma_V / Paras.mu_V) * E3_V;
    P_V = (Paras.alpha + Paras.mu_V) * S_V / (Paras.xi_V * Paras.p_survive);

    ICs = {S_H, E_H, I1_H, I2_H, R_H, S_A, E_A, I_A, P_V, S_V, G_V, E1_V, E2_V, E3_V, I_V};
end
