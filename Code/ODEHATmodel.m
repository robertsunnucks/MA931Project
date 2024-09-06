
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                               %
%   This code solves ordinary differential equations of the Warwick HAT model                                   %
%                                                                                                               %
%   Inputs:                                                                                                     %
%       meff - number denoting effective vector density, calculated by the relation Paras.R0^2 ~ meff           %
%       ICs - cell array containing initial conditions for ODE                                                  %
%       Data - structure containing location-specific historical data                                           %
%       Paras - structure containing location-specific parameters (fixed, fitted and intervention parameters)   %
%       ProjStrat - structure containing parameters associated with future strategy                             %
%                                                                                                               %
%   Outputs:                                                                                                    %
%       Classes - table containing time series of model outputs (e.g. susceptible humans, infectious vectors)   %
%       Aggregate - table containing yearly aggregated outputs (e.g. active/passive stage 1/2 cases, deathss)   %
%                                                                                                               %
%   Note: hosts are (1) low-risk, random participants                                                           %
%                   (2) high-risk, random participants                                                          %
%                   (3) low-risk, non-participants                                                              %
%                   (4) high-risk, non-participants                                                             %
%                   (5) reservoir animals                                                                       %
%                   (6) non-reservoir animals, no dynamics and is ignored                                       %
%                                                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Classes, Aggregate] = ODEHATmodel(meff, ICs, Data, Paras, ProjStrat)
    [S_H, E_H, I1_H, I2_H, R_H, S_A, E_A, I_A, P_V, S_V, G_V, E1_V, E2_V, E3_V, I_V] = ICs{:};
    
    k4 = 1 - Paras.k1 - Paras.k2 - Paras.k3;
    N_A = Data.N_H * Paras.k_A;
    K_V = Data.N_H * Paras.k_V;
    
    k=[Paras.k1 Paras.k2 Paras.k3 k4 Paras.k_A 1];
    s_A=Paras.f_A*((Paras.k1+Paras.k3)+Paras.r*(Paras.k2+k4))*(1+(1-Paras.f_A-Paras.f_H)/(Paras.f_A+Paras.f_H))/(Paras.k_A*(1-Paras.f_A)*(1-Paras.f_A*(1-Paras.f_A-Paras.f_H)/((1-Paras.f_A)*(Paras.f_A+Paras.f_H))));
    s_N=(1-Paras.f_A-Paras.f_H)*((Paras.k1+Paras.k3)+Paras.r*(Paras.k2+k4)+Paras.k_A*s_A)/(Paras.f_A+Paras.f_H);
    bitepref=[1 Paras.r 1 Paras.r s_A s_N];            %biting preference on hosts (given no animal reservoir)  
    f=(bitepref.*k)/sum(bitepref.*k);
    f(6)=[];
    bitepref(6)=[];
    
    NumberScreening = length(Data.ModelScreeningTime);

    % Active screening: will move random participating Infectious (E_H and I1_H and I2_H) to Recovery at the end of the year
    % Default AS strategy (screen Paras.k1 and Paras.k2 only)
    Nparticipant = Data.N_H * (Paras.k1 + Paras.k2);
    TurnOut1 = Data.ModelPeopleScreened / Nparticipant;
    TurnOut2 = TurnOut1; % same turn out rate as Paras.k1
    TurnOut3 = zeros(1, NumberScreening); % non-participating
    TurnOut4 = TurnOut3; % non-participating
    
    % New AS strategy is inctorduced
    if Data.ModelScreeningTime(1) < ProjStrat.NewASyear && ProjStrat.NewASyear < Data.ModelScreeningTime(end)
        Y = find(Data.ModelScreeningTime == ProjStrat.NewASyear);
        
        switch ProjStrat.NewASstrat
        case 'equal' % door-to-door (all populations get equal probability to be screened)
            TurnOut1(Y:end) = Data.ModelPeopleScreened(Y:end) / Data.N_H;
            TurnOut2(Y:end) = TurnOut1(Y:end);
            TurnOut3(Y:end) = TurnOut1(Y:end);
            TurnOut4(Y:end) = TurnOut1(Y:end);
        case 'high' % work place screening (k4 group gets screened first and then equally screened Paras.k1 and Paras.k2)
            TurnOut4(Y:end) = min(Data.ModelPeopleScreened(Y:end) / Data.N_H / k4, 1);
            TurnOut1(Y:end) = max((Data.ModelPeopleScreened(Y:end) / Data.N_H - k4) / (Paras.k1 + Paras.k2), 0);
            TurnOut2(Y:end) = TurnOut1(Y:end);
        end
    end
    TurnOut = [TurnOut1; TurnOut2; TurnOut3; TurnOut4];
    
    % Specificity by screening
    ScreeningSpecificity = repmat(Paras.specificity, 1, NumberScreening);
    % Lower specificity can exist in some years (eg MSF interventions in
    % Ango and Ganga HZ of Orientale former province).
    altSpec = find(Data.ModelScreeningTime <= Paras.Last_year);
    ScreeningSpecificity(altSpec) = ScreeningSpecificity(altSpec) * Paras.b_specificity;
    % Also change sensitivity
    ScreeningSensitivity = repmat(Paras.Sensitivity, 1, NumberScreening);
    ScreeningSensitivity(altSpec) = Paras.SensitivityMSF;
    
    % Passive screening: move Infected (I1 and I2) to Recovery continuously
    if Data.ModelScreeningTime(1) < ProjStrat.RDTyear && ProjStrat.RDTyear < Data.ModelScreeningTime(end)
        Y = find(Data.ModelScreeningTime == ProjStrat.RDTyear);
    else
        Y = length(Data.ModelScreeningTime) + 1;
    end
    yearlyeta_H = [(1 + Paras.eta_H_amp ./ (1 + exp(-Paras.d_steep * (double(Data.ModelScreeningTime(1:Y-1)) - (Paras.d_change+Paras.eta_H_lag))))) * Paras.eta_H,...
                   (1 + ProjStrat.RDTincrease) * (1 + Paras.eta_H_amp ./ (1 + exp(-Paras.d_steep * (double(Data.ModelScreeningTime(Y-1)) - (Paras.d_change+Paras.eta_H_lag))))) * Paras.eta_H * ones(1, NumberScreening - (Y-1))];
    yearlygamma_H = [(1 + Paras.gamma_H_amp ./ (1 + exp(-Paras.d_steep * (double(Data.ModelScreeningTime(1:Y-1)) - Paras.d_change)))) * Paras.gamma_H,...
                   (1 + ProjStrat.RDTincrease) * (1 + Paras.gamma_H_amp ./ (1 + exp(-Paras.d_steep * (double(Data.ModelScreeningTime(Y-1)) - Paras.d_change)))) * Paras.gamma_H * ones(1, NumberScreening - (Y-1))];
    
    if Data.ModelScreeningTime(1) == ProjStrat.RDTyear
        yearlyeta_H = (1 + ProjStrat.RDTincrease) * (1 + Paras.eta_H_amp ./ (1 + exp(-Paras.d_steep * (double(Data.ModelScreeningTime(1) - 1) - (Paras.d_change+Paras.eta_H_lag))))) * Paras.eta_H * ones(1, NumberScreening);
        yearlygamma_H = (1 + ProjStrat.RDTincrease) * (1 + Paras.gamma_H_amp ./ (1 + exp(-Paras.d_steep * (double(Data.ModelScreeningTime(1) - 1) - Paras.d_change)))) * Paras.gamma_H * ones(1, NumberScreening);
    end
    
    % calculate yearly uVector maintaining a constant death rate
    death_rate = (1-Paras.u) * Paras.gamma_H;
    uVector = 1- death_rate ./ yearlygamma_H;
    
    % Vector control
    p_targetdie(1 : NumberScreening+1) = 0;
    TargetFreq(1 : NumberScreening+1) = 1;
    if Paras.VCstart ~= 0
        Y1 = find(Data.ModelScreeningTime >= Paras.VCstart, 1);
        p_targetdie(Y1:end) = Paras.TargetDie;
        TargetFreq(Y1:end) = Paras.TargetFreq;
    end
    if ProjStrat.NewVCyear ~= 0
        Y2 = find(Data.ModelScreeningTime >= ProjStrat.NewVCyear, 1);
        p_targetdie(Y2:end) = ProjStrat.NewTargetDie;
        TargetFreq(Y2:end) = ProjStrat.NewTargetFreq;
    end
    
    
    Pop = [S_H S_A E_H E_A I1_H I_A I2_H 0 R_H 0 P_V S_V G_V E1_V E2_V E3_V I_V];
    T = 0;
    
    Data.ModelScreeningTime = [Data.ModelScreeningTime Data.ModelScreeningTime(end)+1];
    [ActiveM1, ActiveM2, PassiveM1, PassiveM2, DeathsM, PersonYrsM1, PersonYrsM2, NewInfM] = deal(zeros(1, length(Data.Years)));
    [Active1, Active2, Passive1, Passive2, Deaths, PersonYrs1, PersonYrs2, NewInf] = deal(zeros(1, NumberScreening));
    
    Y = length(Data.Years);
    for s = 1 : NumberScreening
        % No transmissons after EOT 
        if s > 1 && Paras.alpha ~= 0 && Data.ModelScreeningTime(s) == round(Data.ModelScreeningTime(s)) && sum(NewInf(floor(Data.ModelScreeningTime)==floor(Data.ModelScreeningTime(s-1)))) < 1 * Paras.PopGrowth .^ double(Data.PopSizeYear - floor(Data.ModelScreeningTime(s-1)))            
            Paras.alpha = 0;
        end
        
        % Expected active screening detections (S1/S2) each time interval
        TruePos = TurnOut(:,s) * ScreeningSensitivity(s);
        FalsePos = TurnOut(:,s) * (1 - ScreeningSpecificity(s));
                
        if ScreeningSpecificity(s) ~= 1 && ...
               ((Paras.year_spec_100pct == 0 && ((E_H(end,:) + I1_H(end,:) + I2_H(end,:)) * (TruePos - FalsePos) * 10000 / Data.N_H < 0.5 && Data.ModelScreeningTime(s) >= 2018)) || ...
                (Paras.year_spec_100pct ~= 0 && Data.ModelScreeningTime(s) >= Paras.year_spec_100pct))
            ScreeningSpecificity(s:end) = 1;
            FalsePos = TurnOut(:,s) * (1 - ScreeningSpecificity(s));
        end
         
        Active1(s) = I1_H(end,:) * TruePos + S_H(end,:) * FalsePos;
        Active2(s) = I2_H(end,:) * TruePos; %assume false positives are detected as stage 1 
                
        % Dynamics
        DandT = TurnOut(:,s)' * ScreeningSensitivity(s) * Paras.Compliance; % proportional change in different HUMAN group
        ICs = [S_H(end,:) Pop(end,5) E_H(end,:).*(1-DandT) Pop(end,10) I1_H(end,:).*(1-DandT) Pop(end,15) I2_H(end,:).*(1-DandT) Pop(end,20) R_H(end,:)+(E_H(end,:)+I1_H(end,:)+I2_H(end,:)).*DandT Pop(end,25:32)];
        
        parameter = Paras;
        parameter.f = f';
        parameter.meff = meff;
        parameter.mu_H = [Paras.mu_H*ones(1,4) Paras.mu_A]';
        parameter.sigma_H = [Paras.sigma_H*ones(1,4) Paras.sigma_A]';
        parameter.phi_H = [Paras.phi_H*ones(1,4) Paras.phi_A]';
        parameter.omega_H = [Paras.omega_H*ones(1,4) Paras.omega_A]';
        parameter.K_V = K_V;
        
        parameter.gamma_H = [yearlygamma_H(s)*ones(1,4) Paras.gamma_A]';
        parameter.eta_H = [yearlyeta_H(s)*ones(1,4) Paras.eta_A]';
        parameter.p_targetdie = p_targetdie(s);
        parameter.TargetFreq = TargetFreq(s);
            
        if (Data.ModelScreeningTime(s) < Paras.VCstart && Paras.VCstart < Data.ModelScreeningTime(s+1)) || (Data.ModelScreeningTime(s) < ProjStrat.NewVCyear && ProjStrat.NewVCyear < Data.ModelScreeningTime(s+1))
            Tbreak = 365 * min(abs(Paras.VCstart - double(Data.ModelScreeningTime(s))), abs(ProjStrat.NewVCyear - double(Data.ModelScreeningTime(s))));
            [t1, pop1] = ode45(@diffHATmodel, [sum(Data.ModelScreeningFreq(1 : s-1)) sum(Data.ModelScreeningFreq(1 : s-1)) + Tbreak], ICs, [], parameter);
            
            %update VC
            parameter.p_targetdie = p_targetdie(s+1);
            parameter.TargetFreq = TargetFreq(s+1);
            [t2, pop2] = ode45(@diffHATmodel, [sum(Data.ModelScreeningFreq(1 : s-1)) + Tbreak sum(Data.ModelScreeningFreq(1 : s))], pop1(end,:), [], parameter);
            
            pop = [pop1; pop2(2:end,:)];
            t = [t1; t2(2:end,:)];
        else
            [t, pop] = ode45(@diffHATmodel, [sum(Data.ModelScreeningFreq(1 : s-1)) sum(Data.ModelScreeningFreq(1 : s))], ICs, [], parameter);
        end
        Pop = [Pop; pop];
        T = [T; t];
       
        S_H = pop(:, 1:4);
        E_H = pop(:, 6:9);
        I1_H = pop(:, 11:14);
        I2_H = pop(:, 16:19);
        R_H = pop(:, 21:24);
        I_V = pop(:, 32);
        dNH = S_H + E_H + I1_H + I2_H + R_H;
        dNH(dNH == 0) = 1;

                
        % Person years infected in S1 and S2
        % screening per year
        PersonYrs1(s) = trapz(t, sum(I1_H,2)) / 365;
        PersonYrs2(s) = trapz(t, sum(I2_H,2)) / 365;
        
        % Passive detections (S1/S2) each time interval
        Passive1(s) = yearlyeta_H(s) * PersonYrs1(s) * 365;
        Passive2(s) = uVector(s) * yearlygamma_H(s) * PersonYrs2(s) * 365;

        % Deaths
        Deaths(s) = (1 - uVector(s)) * yearlygamma_H(s) * PersonYrs2(s) * 365;
        
        % New infections (influx into I_1H)
        FOI = Paras.alpha * meff * bsxfun(@times, f(1:4), S_H./dNH) .* I_V;
        NewInf(s) = sum(trapz(t, FOI));

        
    end


    
% Output
    % All timepoints
    Classes = array2table([double(Data.ModelScreeningTime(1))+T/365 Pop],... 
              'VariableNames', {'Time', 'S_H1', 'S_H2', 'S_H3', 'S_H4', 'S_A', 'E_H1', 'E_H2', 'E_H3', 'E_H4', 'E_A',...
                                'I1_H1', 'I1_H2', 'I1_H3', 'I1_H4', 'I1_A', 'I2_H1', 'I2_H2', 'I2_H3', 'I2_H4', 'I2_A',...
                                'R_H1', 'R_H2', 'R_H3', 'R_H4', 'R_A',... 
                                'P_V', 'S_V', 'G_V', 'E1_V', 'E2_V', 'E3_V', 'I_V'});
    %Y 
    for y = 1 : Y
        s = find(floor(double(Data.ModelScreeningTime)) == Data.Years(y));
        ActiveM1(y) = sum(Active1(s));
        ActiveM2(y) = sum(Active2(s));
        PassiveM1(y) = sum(Passive1(s));
        PassiveM2(y) = sum(Passive2(s));
        DeathsM(y) = sum(Deaths(s));
        PersonYrsM1(y) = sum(PersonYrs1(s));
        PersonYrsM2(y) = sum(PersonYrs2(s));
        NewInfM(y) = sum(NewInf(s));
    end
    Aggregate = table(Data.Years', ActiveM1', ActiveM2', PassiveM1', PassiveM2', DeathsM', PersonYrsM1', PersonYrsM2', NewInfM', ...
                'VariableNames', {'Year', 'ActiveM1', 'ActiveM2', 'PassiveM1', 'PassiveM2', 'DeathsM', 'PersonYrsM1', 'PersonYrsM2', 'NewInfM'});
    
   
% Main ODE code
function dPop = diffHATmodel(t, pop, parameter)

%Compute vector reduction function
if parameter.p_targetdie==0
    f_T = 0;
else
    f_T = parameter.p_targetdie * (1 - 1/(1+exp(-25/365*(mod(t,365/parameter.TargetFreq)-0.35*365))));
end

%Get populations from inputs
S_H = pop(1:5); E_H = pop(6:10); I1_H = pop(11:15); I2_H = pop(16:20); R_H = pop(21:25);
P_V = pop(26); S_V = pop(27); G_V = pop(28); E1_V=pop(29); E2_V=pop(30); E3_V=pop(31); I_V=pop(32);

N_H = S_H + E_H + I1_H + I2_H + R_H;
N_V = S_V + G_V + E1_V + E2_V + E3_V + I_V;

dNH = N_H;
dNH(N_H==0) = 1;

%Human infection dynamics 
dS_H = parameter.mu_H .* N_H + parameter.omega_H .* R_H - I_V * parameter.alpha * parameter.meff .* parameter.f .* S_H ./ dNH - parameter.mu_H .* S_H;
dE_H = I_V * parameter.alpha * parameter.meff .* parameter.f .* S_H ./ dNH - (parameter.sigma_H + parameter.mu_H) .* E_H;
dI1_H = parameter.sigma_H .* E_H - (parameter.eta_H + parameter.phi_H + parameter.mu_H) .* I1_H;
dI2_H = parameter.phi_H .* I1_H -  (parameter.gamma_H + parameter.mu_H) .* I2_H;
dR_H =  parameter.eta_H .* I1_H + parameter.gamma_H .* I2_H - (parameter.omega_H + parameter.mu_H) .* R_H;

%Tsetse Infection dynamics
%Pupa
dP_V = parameter.B_V * N_V - (parameter.xi_V + P_V/parameter.K_V) * P_V;
%Teneral
dS_V = parameter.xi_V * parameter.p_survive * P_V - parameter.alpha * S_V - parameter.mu_V * S_V;
%Non-teneral
dG_V = parameter.alpha * (1 - f_T) * (1 - sum(parameter.f .* (I1_H + I2_H) ./ dNH) * parameter.p_V) * S_V - parameter.alpha * ((1 - f_T) * parameter.epsilon * sum(parameter.f .* (I1_H + I2_H) ./ dNH) * parameter.p_V + f_T) * G_V - parameter.mu_V * G_V;
%Exposed
dE1_V = parameter.alpha * (1 - f_T) * sum(parameter.f .* (I1_H + I2_H) ./ dNH) * parameter.p_V * (S_V + parameter.epsilon * G_V) - 3 * parameter.sigma_V * E1_V - (parameter.mu_V + parameter.alpha * f_T) * E1_V;
dE2_V = 3 * parameter.sigma_V * E1_V - (3 * parameter.sigma_V + parameter.mu_V + parameter.alpha * f_T) * E2_V;
dE3_V = 3 * parameter.sigma_V * E2_V - (3 * parameter.sigma_V + parameter.mu_V + parameter.alpha * f_T) * E3_V;
%Infected
dI_V= 3 * parameter.sigma_V * E3_V - (parameter.mu_V + parameter.alpha * f_T) * I_V;


dPop = [dS_H; dE_H; dI1_H; dI2_H; dR_H; dP_V; dS_V; dG_V; dE1_V; dE2_V; dE3_V; dI_V];
