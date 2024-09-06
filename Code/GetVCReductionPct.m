%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                               %
%   This code computes VC reduction percentage from TargetDie in posterior                                      %
%                                                                                                               %
%   Inputs:                                                                                                     %
%       Paras - structure containing location-specific parameters (fixed, fitted and intervention parameters)   %                                        %
%       p_targetDie -  known value for Target effectiveness 
%       TargetCycle - 
%       TargetFreq - 
%                                        
%                                                                                                               %
%   Outputs:                                                                                                    %
%       ReductionPct - tsetse reduction after time TrapCycle for p_targetdie value                                     %
%                                                                                                               %                                                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function ReductionPct = GetVCReductionPct(p_targetdie, TrapCycle,TargetFreq)
    
    ODEoptions = odeset('NonNegative', 1:3, 'RelTol', 1e-8, 'AbsTol', 1e-8);
                

    Paras = struct('mu_V',1/33,'xi_V',1/27,'p_survive',0.75,'B_V',0.0505,'alpha',1/3);

                
    % rescale initial condition and K to N_V=100
    Paras.k_V = 100 * Paras.mu_V / (Paras.xi_V^2 * Paras.p_survive * (Paras.p_survive * Paras.B_V / Paras.mu_V - 1));
    IC = 100 * [Paras.mu_V / (Paras.xi_V * Paras.p_survive) (Paras.mu_V / (Paras.mu_V + Paras.alpha)) (Paras.alpha/(Paras.mu_V + Paras.alpha))];
    

    Paras.p_targetdie = p_targetdie;
    Paras.TargetFreq = TargetFreq;
            
    [t, pop] = ode45(@TsetseDyn, [0 TrapCycle], IC, ODEoptions, Paras);
    
    %Compute N_V as the sum of S_V and G_V tsetse
    N_V = sum(pop(:,2:3),2);

    ReductionPct = 100 - sum(pop(end,2:3));
    ReductionPct = round(ReductionPct, 2);

    %Plot dynamics
    %plot([-TrapCycle/2; t],[N_V(1); N_V],'LineWidth',2)
    %axis([ -Inf Inf 0 105])
    %ylabel('Percentage of original population')
    %xlabel('Time (days)')
end

function dPop = TsetseDyn(t, pop, parameter)
    if parameter.p_targetdie==0
        f_T = 0;
    else
        f_T = parameter.p_targetdie * (1 - 1/(1+exp(-25/365*(mod(t,365/parameter.TargetFreq)-0.35*365))));
    end
    
    dPop = zeros(3,1);
    P_V = pop(1);
    S_V = pop(2);
    G_V = pop(3);
    N_V = S_V + G_V;

    % Pupa
    dPop(1) = parameter.B_V * N_V - (parameter.xi_V + P_V/parameter.k_V) * P_V;

    % Teneral (feed twice as fast as other adults)
    dPop(2) = parameter.xi_V * parameter.p_survive * P_V - parameter.alpha * S_V - parameter.mu_V * S_V;

    % Non-teneral
    dPop(3)= parameter.alpha * (1 - f_T) * S_V - parameter.alpha * f_T * G_V - parameter.mu_V * G_V;
    
end