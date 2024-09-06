%function to update the prior using gibbs sampling MCMC method

%parameters
% 1) p_targetdie
% 2) g_0
% 3) overdispersion
function [ptd_MC,g0_MC,OD_MC] = UpdatePrior(total_steps, TRAP_TOTALS, NumTraps, Trap_Times, p_td_vals, g0_vals,OD_vals,ptd_pdf,g0_pdf,OD_pdf,PreCalc)   
    %Initialise arrays of MCMC
    ptd_MC = zeros(1,total_steps+1);
    g0_MC = zeros(1,total_steps+1);
    OD_MC = zeros(1,total_steps+1);
    
    Posterior_ptd = ptd_pdf;
    Posterior_g0 = g0_pdf;
    Posterior_OD = OD_pdf;
    
    ptd_MC(1) = drawfrompdf(p_td_vals, Posterior_ptd);
    g0_MC(1) = drawfrompdf(g0_vals, Posterior_g0);
    OD_MC(1) = drawfrompdf(OD_vals, Posterior_OD);
    
    %Loop through the markov chain up to total steps
    for i = 1:total_steps
        %Calculate the likelihood and posterior
        Like = Likelihood_pdf(TRAP_TOTALS,Trap_Times,NumTraps,1,p_td_vals,[g0_MC(i), OD_MC(i)],PreCalc);
        if or(max(Like .* ptd_pdf) == 0, sum(isnan(Like)) > 0)
            disp('Likelihood error')
            Like = ones(size(Like));
        end
        Posterior_ptd =  Like .* ptd_pdf;
        Posterior_ptd = Posterior_ptd ./ trapz(p_td_vals,Posterior_ptd);
        %Draw from the new posterior
        ptd_MC(i+1) = drawfrompdf(p_td_vals, Posterior_ptd);
        
        %Calculate the likelihood and posterior
        Like = Likelihood_pdf(TRAP_TOTALS,Trap_Times,NumTraps,2,g0_vals,[ptd_MC(i+1), OD_MC(i)],PreCalc);
        if or(max(Like) < 0.00001, sum(isnan(Like)) > 0)
            disp('Likelihood error')
            Like = ones(size(Like));
        end
        Posterior_g0 = Like .* g0_pdf;
        Posterior_g0 = Posterior_g0 ./ trapz(g0_vals,Posterior_g0);
        %Draw from the new posterior
        g0_MC(i+1) = drawfrompdf(g0_vals, Posterior_g0);
        
        %Calculate the likelihood and posterior
        Like = Likelihood_pdf(TRAP_TOTALS,Trap_Times,NumTraps,3,OD_vals,[ptd_MC(i+1), g0_MC(i+1)],PreCalc);
        if or(max(Like) < 0.00001, sum(isnan(Like)) > 0)
            disp('Likelihood error')
            Like = ones(size(Like));
        end
        Posterior_OD = Like .* OD_pdf;
        Posterior_OD = Posterior_OD ./ trapz(OD_vals,Posterior_OD);
        %Draw from the new posterior
        OD_MC(i+1) = drawfrompdf(OD_vals, Posterior_OD);
    end
end