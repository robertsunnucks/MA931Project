% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Pre-Calculate ODE solutions
% 
% PreCalcTimes = linspace(0,365,366);
% PreCalc_PTD_RANGE = cat(2,linspace(0,0.3,10000),linspace(0.3001,1,1000),linspace(1.01,50,1000));
% 
% PreCalcPERCENT_REDUCTION = zeros(366,12000);
% 
% for i = 2:366
%     disp(i)
%     for j = 1:12000
%         PreCalcPERCENT_REDUCTION(i,j) = GetVCReductionPct(PreCalc_PTD_RANGE(j),PreCalcTimes(i),2);
%     end
% end
% 
% 
% 
% PreCalc.Times = PreCalcTimes;
% PreCalc.ptd = PreCalc_PTD_RANGE;
% PreCalc.percent = PreCalcPERCENT_REDUCTION;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save("PreCalc.mat","PreCalc")
load("PreCalc.mat")




%Set up prior distributions
reduction_vals = cat(2, zeros(1), linspace(4,100,1000));
density = betapdf(0.01*reduction_vals,12,3);

[p_td_vals, ptd_pdf] = Get_PDF_ptd(reduction_vals, density,PreCalc);

gmax = 1000;
g0_vals = linspace(0,gmax,1000);
g0_pdf = unifpdf(g0_vals,0,gmax);

OD_vals = cat(2,linspace(0,1,100), linspace(1.01,10,50));
OD_pdf = unifpdf(OD_vals,0,10);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%We take this many steps in our MCMC method
total_steps = 2000;


%Calculate synthetic data
NumTraps = [1000 1000 1000 1000 1000];
Trap_Times = [0 90 180 270 360];
TRAP_TOTALS = zeros(size(NumTraps));

artificial_OD = 0.5;
artificial_g0 = 50;
artificial_ptd = 0.075;
for i = 1:length(Trap_Times)
    artificial_reduction = FastGetReduction(artificial_ptd,Trap_Times(i),PreCalc);
    r = NumTraps(i) / artificial_OD;
    p = 1 / (1 + artificial_OD * artificial_g0 * (100-artificial_reduction)/100);
    TRAP_TOTALS(i) = nbinrnd(r,p);
end




%Update our priors
[ptd_MC,g0_MC,OD_MC] = UpdatePrior(total_steps, TRAP_TOTALS, NumTraps, Trap_Times, p_td_vals, g0_vals,OD_vals,ptd_pdf,g0_pdf,OD_pdf,PreCalc);


%Plot the posterior distribution of p target die
figure(1)
clf(1)
hold on
h=histogram(ptd_MC(200:end),200,"Normalization","pdf")
plot(p_td_vals,ptd_pdf, 'DisplayName','Prior p target die')
plot([artificial_ptd, artificial_ptd], [0, 1.1*max(h.Values)],'-x', 'DisplayName','True value')
legend()
hold off

%Plot the posterior distribution of the scale g0
figure(2)
clf(2)
hold on
h=histogram(g0_MC(200:end),200,"Normalization","pdf")
plot(g0_vals,g0_pdf, 'DisplayName','Prior g0 (scale)')
plot([artificial_g0, artificial_g0], [0, 1.1*max(h.Values)],'-x', 'DisplayName','True value')
hold off
legend()

%Plot the posterior distribution of overdispersion
figure(3)
clf(3)
hold on
h=histogram(OD_MC(200:end),200,"Normalization","pdf")
plot(OD_vals,unifpdf(OD_vals,min(OD_vals),max(OD_vals)), 'DisplayName','Prior overdispersion')
plot([artificial_OD, artificial_OD], [0, 1.1*max(h.Values)],'-x', 'DisplayName','True value')
hold off
legend()