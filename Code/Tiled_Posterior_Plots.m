%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Pre-Calculate ODE solutions - can comment this out if already have the
%%precalculated results stored in the PreCalc file 
%
% PreCalcTimes = linspace(0,365,366);
% PreCalc_PTD_RANGE = linspace(0,0.3,1000);
% 
% PreCalcPERCENT_REDUCTION = zeros(366,1000);
% 
% for i = 2:366
%     disp(i)
%     for j = 1:1000
%         PreCalcPERCENT_REDUCTION(i,j) = GetVCReductionPct(PreCalc_PTD_RANGE(j),PreCalcTimes(i),2);
%     end
% end
% 
% PreCalc.Times = PreCalcTimes;
% PreCalc.ptd = PreCalc_PTD_RANGE;
% PreCalc.percent = PreCalcPERCENT_REDUCTION;
% save("PreCalc.mat","PreCalc")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load("PreCalc.mat")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

myfontsize = 14;

%Initialise prior distributions

reduction_vals = linspace(0,100,1000);
density = betapdf(0.01*reduction_vals,12,3);
[p_td_vals, ptd_pdf] = Get_PDF_ptd(reduction_vals, density,PreCalc);

gmax = 100;
g0_vals = linspace(0,gmax,1000);
g0_pdf = unifpdf(g0_vals,0,gmax);

OD_vals = cat(2,linspace(0,1,100), linspace(1.01,10,50));
OD_pdf = unifpdf(OD_vals,0,10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sequential_colour_palette_1 = ["#b3cde3","#8c96c6","#8856a7","#810f7c"];
sequential_colour_palette_2 = ["#fef0d9","#fdcc8a","#fc8d59","#d7301f"];
categorical_colour_palette = ["#a6cee3","#b2df8a","#33a02c","#1f78b4"];

%Initialise parameters
total_steps = 10000;
Possible_Trap_Times = [0 90 180 270 360];
artificial_OD = 0.1;
artificial_g0 = 1;
artificial_ptd = 0.07;

%Create synthetic data using artificial parameter choices
data_points = 200;
synthetic_data = zeros(data_points,length(Possible_Trap_Times));
synthetic_num_traps = 500;
artificial_reduction = zeros(size(Possible_Trap_Times));
for j = 1:length(Possible_Trap_Times)
    artificial_reduction(j) = FastGetReduction(artificial_ptd,Possible_Trap_Times(j),PreCalc);
    if artificial_OD > 0
        r = synthetic_num_traps / artificial_OD;
        p = 1 / (1 + artificial_OD * artificial_g0 * (100-artificial_reduction(j))/100);
        synthetic_data(:,j) = nbinrnd(r,p,data_points,1);
    else
        mu = artificial_g0 * (100-artificial_reduction(j))/100;
        synthetic_data(:,j) = poissrnd(synthetic_num_traps * mu,data_points,1);
    end
end

%Plot our synthetic data with box and whisker
figure(10)
clf(10)
ax = gca; 
ax.FontSize = myfontsize;
hold on
Time_vals_plot = cat(2,linspace(-90,-1,10),linspace(0,360,100));
reduction_plot = zeros(1,110);
for i = 1:100
    reduction_plot(i+10) = FastGetReduction(artificial_ptd,Time_vals_plot(i+10),PreCalc);
end

reduction_plot = synthetic_num_traps * artificial_g0 * (100-reduction_plot)/100;

plot(Time_vals_plot/90 + 2,reduction_plot,'--','LineWidth',3)
Modified_Plotting_Times = cat(2,-90,Possible_Trap_Times);
data_plot = zeros(data_points,length(Possible_Trap_Times)+1);
data_plot(:,1) = synthetic_data(:,1);
data_plot(:,2) = -1000;
data_plot(:,3:end) = synthetic_data(:,2:end);
boxplot(data_plot,Modified_Plotting_Times)
plot([2,2],[0,1.1*max(synthetic_data,[],"all")],'k')
plot([4,4],[0,1.1*max(synthetic_data,[],"all")],'k')
ylim([0,1.1*max(synthetic_data,[],"all")])
ylabel('Amount caught in traps (with 500 traps)','FontSize',myfontsize)
xlabel('Time (days)','FontSize',myfontsize)
title('Synthetic trapping data, with vector control deployment twice per year','FontSize',myfontsize)
hold off



   
figure(1)
clf(1)
ax = gca; 
ax.FontSize = myfontsize;

Individual_Traps = 1000;
height = 0;
min_hist_val = 1;
max_hist_val = 0;

hold on

disp(1)

%Set up which traps are being used
WhichTraps = [1 0 0 0 0];

total = round(sum(WhichTraps));
NumTraps = ones(1,total) .* Individual_Traps;

%Create trapping times
Trap_Times = zeros(1,total);
count = 0;
for i = 1:length(WhichTraps)
    if WhichTraps(i) > 0
        count = count + 1;
        Trap_Times(count) = Possible_Trap_Times(i);
    end
end


%Create trapping data using artificial params
TRAP_TOTALS = zeros(size(Trap_Times));
for i = 1:length(Trap_Times)
    artificial_reduction = FastGetReduction(artificial_ptd,Trap_Times(i),PreCalc);
    if artificial_OD > 0
        r = NumTraps(i) / artificial_OD;
        p = 1 / (1 + artificial_OD * artificial_g0 * (100-artificial_reduction)/100);
        TRAP_TOTALS(i) = nbinrnd(r,p);
    else
        mu = NumTraps(i) * artificial_g0 * (100-artificial_reduction)/100;
        TRAP_TOTALS(i) = poissrnd(mu);
    end
end

%Update the prior and plot
[ptd_MC,g0_MC,OD_MC] = UpdatePrior(total_steps, TRAP_TOTALS, NumTraps, Trap_Times, p_td_vals, g0_vals,OD_vals,ptd_pdf,g0_pdf,OD_pdf,PreCalc);

h=histogram(ptd_MC(200:end),150,"Normalization","pdf", 'DisplayName','Posterior from only trapping before intervention', 'FaceAlpha',0.8)
height = max(height, max(h.Values));
min_hist_val = min(min_hist_val, prctile(ptd_MC(200:end),2));
max_hist_val = max(max_hist_val, prctile(ptd_MC(200:end),98));
    
plot(p_td_vals,ptd_pdf,'k', 'DisplayName','Prior p target die', 'LineWidth',5)
plot([artificial_ptd, artificial_ptd], [0, 1.1*height],'--x', 'DisplayName','True value', 'LineWidth',1)
hold off
xlim([min_hist_val, max_hist_val])
xlabel('p target die','FontSize',myfontsize)
ylabel('probability density','FontSize',myfontsize)
legend('FontSize',myfontsize)
title('Posterior convergence for only trapping before intervention, with 1000 traps','FontSize',myfontsize)




%Reset the total steps. Time taken to run scales linearly with this but the
%higher it is, the smoother the posterior curve will be
total_steps = 1500;
figure(2)
clf(2)
ax = gca; 
ax.FontSize = myfontsize;
Individual_Traps = synthetic_num_traps;
height = 0;
min_hist_val = 1;
max_hist_val = 0;

hold on

%Loop through different monitoring strategies
for k = 1:4
    disp(k+1)
    %Set up which traps are being used
    WhichTraps = [1 0 0 0 0];
    WhichTraps(k+1) = 1;

    total = round(sum(WhichTraps));
    NumTraps = ones(1,total) .* Individual_Traps;
    
    %Calculate trap times
    Trap_Times = zeros(1,total);
    count = 0;
    for i = 1:length(WhichTraps)
        if WhichTraps(i) > 0
            count = count + 1;
            Trap_Times(count) = Possible_Trap_Times(i);
        end
    end


    %Initialise array to hold all posterior samples
    Posterior_Samples = [];
    %Loop through sample data
    for i = 1:data_points
        %Set up trap data
        TRAP_TOTALS = zeros(1,total);
        count = 0;
        disp(string(i * 100 / data_points) + ' percent')
        for j = 1:length(WhichTraps)
            if WhichTraps(j) > 0
                count = count + 1;
                TRAP_TOTALS(count) = synthetic_data(i,j);
            end
        end
        %Run MCMC posterior calculation and add to our array
        [ptd_MC,g0_MC,OD_MC] = UpdatePrior(total_steps, TRAP_TOTALS, NumTraps, Trap_Times, p_td_vals, g0_vals,OD_vals,ptd_pdf,g0_pdf,OD_pdf,PreCalc);
        Posterior_Samples = cat(2,Posterior_Samples,ptd_MC(200:end));
    end
    
    %Plot the averaged out pdf of the posterior
    [hist_sizes,hist_edges] = histcounts(Posterior_Samples,150);
    Posterior_ptd_VALUES = hist_edges(1:length(hist_edges)-1) + (hist_edges(2:end) - hist_edges(1:length(hist_edges)-1)) ./ 2;
    Posterior_ptd_PDF = hist_sizes ./ trapz(Posterior_ptd_VALUES, hist_sizes);
    plot(Posterior_ptd_VALUES,Posterior_ptd_PDF,'DisplayName','Trapping at 0 and ' + string(Trap_Times(2)) + ' days', 'LineWidth',2,'Color',categorical_colour_palette(k))
    height = max(height, max(Posterior_ptd_PDF));
    min_hist_val = min(min_hist_val, prctile(Posterior_Samples,2));
    max_hist_val = max(max_hist_val, prctile(Posterior_Samples,98));
end

%Plot the prior and true value of p target die
plot(p_td_vals,ptd_pdf,'k', 'DisplayName','Prior p target die', 'LineWidth',3)
plot([artificial_ptd, artificial_ptd], [0, 1.1*height],'--xk', 'DisplayName','True value', 'LineWidth',1)
hold off
xlim([min_hist_val, max_hist_val])
xlabel('p target die','FontSize',myfontsize)
ylabel('probability density','FontSize',myfontsize)
legend('FontSize',myfontsize)
title('Comparing posterior convergence for different trap times, with 500 traps each time','FontSize',myfontsize)

      
figure(3)
clf(3)
ax = gca; 
ax.FontSize = myfontsize;
height = 0;
min_hist_val = 1;
max_hist_val = 0;
hold on

%Loop through different monitoring strategies
for k = 1:4
    disp(5+k)
    %Set up which traps are being used
    WhichTraps = [1 0 0 0 0];
    for i = 1:k
        WhichTraps(6 - i) = 1;
    end

    %Calculate trap times
    total = round(sum(WhichTraps));
    Individual_Traps = round(1200/total);
    NumTraps = ones(1,total) .* Individual_Traps;
    Trap_Times = zeros(1,total);

    count = 0;
    for i = 1:length(WhichTraps)
        if WhichTraps(i) > 0
            count = count + 1;
            Trap_Times(count) = Possible_Trap_Times(i);
        end
    end
 
    %Initialise array to hold all posterior samples
    Posterior_Samples = [];
    %Loop through sample data
    for i = 1:data_points
        %Set up trap data
        TRAP_TOTALS = zeros(1,total);
        count = 0;
        disp(string(i * 100 / data_points) + ' percent')
        for j = 1:total
            t_val = Trap_Times(j);
            artificial_reduction = FastGetReduction(artificial_ptd,t_val,PreCalc);
            if artificial_OD > 0
                r = NumTraps(j) / artificial_OD;
                p = 1 / (1 + artificial_OD * artificial_g0 * (100-artificial_reduction)/100);
                TRAP_TOTALS(j) = nbinrnd(r,p);
            else
                mu = NumTraps(j) * artificial_g0 * (100-artificial_reduction)/100;
                TRAP_TOTALS(j) = poissrnd(mu);
            end
        end
        %Run MCMC posterior calculation and add to our array
        [ptd_MC,g0_MC,OD_MC] = UpdatePrior(total_steps, TRAP_TOTALS, NumTraps, Trap_Times, p_td_vals, g0_vals,OD_vals,ptd_pdf,g0_pdf,OD_pdf,PreCalc);
        Posterior_Samples = cat(2,Posterior_Samples,ptd_MC(200:end));
    end

    %Plot the averaged out pdf of the posterior
    label_text = '';
    for i = 1:(length(Trap_Times) - 1)
        label_text = label_text + string(Trap_Times(i)) + ', ';
    end
    label_text = label_text + 'and '+ string(Trap_Times(end));
    label_text = 'Trapping at ' + label_text + ' days';


    [hist_sizes,hist_edges] = histcounts(Posterior_Samples,150);
    Posterior_ptd_VALUES = hist_edges(1:length(hist_edges)-1) + (hist_edges(2:end) - hist_edges(1:length(hist_edges)-1)) ./ 2;
    Posterior_ptd_PDF = hist_sizes ./ trapz(Posterior_ptd_VALUES, hist_sizes);
    plot(Posterior_ptd_VALUES,Posterior_ptd_PDF,'DisplayName',label_text, 'LineWidth',2,'Color',sequential_colour_palette_2(k))
    height = max(height, max(Posterior_ptd_PDF));
    min_hist_val = min(min_hist_val, prctile(Posterior_Samples,2));
    max_hist_val = max(max_hist_val, prctile(Posterior_Samples,98)); 
end

%Plot the prior and true value of p target die
plot(p_td_vals,ptd_pdf,'k', 'DisplayName','Prior p target die', 'LineWidth',3)
plot([artificial_ptd, artificial_ptd], [0, 1.1*height],'--xk', 'DisplayName','True value', 'LineWidth',1)
hold off
xlim([min_hist_val, max_hist_val])
xlabel('p target die','FontSize',myfontsize)
ylabel('probability density','FontSize',myfontsize)
legend('FontSize',myfontsize)
title('Comparing posterior convergence for different frequencies of trapping, with 1200 traps used total','FontSize',myfontsize)






figure(4)
clf(4)
ax = gca; 
ax.FontSize = myfontsize;
height = 0;
min_hist_val = 1;
max_hist_val = 0;
hold on

%Loop through different monitoring strategies
for k = 1:3
    disp(9+k)

    %Set up traps numbers
    Individual_Traps = 6 * power(5,k);
    total = 2;
    NumTraps = ones(1,total) .* Individual_Traps;
    Trap_Times = [0 360];

    
    %Initialise array to hold all posterior samples
    Posterior_Samples = [];
    %Loop through sample data
    for i = 1:data_points
        %set up trap data
        TRAP_TOTALS = zeros(1,total);
        count = 0;
        disp(string(i * 100 / data_points) + ' percent')
        for j = 1:total
            t_val = Trap_Times(j);
            artificial_reduction = FastGetReduction(artificial_ptd,t_val,PreCalc);
            if artificial_OD > 0
                r = NumTraps(j) / artificial_OD;
                p = 1 / (1 + artificial_OD * artificial_g0 * (100-artificial_reduction)/100);
                TRAP_TOTALS(j) = nbinrnd(r,p);
            else
                mu = NumTraps(j) * artificial_g0 * (100-artificial_reduction)/100;
                TRAP_TOTALS(j) = poissrnd(mu);
            end
        end
        %Run MCMC posterior calculation and add to our array
        [ptd_MC,g0_MC,OD_MC] = UpdatePrior(total_steps, TRAP_TOTALS, NumTraps, Trap_Times, p_td_vals, g0_vals,OD_vals,ptd_pdf,g0_pdf,OD_pdf,PreCalc);
        Posterior_Samples = cat(2,Posterior_Samples,ptd_MC(200:end));
    end

    %Plot the averaged out pdf of the posterior
    label_text = 'Trapping with ' + string(Individual_Traps) + ' traps';
    [hist_sizes,hist_edges] = histcounts(Posterior_Samples,150);
    Posterior_ptd_VALUES = hist_edges(1:length(hist_edges)-1) + (hist_edges(2:end) - hist_edges(1:length(hist_edges)-1)) ./ 2;
    Posterior_ptd_PDF = hist_sizes ./ trapz(Posterior_ptd_VALUES, hist_sizes);
    plot(Posterior_ptd_VALUES,Posterior_ptd_PDF,'DisplayName',label_text, 'LineWidth',2,'Color',sequential_colour_palette_1(k))
    height = max(height, max(Posterior_ptd_PDF));
    min_hist_val = min(min_hist_val, prctile(Posterior_Samples,2));
    max_hist_val = max(max_hist_val, prctile(Posterior_Samples,98));
end

%Plot the prior and true value of p target die
plot(p_td_vals,ptd_pdf,'k', 'DisplayName','Prior p target die', 'LineWidth',3)
plot([artificial_ptd, artificial_ptd], [0, 1.1*height],'--xk', 'DisplayName','True value', 'LineWidth',1)
hold off
xlim([min_hist_val, max_hist_val])
xlabel('p target die','FontSize',myfontsize)
ylabel('probability density','FontSize',myfontsize)
legend('FontSize',myfontsize)
title('Comparing posterior convergence for different numbers of traps, being used before and 360 days after intervention','FontSize',myfontsize)