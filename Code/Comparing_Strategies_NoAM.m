
Cloc = 'DRC';
Ploc = 1;
Zloc = 29;
ParaStr = 'DRC101';
RunProjection = 1000; % options: 1-1000 (number of realizations for Projection)
RunSamples = 0; % options: any positive integer (number of samples per realisation, usually 10) or 0 (ODE results only)
RunPlot = 1;

%%% Import raw data
if Zloc ~= 0 % health Zone level simulation
    p = dir(['../Data/', Cloc, '/P', num2str(Ploc), '_*']);
    load(['../Data/', Cloc, '/', p.name, '/Data.mat']);
    z.name = strcat('Z', num2str(Zloc), '_', CCLOC{Zloc});
    names = {z.name, p.name, Cloc};
    location = Zloc;
    ResultDir = ['../Result/', Cloc, '/', p.name, '/', z.name, '/'];
    PosteriorDir = ['../Posteriors/', p.name, '/', z.name, '/'];
    ICsDir = ['../ICs/', p.name, '/', z.name, '/'];
elseif Ploc ~= 0 % province level simulation
    load(['../Data/', Cloc, '/Data.mat']);
    p.name = strcat('P', num2str(Ploc), '_', CCLOC{Ploc});
    names = {p.name, Cloc};
    location = Ploc;
    ResultDir = ['../Result/', Cloc, '/', p.name, '/'];
    PosteriorDir = ['../Posteriors/', p.name, '/'];
    ICsDir = ['../ICs/', p.name, '/'];
else % country level simulation
    load(['../Data/Data.mat']);
    names = {Cloc};
    location = find(strcmp(COUNTRY, Cloc));
    ResultDir = ['../Result/', Cloc, '/'];
    PosteriorDir = ['../Posteriors/'];
    ICsDir = ['../ICs/'];
end
mkdir(ResultDir)
LocStr = LOCSTR{location}; % location info string ending with the name of smallest scale
IDStr = ['_ID', ParaStr];
M = 'M4';
FileStr = ['_', M, '_'];
input1 = load([PosteriorDir, 'Posterior', FileStr, LocStr, IDStr, '.mat']);

if size(input1.Posterior, 1) == 1
    Message = ['No inference performed in the selected location (', LocStr, ')']
    save([ResultDir, 'NoInferencePerformed.mat'], 'Message')
    return;
end

%%% Parameters
load(['Paras_', ParaStr, '.mat']);

% FixedParameters for current location...
x = 0;
a = 0;
while x == 0
    a = a + 1;
    x = strcmp(FixedParameters.Location, names{a});
end
locFixedPar = FixedParameters(strcmp(FixedParameters.Location, names{a}), ...
                              ~strcmp(FixedParameters.Properties.VariableNames,'Location'));
                          
% Intervention parameters
x = 0;
a = 0;
while x == 0
    a = a + 1;
    x = contains(InterventionParameters.Location, names{a});
end
Irow = find(contains(InterventionParameters.Location, names{a}));
    
% FittedParameters for current location...
locFittedPar = cell2table(cell(0,width(FittedParameters)),'VariableNames',FittedParameters.Properties.VariableNames);
% work from most-local to least-local area name
for i = 1:length(names)
    % Keep line if 1) Location matches current level, and
    %              2) parameter hasn't been matched at a more local level.
    choose = (strcmp(FittedParameters.Location, names{i}) .* ...        
              ~ismember(FittedParameters.Notation, char(locFittedPar.Notation))) == 1;
    locFittedPar = [locFittedPar; FittedParameters(choose,:)];
end
locFittedPar = locFittedPar(:,~strcmp(locFittedPar.Properties.VariableNames,'Location')); % drop location
locFittedPar = [locFittedPar(~startsWith(locFittedPar.Notation,'active_neg_'),:);...                       % if fitting neg active tests,
                sortrows(locFittedPar(startsWith(locFittedPar.Notation,'active_neg_'),:), 'Notation')];    % put that chunk, sorted, last.
            
% FittedParameters for current model only
FittedAll = locFittedPar(contains(locFittedPar.Model, [string('All'), M]), {'Notation','Initial'});
fitted_para_names = locFittedPar.Notation(contains(locFittedPar.Model, [string('All'), M]) & ~contains(locFittedPar.Distribution, 'Delta'))';

% Parameter reorganisation - all parameters 
Paras = table2struct([cell2table(num2cell(FittedAll.Initial)', 'VariableNames', FittedAll.Notation'), locFixedPar, InterventionParameters(Irow, 2:end)]); % input parameters for main functions

Paras.eta_H_amp = Paras.eta_H_amp * (Paras.d_change ~=  0);     % set eta_H_amp   = 0 if d_change = 0
Paras.gamma_H_amp = Paras.gamma_H_amp * (Paras.d_change ~=  0); % set gamma_H_amp = 0 if d_change = 0
    
% Final year in which the active screening specifity can be lower (eg for MSF screenings)
Paras.Last_year = locFittedPar.Last_year(strcmp(locFittedPar.Notation, 'b_specificity'));
                                  
% localise Strategy
x = 0;
a = 0;
while x == 0
    a = a + 1;
    x = strcmp(Strategy.Location, names{a});
end
locStrategy = Strategy(strcmp(Strategy.Location, names{a}), ...
                       ~strcmp(Strategy.Properties.VariableNames,'Location'));
%NumStrat = size(locStrategy,1) - 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Select number of strategies
NumStrat = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Distribution of tiny target effectiveness
load("ptd_prior.mat")
load("ptd_posterior.mat")


%%% Data processing
N_H = PopSize(location);
Data = struct('Years', YEAR, 'N_H', N_H, 'PopSizeYear', PopSizeYear,...
              'MeanPeopleScreened', MeanPeopleScreened(location), 'MaxPeopleScreened', MaxPeopleScreened(location),...
              'LocStr', LocStr, 'IDStr', IDStr); % input data for main functions 'DirStr', Dir
Data.FileStr = ['_', M, '_']; % for MCMC, Fitted dynamics and ReactInfo
FileStr = ['_', M, '_React0_'];



   
%%% Run simulation
NumPosterior = RunProjection; % number of realizations used in Porjection
samples = RunSamples;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Change Years to start at 2020
%Years = YEAR(end) + 1 : locStrategy{'Strat1', 'SIMyear'};
Years = 2020 : locStrategy{'Strat1', 'SIMyear'};

NumYear = length(Years); % fixed value for all strategies
MtoAbsScaling = Paras.PopGrowth .^ double(Years - Data.PopSizeYear);
        
[Active1All, Active2All, Passive1All, Passive2All, DeathsAll, PersonYrs1All, PersonYrs2All, NewInfAll] = deal(zeros(NumPosterior * samples, NumYear, NumStrat));
[SampledActive1All, SampledActive2All, SampledPassive1All, SampledPassive2All, SampledDeathsAll] = deal(zeros(NumPosterior * samples, NumYear, NumStrat));
[YEPHP, YEOT, SampledYEPHP] = deal(zeros(NumPosterior * samples, NumStrat));
[PEPHP, PEOT, SampledPEPHP] = deal(zeros(NumYear, NumStrat));
   
input2 = load([ICsDir, 'ProjectionICs', Data.FileStr, Data.LocStr, IDStr, '.mat']);
RowID = datasample(1:1000, NumPosterior, 'Replace', false); % randomly select realizations from PostID in ICs
PostID = input2.ProjectionICs.PostID(RowID);
SampPostID = reshape(repmat(PostID', samples, 1), [], 1);


%NumPosterior = 1;
N_V = zeros(NumPosterior, NumStrat, 2000);
EFF = zeros(NumPosterior, NumStrat, 2000);
Times = zeros(NumPosterior, NumStrat, 2000);
lengths = zeros(NumPosterior, NumStrat);
YEOTs = zeros(NumPosterior, NumStrat);
PYs = zeros(NumPosterior, NumStrat);

disp('Starting Simulations')
for p = 1 : NumPosterior
    ['Posterior', num2str(PostID(p))];
    % Replace values of fitted parameters in Paras by the values from Posterior
    for i = 1 : length(fitted_para_names)
        Paras.(fitted_para_names{i}) = input1.Posterior{PostID(p), i};
    end
    ICs = table2struct(input2.ProjectionICs(RowID(p), :));

    myICS = {[ICs.S_H1, ICs.S_H2, ICs.S_H3, ICs.S_H4], [ICs.E_H1, ICs.E_H2, ICs.E_H3, ICs.E_H4], [ICs.I1_H1, ICs.I1_H2, ICs.I1_H3, ICs.I1_H4],...
           [ICs.I2_H1, ICs.I2_H2, ICs.I2_H3, ICs.I2_H4], [ICs.R_H1, ICs.R_H2, ICs.R_H3, ICs.R_H4],...
            ICs.S_A, ICs.E_A, ICs.I1_A, ICs.P_V, ICs.S_V, ICs.G_V, ICs.E1_V, ICs.E2_V, ICs.E3_V, ICs.I_V};

    [meff, ~] = GetEndemicEq(Data.N_H, Paras); 

    Data.Years = Years;

    %Set strategy
    for s = 1:NumStrat
        S = ['Strat', num2str(7)];
            
        % Strategy Parameters
        ProjStrat = table2struct(Strategy(S,:));

        %Edit Strategy
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %HERE WE EDIT THE STRATEGY DEPENDING ON WHAT WE WANT
        ProjStrat.NewVCyear = 2021;
        if s <= 1
            Paras.TargetFreq = 2;
            Paras.TargetDie = 0.08;
            Paras.VCstart = 1;
            ProjStrat.NewTargetFreq = 0;
            ProjStrat.NewTargetDie = 0;
        else
            ProjStrat.NewTargetFreq = 2;
            ProjStrat.NewTargetDie = 0.08;
            Paras.TargetDie = 0;
            Paras.TargetFreq = 0;
            Paras.VCstart = 1;
        end

        switch ProjStrat.NewASnum
            case 'mean'
                Data.PeopleScreened = Data.MeanPeopleScreened * ones(1, length(Data.Years));
            case 'max'
                Data.PeopleScreened = Data.MaxPeopleScreened * ones(1, length(Data.Years));
            otherwise % percentage
                Data.PeopleScreened = round(0.01 * str2num(ProjStrat.NewASnum(1 : end - 1)) * Data.N_H * Paras.PopGrowth^double(Data.Years(1) - 3 - Data.PopSizeYear)) * ones(1, length(Data.Years));
        end
        ScaledPeopleScreened = Data.PeopleScreened ./ MtoAbsScaling;
        ModelScreeningFreq = [];
        ModelScreeningTime = [];
        ModelPeopleScreened = [];
        ExpectedScreeningTimes = ceil(ScaledPeopleScreened / Data.N_H / Paras.k1);
        ExpectedScreeningTimes(ExpectedScreeningTimes == 0) = 1;
        for y = 1 : length(Data.Years)
            freq = 365 / ExpectedScreeningTimes(y);
            time = 1 / ExpectedScreeningTimes(y);
            number = ScaledPeopleScreened(y) / ExpectedScreeningTimes(y);
            
            ModelScreeningFreq = [ModelScreeningFreq freq * ones(1, ExpectedScreeningTimes(y))];
            ModelScreeningTime = [ModelScreeningTime Data.Years(y) + (0 : ExpectedScreeningTimes(y) - 1) * time];
            ModelPeopleScreened = [ModelPeopleScreened number * ones(1, ExpectedScreeningTimes(y))];
        end
    
        Data.ModelScreeningFreq = ModelScreeningFreq;
        Data.ModelScreeningTime = ModelScreeningTime;
        Data.ModelPeopleScreened = ModelPeopleScreened;
        
    
        [Classes, Aggregate] = ODEHATmodel(meff, myICS, Data, Paras, ProjStrat);
        
        YEOTs(p,s) = max([find(Aggregate.NewInfM' .* MtoAbsScaling >= 1.0, 1, 'last') 0]) + Years(1); % no transmission threshold = 1
        PYs(p,s) = sum(Aggregate.PersonYrsM1) + sum(Aggregate.PersonYrsM2);

        lengths(p,s) = length(Classes.Time);
        N_V(p,s,1:length(Classes.Time)) = Classes.P_V + Classes.S_V + Classes.G_V + Classes.E1_V + Classes.E2_V + Classes.E3_V + Classes.I_V;
        Times(p,s, 1:length(Classes.Time)) = Classes.Time;

        
        time_vals = Classes.Time;
        p_targetdie = zeros(size(time_vals));
        TargetFreq = ones(size(time_vals));
        
        if Paras.VCstart ~= 0
            Y1 = find(time_vals >= Paras.VCstart, 1);
            p_targetdie(Y1:end) = Paras.TargetDie;
            TargetFreq(Y1:end) = Paras.TargetFreq;
        end
        if ProjStrat.NewVCyear ~= 0
            Y2 = find(time_vals >= ProjStrat.NewVCyear, 1);
            p_targetdie(Y2:end) = ProjStrat.NewTargetDie;
            TargetFreq(Y2:end) = ProjStrat.NewTargetFreq;
        end

        eff = tinytargets(time_vals, p_targetdie,TargetFreq);
        EFF(p,s, 1:length(Classes.Time)) = eff;

    end
    percent = 100 * p / NumPosterior;
    if mod(percent,5) < 0.01
        disp(string(percent) + '% complete')
    end
end
disp('Simulations Finished')



desired_time_points = linspace(2020,2030,6*10+1);

index_values = zeros(NumPosterior,NumStrat,length(desired_time_points));
New_N_V = zeros(NumPosterior,NumStrat,length(desired_time_points));
New_eff = zeros(NumPosterior,NumStrat,length(desired_time_points));

disp('Finding time points')
for i = 1:length(desired_time_points)
    if mod(i,10) < 1
        disp(string(i) + ' out of ' + string(length(desired_time_points)))
    end
    tpoint = desired_time_points(i);
    for j = 1:NumPosterior
        for k = 1:NumStrat
            ind = find(Times(j,k,:) >= tpoint, 1);
            index_values(j,k,i) = ind;
            New_N_V(j,k,i) = N_V(j,k,ind);
            New_eff(j,k,i) = EFF(j,k,ind);
        end
    end
end




figure(1)
clf(1)
title('(Average) Vector Population')
hold on
disp('Plotting Vector Population')

for i = 1:NumStrat
    avgvals = mean(squeeze(New_N_V(:,i,:)), 1);
    plot(desired_time_points, avgvals, 'DisplayName','Strategy ' + string(i))
    legend()
end
hold off


figure(2)
clf(2)
title('(Average) Trap Effectiveness')
hold on
disp('Plotting Trap Effectiveness')
for i = 1:NumStrat
    avgvals = mean(squeeze(New_eff(:,i,:)), 1);
    plot(desired_time_points, avgvals, 'DisplayName','Strategy ' + string(i))
    legend()
end
hold off



figure(3)
clf(3)
title('Probability of Elimination of Transmission')
hold on
disp('Plotting Prob of EOT')
for i = 1:NumStrat
    Prob = zeros(1,14);
    for y = 2017 : 2030
        Prob(y-2016) = mean(YEOTs(:,i) <= y);
    end
    plot(2017:2030, Prob, 'DisplayName','Strategy ' + string(i))
    legend()
end
hold off


figure(4)
clf(4)
title('Person Years')
hold on
disp('Plotting Person Years')
meanPYs = zeros(1,NumStrat);
for i = 1:NumStrat
    meanPYs(i) = mean(squeeze(PYs(:,i)), 1);
end
bar(meanPYs)
hold off



figure(5)
clf(5)
title('Costs')
disp('Plotting Costs')

DALY_VAL = 3000; %Needs actual value
DALY_Cost = DALY_VAL * meanPYs;

Yearly_Cost = 100000; %Need actual value
VC_Years = [1, 9];
VisitsPerYear = [2, 2];
VC_Cost = Yearly_Cost .* VisitsPerYear .* VC_Years;

Single_Monitor_Cost = 10000; %Need actual value
Monitoring_Per_Year = [0, 0];
Monitoring_Costs = Single_Monitor_Cost .* Monitoring_Per_Year .* VC_Years;

Total_Costs = DALY_Cost + VC_Cost + Monitoring_Costs;
bar(Total_Costs)



function Output = tinytargets(Time_years, p_target_die, Target_frequency)
    Time_days = 365 * (Time_years);
    t = mod(Time_days,365./Target_frequency);
    Output = p_target_die .* (1 - 1./(1+exp(-25/365*(t-0.35*365))));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NOTE
%SHOULD PROBABLY RUN MORE SIMS THAN JUST 1000 - NOT SURE IF THIS IS ENOUGH
%FOR ACCURATE EXPECTED VALUES

%Also:
%Next step is to be able to have strategies that change at a certain time
%point (e.g. after half a year)
%This could prove complicated to run
%Either need to change ODE code, OR, run the model up to the year of change
%and then run again with ICs updated