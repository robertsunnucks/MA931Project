%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Pre-Calculate ODE solutions
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load("PreCalc.mat")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%Distributions of our variables for trap data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reduction_vals = cat(2,linspace(0,98,100),linspace(98.1,99.99,25));

density = betapdf(0.01*reduction_vals,12,3);

[p_td_vals, ptd_pdf] = Get_PDF_ptd(reduction_vals, density,PreCalc);

gmax = 100;
g0_vals = linspace(0,gmax,1000);
g0_pdf = unifpdf(g0_vals,0,gmax);


OD_vals = cat(2,linspace(0,1,100), linspace(1.01,10,50));
OD_pdf = unifpdf(OD_vals,0,10);

Prior_ptd_PDF = ptd_pdf;
Prior_ptd_VALUES = p_td_vals;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





Cloc = 'DRC';
Ploc = 1;
Zloc = 29;
ParaStr = 'DRC101';
RunProjection = 10; % options: 1-1000 (number of realizations for Projection)
RunSamples = 0; % options: any positive integer (number of samples per realisation, usually 10) or 0 (ODE results only)
RunPlot = 1;
% 
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
NumStrat = 4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%% Data processing
N_H = PopSize(location);
Data = struct('Years', YEAR, 'N_H', N_H, 'PopSizeYear', PopSizeYear,...
              'MeanPeopleScreened', MeanPeopleScreened(location), 'MaxPeopleScreened', MaxPeopleScreened(location),...
              'LocStr', LocStr, 'IDStr', IDStr); % input data for main functions 'DirStr', Dir
Data.FileStr = ['_', M, '_']; % for MCMC, Fitted dynamics and ReactInfo
FileStr = ['_', M, '_React0_'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%% Run simulation
NumPosterior = RunProjection; % number of realizations used in Porjection
samples = RunSamples;

Years = 2020 : 2050;

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

%Initialise variables
N_V = zeros(NumPosterior, NumStrat, 2000);
EFF = zeros(NumPosterior, NumStrat, 2000);
Times = zeros(NumPosterior, NumStrat, 2000);
lengths = zeros(NumPosterior, NumStrat);
YEOTs = zeros(NumPosterior, NumStrat);
DALYS = zeros(NumPosterior, NumStrat);
Screening_Active_Inf = DALYS;
Screening_Passive_Inf = DALYS;
Screening_NotInf = DALYS;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Output = tinytargets(Time_years, p_target_die, Target_frequency)
    Time_days = 365 * (Time_years);
    t = mod(Time_days,365./Target_frequency);
    Output = p_target_die .* (1 - 1./(1+exp(-25/365*(t-0.35*365))));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

p_td_range = p_td_vals;
Total_Costs = zeros(length(p_td_range), NumStrat);

disp('Starting Simulations')

%Loop through all values of p target die

ACTUAL_PTD = 0.05;

for p = 1 : NumPosterior
    disp(p)
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

        % Set up the strategy, using strategies 1,2,3,4
        ProjStrat = table2struct(Strategy(S,:));
        ProjStrat.NewVCyear = 2021;
        ProjStrat.SIMyear = 2050;
        if s <= 2
            Paras.TargetFreq = 2;
            Paras.TargetDie = ACTUAL_PTD;
            Paras.VCstart = 1;
        else
            Paras.TargetDie = 0;
            Paras.TargetFreq = 0;
            Paras.VCstart = 1;
        end
        if mod(s,2) < 0.5
            ProjStrat.NewTargetFreq = 0;
            ProjStrat.NewTargetDie = 0;
        else
            ProjStrat.NewTargetFreq = 2;
            ProjStrat.NewTargetDie = ACTUAL_PTD;
        end



        %We use mean AS
        ProjStrat.NewASnum = 'mean';
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

        %Run the ODE
        [Classes, Aggregate] = ODEHATmodel(meff, myICS, Data, Paras, ProjStrat);

        YEOTs(p,s) = max([find(Aggregate.NewInfM' .* MtoAbsScaling >= 1.0, 1, 'last') 0]) + Years(1); % no transmission threshold = 1

        %Calculate DALYS

        %Need correct DALY weightings
        Inf1_weight = 0.1;
        Inf2_weight = 0.5;

        %Expected years dead is life expetancy as an adult, minus the
        %average age of death from HAT
        Years_Dead = 70 - 28;
        %3 percent yearly discounting
        Discounting = transpose(power(0.97,linspace(0,(2050-2020),(2050-2020+1))));
        DALYS(p,s) = Inf1_weight*sum(Aggregate.PersonYrsM1 .* Discounting) + Inf2_weight*sum(Aggregate.PersonYrsM2 .* Discounting) + Years_Dead * sum(Aggregate.DeathsM .* Discounting);

        %Calculate Active Screening
        Screening_Active_Inf(p,s) = sum(Aggregate.ActiveM1 + Aggregate.ActiveM2);
        Screening_Passive_Inf(p,s) = sum(Aggregate.PassiveM1 + Aggregate.PassiveM2);
        Screening_NotInf(p,s) = sum(Data.ModelPeopleScreened) - Screening_Active_Inf(p,s);

        %Store infection dynamics
        lengths(p,s) = length(Classes.Time);
        N_V(p,s,1:length(Classes.Time)) = Classes.P_V + Classes.S_V + Classes.G_V + Classes.E1_V + Classes.E2_V + Classes.E3_V + Classes.I_V;
        Times(p,s, 1:length(Classes.Time)) = Classes.Time;

        %Calculate tiny target effectiveness over time
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
end

myfontsize = 14;

figure(1)
clf(1)
hold on
time_end = find(squeeze(Times(1,1,:)) >= 2023, 1) - 1;
plot(squeeze(Times(1,1,1:time_end)),squeeze(EFF(1,1,1:time_end)),'k')
xlabel('Time (year)','FontSize',myfontsize)
ylabel('Tiny target effectiveness','FontSize',myfontsize)
ylim([0, ACTUAL_PTD * 1.1])
hold off


Strategy_Colour_Palette = ["#b3cde3","#8c96c6","#8856a7","#810f7c"];

figure(2)
clf(2)
hold on
for i = 1:2
    time_end = find(squeeze(Times(1,i,:)) >= 2023, 1) - 1;
    plot(squeeze(Times(1,i,1:time_end)),squeeze(N_V(1,i,1:time_end)), 'DisplayName','Strategy' + string(i),'LineWidth',2)
end
xlabel('Time (year)','FontSize',myfontsize)
ylabel('Vector Population','FontSize',myfontsize)
legend('FontSize',myfontsize)
hold off