
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                               %
%   This code runs forcasting of all strategies considered in the Warwick HAT model                             %
%                                                                                                               %
%   Inputs:                                                                                                     %
%       Data - structure containing location-specific historical data                                           %
%       Paras - structure containing location-specific parameters (fixed, fitted and intervention parameters)   %
%       ICs - structure containing initial conditions for a single set of parameters                            %
%       Strategy - table containing parameters of all strategies                                                %
%       samples - number denoting sample size from ODE                                                          %
%                                                                                                               %
%   Outputs:                                                                                                    %
%       Outputs - structure containing yearly aggregated outputs                                                %
%                                                                                                               %
%   Functions required: GetEndemicEq & ODEHATmodel & betabinornd                                                %
%                                                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Outputs = Projection(Data, Paras, ICs, Strategy, samples)
    NumStrat = size(Strategy, 1) - 1;
    Years = Data.Years(end) + 1 : Strategy{'Strat1', 'SIMyear'}; % all strategis must have the same SIMyear
    MtoAbsScaling = Paras.PopGrowth .^ double(Years - Data.PopSizeYear);
    Pop = round(Data.N_H * MtoAbsScaling);
    NumYear = length(Years);
    
    [Active1, Active2, Passive1, Passive2, Deaths, PersonYrs1, PersonYrs2, NewInf] = deal(zeros(1, NumYear, NumStrat));
    [YEPHP, YEOT] = deal(zeros(1, NumStrat));
    [SampledActive1, SampledActive2, SampledPassive1, SampledPassive2, SampledDeaths] = deal(zeros(samples, NumYear, NumStrat));
    SampledYEPHP = zeros(samples, NumStrat);

    % Get meff
    [meff, ~] = GetEndemicEq(Data.N_H, Paras); 

    % Get ICs
    ICs = {[ICs.S_H1, ICs.S_H2, ICs.S_H3, ICs.S_H4], [ICs.E_H1, ICs.E_H2, ICs.E_H3, ICs.E_H4], [ICs.I1_H1, ICs.I1_H2, ICs.I1_H3, ICs.I1_H4],...
           [ICs.I2_H1, ICs.I2_H2, ICs.I2_H3, ICs.I2_H4], [ICs.R_H1, ICs.R_H2, ICs.R_H3, ICs.R_H4],...
            ICs.S_A, ICs.E_A, ICs.I1_A, ICs.P_V, ICs.S_V, ICs.G_V, ICs.E1_V, ICs.E2_V, ICs.E3_V, ICs.I_V};
        
    % Update simulation years for ODE
    Data.Years = Years;

    
    for s = 1 : NumStrat
        S = ['Strat', num2str(s)];
        
        % Strategy Parameters
        ProjStrat = table2struct(Strategy(S,:));
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
        
        
        % Run projections
        [Classes, Aggregate] = ODEHATmodel(meff, ICs, Data, Paras, ProjStrat);

        
        % Sampling from projected dynamics
        ActiveS = betabinornd(repmat(Data.PeopleScreened, samples, 1), repmat((Aggregate.ActiveM1' + Aggregate.ActiveM2')./ScaledPeopleScreened, samples, 1), Paras.disp_act);
        ActiveS1 = binornd(ActiveS, repmat(Aggregate.ActiveM1' ./ (Aggregate.ActiveM1' + Aggregate.ActiveM2'), samples, 1));
        PassiveS = betabinornd(repmat(Pop, samples, 1), repmat((Aggregate.PassiveM1' + Aggregate.PassiveM2') / Data.N_H, samples, 1), Paras.disp_pass);
        PassiveS1 = binornd(PassiveS, repmat(Aggregate.PassiveM1' ./ (Aggregate.PassiveM1' + Aggregate.PassiveM2'), samples, 1));
        DeathsS = binornd(repmat(Pop, samples, 1), repmat(Aggregate.DeathsM' / Data.N_H, samples, 1));
    
        SampledActive1(:, :, s) = ActiveS1;
        SampledActive2(:, :, s) = ActiveS - ActiveS1;
        SampledPassive1(:, :, s) = PassiveS1;
        SampledPassive2(:, :, s) = PassiveS - PassiveS1;
        SampledDeaths(:, :, s) = DeathsS;

        % All dynamics
        Active1(:,:,s) = Aggregate.ActiveM1' .* MtoAbsScaling;
        Active2(:,:,s) = Aggregate.ActiveM2' .* MtoAbsScaling;
        Passive1(:,:,s) = Aggregate.PassiveM1' .* MtoAbsScaling;
        Passive2(:,:,s) = Aggregate.PassiveM2' .* MtoAbsScaling;
        Deaths(:,:,s) = Aggregate.DeathsM' .* MtoAbsScaling;
        PersonYrs1(:,:,s) = Aggregate.PersonYrsM1' .* MtoAbsScaling;
        PersonYrs2(:,:,s) = Aggregate.PersonYrsM2' .* MtoAbsScaling;
        NewInf(:,:,s) = Aggregate.NewInfM' .* MtoAbsScaling;
        
        % Elimination years
        YEPHP(s) = max([find(sum(Aggregate{:, {'ActiveM1','ActiveM2','PassiveM1','PassiveM2'}}, 2) * 10000 >= Data.N_H, 1, 'last') 0]) + Years(1);
        YEOT(s) = max([find(Aggregate.NewInfM' .* MtoAbsScaling >= 1.0, 1, 'last') 0]) + Years(1); % no transmission threshold = 1
        SampledCases = SampledActive1(:,:,s) + SampledActive2(:,:,s) + SampledPassive1(:,:,s) + SampledPassive2(:,:,s);
        SampledYEPHP(:, s) = table2array(rowfun(@(SampledCases)(max([find(SampledCases * 10000 >= Pop, 1, 'last') 0])), table(SampledCases))) + double(Years(1)) * ones(samples, 1);
    end
    

    % Output
    Outputs = struct('Years', Years, ...
                     'Active1', Active1, 'Active2', Active2, 'Passive1', Passive1, 'Passive2', Passive2, ...
                     'Deaths', Deaths, 'PersonYrs1', PersonYrs1, 'PersonYrs2', PersonYrs2, 'NewInf', NewInf, ...
                     'YEPHP', YEPHP, 'YEOT', YEOT, 'SampledYEPHP', SampledYEPHP,...
                     'SampledActive1', SampledActive1, 'SampledActive2', SampledActive2, ...
                     'SampledPassive1', SampledPassive1, 'SampledPassive2', SampledPassive2, ...
                     'SampledDeaths', SampledDeaths);

    
    

