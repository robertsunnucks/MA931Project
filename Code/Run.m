
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                     %
%   This code processes the raw data, reorganises the parameters and runs the simulations of the Warwick HAT model.   %
%                                                                                                                     %
%   Inputs:                                                                                                           %
%       Cloc - a string in the format of ALPHA-3 country codes                                                        %
%       Ploc - an integer, provine index                                                                              %
%       Zloc - an integer, health zone index                                                                          %
%       ParaStr - a string of ALPHA-3 country code and 3-digits related to parameter settings                         %
%       RunProjection - an integer denoting the number of realizations used in Projection                             %
%       RunSamples - an integer denoting the number of samples from ODE                                               %
%                                                                                                                     %
%   Main output files:                                                                                                %
%       Projection - matrices containing key outputs (e.g. cases, new infections, deaths) in forcasting period        %
%       Elimination - matrices containing years and probabilities of EPHP and EOT under different strategies          %
%                                                                                                                     %
%   Files required: Data.mat & Paras*.mat & Posterior*.mat                                                            %
%                                                                                                                     %
%   Note: hosts are (1) low-risk, random participants                                                                 %
%                   (2) high-risk, random participants                                                                %
%                   (3) low-risk, non-participants                                                                    %
%                   (4) high-risk, non-participants                                                                   %
%                   (5) reservoir animals                                                                             %
%                   (6) non-reservoir animals, no dynamics and is ignored                                             %
%                                                                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Run(Cloc, Ploc, Zloc, ParaStr, RunProjection, RunSamples, RunPlot)
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
    NumStrat = size(locStrategy,1) - 1;

    
  
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

    Years = YEAR(end) + 1 : locStrategy{'Strat1', 'SIMyear'};
    NumYear = length(Years); % fixed value for all strategies
            
    [Active1All, Active2All, Passive1All, Passive2All, DeathsAll, PersonYrs1All, PersonYrs2All, NewInfAll] = deal(zeros(NumPosterior * samples, NumYear, NumStrat));
    [SampledActive1All, SampledActive2All, SampledPassive1All, SampledPassive2All, SampledDeathsAll] = deal(zeros(NumPosterior * samples, NumYear, NumStrat));
    [YEPHP, YEOT, SampledYEPHP] = deal(zeros(NumPosterior * samples, NumStrat));
    [PEPHP, PEOT, SampledPEPHP] = deal(zeros(NumYear, NumStrat));
       
    input2 = load([ICsDir, 'ProjectionICs', Data.FileStr, Data.LocStr, IDStr, '.mat']);
    RowID = datasample(1:1000, NumPosterior, 'Replace', false); % randomly select realizations from PostID in ICs
    PostID = input2.ProjectionICs.PostID(RowID);
    SampPostID = reshape(repmat(PostID', samples, 1), [], 1);
    
    parfor p = 1 : NumPosterior
        ['Posterior', num2str(PostID(p))];
        Paraz(p) = Paras;
        Dataz(p) = Data;
        % Replace values of fitted parameters in Paras by the values from Posterior
        for i = 1 : length(fitted_para_names)
            Paraz(p).(fitted_para_names{i}) = input1.Posterior{PostID(p), i};
        end
        ICs(p) = table2struct(input2.ProjectionICs(RowID(p), :));
                
        ProjectionOutputs(p) = Projection(Dataz(p), Paraz(p), ICs(p), locStrategy, samples); % single realization and all strategies
    end
    
    
    
    %%% Output
    for p = 1 : NumPosterior
        Outputs = ProjectionOutputs(p);
        rsamples = max([samples 1]);
        range = (p-1) * rsamples + 1 : p * rsamples;
        Active1All(range,:,:) = repmat(Outputs.Active1, rsamples, 1);
        Active2All(range,:,:) = repmat(Outputs.Active2, rsamples, 1);
        Passive1All(range,:,:) = repmat(Outputs.Passive1, rsamples, 1);
        Passive2All(range,:,:) = repmat(Outputs.Passive2, rsamples, 1);
        DeathsAll(range,:,:) = repmat(Outputs.Deaths, rsamples, 1);
        PersonYrs1All(range,:,:) = repmat(Outputs.PersonYrs1, rsamples, 1);
        PersonYrs2All(range,:,:) = repmat(Outputs.PersonYrs2, rsamples, 1);
        NewInfAll(range,:,:) = repmat(Outputs.NewInf, rsamples, 1);
        YEPHP(range,:) = repmat(Outputs.YEPHP, rsamples, 1);
        YEOT(range,:) = repmat(Outputs.YEOT, rsamples, 1);
                
        range = (p-1) * samples + 1 : p * samples;
        SampledActive1All(range,:,:) = Outputs.SampledActive1;
        SampledActive2All(range,:,:) = Outputs.SampledActive2;
        SampledPassive1All(range,:,:) = Outputs.SampledPassive1;
        SampledPassive2All(range,:,:) = Outputs.SampledPassive2;
        SampledDeathsAll(range,:,:) = Outputs.SampledDeaths;
        SampledYEPHP(range,:) = Outputs.SampledYEPHP;
                
    end
    
    % Calculate elimination probabilities
    for Y = Years
        PEPHP(Y-Years(1)+1, :) = mean(YEPHP <= Y);
        PEOT(Y-Years(1)+1, :) = mean(YEOT <= Y);
        SampledPEPHP(Y-Years(1)+1, :) = mean(SampledYEPHP <= Y);
    end
    
    ElimByYears = Years';
    save([ResultDir, 'Elimination', FileStr, Data.LocStr, IDStr, '.mat'], 'PostID', 'SampPostID', 'ElimByYears', 'YEPHP', 'YEOT', 'SampledYEPHP', 'PEPHP', 'PEOT', 'SampledPEPHP');
            
    
    for s = 1 : NumStrat
        S = ['Strat', num2str(s)];
        FileStr = ['_', M, '_', S, '_React0_'];
        
        Active1 = Active1All(:, :, s);
        Active2 = Active2All(:, :, s);
        Passive1 = Passive1All(:, :, s);
        Passive2 = Passive2All(:, :, s);
        Deaths = DeathsAll(:, :, s);
        PersonYrs1 = PersonYrs1All(:, :, s);
        PersonYrs2 = PersonYrs2All(:, :, s);
        NewInf = NewInfAll(:, :, s);

        SampledActive1 = SampledActive1All(:, :, s);
        SampledActive2 = SampledActive2All(:, :, s);
        SampledPassive1 = SampledPassive1All(:, :, s);
        SampledPassive2 = SampledPassive2All(:, :, s);
        SampledDeaths = SampledDeathsAll(:, :, s);
                
        save([ResultDir, 'Projection', FileStr, Data.LocStr, IDStr, '.mat'], 'PostID', 'SampPostID', 'Years', ...
                                                                       'Active1', 'Active2', 'Passive1', 'Passive2', 'Deaths', 'PersonYrs1', 'PersonYrs2', 'NewInf', ...
                                                                       'SampledActive1', 'SampledActive2', 'SampledPassive1', 'SampledPassive2', 'SampledDeaths') 
    end
    
    % Plot figure 1
    if RunPlot == 1
        PlotFigure1(Ploc, Zloc)
    end
end

