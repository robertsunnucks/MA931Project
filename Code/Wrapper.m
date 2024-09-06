
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                               %
%   This is the wrapper to run Warwick HAT model on your local machine                          %
%                                                                                               %
%   Inputs:                                                                                     %
%       Cloc - a string in the format of ALPHA-3 country codes                                  %
%       Ploc - an integer, provine index                                                        %
%       Zloc - an integer, health zone index                                                    %
%       ParaStr - a string of ALPHA-3 country code and 3-digits related to parameter settings   %
%       RunProjection - an integer denoting the number of realizations used in Projection       %
%       RunSamples - an integer denoting the number of samples from ODE                         %
%                                                                                               %
%   Note: a Paras_Clocxxx.mat file stored in the same directory is required                     %
%                                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cloc = 'DRC';
Ploc = 1;
Zloc = 29;
ParaStr = 'DRC101';
RunProjection = 10; % options: 1-1000 (number of realizations for Projection)
RunSamples = 0; % options: any positive integer (number of samples per realisation, usually 10) or 0 (ODE results only)

if Ploc + Zloc ~= 0
    Run(Cloc, Ploc, Zloc, ParaStr, RunProjection, RunSamples, 1)
else
    Pmax = 11;
    Zmax = [52 31 69 44 49 67 35 18 34 83 34];
    
    for p = 1 : Pmax
        for z = 1 : Zmax(p)
            Run(Cloc, p, z, ParaStr, RunProjection, RunSamples, 0)
        end
    end
end









