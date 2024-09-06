
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                         %
%   This code plots multiple strategies projection results of the Warwick HAT model                                       %
%                                                                                                                         %
%   Inputs:                                                                                                               %
%       Cloc - string containing ALPHA-3 country codes, defined in Wrapper                                                %
%       Ploc - number denoting the provine index, defined in Wrapper                                                      %
%       Zloc - number denoting the health zone index, defined in Wrapper                                                  %
%       ParaStr - string containing ALPHA-3 country code and 3-digits related to parameter settings, defined in Wrapper   %
%       StratStr - string containing 4 strategy IDs                                                                       %
%                                                                                                                         %
%   File required: Data.mat & Projection_*.mat & Elimination_*.mat                                                        %
%                                                                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotFigure1(Ploc, Zloc)
Cloc = 'DRC';
ParaStr = 'DRC101';
StratStr = [1 9 3 11]; % strategy IDs

NumStrat = length(StratStr);
Y = 14; % number of years of projections

MyBlue = [67 147 195]/255;
MyRed = [178 24 42]/255;
MyDarkPurple = [118,42,131]/255;
MyDarkGreen = [27,120,55]/255;
MyGreen = [127,191,123]/255;
MyPink = [231 41 138]/255;
MyOrange = [254 153 41]/255;
MyGrey = [0.6 0.6 0.6];
MyDarkColor2 = {MyBlue, MyRed, MyDarkPurple, MyDarkGreen};

p = dir(['../Data/', Cloc, '/P', num2str(Ploc), '_*']);
load(['../Data/', Cloc, '/', p.name, '/Data.mat']);
z.name = strcat('Z', num2str(Zloc), '_', CCLOC{Zloc});
names = {z.name, p.name, Cloc};
location = Zloc;
Dir = ['../Result/', Cloc, '/', p.name, '/', z.name, '/'];
PostDir = ['../Posteriors/', p.name, '/', z.name, '/'];

LocStr = LOCSTR{location};

load([PostDir, 'Posterior_M4_', LocStr, '_ID', ParaStr, '.mat']);
if size(Posterior, 1) == 1
    Message = ['No inference performed in the selected location (', LocStr, ')']
    return;
end

MeanScreened = MeanPeopleScreened(location);
MaxScreened = MaxPeopleScreened(location);

ActiveQ = [];
PassiveQ = [];
NewInfQ = [];
strat = [];
for s = StratStr
    load([Dir, 'Projection_M4_Strat', num2str(s), '_React0_', LocStr, '_ID', ParaStr, '.mat']);
    if isempty(SampPostID) == 0
        ActiveQ = [ActiveQ quantile(SampledActive1(:, 1 : Y) + SampledActive2(:,  1 : Y), [0.025 0.25 0.5 0.5 0.75 0.975])];
        PassiveQ = [PassiveQ quantile(SampledPassive1(:, 1 : Y) + SampledPassive2(:,  1 : Y), [0.025 0.25 0.5 0.5 0.75 0.975])];
    else
        ActiveQ = [ActiveQ quantile(Active1(:, 1 : Y) + Active2(:,  1 : Y), [0.025 0.25 0.5 0.5 0.75 0.975])];
        PassiveQ = [PassiveQ quantile(Passive1(:, 1 : Y) + Passive2(:,  1 : Y), [0.025 0.25 0.5 0.5 0.75 0.975])];
    end
NewInfQ = [NewInfQ quantile(NewInf(:, 1 : Y), [0.025 0.25 0.5 0.5 0.75 0.975])];
strat = [strat repmat({num2str(s)}, 1, Y)];
end
year = repmat({'2017', '2018', '2019', '2020', '2021', '2022', '2023', '2024', '2025', '2026', '2027', '2028', '2029', '2030'}, 1, length(StratStr));
pos = sort([0.17:Y 0.17+0.22:Y 0.17+0.44:Y 0.17+0.66:Y]);

% Figure setting
figure('Name', ['Projections_', LocStr], 'NumberTitle', 'off')    
 
% ===== Active =====
subplot(4, 1, 1)
ActiveBox = boxplot(ActiveQ, {year, strat}, 'symbol', '', 'colors', repmat([MyBlue; MyRed; MyDarkPurple; MyDarkGreen], Y, 1), 'position', pos);
set(ActiveBox, {'linew'}, {0.5})

uw = findobj(ActiveBox, 'tag', 'Upper Whisker');   % get handle to "Upper Whisker" line
uav = findobj(ActiveBox, 'tag', 'Upper Adjacent Value');   %get handle to "Upper Adjacent Value" line
lw = findobj(ActiveBox, 'tag', 'Lower Whisker');   % get handle to "Lower Whisker" line
lav = findobj(ActiveBox, 'tag', 'Lower Adjacent Value');   %get handle to "Lower Adjacent Value" line
m = findobj(ActiveBox, 'tag', 'Median');   %get handle to "Median" line
out = findobj(ActiveBox, 'tag', 'Outliers');   %get handle to outliers
b = findobj(ActiveBox, 'tag', 'Box');   %get handle to box
   
for i = 1 : length(StratStr) * Y
    %Ensure whiskers are at 97.5% and 2.5% give solid whiskers
    q = Y * mod(i - 1, NumStrat) + ceil(i / NumStrat);
    uw(i).YData(:) = [ActiveQ(5,q) ActiveQ(6,q)];
    uw(i).LineStyle = '-';
    uav(i).YData(:) = [ActiveQ(6,q) ActiveQ(6,q)];
    uav(i).XData(:) = mean(uav(i).XData(:)) + 0.05 * [-1 1];
    lw(i).YData(:) = [ActiveQ(1,q) ActiveQ(2,q)];
    lw(i).LineStyle = '-';
    lav(i).YData(:) = [ActiveQ(1,q) ActiveQ(1,q)];

    %Fill box
    patch(get(b(i), 'XData'), get(b(i), 'YData'), [MyDarkColor2{mod(i - 1, NumStrat) + 1}], 'FaceAlpha', 0.3);
    m(i).LineWidth = 1;
end

plt = gca;
set(plt,'children',flipud(get(gca,'children')))
xticks([3 4 8 9 13 14])
xticklabels({'          2020','','          2025','','          2030',''})
plt.XRuler.TickLabelGapOffset = -5;
plt.XAxis.TickDirection = 'out';
plt.XAxis.TickLength = [0.018 1];
plt.YAxis.TickLength = [0.012 1];
plt.XAxis.MinorTick = 'on';
plt.XAxis.MinorTickValues = [0 1 2 5 6 7 10 11 12];
box off;

xlim([0, Y])
ylim([0, 1.15 * max(ActiveQ(6,:))])

hold on
plot([3 3], ylim, 'k', 'LineStyle', '--', 'LineWidth', 1)

ylabel({'Active', 'cases'})
ax2 = axes('Position', get(plt, 'Position'), 'Color', 'None', 'LineWidth', 0.5, 'XAxisLocation', 'top', 'XTick', [], 'YAxisLocation', 'right', 'YTick', []);
hold off


% ===== Passive =====
subplot(4, 1, 2)
PassiveBox = boxplot(PassiveQ, {year, strat}, 'symbol', '', 'colors', repmat([MyBlue; MyRed; MyDarkPurple; MyDarkGreen], Y, 1), 'position', pos);
set(PassiveBox, {'linew'}, {0.5})

uw = findobj(PassiveBox, 'tag', 'Upper Whisker');   % get handle to "Upper Whisker" line
uav = findobj(PassiveBox, 'tag', 'Upper Adjacent Value');   %get handle to "Upper Adjacent Value" line
lw = findobj(PassiveBox, 'tag', 'Lower Whisker');   % get handle to "Lower Whisker" line
lav = findobj(PassiveBox, 'tag', 'Lower Adjacent Value');   %get handle to "Lower Adjacent Value" line
m = findobj(PassiveBox, 'tag', 'Median');   %get handle to "Median" line
out = findobj(PassiveBox, 'tag', 'Outliers');   %get handle to outliers
b = findobj(PassiveBox, 'tag', 'Box');   %get handle to box
   
for i = 1 : length(StratStr) * Y
    %Ensure whiskers are at 97.5% and 2.5% give solid whiskers
    q = Y * mod(i - 1, NumStrat) + ceil(i / NumStrat);
    uw(i).YData(:) = [PassiveQ(5,q) PassiveQ(6,q)];
    uw(i).LineStyle = '-';
    uav(i).YData(:) = [PassiveQ(6,q) PassiveQ(6,q)];
    uav(i).XData(:) = mean(uav(i).XData(:)) + 0.05 * [-1 1];
    lw(i).YData(:) = [PassiveQ(1,q) PassiveQ(2,q)];
    lw(i).LineStyle = '-';
    lav(i).YData(:) = [PassiveQ(1,q) PassiveQ(1,q)];

    %Fill box
    patch(get(b(i), 'XData'), get(b(i), 'YData'), [MyDarkColor2{mod(i - 1, NumStrat) + 1}], 'FaceAlpha', 0.3);
    m(i).LineWidth = 1;
end

plt = gca;
set(plt,'children',flipud(get(gca,'children')))
xticks([3 4 8 9 13 14])
xticklabels({'          2020','','          2025','','          2030',''})
plt.XRuler.TickLabelGapOffset = -5;
plt.XAxis.TickDirection = 'out';
plt.XAxis.TickLength = [0.018 1];
plt.YAxis.TickLength = [0.012 1];
plt.XAxis.MinorTick = 'on';
plt.XAxis.MinorTickValues = [0 1 2 5 6 7 10 11 12];
box off

xlim([0, Y])
ylim([0, 1.15 * max(PassiveQ(6,:))])

hold on
plot([3 3], ylim, 'k', 'LineStyle', '--', 'LineWidth', 1)

ylabel({'Passive', 'cases'})
ax2 = axes('Position', get(plt, 'Position'), 'Color', 'None', 'LineWidth', 0.5, 'XAxisLocation', 'top', 'XTick', [], 'YAxisLocation', 'right', 'YTick', []);
hold off


% ===== New infections =====
subplot(4, 1, 3)
NewInfBox = boxplot(NewInfQ, {year, strat}, 'symbol', '', 'colors', repmat([MyBlue; MyRed; MyDarkPurple; MyDarkGreen], Y, 1), 'position', pos);
set(NewInfBox, {'linew'}, {0.5})

uw = findobj(NewInfBox, 'tag', 'Upper Whisker');   % get handle to "Upper Whisker" line
uav = findobj(NewInfBox, 'tag', 'Upper Adjacent Value');   %get handle to "Upper Adjacent Value" line
lw = findobj(NewInfBox, 'tag', 'Lower Whisker');   % get handle to "Lower Whisker" line
lav = findobj(NewInfBox, 'tag', 'Lower Adjacent Value');   %get handle to "Lower Adjacent Value" line
m = findobj(NewInfBox, 'tag', 'Median');   %get handle to "Median" line
out = findobj(NewInfBox, 'tag', 'Outliers');   %get handle to outliers
b = findobj(NewInfBox, 'tag', 'Box');   %get handle to box
   
for i = 1 : length(StratStr) * Y
    %Ensure whiskers are at 97.5% and 2.5% give solid whiskers
    q = Y * mod(i - 1, NumStrat) + ceil(i / NumStrat);
    uw(i).YData(:) = [NewInfQ(5,q) NewInfQ(6,q)];
    uw(i).LineStyle = '-';
    uav(i).YData(:) = [NewInfQ(6,q) NewInfQ(6,q)];
    uav(i).XData(:) = mean(uav(i).XData(:)) + 0.05 * [-1 1];
    lw(i).YData(:) = [NewInfQ(1,q) NewInfQ(2,q)];
    lw(i).LineStyle = '-';
    lav(i).YData(:) = [NewInfQ(1,q) NewInfQ(1,q)];

    %Fill box
    patch(get(b(i), 'XData'), get(b(i), 'YData'), [MyDarkColor2{mod(i - 1, NumStrat) + 1}], 'FaceAlpha', 0.3);
    m(i).LineWidth = 1;
end

plt = gca;
set(plt,'children',flipud(get(gca,'children')))
xticks([3 4 8 9 13 14])
xticklabels({'          2020','','          2025','','          2030',''})
plt.XRuler.TickLabelGapOffset = -5;
plt.XAxis.TickDirection = 'out';
plt.XAxis.TickLength = [0.018 1];
plt.YAxis.TickLength = [0.012 1];
plt.XAxis.MinorTick = 'on';
plt.XAxis.MinorTickValues = [0 1 2 5 6 7 10 11 12];
box off

xlim([0, Y])
ylim([0, 1.15 * max(NewInfQ(6,:))])

hold on
plot([3 3], ylim, 'k', 'LineStyle', '--', 'LineWidth', 1)

ylabel({'New', 'infections'})
ax2 = axes('Position', get(plt, 'Position'), 'Color', 'None', 'LineWidth', 0.5, 'XAxisLocation', 'top', 'XTick', [], 'YAxisLocation', 'right', 'YTick', []);
hold off


% ===== PEOT =====
load([Dir, 'Elimination_M4_React0_', LocStr, '_ID', ParaStr, '.mat']);
EOT = YEOT(:, StratStr);
Prob = [];
for y = 2017 : 2030
    Prob(y-2016, :) = mean(EOT <= y);
end

subplot(4, 1, 4)
scatter(0.5:13.5, Prob(:, 1), 30, 'filled', 'MarkerEdgeColor', MyBlue, 'MarkerFaceAlpha', 0.3, 'MarkerFaceColor', MyBlue)
hold on
plot(0.5:13.5, Prob(:, 1), '-o', 'LineWidth', 1, 'Color', MyBlue, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'w')
scatter(0.5:13.5, Prob(:, 2), 30, 'filled', 'MarkerEdgeColor', MyRed, 'MarkerFaceAlpha', 0.3, 'MarkerFaceColor', MyRed)
plot(0.5:13.5, Prob(:, 2), '-o', 'LineWidth', 1, 'Color', MyRed, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'w')
scatter(0.5:13.5, Prob(:, 3), 30, 'filled', 'MarkerEdgeColor', MyDarkPurple, 'MarkerFaceAlpha', 0.3, 'MarkerFaceColor', MyDarkPurple)
plot(0.5:13.5, Prob(:, 3), '-o', 'LineWidth', 1, 'Color', MyDarkPurple, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'w')
scatter(0.5:13.5, Prob(:, 4), 30, 'filled', 'MarkerEdgeColor', MyDarkGreen, 'MarkerFaceAlpha', 0.3, 'MarkerFaceColor' ,MyDarkGreen)
plot(0.5:13.5, Prob(:, 4), '-o', 'LineWidth', 1, 'Color', MyDarkGreen, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'w')
    
plt = gca;
set(gca, 'children', flipud(get(gca, 'children')))
xticks([3 4 8 9 13 14])
xticklabels({'          2020','','          2025','','          2030',''})
plt.XAxis.TickLength = [0.018 1];
plt.YAxis.TickLength = [0.012 1];
plt.XAxis.MinorTick = 'on';
plt.XAxis.MinorTickValues = [0 1 2 5 6 7 10 11 12];
box on

xlim([0, Y])
ylim([0, 1])
    
plot([3 3], ylim, 'k', 'LineStyle', '--', 'LineWidth', 1)
    
ylabel('PEOT')
hold off
    

% ===== Add legend manually =====
annotation('textbox', [0.14, 0.084, 0.2, 0.0], 'String', ['N_H = ', num2str(round(PopSize(location)))], 'FontSize', 12, 'HorizontalAlignment', 'Left', 'LineStyle', 'none')
annotation('line', [0.15 0.2], [0.02 0.02], 'LineWidth', 1 , 'Color', 'k', 'LineStyle', '--')
annotation('textbox', [0.2, 0.044, 0.2, 0.0], 'String', 'VC starts', 'FontSize', 12, 'HorizontalAlignment', 'Left', 'LineStyle', 'none')

load(['Paras_', ParaStr, '.mat']);
AS = Strategy{['Strat', num2str(StratStr(1))], 'NewASnum'}{:};
VC = Strategy{['Strat', num2str(StratStr(1))], 'NewVC'};
str = [AS, 'AS+', num2str(VC), '%VC'];
annotation('line', [0.37 0.42], [0.06 0.06], 'LineWidth', 1 , 'Color', MyBlue)
annotation('textbox', [0.42, 0.084, 0.2, 0.0], 'String', str, 'FontSize', 12, 'HorizontalAlignment', 'Left', 'LineStyle', 'none')

AS = Strategy{['Strat', num2str(StratStr(2))], 'NewASnum'}{:};
VC = Strategy{['Strat', num2str(StratStr(2))], 'NewVC'};
str = [AS, 'AS+', num2str(VC), '%VC'];
annotation('line', [0.65 0.7], [0.06 0.06], 'LineWidth', 1 , 'Color', MyRed)
annotation('textbox', [0.7, 0.084, 0.2, 0.0], 'String', str, 'FontSize', 12, 'HorizontalAlignment', 'Left', 'LineStyle', 'none')

AS = Strategy{['Strat', num2str(StratStr(3))], 'NewASnum'}{:};
VC = Strategy{['Strat', num2str(StratStr(3))], 'NewVC'};
str = [AS, 'AS+', num2str(VC), '%VC'];
annotation('line', [0.37 0.42], [0.02 0.02], 'LineWidth', 1 , 'Color', MyDarkPurple)
annotation('textbox', [0.42, 0.044, 0.2, 0.0], 'String', str, 'FontSize', 12, 'HorizontalAlignment', 'Left', 'LineStyle', 'none')

AS = Strategy{['Strat', num2str(StratStr(4))], 'NewASnum'}{:};
VC = Strategy{['Strat', num2str(StratStr(4))], 'NewVC'};
str = [AS, 'AS+', num2str(VC), '%VC'];
annotation('line', [0.65 0.7], [0.02 0.02], 'LineWidth', 1 , 'Color', MyDarkGreen)
annotation('textbox', [0.7, 0.044, 0.2, 0.0], 'String', str, 'FontSize', 12, 'HorizontalAlignment', 'Left', 'LineStyle', 'none')
