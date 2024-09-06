reduction_vals_1 = linspace(0,100,1000);
density = betapdf(0.01*reduction_vals_1,12,3);

figure(1)
clf(1)

plot(reduction_vals_1, density)
xlabel('Percentage reduction after 1 year')
ylabel('Probability Density')


TC = 365;
VCred = 85;

reduction_vals = linspace(1,95,40);
frequencies = [1, 2, 4, 365];
ptds = zeros(size(reduction_vals));

figure(2)
clf(2)
hold on
for i = 1:length(frequencies)
    f_val = frequencies(i);
    for j = 1:length(reduction_vals)
        L = length(reduction_vals);
        disp(string(j + (i-1)*L) + ' out of ' + string(L * length(frequencies)))
        VCred = reduction_vals(j);
        ptds(j) = GetTargetDie(VCred, f_val, TC, 5000);
    end
    plot(reduction_vals, ptds,  'DisplayName', string(f_val) + ' Target Deployments per year')
end
legend()
xlabel('Percentage reduction after 1 year')
ylabel('p_{target die}')
hold off



p_td_vals = ones(size(reduction_vals_1));
for i = 1:length(reduction_vals_1)
    p_td_vals(i) = GetTargetDie(reduction_vals_1(i), 2, 365, 10000);
    disp(string(i) + ' out of ' + string(length(reduction_vals_1)))
end

disp('Problem percentages are:')
invalid_indices = find(p_td_vals < -0.5);
disp(reduction_vals_1(invalid_indices))

figure(3)
clf(3)
hold on
valid_indices = find(p_td_vals > -0.5);

area_under = trapz(p_td_vals(valid_indices), density(valid_indices));
density = density ./ area_under;
plot(p_td_vals(valid_indices), density(valid_indices),'DisplayName', 'Transformed Prior')

xlabel('p_{target die}')
ylabel('Probability Density')
legend()