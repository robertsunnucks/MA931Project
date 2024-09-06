%Uses precalculated ODE solutions and linear interpolation to find the
%percentage reduction in tsetse population for given parameters
function Output = FastGetReduction(p_target_die, Time, PreCalc)
    t_index = find(PreCalc.Times >= Time, 1);
    ptd_index = find(PreCalc.ptd >= p_target_die, 1) - 1;
    if ptd_index < 1
        ptd_index = 1;
    end
    coeff = (p_target_die - PreCalc.ptd(ptd_index)) / (PreCalc.ptd(ptd_index+1) - PreCalc.ptd(ptd_index));
    Output = PreCalc.percent(t_index, ptd_index) + coeff * (PreCalc.percent(t_index, ptd_index + 1) - PreCalc.percent(t_index, ptd_index));
end