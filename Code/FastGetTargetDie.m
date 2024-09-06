%Uses precalculated ODE solutions and linear interpolation to find the
%p target die parameter, given reduction in tsetse population
function Output = FastGetTargetDie(VCreduction, Time, PreCalc)
    t_index = find(PreCalc.Times >= Time, 1);
    reduction_index = find(PreCalc.percent(t_index,:) >= VCreduction, 1) - 1;
    if reduction_index < 1
        reduction_index = 1;
    end
    coeff = (VCreduction - PreCalc.percent(t_index,reduction_index)) / (PreCalc.percent(t_index,reduction_index+1) - PreCalc.percent(t_index,reduction_index));
    Output = PreCalc.ptd(reduction_index) + coeff * (PreCalc.ptd(reduction_index + 1) - PreCalc.ptd(reduction_index));
    if isempty(reduction_index)
        Output = 1/0;
        disp('Error, impossibly high reduction in this time')
    end
end