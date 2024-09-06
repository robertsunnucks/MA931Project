%Takes in a distribution for the percentage reduction and returns out the
%distribution for p target die, using precalculated ODE results
function [p_td_vals, ptd_pdf] = Get_PDF_ptd(reduction_vals, density,PreCalc)
    p_td_vals_test = ones(size(reduction_vals));
    for i = 1:length(reduction_vals)
        p_td_vals_test(i) = FastGetTargetDie(reduction_vals(i), 365, PreCalc);
    end
    
    valid_indices = find(p_td_vals_test > -0.5);
    
    p_td_vals = p_td_vals_test(valid_indices);
    
    area_under = trapz(p_td_vals, density(valid_indices));
    ptd_pdf = density(valid_indices) ./ area_under;
end