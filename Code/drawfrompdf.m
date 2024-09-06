%function that will, given the values and the pdf of a distribution over 
%these values, draw a random value from this distribution
function Output = drawfrompdf(vals,pdf)
    %Find cumulative distribution
    cumulative = cumtrapz(vals,pdf);
    %Normalise it (in case of error from user)
    cumulative = cumulative ./ trapz(vals,pdf);
    %Draw a uniform random number
    r = rand;
    %Use cumulative and random num to pick a value
    i = find(cumulative >= r, 1);
    Output = vals(i-1) + (r - cumulative(i-1)) / (cumulative(i) - cumulative(i-1)) * (vals(i) - vals(i-1));
    if isnan(Output)
        Output = vals(i-1);
    end
end