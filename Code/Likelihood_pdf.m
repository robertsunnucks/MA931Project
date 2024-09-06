%Function to calculate the likelihood function's pdf
function Output = Likelihood_pdf(x,t, N, varied_param,varied_param_range,known_param_vals,PreCalc)
    reductionpct = zeros(size(x));
    gi = zeros(size(x));

    %We store the value of the 2 parameters that are remaining fixed
    count = 1;
    param = ones(1,3);
    for j = 1:3
        if j ~= varied_param
            param(j) = known_param_vals(count);
            count = count + 1;
        else
            param(j) = 0;
        end
    end
    %If we are keeping ptd fixed, we can calculate the percentage reduction
    if varied_param > 1.5
        for j = 1:length(t)
            reductionpct(j) = (100-FastGetReduction(param(1), t(j), PreCalc))/100;
        end
    end
    %initialise log likelihood array
    loglike = zeros(size(varied_param_range));
    %Loop through the parameter we are varying
    for i = 1:length(varied_param_range)
        param(varied_param) = varied_param_range(i);
        %Calculate g_i
        if varied_param == 1
            for j = 1:length(t)
                gi(j) = param(2) * (100-FastGetReduction(param(1), t(j), PreCalc)) / 100;
            end
        else
            for j = 1:length(t)
                gi(j) = param(2) * reductionpct(j);
            end
        end
        %Calculate the log likelihood at this parameter value
        loglike(i) = 0;
        %Loop over the data points, x
        for j = 1:length(x)
            if param(3) > 0
                %We use Poisson-Gamma if we have nonzero overdispersion
                term1 = gammaln(x(j) + N(j) / param(3)) - gammaln(x(j) + 1) - gammaln(N(j) / param(3));
                term2 = x(j) * log(param(3) * gi(j)) - (N(j) / param(3) + x(j)) * log(1+param(3) * gi(j));
                x_pdf = term1 + term2;
            else
                %For 0 overdispersion, we instead use a Poisson
                mu = N(j) * gi(j);
                x_pdf = x(j) * log(mu) - mu - gammaln(x(j) + 1);
            end
            loglike(i) = loglike(i) + x_pdf;
        end
    end
    %Take exponentional of the log likelihood to get back to the likelihood
    Output = exp(loglike - max(loglike));
end