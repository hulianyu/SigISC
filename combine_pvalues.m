function res = combine_pvalues(pvalues, method, weights)
% COMBINE_PVALUES Combine p-values from independent tests.
%
%   RES = COMBINE_PVALUES(PVALUES) combines the p-values using Fisher's method
%   by default.
%
%   RES = COMBINE_PVALUES(PVALUES, METHOD) specifies the method to use. 
%   Valid methods are:
%       'fisher'             - Fisher's method
%       'pearson'            - Pearson's method
%       'mudholkar_george'   - Mudholkar and George's method
%       'tippett'            - Tippett's method
%       'stouffer'           - Stouffer's Z-score method
%
%   RES = COMBINE_PVALUES(PVALUES, METHOD, WEIGHTS) specifies weights for
%   Stouffer's method. WEIGHTS must be the same length as PVALUES.
%
%   The output RES is a struct with fields:
%       RES.statistic  - The combined statistic
%       RES.pvalue     - The combined p-value

    % Set default method to 'fisher' if not provided
    if nargin < 2 || isempty(method)
        method = 'fisher';
    else
        method = lower(method);
    end

    % Initialize the result structure
    res = struct('statistic', NaN, 'pvalue', NaN);

    % Convert pvalues to a column vector for consistency
    pvalues = pvalues(:);

    % Handle empty pvalues
    if isempty(pvalues)
        return;
    end

    % Number of p-values
    n = length(pvalues);

    switch method
        case 'fisher'
            % Fisher's method
%             statistic = -2 * sum(log(pvalues));    %
            safe_p = max(pvalues, realmin);          % Ensure p-values are above the smallest positive value
            sum_logp = sum(log(safe_p));
            statistic = -2 * sum_logp;
            %%
            df = 2 * n;
%             pval = 1 - chi2cdf(statistic, df);     %
            pval = chi2cdf(statistic, df, 'upper');  
        case 'pearson'
            % Pearson's method
            statistic = 2 * sum(log(1 - pvalues)); %
            df = 2 * n;
            pval = chi2cdf(-statistic, df);
            
        case 'mudholkar_george'
            % Mudholkar and George's method
            normalizing_factor = sqrt(3 / n) / pi;
            statistic = -sum(log(pvalues)) + sum(log(1 - pvalues));
            nu = 5 * n + 4;
            approx_factor = sqrt(nu / (nu - 2));
            t_stat = statistic * normalizing_factor * approx_factor;
            pval = tcdf(t_stat, nu, 'upper');  % Survival function
            
        case 'tippett'
            % Tippett's method
            statistic = min(pvalues);
            pval = betacdf(statistic, 1, n);

        case 'stouffer'
            % Stouffer's Z-score method
            if nargin < 3 || isempty(weights)
                weights = ones(n, 1);
            else
                weights = weights(:);
                if length(weights) ~= n
                    error('pvalues and weights must be of the same size.');
                end
            end
            % Calculate Z-scores
%             Zi = norminv(1 - pvalues); % double precision issue -> Inf
            Zi = -norminv(pvalues);
            statistic = (weights' * Zi) / norm(weights);
            pval = normcdf(statistic,'upper');

        otherwise
            error(['Invalid method ''' method '''. Valid methods are ', ...
                   '''fisher'', ''pearson'', ''mudholkar_george'', ', ...
                   '''tippett'', and ''stouffer''.']);
    end

    % Assign results
    res.statistic = statistic;
    res.pvalue = pval;
end
