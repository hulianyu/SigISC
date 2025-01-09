function [pi,consensus_pval] = consensus_chi(X,m,M)
pi = (X(:,m)~=1) + 1;
X_drop = X;
X_drop(:,m) = [];
pvals = [];
for m_drop=1:M-1
    x = X_drop(:,m_drop);
    if length(unique(x)) ~=1
        [oc,ec,chi] = chi_squared_test_table_two_variable(x,pi);
        % Calculate expected frequencies
        N = sum(oc(:));  % Total sample size
        % Check if any expected frequency is less than 5
        if any(ec(:) < 5)
            % Apply Yates correction if any expected frequency is less than 5
            a = oc(1, 1);
            b = oc(1, 2);
            c = oc(2, 1);
            d = oc(2, 2);
            chi = (N * (abs(a * d - b * c) - N / 2)^2) / ((a + b) * (c + d) * (a + c) * (b + d));
        end
        %% each chi-squared pval
        pvals = [pvals simplifiedChi2cdf(chi,1)];
%     else
%         pvals = [pvals 1];
    end
end
consensus_pval = combine_pvalues(pvals, 'stouffer').pvalue;
end