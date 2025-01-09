function [minpval,bpi,bm] = find_best_pattern(X)
X(:,1) = []; % objsID
M = size(X,2);
pv = cell(M,1);
pi = cell(M,1);
for m=1:M
    on_off = length(unique(X(:,m)));
    if on_off>1
        [pi_current, pv_current] = consensus_chi(X,m,M); % add Yates correction
        pi{m,1} = [pi{m,1} pi_current];
        pv{m,1} = [pv{m,1} pv_current];
    else
        pi{m,1} = ones(size(X,1),1);
        pv{m,1} = 1;
    end
end
% find the best m
[minpval, bm] = min(cellfun(@min, pv));
% find the best category according to the recorded discat
bq = find(pv{bm}==minpval,1);
% find the best partition
bpi = pi{bm}(:,bq);
end