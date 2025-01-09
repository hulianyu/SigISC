function [oc,ec,chi] = chi_squared_test_table_two_variable(Xm,pi)
N = length(Xm);
oc = zeros(3,3);% observed counts
var = 2 - pi;
tmp = Xm + var;
tmp2 = var - Xm;
%% oc
oc(1,1) = sum(tmp==2);
oc(1,2) = sum(tmp2==1);
oc(2,1) = sum(tmp2==-1);
oc(2,2) = sum(tmp==0);
%%
oc(3,:) = sum(oc);
oc(:,3) = sum(oc,2);
ec = (oc(end,1:end-1).*oc(1:end-1,end))./N;
oc = oc(1:end-1,1:end-1);
chi = sum(((oc-ec).^2)./ec,'all');% chi-squared statistic
end