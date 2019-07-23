function [CI] = get95CI(data,dim)
%function to detect 95% CI for normal data
%created 10-06-2016, modified 04-09-2018
sem_data=sem(data,dim);
tval_data=tinv(0.975,size(data,dim)-1);
CI=tval_data*sem_data;

end