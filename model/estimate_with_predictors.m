function [alpha,beta,gamma_1,gamma_2,gamma_3,epsilon]=estimate_with_predictors(ydata,exo)
var=rfvar3(ydata,1,[ones(size(ydata,1),1),exo(1:length(ydata),:)],[],0,0);
alpha=var.By;
beta=var.Bx(1);
gamma_1=var.Bx(2);
gamma_2=var.Bx(3);
gamma_3=var.Bx(4);
epsilon=normrnd(0,sqrt(cov(var.u)));
end
