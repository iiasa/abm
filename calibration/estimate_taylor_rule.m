function [rho,r_star,xi_pi,xi_gamma,pi_star] = estimate_taylor_rule(r_bar,pi_EA,gamma_EA)
ydata=r_bar;
exo=[pi_EA,gamma_EA];
var=rfvar3(ydata,1,exo(1:length(ydata),:),[],0,0);
alpha=var.By;
gamma_1=var.Bx(1);
gamma_2=var.Bx(2);

rho=alpha;
xi_pi=gamma_1/(1-rho);
xi_gamma=gamma_2/(1-rho);
pi_star=(0.02+1).^(1/4)-1;
r_star=pi_star*(xi_pi-1);
end