function [nominal_gdp,real_gdp,nominal_gva,real_gva,nominal_household_consumption,real_household_consumption,nominal_government_consumption,real_government_consumption,nominal_capitalformation,real_capitalformation,nominal_fixed_capitalformation,real_fixed_capitalformation,nominal_fixed_capitalformation_dwellings,real_fixed_capitalformation_dwellings,nominal_exports,real_exports,nominal_imports,real_imports,operating_surplus,compensation_employees,wages,taxes_production,nominal_sector_gva,real_sector_gva,euribor,gdp_deflator_growth_ea,real_gdp_ea,E_CB,D_RoW,L_G,D_k,D_i,D_h,E_k,L_i]=abm(G,H_act,H_inact,J,L,tau_INC,tau_FIRM,tau_VAT,tau_SIF,tau_SIW,tau_EXPORT,tau_CF,tau_G,theta_UB,psi,psi_H,theta_DIV,theta,mu,r_G,zeta,zeta_LTV,zeta_b,alpha_bar_i,beta_i,kappa_i,delta_i,w_bar_i,tau_Y_i,tau_K_i,b_CF_g,b_CFH_g,b_HH_g,c_G_g,c_E_g,c_I_g,a_sg,G_i,T,T_prime,T_max,alpha_pi_EA,beta_pi_EA,sigma_pi_EA,alpha_Y_EA,beta_Y_EA,sigma_Y_EA,rho,r_star,xi_pi,xi_gamma,pi_star,alpha_G,beta_G,sigma_G,alpha_E,beta_E,sigma_E,alpha_I,beta_I,sigma_I,P_i,K_i,M_i,S_i,N_i,D_i,L_i,D_h,w_h,K_h,L_G,E_k,E_CB,D_RoW,O_h,sb_inact,sb_other,Y,pi,r_bar,Y_EA,pi_EA,C_G,C_E,Y_I,P_bar,P_bar_g,P_bar_HH,P_bar_CF,Q_d_i,Pi_i,Pi_k,D_k,C)
nominal_gdp=zeros(1,T);
real_gdp=zeros(1,T);
nominal_gva=zeros(1,T);
real_gva=zeros(1,T);
nominal_household_consumption=zeros(1,T);
real_household_consumption=zeros(1,T);
nominal_government_consumption=zeros(1,T);
real_government_consumption=zeros(1,T);
nominal_capitalformation=zeros(1,T);
real_capitalformation=zeros(1,T);
nominal_fixed_capitalformation=zeros(1,T);
real_fixed_capitalformation=zeros(1,T);
nominal_fixed_capitalformation_dwellings=zeros(1,T);
real_fixed_capitalformation_dwellings=zeros(1,T);
nominal_exports=zeros(1,T);
real_exports=zeros(1,T);
nominal_imports=zeros(1,T);
real_imports=zeros(1,T);
operating_surplus=zeros(1,T);
compensation_employees=zeros(1,T);
wages=zeros(1,T);
taxes_production=zeros(1,T);
nominal_sector_gva=zeros(T,G);
real_sector_gva=zeros(T,G);
euribor=zeros(1,T);
gdp_deflator_growth_ea=zeros(1,T);
real_gdp_ea=zeros(1,T);

% This is the main loop of the model
for t=1:T_max
    % Firms? expectations regarding economic growth and inflation are formed 
    % using simple but optimal AR(1) rules (Eqs. (A.6) on p. 2 of the Appendix)
    [alpha_Y,beta_Y,epsilon_Y]=estimate(log(Y(1:T_prime+t-1)));
    Y_e=exp(alpha_Y*log(Y(T_prime+t-1))+beta_Y+epsilon_Y);
    % Equ. (A.7) on p. 2
    gamma_e=Y_e/Y(T_prime+t-1)-1;
    
    % Expectations on inflation (A.10) on p. 3
    [alpha_pi,beta_pi,epsilon_pi]=estimate(pi(1:T_prime+t-1));
    pi_e=exp(alpha_pi*pi(T_prime+t-1)+beta_pi+epsilon_pi)-1;
    
    % Expectations on economic growth and inflation of the euro area (A.70)
    % The epsilon function correlates shocks for exports and imports with
    % euro area shocks
    [epsilon_Y_EA,epsilon_E,epsilon_I]=epsilon(C);
    gamma_EA=exp(alpha_Y_EA*log(Y_EA)+beta_Y_EA+epsilon_Y_EA)/Y_EA-1;
    Y_EA=exp(alpha_Y_EA*log(Y_EA)+beta_Y_EA+epsilon_Y_EA);
    
    epsilon_pi_EA=normrnd(0,sigma_pi_EA);
    pi_EA=exp(alpha_pi_EA*log(1+pi_EA)+beta_pi_EA+epsilon_pi_EA)-1;
    
    % (A.69) on p. 13: this is our version of the Taylor rule
    r_bar=rho*r_bar+(1-rho)*(r_star+pi_star+xi_pi*(pi_EA-pi_star)+xi_gamma*gamma_EA);
    r=r_bar+mu;
    
    % (A.5) on p. 2 - firms use expected growth rate for their supply
    % choice
    Q_s_i=Q_d_i*(1+gamma_e);
    
    % (A.9) on p. 3: Firms set their prices according to expected inflation
    % of costs!
    pi_c_i=(1+tau_SIF).*w_bar_i./alpha_bar_i.*(P_bar_HH./P_i-1)+1./beta_i.*(sum(a_sg(:,G_i).*P_bar_g)./P_i-1)+delta_i./kappa_i.*(P_bar_CF./P_i-1);
    % Cost-push inflation and "built-in inflation" are added, equ. (A.8)
    P_i=P_i.*(1+pi_c_i)*(1+pi_e);
    
    % (A.14) p. 4 desired investment demand is calculated from depreciation
    % rate, supply, and capital productivity. 
    I_d_i=delta_i./kappa_i.*min(Q_s_i,K_i.*kappa_i);
    
    % (A.18) p. 4 desired amount of intermediate inputs
    % Minimum function is there to avoid investment above the productive capacity
    DM_d_i=min(Q_s_i,K_i.*kappa_i)./beta_i;
    
    % (A.22) p. 5 planned amount of employment
    N_d_i=max(1,round(min(Q_s_i,K_i.*kappa_i)./alpha_bar_i));
    
    % (A.28) p. 6 expected profit of firm i
    Pi_e_i=Pi_i*(1+pi_e)*(1+gamma_e);

    % (A.27) p. 6 Expected net cash flow
    DD_e_i=Pi_e_i-theta*L_i-tau_FIRM*max(0,Pi_e_i)-(theta_DIV*(1-tau_FIRM))*max(0,Pi_e_i);
    % (A.29) desired amount of bank loans
    DL_d_i=max(0,-DD_e_i-D_i);
    
    % Expected capital stock of firm i (A.63), using the average price
    % index of investment goods to account capital stock
    K_e_i=P_bar_CF*(1+pi_e)*K_i;
    % expected loans of firm i
    L_e_i=(1-theta)*L_i;
    % Search and matching of loans for firm i
    DL_i=search_and_matching_credit(DL_d_i,K_e_i,L_e_i,E_k,zeta,zeta_LTV);
    
    % Labour market and search and matching
    % This is the number of vacancies, i.e. employees desired by firms
    V_i=N_d_i-N_i;
    % This is search and matching for the labor market
    % New amount of employment per firm N_i, and the employment status
    % variable O_h of who is unemployed or working at which firm
    [N_i,O_h]=search_and_matching_labor(N_i,V_i,O_h);
    
    % (A.26) auf S. 5: overtime and part-time change the wage
    % Amount of hours of all employees of each firm are reduced.
    w_i=w_bar_i.*min(1.5,min(Q_s_i,min(K_i.*kappa_i,M_i.*beta_i))./(N_i.*alpha_bar_i));
    % adapt labour productivity of firm to part- and overtime!
    alpha_i=alpha_bar_i.*min(1.5,min(Q_s_i,min(K_i.*kappa_i,M_i.*beta_i))./(N_i.*alpha_bar_i));
    
    % HERE production happens according to supply choice and Leontief
    % limitational production function (A.13) p. 3
    Y_i=min(Q_s_i,min(N_i.*alpha_i,min(K_i.*kappa_i,M_i.*beta_i)));
    
    I=length(G_i);
    H_W=H_act-I-1;
    % Here, the Status is updated wenn as soon as somebody finds a job,
    % (A.41), wages are allocated to employed persons per firm (sectorally
    % average wages)
    % There exist version of the model where wages are per firm and income
    % distribution of Austria are included (in progress)
    for h=1:H_W
       i=O_h(h);
       if i~=0
          w_h(h)=w_i(i);
       end
    end
    
    % We update the social benefits with the expected growth rate (inflation accounting is later in code), (A.42 and A.43)
    sb_other=sb_other*(1+gamma_e);
    sb_inact=sb_inact*(1+gamma_e);
    
    % Expected profits of bank updated by exp. inflation and growth rates
    % A.45
    Pi_e_k=Pi_k*(1+pi_e)*(1+gamma_e);
    
    % Here, households calcualate their expected income, (A.44) on p. 8
    % Household consume first in the period, at the end they receive their
    % wages
    H=H_act+H_inact;
    Y_e_h=zeros(1,H);
    for h=1:H
        if h<=H_W
            if O_h(h)~=0
                Y_e_h(h)=(w_h(h)*(1-tau_SIW-tau_INC*(1-tau_SIW))+sb_other)*P_bar_HH*(1+pi_e);
            else
                Y_e_h(h)=(theta_UB*w_h(h)+sb_other)*P_bar_HH*(1+pi_e);
            end
        elseif h>H_W && h<=H_W+H_inact
            Y_e_h(h)=(sb_inact+sb_other)*P_bar_HH*(1+pi_e);
        elseif h>H_W+H_inact && h<=H_W+H_inact+I
            i=h-(H_W+H_inact);
            Y_e_h(h)=theta_DIV*(1-tau_INC)*(1-tau_FIRM)*max(0,Pi_e_i(i))+sb_other*P_bar_HH*(1+pi_e);
        elseif h>H_W+H_inact+I && h<=H
            Y_e_h(h)=theta_DIV*(1-tau_INC)*(1-tau_FIRM)*max(0,Pi_e_k)+sb_other*P_bar_HH*(1+pi_e);
        end
    end
    
    % Consumption is a fraction of expected disposable HH income, corrected
    % for VAT, (A.46) p. 8
    C_d_h=(psi*Y_e_h)/(1+tau_VAT);
    
    % Household investment is a fraction of expected disposable HH income, corrected
    % for this investment tax, (A.49) p. 9
    I_d_h=psi_H*Y_e_h/(1+tau_CF);
    
    % Econometric estimations of exogenous variables gov. cons, exports,
    % imports takes place, this is what is different in ABMX!
    % Normally distributed shock
    epsilon_G=normrnd(0,sigma_G);
    % Gov. cons (A.55)
    C_G=exp(alpha_G*log(C_G)+beta_G+epsilon_G);
    C_d_j=C_G/J*ones(1,J)*sum(c_G_g.*P_bar_g)*(1+pi_e);
    
    % Exports in RoW (A.81)
    C_E=exp(alpha_E*log(C_E)+beta_E+epsilon_E);
    C_d_l=C_E/L*ones(1,L)*sum(c_E_g.*P_bar_g)*(1+pi_e);
    
    % Imports (A.76)
    Y_I=exp(alpha_I*log(Y_I)+beta_I+epsilon_I);
    Y_m=c_I_g'*Y_I;
    P_m=P_bar_g'*(1+pi_e);
    
    % This is where the goods markets for consumption, investment,
    % intermediate inputs and exports takes place - this I guess takes the
    % most time in the model!
    [Q_d_i,Q_d_m,P_bar_i,DM_i,P_CF_i,I_i,P_bar_h,C_h,P_bar_CF_h,I_h,P_j,C_j,P_l,C_l]=search_and_matching(P_i,Y_i,S_i,K_i.*kappa_i-Y_i,G_i,P_m,Y_m,a_sg,DM_d_i,b_CF_g,I_d_i,P_bar_g.*b_HH_g/sum(P_bar_g.*b_HH_g),C_d_h,P_bar_g.*b_CFH_g/sum(P_bar_g.*b_CFH_g),I_d_h,P_bar_g.*c_G_g/sum(P_bar_g.*c_G_g),C_d_j,P_bar_g.*c_E_g/sum(P_bar_g.*c_E_g),C_d_l);
    
	% Realized sales domestic Q is limited by production possibilities Y + inventories S and demand Q_d, equ. (A.2)
    Q_i=min(Y_i+S_i,Q_d_i);
    
	% Realized sales from imports (A.79)
    Q_m=min(Y_m,Q_d_m);
    
	% Capital stock law of motion for household investment (A.52)
    K_h=K_h+I_h;
    
	% Inflation rate, equ. (A.11) on p. 3 
    pi(T_prime+t)=log(sum(P_i.*Y_i)/sum(Y_i)/P_bar);
	% Update P_bar to the current time step
    P_bar=sum(P_i.*Y_i)/sum(Y_i);
    
	% Producer index price for goods, equ. (A.32)
    for g=1:G
        P_bar_g(g)=(sum(P_i(G_i==g).*Q_i(G_i==g))+P_m(g)*Q_m(g))/(sum(Q_i(G_i==g))+Q_m(g));
    end
    
	% Price indices for capital formation CF (investment firms) 
    P_bar_CF=sum(b_CF_g.*P_bar_g);
	% (A.31) consumer price index
    P_bar_HH=sum(b_HH_g.*P_bar_g);
    
	% Capital stock law of motion for firms, equ. (A.17)
    K_i=K_i-delta_i./kappa_i.*Y_i+I_i;
    
    % Intermediate inputs law of motion (A.21)
    M_i=M_i-Y_i./beta_i+DM_i;
	
    % Changes in inventories (A.3)
    DS_i=Y_i-Q_i;
    % Inventories updated with change (A.4)
    S_i=S_i+DS_i;
    
    % Firm profits (A.30)
    Pi_i=P_i.*Q_i+P_i.*DS_i-(1+tau_SIF)*w_i.*N_i*P_bar_HH-1./beta_i.*P_bar_i.*Y_i-delta_i./kappa_i.*P_CF_i.*Y_i-tau_Y_i.*P_i.*Y_i-tau_K_i.*P_i.*Y_i-r*(L_i+max(0,-D_i))+r_bar*max(0,D_i);
    
    % Profit of bank, (A.65)
    Pi_k=r*sum(L_i+max(0,-D_i))+r*sum(max(0,-D_h))+r_bar*max(0,D_k)-r_bar*sum(max(0,D_i))-r_bar*sum(max(0,D_h))-r_bar*max(0,-D_k);
    
    % Bank equity updated, (A.67)
    E_k=E_k+Pi_k-theta_DIV*(1-tau_FIRM)*max(0,Pi_k)-tau_FIRM*max(0,Pi_k);
    
    % Income of households, (A.53)
    Y_h=zeros(1,H);
    for h=1:H
        if h<=H_W
            if O_h(h)~=0
                Y_h(h)=(w_h(h)*(1-tau_SIW-tau_INC*(1-tau_SIW))+sb_other)*P_bar_HH;
            else
                Y_h(h)=(theta_UB*w_h(h)+sb_other)*P_bar_HH;
            end
        elseif h>H_W && h<=H_W+H_inact
            Y_h(h)=(sb_inact+sb_other)*P_bar_HH;
        elseif h>H_W+H_inact && h<=H_W+H_inact+I
            i=h-(H_W+H_inact);
            Y_h(h)=theta_DIV*(1-tau_INC)*(1-tau_FIRM)*max(0,Pi_i(i))+sb_other*P_bar_HH;
        elseif h>H_W+H_inact+I && h<=H
            Y_h(h)=theta_DIV*(1-tau_INC)*(1-tau_FIRM)*max(0,Pi_k)+sb_other*P_bar_HH;
        end
    end
    
    % Deposits of households increasae/decreas with savings and/or interest rates (A.54)
    D_h=D_h+Y_h-(1+tau_VAT)*C_h-(1+tau_CF)*I_h+r_bar*max(0,D_h)-r*max(0,-D_h);
    
    % Profit of central bank, (A.73)
    pi_CB=r_G*L_G-r_bar*D_k;
    
    % Governnment accounting
    % Government revenues, (A.59)
    Y_G=(tau_SIF+tau_SIW)*sum(w_h(O_h~=0))*P_bar_HH+tau_INC*(1-tau_SIW)*P_bar_HH*sum(w_h(O_h~=0))+tau_VAT*sum(C_h)+tau_INC*(1-tau_FIRM)*theta_DIV*(sum(max(0,Pi_i))+max(0,Pi_k))+tau_FIRM*(sum(max(0,Pi_i))+max(0,Pi_k))+tau_CF*sum(I_h)+sum(tau_Y_i.*P_i.*Y_i)+sum(tau_K_i.*P_i.*Y_i)+tau_EXPORT*C_l;
    % Government deficit (or surplus) (A.60)
    Pi_G=C_j+r_G*L_G+H_inact*sb_inact*P_bar_HH+theta_UB*sum(w_h(O_h==0))*P_bar_HH+H*sb_other*P_bar_HH-Y_G;
    % Update of government debt (A.61)
    L_G=L_G+Pi_G;
    
    % Firm accounting
    % Realized cash flow of firm (A.33)
    DD_i=P_i.*Q_i-(1+tau_SIF)*w_i.*N_i*P_bar_HH-DM_i.*P_bar_i-P_CF_i.*I_i-tau_Y_i.*P_i.*Y_i-tau_K_i.*P_i.*Y_i-r*(L_i+max(0,-D_i))+r_bar*max(0,D_i)+DL_i-theta*L_i-tau_FIRM*max(0,Pi_i)-theta_DIV*(1-tau_FIRM)*max(0,Pi_i);
    % Deposits are updated (A.34)
    D_i=D_i+DD_i;
    % Loans are updated (A.35)
    L_i=(1-theta)*L_i+DL_i;
    % Equity of firm updated mark-to-market (A.36)
    E_i=D_i+M_i.*sum(a_sg(:,G_i).*P_bar_g)+P_i.*S_i+P_bar_CF*K_i-L_i;
    
    % Accounting RoW
    % Equity "OeNB" (A.73)
    E_CB=E_CB+pi_CB;
    % (A.75) Capital account is added to net foreign assets
    D_RoW=D_RoW-(1+tau_EXPORT)*C_l+sum(P_m.*Q_m);
    % Reserves at central bank (A.68)
    D_k=sum(D_i)+sum(D_h)+E_k-sum(L_i);
    
    % Times series of aggregate output is updated
    Y(T_prime+t)=sum(Y_i);
    
    % Macro time series are updated (which are the output of the function),
    % see Appendix A.7. Macroeconomic aggregates on p. 16
    nominal_gdp(t)=sum(tau_Y_i.*Y_i.*P_i)+tau_VAT*sum(C_h)+tau_CF*sum(I_h)+tau_G*C_j+tau_EXPORT*C_l+sum((1-tau_Y_i).*P_i.*Y_i)-sum(1./beta_i.*P_bar_i.*Y_i);
    real_gdp(t)=sum(Y_i.*((1-tau_Y_i)-1./beta_i))+sum(tau_Y_i.*Y_i)+tau_VAT*sum(C_h)/P_bar_h+tau_CF*sum(I_h)/P_bar_CF_h+tau_G*C_j/P_j+tau_EXPORT*C_l/P_l;
    nominal_gva(t)=sum((1-tau_Y_i).*P_i.*Y_i)-sum(1./beta_i.*P_bar_i.*Y_i);
    real_gva(t)=sum(Y_i.*((1-tau_Y_i)-1./beta_i));
    nominal_household_consumption(t)=(1+tau_VAT)*sum(C_h);
    real_household_consumption(t)=(1+tau_VAT)*sum(C_h)/P_bar_h;
    nominal_government_consumption(t)=(1+tau_G)*C_j;
    real_government_consumption(t)=(1+tau_G)*C_j/P_j;
    nominal_capitalformation(t)=sum(P_CF_i.*I_i)+(1+tau_CF)*sum(I_h)+sum(DS_i.*P_i)+sum(DM_i.*P_bar_i-1./beta_i.*P_bar_i.*Y_i);
    real_capitalformation(t)=sum(I_i)+(1+tau_CF)*sum(I_h)/P_bar_CF_h+sum(DM_i-Y_i./beta_i)+sum(DS_i);
    nominal_fixed_capitalformation(t)=sum(P_CF_i.*I_i)+(1+tau_CF)*sum(I_h);
    real_fixed_capitalformation(t)=sum(I_i)+(1+tau_CF)*sum(I_h)/P_bar_CF_h;
    nominal_fixed_capitalformation_dwellings(t)=(1+tau_CF)*sum(I_h);
    real_fixed_capitalformation_dwellings(t)=(1+tau_CF)*sum(I_h)/P_bar_CF_h;
    nominal_exports(t)=(1+tau_EXPORT)*C_l;
    real_exports(t)=(1+tau_EXPORT)*C_l/P_l;
    nominal_imports(t)=sum(P_m.*Q_m);
    real_imports(t)=sum(Q_m);
    operating_surplus(t)=sum(P_i.*Q_i+P_i.*DS_i-(1+tau_SIF)*w_i.*N_i*P_bar_HH-1./beta_i.*P_bar_i.*Y_i-tau_Y_i.*P_i.*Y_i-tau_K_i.*P_i.*Y_i);
    compensation_employees(t)=(1+tau_SIF)*sum(w_i.*N_i)*P_bar_HH;
    wages(t)=sum(w_i.*N_i)*P_bar_HH;
    taxes_production(t)=sum(tau_K_i.*Y_i.*P_i);
    
    % Sectoral Gross Value added per product is calculated
    for g=1:G
        nominal_sector_gva(t,g)=sum((1-tau_Y_i(G_i==g)).*P_i(G_i==g).*Y_i(G_i==g))-sum(1./beta_i(G_i==g).*P_bar_i(G_i==g).*Y_i(G_i==g));
        real_sector_gva(t,g)=sum(Y_i(G_i==g).*((1-tau_Y_i(G_i==g))-1./beta_i(G_i==g)));
    end
    
    euribor(t)=r_bar;
    gdp_deflator_growth_ea(t)=pi_EA;
    real_gdp_ea(t)=Y_EA;
    
    % Bankrupt firms are found and their net debtor position is assessed
    % See Appendix A.1.9, p. 7, (A.38) - (A.40)
    insolvent=find(D_i<0 & E_i<0);
    for q=1:length(insolvent)
        i=insolvent(q);
        E_k=E_k-(L_i(i)-D_i(i)-zeta_b*P_bar_CF*K_i(i));
        E_i(i)=E_i(i)+(L_i(i)-D_i(i)-zeta_b*P_bar_CF*K_i(i));
        L_i(i)=zeta_b*P_bar_CF*K_i(i);
        D_i(i)=0;
    end
end

end
