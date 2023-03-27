function [nominal_gdp,real_gdp,nominal_gva,real_gva,nominal_household_consumption,real_household_consumption,nominal_government_consumption,real_government_consumption,nominal_capitalformation,real_capitalformation,nominal_fixed_capitalformation,real_fixed_capitalformation,nominal_fixed_capitalformation_dwellings,real_fixed_capitalformation_dwellings,nominal_exports,real_exports,nominal_imports,real_imports,operating_surplus,compensation_employees,wages,taxes_production,nominal_sector_gva,real_sector_gva,euribor,E_CB,D_RoW,L_G,D_k,D_i,D_h,E_k,L_i]=abmx(G,H_act,H_inact,J,L,tau_INC,tau_FIRM,tau_VAT,tau_SIF,tau_SIW,tau_EXPORT,tau_CF,tau_G,theta_UB,psi,psi_H,theta_DIV,theta,mu,r_G,zeta,zeta_LTV,zeta_b,alpha_bar_i,beta_i,kappa_i,delta_i,w_bar_i,tau_Y_i,tau_K_i,b_CF_g,b_CFH_g,b_HH_g,c_G_g,c_E_g,c_I_g,a_sg,G_i,T,T_prime,T_max,P_i,K_i,M_i,S_i,N_i,D_i,L_i,D_h,w_h,K_h,L_G,E_k,E_CB,D_RoW,O_h,sb_inact,sb_other,Y,pi,r_bar,C_G,C_E,Y_I,P_bar,P_bar_g,P_bar_HH,P_bar_CF,Q_d_i,Pi_i,Pi_k,D_k)
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

for t=1:T_max
    [alpha_Y,beta_Y,gamma_G,gamma_E,gamma_I,epsilon_Y]=estimate_with_predictors(log(Y(1:T_prime+t-1)),[log(C_G(1:T_prime+t-1)),log(C_E(1:T_prime+t-1)),log(Y_I(1:T_prime+t-1))]);
    Y_e=exp(alpha_Y*log(Y(T_prime+t-1))+beta_Y+gamma_G*log(C_G(T_prime+t))+gamma_E*log(C_E(T_prime+t))+gamma_I*log(Y_I(T_prime+t))+epsilon_Y);
    gamma_e=Y_e/Y(T_prime+t-1)-1;
    
    [alpha_pi,beta_pi,gamma_G,gamma_E,gamma_I,epsilon_pi]=estimate_with_predictors(pi(1:T_prime+t-1),[log(C_G(1:T_prime+t-1)),log(C_E(1:T_prime+t-1)),log(Y_I(1:T_prime+t-1))]);
    pi_e=exp(alpha_pi*pi(T_prime+t-1)+beta_pi+gamma_G*log(C_G(T_prime+t))+gamma_E*log(C_E(T_prime+t))+gamma_I*log(Y_I(T_prime+t))+epsilon_pi)-1;
    
    r=r_bar+mu;
    
    Q_s_i=Q_d_i*(1+gamma_e);
    
    pi_c_i=(1+tau_SIF).*w_bar_i./alpha_bar_i.*(P_bar_HH./P_i-1)+1./beta_i.*(sum(a_sg(:,G_i).*P_bar_g)./P_i-1)+delta_i./kappa_i.*(P_bar_CF./P_i-1);
    P_i=P_i.*(1+pi_c_i)*(1+pi_e);
    
    I_d_i=delta_i./kappa_i.*min(Q_s_i,K_i.*kappa_i);
    
    DM_d_i=min(Q_s_i,K_i.*kappa_i)./beta_i;
    
    N_d_i=max(1,round(min(Q_s_i,K_i.*kappa_i)./alpha_bar_i));
    
    Pi_e_i=Pi_i*(1+pi_e)*(1+gamma_e);
    DD_e_i=Pi_e_i-theta*L_i-tau_FIRM*max(0,Pi_e_i)-(theta_DIV*(1-tau_FIRM))*max(0,Pi_e_i);
    DL_d_i=max(0,-DD_e_i-D_i);
    
    K_e_i=P_bar_CF*(1+pi_e)*K_i;
    L_e_i=(1-theta)*L_i;
    DL_i=search_and_matching_credit(DL_d_i,K_e_i,L_e_i,E_k,zeta,zeta_LTV);
    
    V_i=N_d_i-N_i;
    [N_i,O_h]=search_and_matching_labor(N_i,V_i,O_h);
    
    w_i=w_bar_i.*min(1.5,min(Q_s_i,min(K_i.*kappa_i,M_i.*beta_i))./(N_i.*alpha_bar_i));
    alpha_i=alpha_bar_i.*min(1.5,min(Q_s_i,min(K_i.*kappa_i,M_i.*beta_i))./(N_i.*alpha_bar_i));
    
    Y_i=min(Q_s_i,min(N_i.*alpha_i,min(K_i.*kappa_i,M_i.*beta_i)));
    
    I=length(G_i);
    H_W=H_act-I-1;
    for h=1:H_W
       i=O_h(h);
       if i~=0
          w_h(h)=w_i(i);
       end
    end
    
    sb_other=sb_other*(1+gamma_e);
    sb_inact=sb_inact*(1+gamma_e);
    
    Pi_e_k=Pi_k*(1+pi_e)*(1+gamma_e);
    
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
    
    C_d_h=(psi*Y_e_h)/(1+tau_VAT);
    
    I_d_h=psi_H*Y_e_h/(1+tau_CF);
    
    C_d_j=C_G(T_prime+t)/J*ones(1,J)*sum(c_G_g.*P_bar_g)*(1+pi_e);
    
    C_d_l=C_E(T_prime+t)/L*ones(1,L)*sum(c_E_g.*P_bar_g)*(1+pi_e);
    
    Y_m=c_I_g'*Y_I(T_prime+t);
    P_m=P_bar_g'*(1+pi_e);
    
    [Q_d_i,Q_d_m,P_bar_i,DM_i,P_CF_i,I_i,P_bar_h,C_h,P_bar_CF_h,I_h,P_j,C_j,P_l,C_l]=search_and_matching(P_i,Y_i,S_i,K_i.*kappa_i-Y_i,G_i,P_m,Y_m,a_sg,DM_d_i,b_CF_g,I_d_i,P_bar_g.*b_HH_g/sum(P_bar_g.*b_HH_g),C_d_h,P_bar_g.*b_CFH_g/sum(P_bar_g.*b_CFH_g),I_d_h,P_bar_g.*c_G_g/sum(P_bar_g.*c_G_g),C_d_j,P_bar_g.*c_E_g/sum(P_bar_g.*c_E_g),C_d_l);
    
    Q_i=min(Y_i+S_i,Q_d_i);
    
    Q_m=min(Y_m,Q_d_m);
    
    K_h=K_h+I_h;
    
    pi(T_prime+t)=log(sum(P_i.*Y_i)/sum(Y_i)/P_bar);
    P_bar=sum(P_i.*Y_i)/sum(Y_i);
    
    for g=1:G
        P_bar_g(g)=(sum(P_i(G_i==g).*Q_i(G_i==g))+P_m(g)*Q_m(g))/(sum(Q_i(G_i==g))+Q_m(g));
    end
    
    P_bar_CF=sum(b_CF_g.*P_bar_g);
    P_bar_HH=sum(b_HH_g.*P_bar_g);
    
    K_i=K_i-delta_i./kappa_i.*Y_i+I_i;
    
    M_i=M_i-Y_i./beta_i+DM_i;
	
    DS_i=Y_i-Q_i;
    S_i=S_i+DS_i;
    
    Pi_i=P_i.*Q_i+P_i.*DS_i-(1+tau_SIF)*w_i.*N_i*P_bar_HH-1./beta_i.*P_bar_i.*Y_i-delta_i./kappa_i.*P_CF_i.*Y_i-tau_Y_i.*P_i.*Y_i-tau_K_i.*P_i.*Y_i-r*(L_i+max(0,-D_i))+r_bar*max(0,D_i);
    
    Pi_k=r*sum(L_i+max(0,-D_i))+r*sum(max(0,-D_h))+r_bar*max(0,D_k)-r_bar*sum(max(0,D_i))-r_bar*sum(max(0,D_h))-r_bar*max(0,-D_k);
    
    E_k=E_k+Pi_k-theta_DIV*(1-tau_FIRM)*max(0,Pi_k)-tau_FIRM*max(0,Pi_k);
    
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
    
    D_h=D_h+Y_h-(1+tau_VAT)*C_h-(1+tau_CF)*I_h+r_bar*max(0,D_h)-r*max(0,-D_h);
    
    pi_CB=r_G*L_G-r_bar*D_k;
    
    Y_G=(tau_SIF+tau_SIW)*sum(w_h(O_h~=0))*P_bar_HH+tau_INC*(1-tau_SIW)*P_bar_HH*sum(w_h(O_h~=0))+tau_VAT*sum(C_h)+tau_INC*(1-tau_FIRM)*theta_DIV*(sum(max(0,Pi_i))+max(0,Pi_k))+tau_FIRM*(sum(max(0,Pi_i))+max(0,Pi_k))+tau_CF*sum(I_h)+sum(tau_Y_i.*P_i.*Y_i)+sum(tau_K_i.*P_i.*Y_i)+tau_EXPORT*C_l;
    Pi_G=C_j+r_G*L_G+H_inact*sb_inact*P_bar_HH+theta_UB*sum(w_h(O_h==0))*P_bar_HH+H*sb_other*P_bar_HH-Y_G;
    L_G=L_G+Pi_G;
    
    DD_i=P_i.*Q_i-(1+tau_SIF)*w_i.*N_i*P_bar_HH-DM_i.*P_bar_i-P_CF_i.*I_i-tau_Y_i.*P_i.*Y_i-tau_K_i.*P_i.*Y_i-r*(L_i+max(0,-D_i))+r_bar*max(0,D_i)+DL_i-theta*L_i-tau_FIRM*max(0,Pi_i)-theta_DIV*(1-tau_FIRM)*max(0,Pi_i);
    D_i=D_i+DD_i;
    L_i=(1-theta)*L_i+DL_i;
    E_i=D_i+M_i.*sum(a_sg(:,G_i).*P_bar_g)+P_i.*S_i+P_bar_CF*K_i-L_i;
    
    E_CB=E_CB+pi_CB;
    D_RoW=D_RoW-(1+tau_EXPORT)*C_l+sum(P_m.*Q_m);
    D_k=sum(D_i)+sum(D_h)+E_k-sum(L_i);
    
    Y(T_prime+t)=sum(Y_i);
    
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
    
    for g=1:G
        nominal_sector_gva(t,g)=sum((1-tau_Y_i(G_i==g)).*P_i(G_i==g).*Y_i(G_i==g))-sum(1./beta_i(G_i==g).*P_bar_i(G_i==g).*Y_i(G_i==g));
        real_sector_gva(t,g)=sum(Y_i(G_i==g).*((1-tau_Y_i(G_i==g))-1./beta_i(G_i==g)));
    end
    
    euribor(t)=r_bar;
    
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
