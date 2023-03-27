function DL_i=search_and_matching_credit(DL_d_i,K_e_i,L_e_i,E_k,zeta,zeta_LTV)
    DL_i=zeros(size(DL_d_i));
    % Firms that need credit are found
    I_FG=find(DL_d_i>0);
    % Are sorted randomly so indices don't matter
    I_FG=I_FG(randperm(length(I_FG)));
    % Equ. (A.64), referring to equations (A.62 and A.63) on p. 12: firms gets desired loand IF the borrower's
    % loan to value ration AND the capital requirement of bank are
    % fulfilled
    % Here all firms are iterated in a random order
    for f=1:length(I_FG)
        i=I_FG(f);
        % Max. function is here to avoid that a firm does not get a credit
        % if its credit constraint is binding, analogously for the bank
        DL_i(i)=max(0,min(min(DL_d_i(i),zeta_LTV*K_e_i(i)-L_e_i(i)),E_k/zeta-sum(L_e_i+DL_i)));
    end
end

