function [N_i,O_h]=search_and_matching_labor(N_i,V_i,O_h)
% Look for people that are employed, O_h==0 if people are unemployed,
% otherwise the index of the firm.
H_E=find(O_h>0);
% Random re-sorting
H_E=H_E(randperm(length(H_E)));
% See ch. A.1.6 on pp. 5
% FIRING
% Loop over all employed persons in model
for e=1:length(H_E)
    h=H_E(e);
    i=O_h(h);
    % THIS is the "firing function". Bei jeder Firma, die zu viel Personal
    % hat, wird ein zufälliger "Mensch" (Agent) entlassen
    if V_i(i)<0
        % This person becomes unemployed
        O_h(h)=0;
        % Employment of firm i is reduced by one person
        N_i(i)=N_i(i)-1;
        % now the vacancy is = 0, i.e. the the firm has fired all they want
        V_i(i)=V_i(i)+1;
    end
end

% HIRING
% Now we look for the unemployed
H_U=find(O_h==0);
% Firms with vacancies
I_V=find(V_i>0);
% This is the "hiring function"
% "efficient" labour market, i.e. as long as there are vacancies and
% unemployed, search continues. Here, you could include frictional
% unemployment if you wanted to.
while ~isempty(H_U) && ~isempty(I_V)
    I_V=I_V(randperm(length(I_V)));
    for f=1:length(I_V)
        i=I_V(f);
        % Matlab function that returns an integer, uniformly distributed,
        % one unemployed person is picked out of all unemployed
        e=randi(length(H_U));
        h=H_U(e);
        % Is employed at firm i
        O_h(h)=i;
        % Here people are employed, i.e. stock of employees is increased,
        % and vacancies are decreased
        N_i(i)=N_i(i)+1;
        V_i(i)=V_i(i)-1;
        % This person is removed from the vector of unemployed
        H_U(e)=[];
        % As soon as nobody is unemployed (or all vacancies are filled!), loop is finished
        if isempty(H_U)
            break
        end
    end
    % End loop if all vacancies are filled
    I_V=find(V_i>0);
end
end

