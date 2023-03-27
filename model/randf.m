function f=randf(pr_cum_f)
% f=sum(pr_cum_f<rand);
[~,f]=histc(rand,pr_cum_f);
end

