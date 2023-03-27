function x=randpl(n,alpha,N)
% RANDPL generates n observations distributed as powerlaw.
x=(1-rand(n,1)).^(-1/(alpha-1));
x=round(x/sum(x)*N);
x(x<1)=1;
dx=sum(x)-N;
while dx~=0
    if dx<0
        id=randperm(length(x),abs(dx));
        x(id)=x(id)+1;
    elseif dx>0
        id=find(x>1);
        id=id(randperm(length(id),min(abs(dx),length(id))));
        x(id)=x(id)-1;
    end
    dx=sum(x)-N;
end
end

