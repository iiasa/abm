function fts = toannual(ftsa)
m=4;
fts=zeros(size(ftsa,1),size(ftsa,2)/m);
if size(ftsa,1)>1
    for i=1:size(ftsa,2)/m
        fts(:,i)=sum(ftsa(:,i*m-m+1:i*m)');
    end
else
    for i=1:length(ftsa)/m
        fts(i)=sum(ftsa(i*m-m+1:i*m));
    end
end
end
