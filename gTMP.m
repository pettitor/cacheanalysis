
N=1000;
par.alpha = 0.99;
active = 1:1000;

a=exp(-par.alpha .* log(1:N));
zipfcdf = cumsum([0 a]);
par.zipfcdf = zipfcdf/zipfcdf(end);

aa = 0:0.1:0.5;
request = nan(length(aa),100000);
for j=1:length(aa)
alpha = aa(j);
active = 1:1000;
maxid = N;
for i=1:100000
    if rand()<alpha
        maxid = maxid + 1;
        vid = maxid;
        active(randi(N)) = vid;
    else
        vid=active(find(par.zipfcdf>rand(),1,'first')-1);
    end
    request(j,i) = vid;
    
end
end
%%
figure(1);clf;hold all;
for j=1:length(aa)
    count = hist(request(j,:),1:max(request(j,:)));
    plot(sort(count,'descend'))
end
set(gca,'xscale','log','yscale','log')

plot(1:1e4,500*((1:1e4)).^-(0.6))
plot(1:1e4,10000*((1:1e4)).^-(0.9))
%%
count = hist(request(end,:),1:max(request(end,:)));
firstocc = nan(1,length(count));
lastocc = nan(1,length(count));
for i=1:length(count)
    first = find(request(end,:)==i,1,'first');
    if first
        firstocc(i) = first;
    end
    last = find(request(end,:)==i,1,'last');
    if last
        lastocc(i) = last;
    end
end