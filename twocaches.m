
N = 100000;
nuser = 100000;
CISP = 1000;

alpha = 0.99;

pshare = 10.^(-5:0.5:0);

%pshare = 10.^-3;

eps = 1e-4;

pview = (1:N).^(-alpha);
q = pview/sum(pview);

tot = nan(1,length(pshare));

for j=1:length(pshare)
    Nleafs = round(pshare(j)*nuser);

    Cleafs = 4;
    Cleafs = Nleafs*Cleafs;
    Nleafs = 1;
if Cleafs < N
l = nan(Nleafs+1,N);
pin = nan(Nleafs+1,N);
phit = nan(Nleafs+1,N);
tC = nan(Nleafs+1,1);

l(1,:) = q;
[hitrate, pin(1,:), tC(1)]=hitrateLRU(l(1,:),Cleafs,eps);
phit(1,:) = pin(1,:);

l(2,:) = l(1,:).*(1-phit(1,:));
%l(Nleafs+1,:) = l(Nleafs+1,:)/sum(l(Nleafs+1,:));
[hitrate, pin(2,:), tC(2)]=hitrateLRU(l(2,:),CISP,eps);
phit(2,:) = 1-exp(-l(2,:)*max(0,tC(2)-tC(1)));

hitrate = (l./(sum(l,2)*ones(1,length(l))))*phit';
hrs(:,j) = diag(hitrate);

hit1(j) = l(1,:)*phit(1,:)';
hit2(j) = l(2,:)*phit(2,:)';

tothr = (1-cumprod(1-hrs(:,j))); %TODO
tot(j) = tothr(end);
tot2(j) = l(Nleafs,:)*(phit(Nleafs,:)+phit(Nleafs,:))';
else
    hrs(:,j) = [1;0];
    tot(j) = 1; tot2(j) = 1;
end

end