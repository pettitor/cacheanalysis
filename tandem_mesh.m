%% tandem

N = 100000;
nuser = 100000;
CISP = 1000;

alpha = 0.99;

Cleafs = 4;

pshare = 10.^(-5:0.5:-2);

eps = 1e-4;

pview = (1:N).^(-alpha);
q = pview/sum(pview);

tot = nan(1,length(pshare));

for j=1:length(pshare)
Nleafs = round(pshare(j)*nuser);

% Cleafs = Nleafs*Cleafs;
% Nleafs = 1;

l = nan(Nleafs+1,N);
pin = nan(Nleafs+1,N);
phit = nan(Nleafs+1,N);
tC = nan(Nleafs+1,1);

l(1,:) = q;
[hitrate, pin(1,:), tC(1)]=hitrateLRU(l(1,:),Cleafs,eps);
phit(1,:) = pin(1,:);

for i=2:Nleafs 
l(i,:) = l(i-1,:).*(1-phit(i-1,:));
l(i,:) = l(i,:)/sum(l(i,:));
[hitrate, pin(i,:), tC(i)]=hitrateLRU(l(i,:),Cleafs,eps);
phit(i,:) = 1-exp(-l(i,:)*max(0,tC(i)-tC(i-1)));
%sum(l')
end

l(Nleafs+1,:) = l(Nleafs,:).*(1-phit(Nleafs,:));
l(Nleafs+1,:) = l(Nleafs+1,:)/sum(l(Nleafs+1,:));
[hitrate, pin(Nleafs+1,:), tC(Nleafs+1)]=hitrateLRU(l(Nleafs+1,:),CISP,eps);
phit(Nleafs+1,:) = 1-exp(-l(Nleafs+1,:)*(tC(Nleafs+1)-tC(Nleafs)));

hitrate = (l./(sum(l,2)*ones(1,length(l))))*phit';
hrs = diag(hitrate);

tothr = (1-cumprod(1-hrs)); %TODO
tot(j) = tothr(end);
end
%% mesh

for j=1:length(pshare)
Nleafs = round(pshare(j)*nuser);

lhr = q;%+(Nleafs-1)*l*(1-phithr)/Nleafs;
%lhr = (Nleafs-1)*lhr*(1-phithr)/Nleafs+l/Nleafs

%Nleafs*l*(1-phithr)/Nleafs;

tChr=1;tChrtmp = 0;
pinhr = 1-exp(-q*tChr);
while abs(tChr-tChrtmp) > eps
    
   tChrtmp = tChr;
   tChr = tChr - (sum(1-exp(-lhr*tChr))-Cleafs);

   pinhr = 1-exp(-lhr*tChr);
   
   Ahr = q*tChr + (1-Nleafs)*lhr.*(1-pinhr)*tChr/Nleafs;

%TODO decondition
phithr = 1-exp(-Ahr);


   lhr = q+(Nleafs-1)*lhr.*(1-phithr)/Nleafs;
   lhr = lhr/sum(lhr);
   
end
hitratehr = lhr*phithr';

lisp = Nleafs*lhr.*(1-phithr)/Nleafs;
lisp = lisp/sum(lisp);

[hitrateISP, pinisp, tCisp]=hitrateLRU(lisp,CISP,eps);

Aisp = lhr.*(1-pinhr)*max(0,tCisp-tChr)/Nleafs + ...
    (1-Nleafs)*lhr.*(1-pinhr)*tCisp/Nleafs;

phitisp = 1-exp(-Aisp);
hitrateISP = lisp*phitisp';

hrhr(j) = hitratehr;
hrisp(j) = hitrateISP;
hrtot(j) = 1-(1-hitrateISP)*(1-hitratehr);
end
%% tandem LCD

N = 100000;
nuser = 100000;
CISP = 1000;

alpha = 0.99;

Cleafs = 4;

pshare = 10.^(-5:0.5:-2);

eps = 1e-4;

pview = (1:N).^(-alpha);
q = pview/sum(pview);

for j=1:length(pshare)
Nleafs = round(pshare(j)*nuser);

l = nan(Nleafs+1,N);
pin = nan(Nleafs+1,N);
phit = nan(Nleafs+1,N);
tC = nan(Nleafs+1,1);

l(1,:) = q;
[hitrate, pin(1,:), tC(1)]=hitrateLRU(l(1,:),Cleafs,eps);
phit(1,:) = pin(1,:);

for i=2:Nleafs 
l(i,:) = l(i-1,:).*(1-phit(i-1,:));
l(i,:) = l(i,:)/sum(l(i,:));
[hitrate, pin(i,:), tC(i)]=hitrateLRU(l(i,:),Cleafs,eps);
phit(i,:) = 1-exp(-l(i,:)*max(0,tC(i)-tC(i-1)));
%sum(l')
end

l(Nleafs+1,:) = l(Nleafs,:).*(1-phit(Nleafs,:));
l(Nleafs+1,:) = l(Nleafs+1,:)/sum(l(Nleafs+1,:));
[hitrate, pin(Nleafs+1,:), tC(Nleafs+1)]=hitrateLRU(l(Nleafs+1,:),CISP,eps);
phit(Nleafs+1,:) = 1-exp(-l(Nleafs+1,:)*(tC(Nleafs+1)-tC(Nleafs)));

hitrate = (l./(sum(l,2)*ones(1,length(l))))*phit';
hrs = diag(hitrate);

tothr = (1-cumprod(1-hrs)); %TODO
tot(j) = tothr(end);
end