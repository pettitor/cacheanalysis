load('/home/valli/Studis/Sebastian Stauch/asvalues.mat')

%%

N = 100000;

Cleafs = 4;
Nleafs = 100;

alpha = 0.99;
pview = (1:N).^(-alpha);
q = pview/sum(pview);

C = 0.01*N;

par.C = (10.^(-4:.1:-0.1))*N;
hitrate = NaN(1,length(par.C));

eps = 1e-4;

%parD24;

par.pshare = 10.^(-4:0.5:-1);
%par.ASnuser = 10.^(1:6);

%nuser = par.ASnuser;
nuser = [1e4 1e5 1e6];

par.pshare = 10.^(-log10(nuser(k)):0.5:-1);

%hitrate = NaN(length(nuser),length(par.pshare));

hitrate = NaN(1,length(par.pshare));

%for j=1:length(par.C)
%C = par.C(j);
for j=1:length(par.pshare)
%for k=1:length(nuser)
Nleafs = nuser(k)*par.pshare(j);

[hitrateRef, phit, tC] = hitrateLRU(q,C,eps);

% phit UNaDa
[hitrateLeaf, phit, tC2] = hitrateLRU(q,Cleafs,eps);

%improve c.f. tandem

lmu = Nleafs*(q.*(1-phit));
lmu = q.*(1-phit).^Nleafs;
lmu = lmu/sum(lmu);

% phit rest of unadas each item stored only once per AS
% [hitrateLeaf, phit] = hitrateLRU(lmu,(Nleafs-1)*Cleafs,eps);
% 
% lm = lmu.*(1-phit);
% lm = lm/sum(lm);

[~, pin1, tC1] = hitrateLRU(lmu,C,eps);

A = q.*(1-phit)*max(0,tC1 - tC2) + (Nleafs-1)*q.*(1-phit)*tC1;
phitc = 1-exp(-A/Nleafs);

% TODO wrong
Ao = q.*(1-phit).^Nleafs*max(0,tC1 - tC2) + (Nleafs-1)*q.*(1-phit).^Nleafs*tC1;
%Ao = q.*(1-phit).^Nleafs*max(0,tC1 - tC2) + (Nleafs-1)*q.*(1-phit).^Nleafs*tC1;
%Ao = Nleafs*q.*(1-phit)*max(0,tC1 - tC2);
%Ao = lmu.*(1-phit)*max(0,tC1 - tC2) + (Nleafs-1)*lmu.*(1-phit)*tC1;
phito = 1-exp(-Ao/Nleafs);


% hitrateR(j) = hitrateRef;
% hitrate(j) = sum(lmu.*phitc);
% hitrate(j) = sum(lmu.*pin1);
 hitrate(j) = sum(lmu.*phito);
%hitrate(k,j) = sum(lmu.*phito);

end
%end
%%

figure(1);hold all;box on;
%color = copper(size(hitrate,1));
color = copper(3);
for i=1:1%size(hitrate,1)
plot(par.pshare,hitrate(i,:),'-+','linewidth',2,'color',color(3,:))
%plot(par.pshare,hitrateRef,'-+','linewidth',2,'color',color(i,:))
end
%plot(par.pshare,(par.ASnuser/sum(par.ASnuser))*(hitrate),'--','color','black','linewidth',2)
%plot(par.pshare,mean(hitrate),'--','color','black','linewidth',2)
set(gca,'xscale','log')
xlabel('cache size C / catalouge size')
ylabel('hitrate p_{hit}')
%% as rank vs hitrate

figure(2);clf;hold all;box on;
color = copper(size(hitrate,1));
plot(1:5,hitrate(1:5,par.pshare==1e-4),'-+','linewidth',2,'color',color(i,:))
%plot(par.pshare,(par.ASnuser/sum(par.ASnuser))*(hitrate),'--','color','black','linewidth',2)
set(gca,'xscale','log')
xlabel('cache size C / catalouge size')
ylabel('hitrate p_{hit}')