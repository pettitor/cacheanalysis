addpath('randraw')

%%
N = 1.6e6;
alpha = 0.9;
N = 1.e3;
alpha = 1.5;
pview = (1:N).^(-alpha);
q = pview/sum(pview);

C = 20;

par.C = (10.^(-4:.1:0))*N;
hitrate = NaN(1,length(par.C));

for j=1:length(par.C)
C = par.C(j);
    
t = 1;
maxit = 1e2;
for i=1:maxit
   t = t - (sum(1-exp(-q*t))-C);
end
tC = t;
phit = 1-exp(-q*tC);

hitrate(j) = sum(q.*phit);

end
figure(2);
loglog(par.C/N,hitrate,'--','color','black','linewidth',2)

%% validate q
figure(22);clf;box on;hold all;
b = [1/8 2/8 4/8 8/8];
j = 1;
alpha = 1+b(j);
figure(22);clf;hold all;
plot(sort(histc(randraw('zeta', alpha, 1e5),1:1e5)/1e5,'descend'))
pview = (1:N).^(-alpha);
q = pview/sum(pview);
plot(q);
set(gca,'xscale','log','yscale','log')