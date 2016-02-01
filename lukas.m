N = 1e4;
alpha = 0.99;
pview = (1:N).^(-alpha);
q = pview/sum(pview);

C = 20;

par.C = (10.^(-4:.1:-0.1))*N;
hitrate = NaN(1,length(par.C));

for j=1:length(par.C)
C = par.C(j);


% Ctmp = 0;
% tC = 0;
% eps = 1e-4;
% while abs(C-Ctmp) > eps
%    Ctmp = sum(1-exp(-q*tC));
%    tC = tC - (Ctmp-C);
% end
% phit = 1-exp(-q*tC);
%C = 1;

hitrate(j) = hitrateLRU(q,C,1e-4);

end
figure(3);clf;box on;
plot(par.C,hitrate,'--','color','black','linewidth',2)
set(gca,'xscale','log')
xlabel('cache size C / catalouge size')
ylabel('hitrate p_{hit}')