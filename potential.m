N = 100000;
nuser = 100000;
CISP = 1000;

alpha = 0.99;

Cleafs = 4;

pshare = 10.^(-5:0.5:0);

addC = pshare*Cleafs*nuser;

eps = 1e-4;

pview = (1:N).^(-alpha);
q = pview/sum(pview);

hr = nan(4,length(pshare));

% case solo ISP cache
[hitrate, pin, tC]=hitrateLRU(q,CISP,eps);
hr(1,:) = hitrate*ones(1,length(pshare));

[hitrateLeaf, pin2, tC2] = hitrateLRU(q,Cleafs,eps);
for i=1:length(pshare)

    % fully organized
    addC = round(pshare(i)*Cleafs*nuser);
    if CISP +addC < N
        [hitrateM, pinM, tCM]=hitrateLRU(q,CISP+addC,eps);
    else
        hitrateM = 1;
    end
hr(2,i) = hitrateM;
    
    if addC < N
        [hitrateA, pinA, tCA]=hitrateLRU(q,addC,eps);
    else
        hitrateA = 1;
    end

    %case overlay
    Nleafs = round(pshare(i)*nuser);

    lmu = NaN(3,length(q));
    pin1 = NaN(3,length(q));
    tC1 = NaN(3,1);
    nl = [1 Nleafs/2 Nleafs];
for j=1:length(nl);
lmu(j,:) = q.*(1-pin2).^nl(j);
%lmu(j,:) = q.*(1-pinA);
lmu(j,:) = lmu(j,:)/sum(lmu(j,:));
[~, pin1(j,:), tC1(j)] = hitrateLRU(lmu(j,:),CISP,eps);
end

%lm
% phit UNaDa
Ao = q.*(1-pin2).^Nleafs*max(0,tC1(3) - tC2) + (Nleafs-1)*q.*(1-pin2).^Nleafs.*tC1(3);
Ao = Nleafs*q.*(1-pin2).^Nleafs*max(0,tC1(3) - tC2);% + (Nleafs-1)*q.*(1-pin2).^Nleafs.*tC1(3);
%Ao = q.*(1-phit).^Nleafs*max(0,tC1 - tC2) + (Nleafs-1)*q.*(1-phit).^Nleafs*tC1;
%Ao = Nleafs*q.*(1-phit)*max(0,tC1 - tC2);
%Ao = lmu.*(1-phit)*max(0,tC1 - tC2) + (Nleafs-1)*lmu.*(1-phit)*tC1;
Ao = q.*(1-pin2)*max(0,tC1(3) - tC2) + (Nleafs-1)*q.*(1-pin2).*tC1(3);
phito = 1-exp(-Ao/Nleafs);


A1 = q.*(1-pin2).^1*max(0,tC1(1) - tC2) + (Nleafs-1)*q.*(1-pin2).^1*tC1(1);
phit1 = 1-exp(-A1/Nleafs);

A2 = q.*(1-pin2).^(Nleafs/2)*max(0,tC1(2) - tC2) + (Nleafs-1)*q.*(1-pin2).^(Nleafs/2)*tC1(2);
phit2 = 1-exp(-A2/Nleafs);

%hr(3,i) = sum(q.*(1-(1-pin2).^Nleafs))+sum(lmu(3,:).*phito);%-sum(lmu(3,:).*(1-(1-pin2).^Nleafs).*phito);
hrc(i) = sum(q.*(1-(1-pin2).^Nleafs))
cc(i) = sum(lmu(3,:).*phito)
%hr(4,i) = sum(q.*(1-(1-pin2).^1))+sum(lmu(1,:).*phit1);
%hr(5,i) = sum(q.*(1-(1-pin2).^(Nleafs/2)))+sum(lmu(2,:).*phit2);
end
%%
plot(pshare,hr)
set(gca,'xscale','log')

%%
N=1e6;

bw_up = 1/8;

C_isp = 1e4;

n_hr = 1e5;

C_hr = 5;

rep_f = 1;

alpha = 0.99;
pview = (1:N).^(-alpha);
q = pview/sum(pview);

c1 = 1;

request_rate = q*c1;

serve_rate = bw_up*rep_f*
