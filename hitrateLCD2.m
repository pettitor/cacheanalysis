function [phit, pin, tC]=hitrateLCD2(lm1,C1,C2,eps)

tC=[1 1];C1tmp=0;C2tmp=0;
pin = [ones(1,length(lm1))*C1/length(lm1);ones(1,length(lm1))*C2/length(lm1)];
phit = pin;


while abs(C1-C1tmp) > eps || abs(C2-C2tmp) > eps
   phit(1,:) = ((1-pin(1,:)).*phit(2,:)+pin(1,:)).*(1-exp(-lm1*tC(1)));
   pin(1,:) = phit(1,:);
   C1tmp = sum(pin(1,:));
   tC(1) = tC(1) * C1/C1tmp;
    
   lm2 = (1-pin(1,:)).*lm1;
   pin(2,:) = (1-exp(-lm2*tC(2)));
   C2tmp = sum(pin(2,:));
   tC(2) = tC(2) * C2/C2tmp;
   
   phit(2,:) = pin(2,:); %TODO correct?
%    if tC(2) > tC(1)
%        phit(2,:) = (1-exp(-lm2.*(tC(2)-tC(1))))+(1-phit(2,:)).*(1-exp(-lm1*tC(1)));
%    else
%        phit(2,:) = (1-phit(2,:)).*(1-exp(-lm1*tC(2)));
%    end
%    phit(2,phit(2,:)>1) = 1;
end

%hitrate = sum(lm/sum(lm).*pin);

end