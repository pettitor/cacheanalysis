function [hitrate, pin, tC]=hitrateLRU(lm,C,eps)

tC=1;tCtmp = 0;

while abs(tC-tCtmp) > eps
   tCtmp = tC;
   tC = tC - (sum(1-exp(-lm*tC))-C);
end
pin = 1-exp(-lm*tC);

hitrate = sum(lm.*pin);

end