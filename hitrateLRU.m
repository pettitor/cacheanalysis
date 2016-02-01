function [hitrate, pin, tC]=hitrateLRU(lm,C,eps)

Ctmp = 0;
tC = 1;

while abs(C-Ctmp) > eps
   Ctmp = sum(1-exp(-lm*tC));
   tC = tC * C/Ctmp;
end
pin = 1-exp(-lm*tC);

hitrate = sum(lm/sum(lm).*pin);

end