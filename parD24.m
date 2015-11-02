par.nuser = 1829425;

%%% Resource size and distribution (CDN, caches, end-devices)

% model from Internet Cencus, c.f. Sebastian
par.ASn = 100;

par.alpha = 1.5;

a=exp(-par.alpha .* log(1:par.ASn));
zipfcdf = cumsum([a]);
zipfcdf = zipfcdf/zipfcdf(end);

par.AS = nan(1, par.nuser);
for i=1:par.nuser
    par.AS(i) = sum(rand()>zipfcdf)+1;
end

par.ASnuser = hist(par.AS, 1:par.ASn);