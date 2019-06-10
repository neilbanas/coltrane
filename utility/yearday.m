function tc = yearday(t);

%[yr,mo,dy,hr,mn,sc] = datevec(t);
%tc = datenum(0.*yr,mo,dy,hr,mn,sc);

[yr,~] = datevec(t);
t0 = datenum(yr,0,0);
tc = t - reshape(t0,size(t));

