function tc = yearday(t);

[yr,mo,dy,hr,mn,sc] = datevec(t);
tc = datenum(0.*yr,mo,dy,hr,mn,sc);