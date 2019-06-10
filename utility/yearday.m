function tc = yearday(t);

% tc = yearday(t);

[yr,~] = datevec(t);
t0 = datenum(yr,0,0);
tc = t - reshape(t0,size(t));
