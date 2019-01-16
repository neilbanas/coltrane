function [t0,s] = timingCombinations(forcing,p);

% [t0,s] = timingCombinations(forcing,p);


% spawning date
t0 = forcing.t(1) : p.dt_spawn : (forcing.t(end) - 365);

% yearday of diapause exit
tdia_exit = p.tdia_exit;
if isempty(tdia_exit)
	tdia_exit = 0 : p.dt_dia : 365/2;
end

% yearday of diapause entry
tdia_enter = p.tdia_enter;
if isempty(tdia_enter)
	tdia_enter = (max(tdia_exit) + p.dt_dia) : p.dt_dia : 365;
end

% the date that egg production begins relative to t0
dtegg = p.dtegg;
if isempty(dtegg)
	dteggmin = (p.min_genlength_years - 0.5) .* 365;
	dteggmin = max(dteggmin, p.dt_spawn);
	dteggmax = (p.max_genlength_years + 0.5) .* 365;
	dteggmax = min(dteggmax, forcing.t(end) - t0(end));
	dtegg = dteggmin : p.dt_spawn : dteggmax;
end

% strategy vector _s_ (conceptually a vector, but in practice a structure)
[s.tdia_exit, s.tdia_enter, s.dtegg] = ndgrid(tdia_exit, tdia_enter, dtegg);
