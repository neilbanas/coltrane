function pot = coltranePostproc_phi(v,p,forcing);

% postprocessing for the Coltrane potential model
% (see Banas et al., Front. Mar. Sci., 2016)


Ny = p.Nyears+1; % Ny = record length, p.Nyears = time allowed for the 
				 % latest-spawned generation to develop
t = v.t;
pot.Wa_theo = v.Wa_theo;
pot.We_theo = v.We_theo;
F0 = v.phi ./ v.We_theo .* exp(v.lnN) .* double(v.D >= 1);
	% F0(i,j) is potential eggs ("fitness") at time _forcing.t(i,j)_ for a 
	% cohort spawned on _v.tspawn(j)_
pot.F0 = F0;
pot.t = t;
pot.tspawn = v.tspawn;


% rejigger F0 into a square matrix, covering Nyears+1 x Nyears+1,
% with timebase _tspawn_ on both axes
dt_per_dtspawn = (v.tspawn(2)-v.tspawn(1)) / (t(2)-t(1));
N1y = 365 / (v.tspawn(2)-v.tspawn(1)); % = length(v.tspawn)-1, normally
t1 = t(1:dt_per_dtspawn:end);
F1 = F0(1:dt_per_dtspawn:end,:);
Fn = [];
for i=1:p.Nyears+1
	Fn = [Fn [zeros(N1y*(i-1),N1y); F1(1:end-N1y*(i-1),1:N1y)]];
end
Fn(isnan(Fn)) = 0;
pot.Fn = Fn;
pot.t1 = t1;


for j=1:length(v.tspawn)
	pot.Fmax(j) = max(F0(:,j)); % max fitness of a 1-generation strategy
	f = find(mod(t-v.tspawn(j),365)==0 & t>v.tspawn(j));
	f = f(1:p.Nyears);
	pot.Fmult(:,j) = F0(f,j); % fitness at integer # years from spawning
	
	% find optimal sequence of generations
	fs_new = [find(t1==v.tspawn(j)) find(t1==v.tspawn(j)+p.Nyears*365)];
	fs = [];
	while length(fs_new) > length(fs)
		fs = fs_new;
		fs_new = [];
		for k=1:length(fs)-1
			fs_new = [fs_new fs(k)];
			Fkold = Fn(fs(k+1),fs(k));
			% try to insert a new generation in between fs(k) and fs(k+1),
			% see if it increases fitness
			Fknew = [];
			for m = 1 : (fs(k+1)-fs(k)-1)
				Fknew(m) = Fn(fs(k)+m,fs(k)) * Fn(fs(k+1),fs(k)+m);
			end
			if max(Fknew) > Fkold
				mmax = find(Fknew==max(Fknew));
				mmax = mmax(1);
				fs_new = [fs_new fs(k)+mmax];
			end
		end
		fs_new = [fs_new fs(end)];
	end
	Foptj = Fn(fs(2),fs(1));
	for k=2:length(fs)-1
		Foptj = Foptj * Fn(fs(k+1),fs(k));
	end
	pot.Fopt(j) = Foptj; % fitness of optimal sequence of generations
	pot.tgen{j} = t1(fs(1:end-1)); % start times of each generation
	pot.ngen(j) = (length(fs)-1)/p.Nyears; % generations per year
end


% some scalar stats
pot.Fpot = max(pot.Fopt); % total eggs per egg over the model run
f = find(pot.Fopt==max(pot.Fopt));
f = round(mean(f));
pot.tFpot = v.tspawn(f);
pot.Sapot= max(v.S(:,f));
pot.ngenpot = pot.ngen(f); % generations per year for optimal strategy
pot.Fpot_per_yr = pot.Fpot ^ (1/p.Nyears); % eggs per egg per year
pot.Fpot_per_gen = pot.Fpot ^ (1/p.Nyears/pot.ngenpot);
	% eggs per egg per generation, averaged over generations in the
	% optimal strategy
pot.tgen1pot = min(yearday([pot.tgen{f}]));
	% first spawning yearday in the optimal strategy
pot.afracpot = sum(v.a(:,f)==1)/sum(~isnan(v.a(:,f)));
	% amount of time spent not in diapause in the optimal strategy
	
	
% stats for strictly 1 and 2 yrs/gen strategies
pot.Fpot_at1yr = max(pot.Fmult(1,:)); 
pot.tFpot_at1yr = nan;
pot.Sapot_at1yr = nan;
if ~isnan(pot.Fpot_at1yr)
	f = find(pot.Fmult(1,:)==max(pot.Fmult(1,:)));
	f = round(mean(f));
	pot.tFpot_at1yr = v.tspawn(f);
	pot.Sapot_at1yr = max(v.S(:,f));
end	
pot.Fpot_at2yr = nan;
pot.tFpot_at2yr = nan;
pot.Sapot_at2yr = nan;
if p.Nyears > 1
	pot.Fpot_at2yr = max(pot.Fmult(2,:));
	if ~isnan(pot.Fpot_at2yr)
		f = find(pot.Fmult(2,:)==max(pot.Fmult(2,:)));
		f = round(mean(f));
		pot.tFpot_at2yr = v.tspawn(f);
		pot.Sapot_at2yr = max(v.S(:,f));
	end
end
