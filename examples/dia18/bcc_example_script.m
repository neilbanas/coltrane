% bcc = "Bering-Chukchi Coltrane"


% ------------------------------------------------------------------------------
% make a set of one-year return maps from BIOMAS
if 0
	for year = 1979:2016
		map = bcc_returnMap(year);
			% each of these calls takes a couple of hours on my laptop
		save(['returnMaps/bcc_returnMap_' num2str(year)]);
	end
end

% make a set of 3-year trajectories
if 0
	load returnMaps/bcc_returnMap_2008 map
	map08 = map;
	load returnMaps/bcc_returnMap_2009 map
	map09 = map;
	load returnMaps/bcc_returnMap_2010 map
	map10 = map;
	P = bcc_3yr_particles(map08,map09,map10);
		% this takes about 10 min on my laptop
		% for some reason the timebase is one year later than I would expect
	save output/bcc_08_09_10_particles P
end


% ------------------------------------------------------------------------------
% run the separable version of Coltrane along those trajectories
if 1
	% pick some particles to run
	load output/bcc_08_09_10_particles P
	ind = find(P.y(1,:) < 65 & P.y(end,:) > 67);
		% particles that start in the Bering Sea and end in the Chukchi Sea
		% 3 yrs later. This is just an example; this doesn't properly select
		% _animals_ that cross from the Bering to the Chukchi, since the model
		% cohorts are spawned across the first whole year of transport. The
		% way to filter animals in that way si to run a broader set of
		% particles and then filter v.x0,v.y0,v.x1yr,v.y1yr at the very end.

	% big P is particles, little p is parameters
	p = coltraneParams_dia18('m0',0.08,'Ks',3,'u0',0.008,'iceToSat',0.7);

	parallel = 0;	
	if ~parallel % i.e., serial
		vfull = coltraneModel_dia18(P,p,ind); % <------
		v = filterCohorts(vfull,'bcc');
		% there are a lot of somewhat arbitrary choices in getting from vfull
		% down to v
	else
		% run the model in a bunch of parallel slices
		Nsplit = round(size(P.x,2)/500);
		clear vn
		parfor n=1:Nsplit
			disp(n);
			ind = n : Nsplit : size(P.x,2);
			vfull = coltraneModel_dia18(P,p,ind);
			vn(n) = filterCohorts(vfull,'bcc');
		end
		% reassemble them
		v = vn(1);
		fields = fieldnames(v);
		for n=2:Nsplit
			for i=1:length(fields)
				v.(fields{i}) = cat(2,v.(fields{i}),vn(n).(fields{i}));
			end
		end
	end
	
	save -v7.3 output/bcc_08_09_10_cohorts v p
		% the -v7.3 flag is necessary if v is huge
end