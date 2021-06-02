function v = preySaturation(v0,p);

% v = preySaturation(v0,p);
%
% calculate prey saturation based on some settings within p, which is the
% structure that comes out of coltraneParams.m. In Coltrane 1.0, all the
% complexities of this were handled in coltraneForcing.m, before the model was
% run, but doing things in this order makes it possible to vary assumptions
% about the forcing as part of a big parameterisation experiment.

v = v0;

if strcmpi(p.preySatVersion,'biomas_dia21')
	% pelagic phytoplankton as before
	if ~isfield(v,'Ptot')
		v.Ptot = v.flagel + v.diatom;
	end
	% BIOMAS contains two phytoplankton classes. For now just add them
	% together (and ignore the third class of pelagic prey, microzooplankton)

	% prey saturation considering water-column prey only
	v.satWC = v.Ptot ./ (p.Ks + v.Ptot);

	% ice algae index mimicking Castellani et al. 2017, a function of yearday
	% and latitude only. Multiplied by ice cover
	t_init = max(45, 2.78 .* v.y - 117);
	t_max = max(45, 2.08 .* v.y - 30); % linear fits to Castellani et al. Table 5
	t_end = 200; % mid July
	IAind = zeros(size(v.y));
	f = find(v.yday >= t_init & v.yday <= t_max);
	IAind(f)  = v.ice(f) .* (v.yday(f) - t_init(f)) ./ (t_max(f) - t_init(f));
	f = find(v.yday >= t_max & v.yday <= t_end);
	IAind(f) = v.ice(f) .* (1 - (v.yday(f) - t_max(f)) ./ (t_end - t_max(f)));
	
	% prey saturation considering ice-algae only, based on an effective
	% half-saturation (in index units, not chl units). I briefly tried an adjustable
	% iceToSat multiplying this expression but tuning suggested that 1 is about right--
	% whereas KsIA << 1
	v.satIA = IAind ./ (p.KsIA + IAind);
		
	v.sat = max(v.satWC, v.satIA); % overall prey saturation
	
	

elseif strcmpi(p.preySatVersion,'biomas_dia19a')
	% same as biomas_dia19 below, but tIA is a function of latitude chosen so that
	% the time point 1/3 of the way through the dtIA interval aligns with the
	% bloom max a function of latitude from Castellani et al. 2017
	%
	% so dtIA and iceToSat are the free parameters
	
	if ~isfield(v,'Ptot')
		v.Ptot = v.flagel + v.diatom;
	end
	v.satWC = v.Ptot ./ (p.Ks + v.Ptot);

	tIA = max(45, 2.08 .* v.y - 30) - p.dtIA/3; % the new bit
	v.satIA = (v.ice > 0.15) ...
		  .* (v.yday > tIA & v.yday < 200) ...
		  .* p.iceToSat;
	yr = ceil((v.t(:,1)-v.t(1))./365);
	for n=1:max(yr)
		satn = zeros(size(v.satIA));
		satn(yr==n) = v.satIA(yr==n);
		satn(cumsum(double(satn)) > p.dtIA/(v.t(2)-v.t(1))) = 0;
		v.satIA(yr==n) = satn(yr==n);
	end  
	v.sat = max(v.satWC, v.satIA); % overall prey saturation


elseif strcmpi(p.preySatVersion,'biomas_dia19')
	if ~isfield(v,'Ptot')
		v.Ptot = v.flagel + v.diatom;
	end
	% BIOMAS contains two phytoplankton classes. For now just add them
	% together (and ignore the third class of pelagic prey, microzooplankton)

	% prey saturation considering water-column prey only
	v.satWC = v.Ptot ./ (p.Ks + v.Ptot);

	% prey saturation considering (a guess at) ice algae only.
	% days that ice cover is > 15%, after yearday _tIA_ and before yearday
	% 200, for a maximum of _dtIA_ days, all weighted by a highly uncertain weighting
	% factor iceToSat
	v.satIA = (v.ice > 0.15) ...
		  .* (v.yday > p.tIA & v.yday < 200) ...
		  .* p.iceToSat;
	yr = ceil((v.t(:,1)-v.t(1))./365);
	for n=1:max(yr)
		satn = zeros(size(v.satIA));
		satn(yr==n) = v.satIA(yr==n);
		satn(cumsum(double(satn)) > p.dtIA/(v.t(2)-v.t(1))) = 0;
		v.satIA(yr==n) = satn(yr==n);
	end
		  
	v.sat = max(v.satWC, v.satIA); % overall prey saturation
	
	
elseif strcmpi(p.preySatVersion,'satellite_dia18')
	% prey saturation considering water-column prey only
	v.chl(isnan(v.chl) & v.ice > 0.15) = p.chlUnderIce;
	v.chl(isnan(v.chl)) = p.chlUnderPersistentCloud;
	v.satWC = v.chl ./ (p.Ks + v.chl);

	% prey saturation considering (a guess at) ice algae only.
	% days that ice cover is > 15%, after yearday _tIA_ and before yearday
	% (365-tIA), all weighted by a highly uncertain weighting factor iceToSat
	v.satIA = (v.ice > 0.15) ...
		  .* (v.yday > p.tIA & v.yday < 270) ...
		  .* p.iceToSat;
		  
	v.sat = max(v.satWC, v.satIA); % overall prey saturation


else
	% generic version, as in Coltrane 1.0
	v.sat = v.P ./ (p.Ks + v.P);
end