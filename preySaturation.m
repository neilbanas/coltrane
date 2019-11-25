function v = preySaturation(v0,p);

% v = preySaturation(v0,p);
%
% calculate prey saturation based on some settings within p, which is the
% structure that comes out of coltraneParams.m. In Coltrane 1.0, all the
% complexities of this were handled in coltraneForcing.m, before the model was
% run, but doing things in this order makes it possible to vary assumptions
% about the forcing as part of a big parameterisation experiment.

v = v0;

if strcmpi(p.preySatVersion,'biomas_dia19')
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