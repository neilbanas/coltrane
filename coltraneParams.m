function p = coltraneParams(varargin);

% p = coltraneParams('param1',val1,'param2',val2,...);
%
% returns a complete set of parameters for coltraneModel.m. Specify whichever
% non-default values you like and this will fill in the rest.

p = struct(varargin{:});

% the reason for the setDefault() syntax is to allow user-specified parameters
% to be inserted into the sequence at the same place they would be declared by
% default; this makes the dependencies among parameters a little more
% predictable. It may not be necessary at all.
 
p=setDefault(p,'dt_spawn',10); % resolution (d) of spawning date cases
p=setDefault(p,'tdia_exit',[]); % set of diapause exit dates to consider
	% if this is empty, constructs a set using dt_dia
p=setDefault(p,'tdia_enter',[]); % set of diapause entry dates to consider
	% if this is empty, constructs a set using dt_dia
p=setDefault(p,'dt_dia',20);
p=setDefault(p,'dtegg',[]);
	% if this is empty, constructs a set using min_ and max_genlength_years
p=setDefault(p,'min_genlength_years',0);
p=setDefault(p,'max_genlength_years',Inf);
	% range of generation lengths to evaluate (in integer years)


p=setDefault(p,'preySatVersion','biomas_dia21');
p=setDefault(p,'KsIA',0.2);
p=setDefault(p,'iceToSat',1);


p=setDefault(p,'r_ea',0.013); % egg:adult weight ratio if exp_ea=1
p=setDefault(p,'exp_ea',0.62); % egg weight = r_ea * adult weight^exp_ea
	% Ki¿rboe and Sabatini 1995 (Table 1, Fig 1):
	% sac spawners: r_ea = 0.014, exp_ea = 1
	% broadcast spawners: r_ea = 0.013, exp_ea = 0.62
p=setDefault(p,'theta',0.7); % metabolic scaling exponent
p=setDefault(p,'u0',0.008); % food-saturated development rate at T = 0 
p=setDefault(p,'I0',0.37); % food-saturated ingestion rate at S=1 µgC, T=0
	% This was 0.4 in Banas et al. 2016 and Hobbs et al. 2020,
	% but tuning in June 2021 showed that lower I0 was better for the Wa-u0 relationship
p=setDefault(p,'GGE_nominal',0.33); % for relating adult size to I0 and u0.
p=setDefault(p,'Q10g',2.5); % Q10 for growth and ingestion
p=setDefault(p,'Q10d',3); % Q10 for development
p=setDefault(p,'Ks',1.2); % prey half-saturation; same units as P
	% from tuning, June 2021
p=setDefault(p,'Df',0.10); % age of first feeding (0.10 = start of N3)
p=setDefault(p,'Ds',0.35); % age at which to start storing lipids
						   % (0.35 = start of C1)
p=setDefault(p,'Ddia',0.6); % minimum diapause-capable stage (0.6 = start of C4;
							% (alternatively this can be set to 0 or Ds and 
							% evaluated in postprocessing)
p=setDefault(p,'requireActiveSpawning',0);
	% omit cases in which t0 falls during the diapause period, and set E=0
	% during diapause as well
p=setDefault(p,'allowGainAfterD1',0);
	% is it allowed to put on mass after D=1 but before dtegg is reached?
	% (Nov 2020: turning this on doesn't help survival and generates huge R/W values)
p=setDefault(p,'r_assim',0.67); % assimilation efficiency of ingestion
p=setDefault(p,'rm',0.8 * 0.17);
	% active metabolism as fraction of max assimilation
	% 0.8*0.17 means that GGE = 0 at P = 0.25 Ks
p=setDefault(p,'rb',0.25); % metabolism at a=0 as fraction of metabolism at a=1
p=setDefault(p,'rstarv',0.1); % starvation tolerance

% set mortality (m0 = mortality at W = 1 µgC, T = 0)
p=setDefault(p,'m0_over_GGE_I0',0.67);
p=setDefault(p,'m0',p.m0_over_GGE_I0 * p.GGE_nominal * p.I0);
p.m0_over_GGE_I0 = p.m0./p.GGE_nominal./p.I0;
% this is probably overcomplicated at this point. Tuning June 2021 says m0 = 0.065


% ------------------------------------------------------------------------------
function p=setDefault(p0,name,val);

% if (name) is not a field in the param structure p0, fills in the default value 
% (val). otherwise, leaves it alone
p = p0;
if ~isfield(p,name) | isempty(p.(name))
	p.(name) = val;
end