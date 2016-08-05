function p = coltraneParams(varargin);

% p = coltraneParams('param1',val1,'param2',val2,...);
%
% returns a complete set of parameters for coltraneModel.m. Specify whichever
% non-default values you like and this will fill in the rest.
%
% For full model description, see Banas et al (2016) Front. Mar. Res., submitted
% (http://neilbanas.com/projects/coltrane)


p = struct(varargin{:});

% the reason for the setDefault() syntax is to allow user-specified parameters
% to be inserted into the sequence at the same place they would be declared by
% default; this makes the dependencies among parameters a little more
% predictable. It may not be necessary at all.
 

% simulation setup
p=setDefault(p,'modelVersion','latest');
p=setDefault(p,'Nyears',3); % number of years in the simulation
p=setDefault(p,'dt',0.5); % integration timestep
p=setDefault(p,'dt_tspawn',5); % spacing of initial spawning dates to consider
		
% parameters for one cohort
% all must be either scalars or a row vector of length NC
p=setDefault(p,'r_ea',0.013); % egg:adult weight ratio if exp_ea=1
p=setDefault(p,'exp_ea',0.62); % egg weight = r_ea * adult weight^exp_ea
	% Ki¿rboe and Sabatini 1995 (Table 1, Fig 1):
	% sac spawners: r_ea = 0.014, exp_ea = 1
	% broadcast spawners: r_ea = 0.013, exp_ea = 0.62
p=setDefault(p,'theta',0.7); % metabolic scaling exponent
p=setDefault(p,'u0',0.009); % food-saturated development rate at T = 0 
p=setDefault(p,'I0',0.4); % food-saturated ingestion rate at S=1 µgC, T=0
p=setDefault(p,'Wa0',150); % estimated adult size
p=setDefault(p,'GGE_nominal',0.33); % for relating adult size to I0 and u0.
p=setDefault(p,'T_nominal',8);
p=setDefault(p,'infer_what','Wa0');
	% u0, I0, Wa0 are related to each other (via GGE_nominal, T_nominal, Q10d, 
	% Q10g). One of them needs to be adjusted to match the others. This switch 
	% says which.
p=setDefault(p,'Q10g',2.5); % Q10 for growth and ingestion
p=setDefault(p,'Q10d',3); % Q10 for development
p=setDefault(p,'Ks',3); % prey saturation half-saturation; same units as P
p=setDefault(p,'alpha_s',1);
	% shape of S/R allocation curve (>1: more investment in R)
	% note: not used in Coltrane 1.0
p=setDefault(p,'Df',0.10); % age of first feeding (0.10 = start of N3)
p=setDefault(p,'Ds',0.35); % age at which to start storing lipids
p=setDefault(p,'Ddia',0.49);
	% minimum age of diapause (0.73 = C5; 0 = full flexibility)
p=setDefault(p,'tegg',0);
	% time at which to start reproducing, maturation permitting
p=setDefault(p,'fcap',1); % rate of capital egg prod relative to max possible
p=setDefault(p,'finc',1); % rate of income egg prod relative to max possible
p=setDefault(p,'r_assim',0.67); % assimilation as fraction of ingestion
p=setDefault(p,'rm',0.8 * 0.17);
	% active metabolism as fraction of max assimilation
	% 0.8*0.17 means that GGE = 0 at P = 0.25 Ks
	% 0.17 means that GGE = 0 at sat = 0.25
p=setDefault(p,'rb',0.25); % metabolism at a=0 as fraction of metabolism at a=1
p=setDefault(p,'r_starv',0.1); % if R < -r_starv S, starvation
p=setDefault(p,'r_phi_max',1.5);

% harmonize Wa0, u0, I0
% if you ignore egg mass (as in Banas and Campbell, MEPS, submitted)
% Wa ^ (1-theta) ~ (1-theta) GGE (1-Df) (I0 / u0) (qg/qd)
qoverq = (p.Q10g./p.Q10d).^(p.T_nominal./10);
co = (1-p.theta) .* p.GGE_nominal .* (1-p.Df) .* qoverq;
if strcmpi(p.infer_what, 'Wa0')
	p.Wa0 = (co .* p.I0 ./ p.u0) .^ (1./(1-p.theta));
elseif strcmpi(p.infer_what, 'I0')
	p.I0 = p.u0 .* p.Wa0 ^ (1-p.theta) ./ co;
elseif strcmpi(p.infer_what, 'u0')
	p.u0 = co .* p.I0 ./ p.Wa0.^(1-p.theta);
elseif strcmpi(p.infer_what,'GGE_nominal')
	p.GGE_nominal = p.u0 ./ p.I0 .* p.Wa0 .^ (1-p.theta) ./ ...
		( (1-p.Df) .* qoverq .* (1-p.theta) );
end

% estimated egg size from estimated adult size
p=setDefault(p,'We0', p.r_ea .* p.Wa0.^p.exp_ea);

% set mortality
p=setDefault(p,'m0',0.67 * p.GGE_nominal * p.I0);
	% mortality at S = 1 µgC, T = 0




% ------------------------------------------------------------------------------
function p=setDefault(p0,name,val);

% if (name) is not a field in the param structure p0, fills in the default value 
% (val). otherwise, leaves it alone
p = p0;
if ~isfield(p,name) | isempty(p.(name))
	p.(name) = val;
end