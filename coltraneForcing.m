function [forcing0,pf] = coltraneForcing(varargin);

% [forcing0,pf] = coltraneForcing(pf);
%           ... = coltraneForcing(scenarioName);
%           ... = coltraneForcing(scenarioName,'param1',val1,'param2',val2,...);
%
% defines the forcing for a coltrane run, filling in defaults for the 
% controlling params as needed in one of several pre-defined scenarios.
%
% pass in either the name of scenario or a structure containing the scenario 
% name + some or all of the associated parameter values. This will return the 
% complete param structure and a structure of forcing time series.
%
% For details and context, see Banas et al (2016) Front. Mar. Res., submitted
% (http://neilbanas.com/projects/coltrane)

if isstruct(varargin)
	pf = varargin;
else
	pf = struct('scenario',varargin{:});
end

% initialize the forcing time-series structure
forcing0.t = (0:365)';
	% the last and the first need to represent the same yearday
forcing0.T = nan.*forcing0.t;
forcing0.Tb = nan.*forcing0.t;
forcing0.P = nan.*forcing0.t;
N = length(forcing0.t);




% scenarios --------------------------------------------------------------------


% ------------------------------------------------------------------------------
% simple idealized scenarios: gaussian prey availability, constant T
elseif strcmpi(pf.scenario,'simple-g')
	pf = setDefault(pf,'T0',8);
	pf = setDefault(pf,'Tb0_over_T0',0.4); % based on WOA13 analysis
	pf = setDefault(pf,'tPmax',365/2);
	pf = setDefault(pf,'dtP',90);
	pf = setDefault(pf,'P0',10);

	forcing0.T = pf.T0 .* ones(N,1);
	forcing0.Tb = pf.Tb0_over_T0 .* pf.T0 .* ones(N,1);
	forcing0.P = pf.P0 .* exp(-(forcing0.t-365/2).^2./pf.dtP.^2);



% ------------------------------------------------------------------------------
% Bering Sea; valid across M2, M4, M5, M8.
% T, Tb from ad-hoc fits to the long BESTMAS hindcast in Banas et al., JGR, 2016
% P from stats in Sigler et al. 2014
elseif strcmpi(pf.scenario,'bering')
	pf = setDefault(pf,'tice',120); % date of ice melt
	if pf.tice > 75 				% spring bloom date based on tice
		tspr = pf.tice + 10;
	else
		tspr = 150;
	end
	pf = setDefault(pf,'tspr',tspr);
	pf = setDefault(pf,'tfall',265); % constant fall bloom date
	pf = setDefault(pf,'tIA',45); % first appearance of ice algae
	
	% default magnitudes of each part of the P cycle
	pf = setDefault(pf,'P0win',0.1);
	pf = setDefault(pf,'P0IA',0);
	pf = setDefault(pf,'P0spr',16);
	pf = setDefault(pf,'P0sum',pf.P0spr*0.1);
	pf = setDefault(pf,'P0fall',pf.P0spr*0.5);
	pf = setDefault(pf,'dtbloom',15); % bloom duration
	
	T0mean = min(5.3, 5.3 - 0.044*(pf.tice-75));
		% mean SST from regression to tice
	pf = setDefault(pf,'T0mean',T0mean);
	pf = setDefault(pf,'Tbmean',0.51 * pf.T0mean - 0.97);
		% mean bottom T from T0mean
	
	% define surface temp cycle based on T0mean
	t_Tmax = 245; % date of max temp
	T0spr = 1.5 * pf.T0mean - 3.1; % temp 365/4 days before t_Tmax
	dTsprmax = -0.62 * pf.T0mean + 9.6; % winter-spring temp diff
	dTwinspr = 0.76 * pf.T0mean + 0.2; % spring-summer max temp diff
	% sinusoid defining winter
	forcing0.T = T0spr + dTwinspr .* cos(2.*pi./365.*(forcing0.t-t_Tmax));
	% overwrite half of it with another sinusoid defining summer
	f = find(forcing0.t > t_Tmax-91 & forcing0.t < t_Tmax+91);
	forcing0.T(f) = T0spr + dTsprmax .* cos(2.*pi./365.*(forcing0.t(f)-t_Tmax));
	forcing0.T = max(-1.8, forcing0.T);

	% bottom temp, also based on T0mean (if the default Tbmean - T0mean relation 
	% is used)
	tp = 2.*forcing0.t./365 - 1;
	forcing0.Tb = pf.Tbmean + 1.75.*(tp-tp.^5);
	forcing0.Tb = max(-1.8, forcing0.Tb);
	
	% P cycle
	forcing0.P = pf.P0spr .* exp(-((forcing0.t-pf.tspr)./pf.dtbloom).^2);
	forcing0.P = max(forcing0.P, ...
	                 pf.P0fall .* exp(-((forcing0.t-pf.tfall)./pf.dtbloom).^2));
	f = forcing0.t > pf.tspr & forcing0.t < pf.tfall;
	forcing0.P(f) = max(forcing0.P(f), pf.P0sum);
	forcing0.P(~f) = max(forcing0.P(~f), pf.P0win);
	f = forcing0.t > pf.tIA & forcing0.t < pf.tspr & forcing0.t < pf.tice+10;
	forcing0.P(f) = max(forcing0.P(f), pf.P0IA);


% ------------------------------------------------------------------------------
% Disko Bay, 1996--analytical fit
elseif strcmpi(pf.scenario,'disko')
	% temperature ---------------------------------------------------
	pf=setDefault(pf,'Tjan1',-0.7); % T on Jan 1
	pf=setDefault(pf,'Tmin',-1.8); % winter min T
	pf=setDefault(pf,'Tmax',3.7); % summer max T
	pf=setDefault(pf,'tTmin',105); % time when spring T increase starts
	pf=setDefault(pf,'tTmax',250); % time of T max
	nn = find(forcing0.t >=0 & forcing0.t <= pf.tTmin);
	forcing0.T(nn) = linspace(pf.Tjan1,pf.Tmin,length(nn));
	nn = find(forcing0.t >=pf.tTmin & forcing0.t <= pf.tTmax);
	forcing0.T(nn) = linspace(pf.Tmin,pf.Tmax,length(nn));
	nn = find(forcing0.t >=pf.tTmax & forcing0.t <= 365);
	forcing0.T(nn) = linspace(pf.Tmax,pf.Tjan1,length(nn));
	
	pf=setDefault(pf,'Tb0',1); % deep T year-round
	forcing0.Tb = pf.Tb0 .* ones(size(forcing0.t));
	
	pf=setDefault(pf,'dT',0); % overall T correction
	forcing0.T = forcing0.T + pf.dT;
	forcing0.Tb = forcing0.Tb + pf.dT;
	
	% chlorophyll ---------------------------------------------------
	pf = setDefault(pf,'P0win',0.1);
	pf = setDefault(pf,'P0spr',13);
	pf = setDefault(pf,'P0sum',0.05*pf.P0spr);
	pf = setDefault(pf,'P0aut',5);
	pf = setDefault(pf,'dtPspr',15); % bloom duration
	pf = setDefault(pf,'dtPaut',30);
	pf = setDefault(pf,'tPspr',150);
	pf = setDefault(pf,'tPaut',225);
	
	forcing0.P = pf.P0win .* ones(size(forcing0.t));
	forcing0.P = max(forcing0.P, ...
					 pf.P0spr .* exp(-((forcing0.t-pf.tPspr)./pf.dtPspr).^2));
	forcing0.P = max(forcing0.P, ...
	                 pf.P0aut .* exp(-((forcing0.t-pf.tPaut)./pf.dtPaut).^2));
	fsum = forcing0.t > pf.tPspr & forcing0.t < pf.tPaut;
	forcing0.P(fsum) = max(forcing0.P(fsum),pf.P0sum);

% ------------------------------------------------------------------------------
else
	disp(['what is ' pf.scenario '?']);
end
	
	
	

% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
function pf = setDefault(pf0,name,val);

% if (name) is not a field in the param structure pf0, fills in the default value (val).
% otherwise, leaves it alone
pf = pf0;
if ~isfield(pf,name) | isempty(pf.(name))
	pf.(name) = val;
end
