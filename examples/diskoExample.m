% Approximating the Disko Bay experiment from Banas et al. 2016 in Coltrane 2.0
%
% Neil B, Oct 2022


% create forcing year from idealised curves -------------------------------------------

forcing.t = (0:365)'; % start with one year
forcing.T0 = nan.*forcing.t;
forcing.Td = nan.*forcing.t;
forcing.P = nan.*forcing.t;

% surface temperature -------------------------
Tjan1=-0.7; % T on Jan 1
Tmin=-1.8; % winter min T
Tmax=3.7; % summer max T
tTmin=105; % time when spring T increase starts
tTmax=250; % time of T max
nn = find(forcing.t >=0 & forcing.t <= tTmin);
forcing.T0(nn) = linspace(Tjan1,Tmin,length(nn));
nn = find(forcing.t >=tTmin & forcing.t <= tTmax);
forcing.T0(nn) = linspace(Tmin,Tmax,length(nn));
nn = find(forcing.t >=tTmax & forcing.t <= 365);
forcing.T0(nn) = linspace(Tmax,Tjan1,length(nn));

% deep temperature ----------------------------
Td0=1; % deep T year-round
forcing.Td = repmat(Td0, size(forcing.t));

% chlorophyll ---------------------------------
P0win=0.1;
P0spr=13;
P0sum=0.05*P0spr;
P0aut=5;
dtPspr=15; % bloom duration
dtPaut=30;
tPspr=150;
tPaut=225;
forcing.P = P0win .* ones(size(forcing.t));
forcing.P = max(forcing.P, ...
				 P0spr .* exp(-((forcing.t-tPspr)./dtPspr).^2));
forcing.P = max(forcing.P, ...
				 P0aut .* exp(-((forcing.t-tPaut)./dtPaut).^2));
fsum = forcing.t > tPspr & forcing.t < tPaut;
forcing.P(fsum) = max(forcing.P(fsum),P0sum);

% repeat this cycle for 7 years to resolve 3 year life cycles ----------------
Nyears = 7;
ind = [repmat((1:365)',[Nyears 1]); 1];
forcing.t = (0:(Nyears*365))';
forcing.T0 = forcing.T0(ind);
forcing.Td = forcing.Td(ind);
forcing.P = forcing.P(ind);


% main coltrane experiment ----------------------------------------------------------

dt = 20; % 15 is better, 30 is quicker
p00 = coltraneParams('Nyears',3,...
				   'requireActiveSpawning',0,...
				   'tdia_exit',[80 : dt : 120],...
				   'tdia_enter',[260 : dt : 300],...
				   'min_genlength_years',1,...
				   'max_genlength_years',3,...
				   'dt_spawn',15,...
				   'preySatVersion','default',...
				   'I0',0.36,...
				   'Ks',1,...
				   'Ddia',0.5*(stage2D('C2')+stage2D('C3')),... % start of C3
				   'allowGainAfterD1',0);
traits.u0 = 0.005 : 0.002 : 0.009;
coltraneCommunity('diskoEx',forcing,p00,traits);


