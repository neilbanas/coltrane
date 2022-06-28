function coltraneEnsemble(outfileBasename,P,expts,p0);


% set up for saving output ------------
dirname = [outfileBasename '_output/'];
if ~exist(dirname,'dir'), mkdir(dirname); end
summaryFile = [dirname 'summary.mat'];
interimFile = [dirname 'summary-interim.mat'];
allFile = [dirname 'allCohortsAndStrategies.mat'];
% initialise structures to hold summaries
V1 = initialiseSummary(expts);
V2 = initialiseSummary(expts);
V1yr = initialiseSummary(expts);
V2yr = initialiseSummary(expts);


% -------------------------------------
exptFields = fieldnames(expts);
Nexpts = numel(expts.(exptFields{1}));
for i=1:Nexpts
	p = p0;
	for k=1:length(exptFields)
		p.(exptFields{k}) = expts.(exptFields{k})(i);
	end
	
	clear forcing
	forcing.t = P.t;
	forcing.yday = yearday(forcing.t);
	forcing.x = P.x(:,p.trajn);
	forcing.y = P.y(:,p.trajn);
	forcing.Ptot = P.chlmax(:,p.trajn);
	forcing.T0 = P.T0(:,p.trajn);
	forcing.Td = P.Td(:,p.trajn);
	forcing.ice = P.ice(:,p.trajn);

	preySat = preySaturation(forcing,p);
	forcing.satWC(:,i) = preySat.satWC;
	forcing.satIA(:,i) = preySat.satIA;
	forcing.sat(:,i) = preySat.sat;
	
	disp([num2str(i) ' / ' num2str(Nexpts)]);
	vi = coltraneModel(forcing,p,'scalars and fitness');
	vi = rmfield(vi,{'dF1','t'});
		% we save the full fitness landscape in coltrane_integrate.m so that 
		% coltraneModel.m can calculate the 2-generation fitness--
		% but after that it's unnecessary
	
	% A.* = all cohorts/strategies
%	fields = fieldnames(vi); % could make this a reduced list of essentials only
	fields = {'F1','F2','Wa','capfrac','tEcen','D_winter','t0','level'};
	for k=1:length(fields)
		vifk = squeeze(vi.(fields{k}));
		vifk = reshape(vifk,[1 size(vifk)]);
		if i==1
			A.(fields{k}) = repmat(nan.*vifk,[Nexpts 1 1]);
		end
		A.(fields{k})(i,:,:) = vifk;
	end

	% V1.* = optimal strategy/cohort, where optimal = max F1
	firstYear = vi.t0 < forcing.t(1) + 365;
	Fmax = max(vi.F1(firstYear));
	f = find(vi.F1==Fmax & firstYear);
	V1 = summarise(V1,i,vi,f,forcing);
	% where optimal = max F2
	Fmax = max(vi.F2(firstYear));
	f = find(vi.F2==Fmax & firstYear);
	V2 = summarise(V2,i,vi,f,forcing);
	% optimal (F1) for 1-yr, 2-yr generation lengths
	gl = (vi.tEcen - vi.t0)./365;
	Fmax = max(vi.F1(firstYear & round(gl)<=1));
	if ~isempty(Fmax)
		f = find(vi.F1==Fmax & firstYear & round(gl)<=1);
		V1yr = summarise(V1yr,i,vi,f,forcing);
	end
	Fmax = max(vi.F1(firstYear & round(gl)>=2));
	if ~isempty(Fmax)
		f = find(vi.F1==Fmax & firstYear & round(gl)>=2);
		V2yr = summarise(V2yr,i,vi,f,forcing);
	end
			
	if mod(i,100)==0
		save('-v7.3', interimFile, 'V*','expts','p0','P');
	end
end


% clean up V and save -----------------
V1 = cleanUpSummary(V1,expts);
V2 = cleanUpSummary(V2,expts);
V1yr = cleanUpSummary(V1yr,expts);
V2yr = cleanUpSummary(V2yr,expts);
save('-v7.3', summaryFile, 'V*','expts','p0','P');


% clean up A and save -----------------
try
fields = fieldnames(A);
	for k=1:length(fields)
		A.(fields{k}) = reshape(A.(fields{k}),[size(expts.(exptFields{1})) ... 	
											size(A.(fields{k}),2) size(A.(fields{k}),3)]);
	end
catch
	warning('reshaping problem with A.');
end
save('-v7.3',allFile, '-struct', 'A');





% ------
function V = initialiseSummary(expts)
exptFields = fieldnames(expts);
V.F1 = zeros(size(expts.(exptFields{1})));
V.F2 = V.F1;
V.t0 = nan .* V.F1;
suffixes = {'0','Ecen','Gain','Yield','YieldR','YieldC56'};
for k=1:length(suffixes)
	V.(['t' suffixes{k}]) = nan .* V.F1;
end

% ------
function V = summarise(V0,i,vi,f,forcing)
V = V0;
V.level(i) = max(vi.level(:));
if ~isempty(f) & V.level(i) >= 2
	fields = fieldnames(vi);
	for k=1:length(fields)
		V.(fields{k})(i) = vi.(fields{k})(f(1));
	end
end
suffixes = {'0','Ecen','Gain','Yield','YieldR','YieldC56'};
for k=1:length(suffixes)
	V.(['x' suffixes{k}])(i) = interp1(forcing.t,forcing.x,V.(['t' suffixes{k}])(i));
	V.(['y' suffixes{k}])(i) = interp1(forcing.t,forcing.y,V.(['t' suffixes{k}])(i));
end

% ------
function V = cleanUpSummary(V0,expts)
V = V0;
V.tEcen(V.tEcen==0) = nan;
V.gl = (V.tEcen-V.t0)./365;
exptFields = fieldnames(expts);
try
	fields = fieldnames(V);
	for k=1:length(fields)
		% pad the end of the fields of V if the last 1+ cases are terminated early
		dn = numel(expts.(exptFields{1})) - numel(V.(fields{k}));
		if dn>0
			V.(fields{k}) = cat(1,V.(fields{k})(:),nan.*ones(dn,1));
		end
		% make dimensions of V match the dimensions in the array of experiments
		V.(fields{k}) = reshape(V.(fields{k}),size(expts.(exptFields{1})));
	end
catch
	warning('reshaping problem with V.');
end
