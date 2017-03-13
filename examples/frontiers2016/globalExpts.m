% Banas et al., Front. Mar. Res., 2016
% Section 3.2


% filename = 'wherever/you/want/to/keep/full/model/output/';
filename = ''; % if blank, saves summary stats only

fullResolution = 0; % set to 0 if you want something that runs in just a minute
					% or two, 1 to actually replicate results in the paper 

if fullResolution
	forcingCases = {'dtP',5:5:150,...
					'T0',-2:16,...
					'Tb0_over_T0',[0.4 1]};
	paramCases = {'Nyears',4,...
				  'tegg',0,...
				  'u0',0.005:0.001:0.01,...
				  'Ks',1,...
				  'm0',0.08};
else
	forcingCases = {'dtP',20:20:150,...
					'T0',-2:4:16,...
					'Tb0_over_T0',[0.4 1]};
	paramCases = {'Nyears',4,...
				  'tegg',0,...
				  'u0',0.005:0.001:0.01,...
				  'Ks',1,...
				  'm0',0.08};
end


[cases,pot,pop,p1,forcing1] = ...
	coltraneEnsemble('simple-g',forcingCases,paramCases,filename);

save output/runs_global_stats pot cases p1 forcing1
% note: since we didn't vary tegg or choose it carefully, the results from the
% ER model in _pop_ aren't considered meaningful and aren't saved


% quick sketch of Fig. 6 in the paper--but note the difference between the
% input parameter dtP (which is what we have quick access to here) and the 
% analysis parameter dt (see Eq. 38: if you want dt, recreate the forcing for 
% each model case and calculate it from the formula)
figure
dtP = unique(cases.dtP);
T0 = unique(cases.T0);
glacialis = find(unique(cases.u0)==0.008);
contourf(dtP, T0, pot.ngenpot(:,:,1,glacialis)');
