% filename = 'wherever/you/want/to/keep/full/model/output/';
filename = ''; % if blank, saves summary stats only


Nyears = 5;
paramCases = {'Nyears',Nyears,...
			  'tegg',0:5:365*4,...
			  'u0',0.004:0.0005:0.01,...
			  'Ks',1,...
			  'm0',0.06};
forcingCases = {'dT',0};


[cases,pot,pop,p1,forcing1] = ...
	coltraneEnsemble('disko',forcingCases,paramCases,filename);
	
save runs_disko_stats pot pop cases p1 forcing1
