function map = bcc_returnMap(year);

% map = bcc_returnMap(year);
% 
% return maps for the Bering-Chukchi region in one year of BIOMAS output.
% uses the "particulator" package.

% load a biomas model run
modelDir = '/Volumes/sambal/biomas/1979-2016/';
gridDir = '/Volumes/sambal/biomas/info/';
depthRange = [-200 0];
run = modelRun_biomas2d(modelDir,year,gridDir,depthRange);

x00 = 120 : 0.25 : 240;
y00 = 53 : 0.125 : 75;
t0 = run.t(1) + (0:10:360);
KH = 200; % horizontal diffusivity

% make a blank return map in the Bering-Chukchi region
[x0,y0] = meshgrid(x00,y00);
m = run.interpMask(x0,y0) > 0.5;
x0 = x0(m);
y0 = y0(m);
map = returnMap(x0,y0,t0);

% fill in the return map using the model run
map.generate(run,'KH',0,'Nreplicates',1,...
			 'Ninternal',4,'verbose',1,'diffusive',0,...
			 'tracers',{'ice','temp','diatom','flagel'});
map.addDiffusion(KH,30);
