function out = ...
	coltraneEnsemble_rerunOne(forcingCase,forcingOptions,paramOptions,n);

% out = coltraneEnsemble_rerunOne(forcingCase,forcingOptions,paramOptions,n);
%
% reruns one case out of an ensemble to give full output rather than just
% summary stats. This could be done with coltraneModel.m, but this wrapper
% function reduces the need to think too hard about the inputs.
%
% n can be a linear index or a vector of subscripts.
%
% e.g. if you previous ran
% [cases,pot,pop] = coltraneEnsemble(forcingCase,forcingOptions,paramOptions);
% and you decide that you particularly like case (i,j,k) out of the output from
% that ensemble, you can call
% out = coltraneEnsemble(forcingCase,forcingOptions,paramOptions,[i j k]);
% to reproduce what you would have gotten had you run coltraneModel for that
% case by itself.


% generate all combos of params as in coltraneEnsemble
NF = length(forcingOptions)/2;
NP = length(paramOptions)/2;
outputList = '[';
for i=1:NF
	N(i) = length(forcingOptions{i*2});
	outputList = [outputList 'cases.' forcingOptions{i*2-1} ','];
end
for i=1:NP
	N(NF+i) = length(paramOptions{i*2});
	outputList = [outputList 'cases.' paramOptions{i*2-1} ','];
end
outputList(end) = ']';
eval([outputList ' = ndgrid(forcingOptions{2:2:end},paramOptions{2:2:end});']);
N1 = [N(N>1) 1 1]; % with singleton dimensions removed
Ntot = prod(N);

% make n into a linear index if it isn't already
if length(n) > 1
	fields = fieldnames(cases);
	sz = size(cases.(fields{1}));
	n = sub2ind(sz,n);
end

% set up forcing for case n
opt = forcingOptions;
for i=1:NF
	opt{i*2} = cases.(forcingOptions{i*2-1})(n);
end
forcing0 = coltraneForcing(forcingCase,opt{:});

% set up params for case n
opt = paramOptions;
for i=1:NP
	opt{i*2} = cases.(paramOptions{i*2-1})(n);
end
p = coltraneParams(opt{:});

% run it
out = coltraneModel(forcing0,p);


