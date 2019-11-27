function fl = fitnessLandscape(dF1,t,t0,s,ng)

% fl = fitnessLandscape(dF1,t,t0,s,ng)
%
% Function always returns 1- and 2-generation fitness. If number of
% generations ng>1 then the function returns all fitnesses up to the
% (ng+1)-generation. This lets us know when the offspring of the
% ng-generation can spawn.

fl.dF1 = dF1;
fl.t = t;

[~,NC,NS] = size(dF1);
	% NC should match length of t0
	% NS should match length of fields of s

fields = fieldnames(s);
sz = size(s.(fields{1})); % prod(sz) = NS
fl.t0 = repmat(t0(:),[1 sz]);
for k=1:length(fields)
	fl.(fields{k}) = repmat(reshape(s.(fields{k}),[1 sz]),[NC 1]);
end


fl.F1 = sum(dF1); % lifetime egg production for each (t0,s)
fl.F1expected = max(fl.F1,[],3); % expected LEP for each t0, assuming that the 
						   % offspring will take the optimal strategy
F1ex_ = interp1(t0,fl.F1expected,t); % interpolate over all possible t0
F1ex_(isnan(F1ex_)) = 0;
F1ex_ = repmat(F1ex_(:),[1 NC NS]);

fl.dF2 = dF1 .* F1ex_; % contribution to two-generation fitness at each (t,t0,s)
fl.F2 = sum(fl.dF2); % two-generation fitness at each (t0,s)
% fl.F2 = sqrt(sum(fl.dF2)); % two-generation fitness at each (t0,s)
fl.F2expected = max(fl.F2,[],3); % at each t0, assuming that the offspring's
						   % offspring take the optimal strategy
                           
                           % could iterate like this...
                           
                           
% Aidan's additions for dealing with multiple generations per year
if ng > 1
    for i=2:ng
        j=i+1;
        dFj = ['dF' num2str(j)];
        Fj = ['F' num2str(j)];
        
        Fex_ = interp1(t0,fl.(['F' num2str(i) 'expected']),t);
        Fex_(isnan(Fex_)) = 0;
        Fex_ = repmat(Fex_(:),[1 NC NS]);
        
        fl.(dFj) = fl.dF1 .* Fex_;
        fl.(Fj) = sum(fl.(dFj));
%        fl.(Fj) = sum(fl.(dFj)) .^ (1/(j));
        fl.(['F' num2str(j) 'expected']) = max(fl.(Fj),[],3);
    end
end
                           
for i=2:(ng+1)
    fl.(['F' num2str(i)]) = fl.(['F' num2str(i)]) .^ (1/i);
end
    
    
    
