function [ngen, whichGen] = gensPerYear(fl)

% [ngen, whichGen] = gensPerYear(fl)
%
% Returns information on possible number of generations per year. Each
% valid timing strategy must form part of a annual cycle. For each timing
% strategy, function returns total number of generations in the annual
% cycle and which of those generations the strategy belongs to.
% Only applicable to 1 or more gens/year, will not work with fractional
% gens/year.

[NT,NC,NS] = size(fl.dF1);

dt = fl.t(2)-fl.t(1);
t = repmat(fl.t,[1 NC NS]);
t0 = reshape(fl.t0, [1, NC, NS]);

ngen = zeros(1,NC,NS); % each timing strategy falls within a cycle of ngen generations/year
whichGen = nan(1,NC,NS); % which generation in cycle each timing strategy belongs to

S1 = fl.dF1 > 0; % which cohorts can spawn when

t01 = t0 < fl.t(round(365/dt)); % strategies hatched in 1st year

x = diff(S1); % indexes start of spawning for each strategy
x(x==-1) = 0;
x = [zeros(1,NC,NS); x]; % recover correct dimension lost after differencing
tx = t .* x; % times at start of spawning

ta1 = any(tx > 0 & tx < fl.t(round(365/dt))); % strategies spawning in 1st year
ta2 = any(tx > fl.t(round(365/dt)) & tx < fl.t(round(730/dt))); % strategies spawning in 2nd year

% fv lists [NC NS] linear indices of all valid strategies given various
% gen/year cycles, stored in nested structure of total gens/year then
% specific generation
fv.ngen1.gen1 = find(fl.F2 > fl.p.fitnessFilter & t01 & ta2);
minGen = nan; maxGen = nan; % note potential min and max gen/year cycles
if ~isempty(fv.ngen1.gen1)
    minGen = 1;
    maxGen = 1;
end
for jj=2:fl.p.maxGenPerYear % loop over possible gens/year
    lab = ['ngen' num2str(jj)];
    for kk=1:jj % loop over individual generations
        lab2 = ['gen' num2str(kk)];
        if kk==1
            fv.(lab).(lab2) = ...
                find(fl.(['F' num2str(jj-kk+2)]) > fl.p.fitnessFilter & t01 & ta1);
        elseif kk~=jj
            fv.(lab).(lab2) = ...
                find(fl.(['F' num2str(jj-kk+2)]) > fl.p.fitnessFilter & t01 & ta1 & t0 >= tamin);
        else
            fv.(lab).(lab2) = ...
                find(fl.(['F' num2str(jj-kk+2)]) > fl.p.fitnessFilter & t01 & ta2 & t0 >= tamin);
        end                        
        tamin = tx(:,fv.(lab).(lab2)); % earliest spawning date of gen i is earliest possible hatch date of gen i+1,
        tamin = min(tamin(tamin>0));   % used to reduce size of internal nested loops
        if kk==1
            possible = ~isempty(fv.(lab).(lab2));
        else
            possible = possible & ~isempty(fv.(lab).(lab2));
        end
    end
    if possible
        maxGen = jj;
        if isnan(minGen)
            minGen = jj;
        end
    end
end

Gens = minGen:maxGen; % potential numbers of generations per year

if ~isnan(maxGen)    
    for j=1:length(Gens) % loop through gen/year possibilities
        jj=Gens(j);
        lab = ['ngen' num2str(jj)];
        for kk=1:jj % loop through specific generations
            lab2 = ['gen' num2str(kk)];
            t0v.(lab).(lab2) = t0(fv.(lab).(lab2))'; % all valid hatch dates
            [I,~] = find(t(:,fv.(lab).(lab2)) == repmat(t0v.(lab).(lab2),[NT 1]));
            t0v_ind.(lab).(lab2) = I; % timing index of valid hatch dates
            t0v_indu.(lab).(lab2) = unique(I); % unique hatch dates - define loop index sizes
            % S stores valid spawning times given total gens/year and
            % specific gen within the annual cycle - an i gen/year cycle
            % requires the 1st generation to have positive F(i+1) fitness
            % and the ith gen to have positive F2 fitness
            S.(lab).(['S' num2str(jj-kk+2)]) = S1(:,fv.(lab).(lab2));
        end
    end
    % call genCount function for each possible gen/year cycle
    for ijk=1:length(Gens)
        kk=Gens(ijk);        
        count=genCount(kk,t0v_ind,t0v_indu,fv,S,dt,fl.dF1);        
        validCycle = count.ngen > 0;
        ngen(validCycle) = count.ngen(validCycle);
        whichGen(validCycle) = count.whichGen(validCycle);
    end    
end


function out = genCount(ng,t0,t0u,fvalid,S,dt,dF1)

[~,NC,NS] = size(dF1);

out.ngen = zeros(1,NC,NS); % generations per year given timing strategy
out.whichGen = nan(1,NC,NS); % which generation in cycle each timing strategy belongs to

% depth of nested loops depends on number of generations ng

if ng == 1 % 1 gen/year cycles
    % S2 = spawning times of cohorts with positive F2 fitness
    for i=1:length(t0u.ngen1.gen1)
        t0i = t0u.ngen1.gen1(i);
        f1 = t0.ngen1.gen1 == t0i & S.ngen1.S2(t0i+365/dt,:)'; % valid strategies hatched on t0i and spawning 1 year later
        out.ngen(1,fvalid.ngen1.gen1(f1)) = ng;
        out.whichGen(1,fvalid.ngen1.gen1(f1)) = 1;
    end
    
elseif ng == 2 % 2 gen/year cycles
    % S3 = spawning times of cohorts with positive F3 fitness
    for i=1:length(t0u.ngen2.gen1)
        t0i = t0u.ngen2.gen1(i);
        jpos = find(t0u.ngen2.gen2 > t0i);
        for j=min(jpos):max(jpos)
            t0j = t0u.ngen2.gen2(j);
            f1 = t0.ngen2.gen1 == t0i & S.ngen2.S3(t0j,:)'; % hatched at t0i, spawning at t0j
            f2 = t0.ngen2.gen2 == t0j & S.ngen2.S2(t0i+365/dt,:)'; % hatched at t0j, spawning at t0i +365
            if any(f1) && any(f2) % do combined hatch/spawning dates complete a cycle?
                out.ngen(1,[fvalid.ngen2.gen1(f1); fvalid.ngen2.gen2(f2)]) = ng;
                out.whichGen(1,fvalid.ngen2.gen1(f1)) = 1;
                out.whichGen(1,fvalid.ngen2.gen2(f2)) = 2;
            end
        end
    end
    
elseif ng == 3 % 3 gen/year cycles
    % S4 = spawning times of cohorts with positive F4 fitness
    for i=1:length(t0u.ngen3.gen1)
        t0i = t0u.ngen3.gen1(i);
        jpos = find(t0u.ngen3.gen2 > t0i);
        for j=min(jpos):max(jpos)
            t0j = t0u.ngen3.gen2(j);
            kpos = find(t0u.ngen3.gen3 > t0j);
            for k=min(kpos):max(kpos)
                t0k = t0u.ngen3.gen3(k);
                f1 = t0.ngen3.gen1 == t0i & S.ngen3.S4(t0j,:)'; % hatched at t0i, spawning at t0j
                f2 = t0.ngen3.gen2 == t0j & S.ngen3.S3(t0k,:)'; % hatched at t0j, spawning at t0k
                f3 = t0.ngen3.gen3 == t0k & S.ngen3.S2(t0i+365/dt,:)'; % hatched at t0k, spawning at t0i +365
                if any(f1) && any(f2) && any(f3) % do combined hatch/spawning dates complete a cycle?
                    out.ngen(1,[fvalid.ngen3.gen1(f1); fvalid.ngen3.gen2(f2); fvalid.ngen3.gen3(f3)]) = ng;
                    out.whichGen(1,fvalid.ngen3.gen1(f1)) = 1;
                    out.whichGen(1,fvalid.ngen3.gen2(f2)) = 2;
                    out.whichGen(1,fvalid.ngen3.gen3(f3)) = 3;
                end
            end
        end
    end
    
elseif ng == 4 % 4 gen/year cycles
    % S5 = spawning times of cohorts with positive F5 fitness
    for i=1:length(t0u.ngen4.gen1)
        t0i = t0u.ngen4.gen1(i);
        jpos = find(t0u.ngen4.gen2 > t0i);
        for j=min(jpos):max(jpos)
            t0j = t0u.ngen4.gen2(j);
            kpos = find(t0u.ngen4.gen3 > t0j);
            for k=min(kpos):max(kpos)
                t0k = t0u.ngen4.gen3(k);
                lpos = find(t0u.ngen4.gen4 > t0k);
                for l=min(lpos):max(lpos)
                    t0l = t0u.ngen4.gen4(l);
                    f1 = t0.ngen4.gen1 == t0i & S.ngen4.S5(t0j,:)'; % hatched at t0i, spawning at t0j
                    f2 = t0.ngen4.gen2 == t0j & S.ngen4.S4(t0k,:)'; % hatched at t0j, spawning at t0k
                    f3 = t0.ngen4.gen3 == t0k & S.ngen4.S3(t0l,:)'; % hatched at t0k, spawning at t0l
                    f4 = t0.ngen4.gen4 == t0l & S.ngen4.S2(t0i+365/dt,:)'; % hatched at t0l, spawning at t0i +365
                    if any(f1) && any(f2) && any(f3) && any(f4) % do combined hatch/spawning dates complete a cycle?
                        out.ngen(1,[fvalid.ngen4.gen1(f1); fvalid.ngen4.gen2(f2); fvalid.ngen4.gen3(f3); fvalid.ngen4.gen4(f4)]) = ng;
                        out.whichGen(1,fvalid.ngen4.gen1(f1)) = 1;
                        out.whichGen(1,fvalid.ngen4.gen2(f2)) = 2;
                        out.whichGen(1,fvalid.ngen4.gen3(f3)) = 3;
                        out.whichGen(1,fvalid.ngen4.gen4(f4)) = 4;
                    end
                end
            end
        end
    end
    
elseif ng == 5 % 5 gen/year cycles    
    % S6 = spawning times of cohorts with positive F6 fitness    
    for i=1:length(t0u.ngen5.gen1)
        t0i = t0u.ngen5.gen1(i);
        jpos = find(t0u.ngen5.gen2 > t0i);
        for j=min(jpos):max(jpos)
            t0j = t0u.ngen5.gen2(j);
            kpos = find(t0u.ngen5.gen3 > t0j);
            for k=min(kpos):max(kpos)
                t0k = t0u.ngen5.gen3(k);
                lpos = find(t0u.ngen5.gen4 > t0k);
                for l=min(lpos):max(lpos)
                    t0l = t0u.ngen5.gen4(l);
                    mpos = find(t0u.ngen5.gen5 > t0l);
                    for m=min(mpos):max(mpos)
                        t0m = t0u.ngen5.gen5(m);
                        f1 = t0.ngen5.gen1 == t0i & S.ngen5.S6(t0j,:)'; % hatched at t0i, spawning at t0j
                        f2 = t0.ngen5.gen2 == t0j & S.ngen5.S5(t0k,:)'; % hatched at t0j, spawning at t0k
                        f3 = t0.ngen5.gen3 == t0k & S.ngen5.S4(t0l,:)'; % hatched at t0k, spawning at t0l
                        f4 = t0.ngen5.gen4 == t0l & S.ngen5.S3(t0m,:)'; % hatched at t0l, spawning at t0m
                        f5 = t0.ngen5.gen5 == t0m & S.ngen5.S2(t0i+365/dt,:)'; % hatched at t0m, spawning at t0i +365
                        if any(f1) && any(f2) && any(f3) && any(f4) && any(f5) % do combined hatch/spawning dates complete a cycle?
                            out.ngen(1,[fvalid.ngen5.gen1(f1); fvalid.ngen5.gen2(f2); fvalid.ngen5.gen3(f3); fvalid.ngen5.gen4(f4); fvalid.ngen5.gen5(f5)]) = ng;
                            out.whichGen(1,fvalid.ngen5.gen1(f1)) = 1;
                            out.whichGen(1,fvalid.ngen5.gen2(f2)) = 2;
                            out.whichGen(1,fvalid.ngen5.gen3(f3)) = 3;
                            out.whichGen(1,fvalid.ngen5.gen4(f4)) = 4;
                            out.whichGen(1,fvalid.ngen5.gen5(f5)) = 5;
                        end
                    end
                end
            end
        end
    end

    
elseif ng == 6 % 6 gen/year cycles
    % S7 spawning times of cohorts with positive F7 fitness
    for i=1:length(t0u.ngen6.gen1)
        t0i = t0u.ngen6.gen1(i);
        jpos = find(t0u.ngen6.gen2 > t0i);
        for j=min(jpos):max(jpos)
            t0j = t0u.ngen6.gen2(j);
            kpos = find(t0u.ngen6.gen3 > t0j);
            for k=min(kpos):max(kpos)
                t0k = t0u.ngen6.gen3(k);
                lpos = find(t0u.ngen6.gen4 > t0k);
                for l=min(lpos):max(lpos)
                    t0l = t0u.ngen6.gen4(l);
                    mpos = find(t0u.ngen6.gen5 > t0l);
                    for m=min(mpos):max(mpos)
                        t0m = t0u.ngen6.gen5(m);
                        npos = find(t0u.ngen6.gen6 > t0m);
                        for n=min(npos):max(npos)
                            t0n = t0u.ngen6.gen6(n);
                            f1 = t0.ngen6.gen1 == t0i & S.ngen6.S7(t0j,:)'; % hatched at t0i, spawning at t0j
                            f2 = t0.ngen6.gen2 == t0j & S.ngen6.S6(t0k,:)'; % hatched at t0j, spawning at t0k
                            f3 = t0.ngen6.gen3 == t0k & S.ngen6.S5(t0l,:)'; % hatched at t0k, spawning at t0l
                            f4 = t0.ngen6.gen4 == t0l & S.ngen6.S4(t0m,:)'; % hatched at t0l, spawning at t0m
                            f5 = t0.ngen6.gen5 == t0m & S.ngen6.S3(t0n,:)'; % hatched at t0m, spawning at t0n
                            f6 = t0.ngen6.gen6 == t0n & S.ngen6.S2(t0i+365/dt,:)'; % hatched at t0n, spawning at t0i +365
                            if any(f1) && any(f2) && any(f3) && any(f4) && any(f5) && any(f6) % do combined hatch/spawning dates complete a cycle?
                                out.ngen(1,[fvalid.ngen6.gen1(f1); fvalid.ngen6.gen2(f2); fvalid.ngen6.gen3(f3); fvalid.ngen6.gen4(f4); fvalid.ngen6.gen5(f5); fvalid.ngen6.gen6(f6)]) = ng;
                                out.whichGen(1,fvalid.ngen6.gen1(f1)) = 1;
                                out.whichGen(1,fvalid.ngen6.gen2(f2)) = 2;
                                out.whichGen(1,fvalid.ngen6.gen3(f3)) = 3;
                                out.whichGen(1,fvalid.ngen6.gen4(f4)) = 4;
                                out.whichGen(1,fvalid.ngen6.gen5(f5)) = 5;
                                out.whichGen(1,fvalid.ngen6.gen6(f6)) = 6;
                            end
                        end
                    end
                end
            end
        end
    end
end   

