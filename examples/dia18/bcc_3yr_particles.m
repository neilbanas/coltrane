function P = bcc_3yr_particles(map1,map2,map3);

% assemble 3-yr particle tracks from 3 individual years of flow

reps = 30; % much smaller is almost as good
x00 = 120 : 0.25 : 240;
y00 = 53 : 0.125 : 75;
[x00,y00] = meshgrid(x00,y00);
f = map1.siMask(x00(:),y00(:)) > 0.5;

% note: assumes that map1,map2,map3 are all on the same grid with the same
% timestep

% first year
P = map1.integrate(repmat(x00(f),[reps 1]),repmat(y00(f),[reps 1]));

% second year: concatenate
P2 = map2.integrate(P.x(end,:)',P.y(end,:)');
fields = fieldnames(P);
for i=1:length(fields)
	P.(fields{i}) = [P.(fields{i}); P2.(fields{i})(2:end,:)];
end

% third year: concatenate
P2 = map3.integrate(P.x(end,:)',P.y(end,:)');
fields = fieldnames(P);
for i=1:length(fields)
	P.(fields{i}) = [P.(fields{i}); P2.(fields{i})(2:end,:)];
end

% fill in some other useful fields
t1 = map1.t(1) + (map1.t(2)-map1.t(1)) .* (0:size(P.x,1)-1);
P.t = repmat(t1(:),[1 size(P.x,2)]);
P.H = griddata(map1.x,map1.y,map1.H,P.x,P.y);

% reduce to a useful subset of trajectories
P.mask = map1.siMask(P.x,P.y);
f = find(all(P.mask==1) & ...
		 any(P.x>180 & P.y<210 & P.y>55 & P.y<75 & ...
			 P.y-P.x>60-210 & ...
			 P.H<200));
fields = fieldnames(P);
for i=1:length(fields)
	P.(fields{i}) = P.(fields{i})(:,f);
end

% winnow trajectories that occupy the same cell halfway through year 2.
% this makes the sampling more like a set that fills space at the midpoint
% and was tracked backwards and forwards
n = round(1.5 * length(map1.t));
xy = P.x(n,:) + sqrt(-1) .* P.y(n,:);
xyu = unique(xy);
keep = zeros(1,size(P.x,2));
for i=1:length(xyu)
	f = find(xy==xyu(i));
	keep(f(randi(length(f),1))) = 1;
end
f = find(keep);
for i=1:length(fields)
	P.(fields{i}) = P.(fields{i})(:,f);
end



