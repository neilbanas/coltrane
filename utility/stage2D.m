function D = stage2D(stage);

% D = stage2D(stage);
%
% quick lookup for middle of each developmental stage as a fraction of total development,
% based on schedule in Campbell et al. 2001.

stages = {'E','N1','N2','N3','N4','N5','N6','C1','C2','C3','C4','C5'};
Bele_a = [0 595 983 1564 2951 3710 4426 5267 6233 7370 8798 10964 15047];  
n = strmatch(stage,stages);
D = 0.5.*(Bele_a(n+1) + Bele_a(n)) ./ Bele_a(end);
