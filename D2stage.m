function stage = D2stage(D);

% stage = D2stage(D);
%
% convert D to stage, based on schedule in Campbell et al. 2001
% (Maps et al. 2014 is similar)

Bele_a = [0 595 983 1564 2951 3710 4426 5267 6233 7370 8798 10964 15047];  
Dstage = Bele_a./Bele_a(end); % age of entry into each stage
stage = nan.*D;
for i=1:length(Dstage)-1
	stage(D >= Dstage(i) & D < Dstage(i+1)) = i;
end
stage(D>=1) = length(Dstage);
