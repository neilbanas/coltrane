function fname = num2filename(basedir, n);

% fname = num2filename(basedir, n);
%
% convention for converting run numbers to filenames in a nested structure that keeps
% directories from filling up with tens of thousands of files.
%
% num2filename('my/dir/',345) gives
% my/dir/0/345.mat
% num2filename('my/dir/',12345) gives
% my/dir/12/345.mat

nstr = num2str(n);
if length(nstr)<4
	nstr = ['0000' nstr];
	nstr = nstr(end-3:end);
end
fname = [basedir nstr];
fname = [fname(1:end-3) '/' fname(end-2:end)];

  