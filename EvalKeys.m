function EvalKeys(fd)

% Finds a keys file in directory fd (default '.')
% and evaluates it in the caller's workspace
% errors if no keys file found

if nargin > 0
	pushdir(fd);
end

clear ExpKeys
[fd fn fe] = fileparts(FindFile('*keys.m'));
evalin('caller', fn);

if nargin > 0
	popdir;
end
