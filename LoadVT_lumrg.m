function [A,B,C] = LoadVT_lumrg(fn)

% [TS] = LoadVT(fn) : ts-obj
% [X,Y] = LoadVT(fn) : 2 tsd-obj
% [X,Y,PHI] = LoadVT(fn): 3 tsd-obj
%
%
% X,Y returned in camera coordinates
% PHI returned in radians
% NOTE: timestamps returned are in seconds!
%
% Status: PROMOTED
% ADR: 2001
% Version v1.0

A = []; B = []; C = [];

switch nargout
case 1 % TS only
	TS = LoadVT0_lumrg(fn);
	A = ts(TS'/1e6);
case 2 % X,Y only
	[TS,X,Y] = LoadVT0_lumrg(fn);
	A = tsd(TS'/1e6,X');
	B = tsd(TS'/1e6,Y');
case 3 % X,Y,PHI
	[TS,X,Y,PHI] = LoadVT0_lumrg(fn);
	A = tsd(TS'/1e6,X');
	B = tsd(TS'/1e6,Y');
	C = tsd(TS'/1e6,PHI');
otherwise
	error('Invalid function call.');
end
