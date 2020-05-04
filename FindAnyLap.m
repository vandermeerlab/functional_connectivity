function [Lstart, Lend] = FindAnyLap(time,x, varargin)
% 2011-04-13 AndyP  This function finds laps on the DD task based on
% restrictions in the VT data.  Function is based on a similar version
% written for the MT task.
% x input are the x coordinates of VT2
% Lstart, Lend are the output 

splitTime = 2; % sec.  Gives the minimum time between the last point on Lap X and the first point on lap X+1
RemoveNans = true;
process_varargin(varargin);

if RemoveNans
	NaNsCut = isnan(x);
	if any(NaNsCut)
		%disp('Removing NaNs from VT2');
		time(isnan(x)) = [];
		x(isnan(x))=[]; %#ok<NASGU>
	end
end

d = diff(time); %takes the difference between consecutive datapoints
L = find(d > splitTime);
L = [time(1); time(L); time(end)]; %#ok<FNDSB>

nLaps = length(L)-1;
for iL = 1:nLaps
	xr0 = time(time > L(iL));
	Lstart(iL) = xr0(1);
	[~,ix] = max(time(time <= L(iL+1)));
	Lend(iL) = time(ix);
end
