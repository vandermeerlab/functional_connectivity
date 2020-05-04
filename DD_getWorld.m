function [DelayZone, DelayUpStep, DelayDnStep] = DD_getWorld(World)
%2011-05-07 AndyP Simple function to get the delayed zone from any DD.m
%file.
% 2012-07-25 AndyP Renamed from DDgetWorld
% calculate the world


if exist('World', 'var');
	
	if (World.incrLgoL > World.incrLgoR) % L-delay world
		DelayZone = 3;
		DelayDnStep = World.incrLgoR;
		DelayUpStep = World.incrLgoL;
	elseif (World.incrRgoR > World.incrRgoL) % R-delay world
		DelayZone = 4;
		DelayDnStep = World.incrRgoL;
		DelayUpStep = World.incrRgoR;
	else
		error('Unknown world.');
	end
	
else
	error('Unknown world.');
end

end

