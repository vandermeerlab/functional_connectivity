function [OUT] = GetSSN(scope)
%Simple function to grab the Rat SSN:  R###-YYYY-MM-DD
%2011-05-10 AndyP Added scope.

if nargin==0
	D = cd;
	strRat = regexp(D, 'R\d\d\d-\d\d\d\d-\d\d-\d\d', 'once');
	if isempty(strRat);
		scope = 'MultipleSessions';
	else
		scope = 'SingleSession';
	end	
end

switch scope
	case 'MultipleSessions'
		
		D=dir;
		skip = 0;
		for iC=1:length(D);
			IsSession = regexp(D(iC,1).name, 'R\d\d\d-\d\d\d\d-\d\d-\d\d');
			[~,~,ext] = fileparts(D(iC,1).name);
			if strcmp(ext, '')~=1;
				IsSession = 0;
			end
			if IsSession ==1;
				SessioniD = D(iC,1).name;
				iD = iC-skip;
				Session{iD,1}= SessioniD; %#ok<AGROW>
			else
				skip = skip +1;
			end
		end
		try
			OUT = Session;
		catch %#ok<CTCH>
			error('No Sessions found in this directory.');
		end
	case 'SingleSession'
		
		if nargin~=0
			D = cd;
			strRat = regexp(D, 'R\d\d\d-\d\d\d\d-\d\d-\d\d');
		end
		SSNstr = D(strRat:end);
		if ~isempty(SSNstr)
			OUT = SSNstr;
		else
			error('Not in the correct directory to run SingleSession analysis.');
		end
end





end

