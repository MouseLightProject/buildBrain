function [Gfilt,subsfilt] = filterGraph(Gin,subsin,filtarg);
%FILTERGRAPH wrapper for graph filtering functions
%
% [OUTPUTARGS] = FILTERGRAPH(INPUTARGS) Explain usage here
%
% Inputs:
%
% Outputs:
%
% Examples:
%     tic
%     filt{1}{1} = 'filtGsize';
%     filt{1}{2} = {sizeThr};
%     [Gfilt,subsfilt] = filterGraph(Gin,subsin,filt);
% runs size based filtering for graphs 
%
% See also: List related files here

% $Author: base $	$Date: 2017/03/29 17:29:31 $	$Revision: 0.1 $
% Copyright: HHMI 2017

%% filter based on size
for ifilt = 1:length(filtarg)
    [Gfilt,subsfilt] = feval(filtarg{ifilt}{1},Gin,subsin,filtarg{ifilt}{2});
end



end
