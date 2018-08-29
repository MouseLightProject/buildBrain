function [Gfilt,subsfilt] = pruneGraph(Ain,subsin,sizeThr,opt)
%PRUNEGRAPH Summary of this function goes here
%
% [OUTPUTARGS] = PRUNEGRAPH(INPUTARGS) Explain usage here
%
% Inputs:
%
% Outputs:
%
% Examples:
%
% Provide sample usage code here
%
% See also: List related files here

% $Author: base $	$Date: 2017/03/29 16:19:40 $	$Revision: 0.1 $
% Copyright: HHMI 2017
% distThr = 10;
% sizeThr = 10;
Gin = graph(max(Ain,Ain'));
if 1
    if 1
        tic
        filt{1}{1} = 'graphfuncs.filtGsize';
        filt{1}{2} = {2*sizeThr};
        [Gfilt,subsfilt] = filterGraph(Gin,subsin,filt);
        %     [Gfilt,subsfilt] = filterGraph(Gin,subsin,{'filtGsize'},{sizeThr});
        sprintf('FILTER DONE in %d sec',round(toc))
    else
        tic
        filt{1}{1} = 'graphfuncs.filtSsize';
        filt{1}{2} = {sizeThr,opt};
        [Gfilt,subsfilt] = filterGraph(Gin,subsin,filt);
        sprintf('FILTER DONE in %d sec',round(toc))
    end
    
    tic
    [Gfilt,subsfilt] = graphfuncs.deleteLoops(Gfilt,subsfilt);
    sprintf('del loop DONE in %d sec',round(toc))
end
tic
iter=0;
G = Gfilt;
subs=subsfilt;
pruned = 1;
while pruned
    iter=iter+1;
    % [iter length(leafnodes)]
    [G,subs,pruned] = graphfuncs.pruneleafbranches(G,subs,sizeThr);
    [G,subs,totdel] = graphfuncs.deleteLoops(G,subs);
end
[iter pruned]
sprintf('PRUNE DONE in %d sec in %d iter',round(toc),iter)
Gfilt = G;
subsfilt = subs;

end
