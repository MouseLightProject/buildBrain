function [A,subs] = filterEdges(A,subs,junk,params)
if isempty(junk)
    s = 16;
    for cr = 1:s:100;
        clustrange = cr:cr+s-1;
        vizClustersforJunks(A,subs,params,clustrange)
    end
    close all
end
%%
73336.7, 14360.0, 34718.8
if ~isempty(junk)
    [A,subs] = maskEdges(A,subs,junk);
end
