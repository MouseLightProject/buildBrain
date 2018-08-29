function hits = maskWithInputLocations(subs,querrylocs,distThr)
if nargin<3
    distThr = 5;
end
% search eps-neighborhood of inputlocs in subs
subs_in_range = find(all(subs>=min(querrylocs) & subs<=max(querrylocs),2));
subs_ = subs(subs_in_range,:);
Idx = rangesearch(subs_,querrylocs,distThr);
Idx_ = unique(cat(2,Idx{:}));
hits = subs_in_range(Idx_);
