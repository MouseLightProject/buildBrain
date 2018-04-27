function hits = maskWithInputLocations(subs,querrylocs)
% search eps-neighborhood of inputlocs in subs
subs_in_range = find(all(subs>=min(querrylocs) & subs<=max(querrylocs),2));
subs_ = subs(subs_in_range,:);
Idx = rangesearch(subs_,querrylocs,15);
Idx_ = unique(cat(2,Idx{:}));
hits = subs_in_range(Idx_);
