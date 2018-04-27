function vizClustersforJunks(A,subs,params,clustrange)
[S,Comps] = graphconncomp(A,'DIRECTED',false);
Y = histcounts(Comps,1:S+1);
% filter out results based on manual mask - or allen registration or some other criteria
[~,mY] = max(Y);
[ia,ib]=sort(Y,'descend');
% close all
%%
figure12

for ic = clustrange
    ic_ids=Comps==ib(ic);
    subs_ = subs(ic_ids,:);
    A_ = A(ic_ids,:); A_ = A_(:,ic_ids);
    subplot(4,4,ic-clustrange(1)+1);
    gplot3(A_,subs_,'-');
    tit = sprintf('%d:\n [%d %d %d]',round([ic pix2um(params,median(subs_))]));
    title(tit);
    disp(tit)
end

