function frag_assignment = assign_frags_to_GT(branches,subs,GT)
%%
n_branches = length(branches);
n_sources = length(GT.swcpixlocs);

%%
if 1
    % create an accumulator for GT
    [swc_locs, labels] = deal(cell(1,n_sources));
    parfor is = 1:n_sources
        swc_locs{is} = GT.swcpixlocs{is};
        labels{is} = is*ones(1,size(swc_locs{is},1));
    end
    swc_locs = cat(1,swc_locs{:});
    labels = cat(2,labels{:});
    [initass,initR] = knnsearch(swc_locs,subs,'k',1);
    % 5331
    %%
    frag_assignment = cell(1,n_sources);
    for ibr = 1:n_branches%[5719        5747        5748]%
        %%
        assign = false;
        assign_from_maxima = false;
        labs = labels(initass(branches(ibr).inds,:));
        hits = double(initR(branches(ibr).inds,:)<5);
        localmaximas = find(imregionalmax(bwdist(~hits)).*hits);
        if labs(1)==labs(end) & median(initR(branches(ibr).inds,:))<15
            assign = true;
        elseif median(initR(branches(ibr).inds,:))<5
            assign = true;
        elseif length(localmaximas)>1 & all(localmaximas~=[1 length(hits)])
            assign = false;
            assign_from_maxima = true;
        end
        
        if assign
            lab = mode(labels(initass(branches(ibr).inds,:)));
            if round(lab)==lab
                frag_assignment{lab}(end+1) = ibr;
            end
        elseif assign_from_maxima
            lab = labels(initass(branches(ibr).inds(localmaximas(1))));
            frag_assignment{lab}(end+1) = ibr;
        end
    end
else
    %%
    br_stat = zeros(n_branches,n_sources,5);
    try parfor_progress(0);catch;end
    parfor_progress(n_branches)
    parfor ibr = 1:n_branches
        parfor_progress()
        br_subs = branches(ibr).subs;
        for is = 1:n_sources
            pd = min(pdist2(br_subs,GT.swcpixlocs{is}),[],2);
            br_stat(ibr,is,:) = [min(pd) max(pd) median(pd) pd(1) pd(2)];
        end
    end
    parfor_progress(0)
    %%
end
% %%
% is = 2
% figure(32),cla
% hold on
% for ii = 1:length(frag_assignment{is})
%     myplot3(branches(frag_assignment{is}(ii)).subs,'-')
% end
% gplot3(GT.connMatrix{is},GT.swcpixlocs{is},'o')

%%
% %%
% valid_ones = find(initR<sqrt(12));
%
% %%
% % only keep nodes that have single id
% tr=[];
% for it = 1:length(initass)
%     if length(nodeBrid{initass(it)})==1
%         tr{end+1} = nodeBrid{initass(it)};
%     end
% end
% tr = unique([tr{:}]);
% tmp{is} = tr;















