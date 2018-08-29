function workflow2(opt,tempfold,tag,G,subs,branches,nodeBrid)
%%
params=opt.params;
%% dump each branch as a seperate swc file
fold_path = '/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/2017-09-25/striatum/'
branch_swcoutfolder = fullfile(fold_path,'workflow2-updated-full')
mkdir(branch_swcoutfolder)
if exist('branch_swcoutfolder','var') & 1
    branch2swc(params,branches,branch_swcoutfolder,1)
end
%%
if 0
    workflow1(G,subs,opt)
    return
end

if isempty(tag)
    tag = [];%'051617'
end

%%
%%
% RUN graph simalirity functions
if 1
    querdist = 50;
    branch_pair_distance = graphfuncs.branchConn(branches,subs,nodeBrid,querdist); % pwdist of branches
    [aa,bb,cc] = find(branch_pair_distance);
    branch_pair_connectivity = sparse(aa,bb,cc<=querdist,size(branch_pair_distance,1),size(branch_pair_distance,1)); %pwcost of branches
    branch_pair_dissimilarity = graphfuncs.calcDists(branches,branch_pair_connectivity);
    %%
    save(fullfile(tempfold,['branchConn',tag]),'branch_pair_distance','branch_pair_dissimilarity')
    save(fullfile(tempfold,['branchConn']),'branch_pair_distance','branch_pair_dissimilarity')
else
    load(fullfile(tempfold,['branchConn',tag]),'branch_pair_distance','branch_pair_dissimilarity')
end
%%
% load concensus for validation if exists
if 1
    %%
    clear GT
    swcfold = '/nrs/mouselight/seggui/swcfiles/GT/2017-09-25_striatum_neurons_CA'
    swcfold = '/nrs/mouselight/seggui/swcfiles/GT/2017-09-25_striatum_neurons_temp'
    myfiles=dir(fullfile(swcfold,'*.swc'));
    swcgts={myfiles.name};
    swcgts = swcgts(~cellfun(@(x) contains(x,'mask'),swcgts));
    n_sources = length(swcgts);
    % sort based on file indicies <tag>_L/R_idx.swc
    fileinds = zeros(1,n_sources);
    for ii=1:n_sources
        str_ = strsplit(swcgts{ii},'_');
        str_ = strsplit(str_{end},'.');
        fileinds(ii) = str2double(str_{1});
    end
    [fileinds,sortedinds] = sort(fileinds);
    swcgts = swcgts(sortedinds);
    
    swcgts=cellfun(@(x) fullfile(swcfold,x),swcgts, 'UniformOutput',false);
    index = 1:n_sources;
    index = fileinds;
    [swcpixlocs,connMatrix,swcfilepath] = loadGT(swcgts,params,index);
    GT.swcpixlocs = swcpixlocs;
    GT.connMatrix = connMatrix;
    GT.swcfilepath = swcfilepath;
    %% make an initial assignment of branches wrto GT
    % for every branch get head-mid-tail, then get acc of GT locations,
    % make assignment of branches to GT
    if 1
        %%
        frag_assignment = assign_frags_to_GT(branches,subs,GT);
        GT.frag_assignment = frag_assignment;
        
        %%
        is = 25
        figure(32),cla
        hold on
        for ii = 1:length(frag_assignment{is})
            myplot3(branches(frag_assignment{is}(ii)).subs,{'-','LineWidth',6});
        end
%         myplot3(branches(5748).subs,{'-','LineWidth',6});
        gplot3(GT.connMatrix{is},GT.swcpixlocs{is},'-','MarkerSize',4);
        
    end
%%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    %% create init from GT:
    1
    init = GT;
    for is=1:length(swcgts)
        %%
        gt_A = max(GT.connMatrix{is},GT.connMatrix{is}');
        gt_subs = GT.swcpixlocs{is};
        [L,list] = getBranches(gt_A,1);
        % use the closest 3 branches for initialization
        %%
        parBr = zeros(1,length(L));
        for il = 1:length(L)
            parBr(il) = L(il).parentbranch;
        end
        theseBranches = find(parBr<3);
        initsubs_il = [];
        for il = theseBranches
            initsubs_il{il} = L(il).set;
        end
        initsubs_il = unique([1 cat(2,initsubs_il{:})]);
        init.swcpixlocs{is} = gt_subs(initsubs_il,:);
        init.connMatrix{is} = gt_A(initsubs_il,:);
        init.connMatrix{is} = init.connMatrix{is}(:,initsubs_il);
        %%
        % find the frags that are within initsub
        these_branches=[];
        for il = GT.frag_assignment{is}
            pd=min(pdist2(branches(il).subs,gt_subs(initsubs_il,:)),[],2);
            if median(pd(:))<5
                these_branches(end+1) = il;
            end
        end
        if isempty(these_branches)
            init.frag_assignment{is} = GT.frag_assignment{is}(randperm(length(GT.frag_assignment{is}),5));
        else
            init.frag_assignment{is} = these_branches;
        end
    end
    %%
    % make sure inits do not intersect
    clc
    for is1=1:length(swcgts)
        is_frag1 = init.frag_assignment{is1};
        for is2=is1+1:length(swcgts)
            is_frag2 = init.frag_assignment{is2};
            int=intersect(is_frag1,is_frag2);
            if int
                [is1 is2]
                warning('matching branches')
            else
                init.frag_assignment{is1} = ...
                    setdiff(init.frag_assignment{is1},int);
                init.frag_assignment{is2} = ...
                    setdiff(init.frag_assignment{is2},int);
            end
        end
    end
end

%%
if 0
    clear init
    % load dendrite/main axon for validation
    swcfold = '/groups/mousebrainmicro/mousebrainmicro/erhan_dm11/annotations_2015-06-19/GN_initialization/';
    myfiles=dir([swcfold,'*.swc']);
    swcinits={myfiles.name};
    swcinits=cellfun(@(x) fullfile(swcfold,x),swcinits, 'UniformOutput',false);
    index = [1 3 4 5 6 7 8 9 10 11 12 13];
    [swcpixlocs,connMatrix,swcfilepath] = loadGT(swcinits,params,index);
    init.swcpixlocs = swcpixlocs;
    init.connMatrix = connMatrix;
    init.swcfilepath = swcfilepath;
end
%%

iter = 1;
clear neuron
fnames = fieldnames(init)';
for idx=index(:)'
    % initialization
    for fn = fnames
        neuron(iter).init.(fn{1}) = init.(fn{1}){idx};
    end
    neuron(iter).init.somaloc = neuron(iter).init.swcpixlocs(1,:);
    if idx<=length(GT.(fnames{1})) & ~isempty(GT.(fnames{1}){idx})
        % check if GT is avaliable
        neuron(iter).GTtrue = 1;
        for fn = fnames
            neuron(iter).GT.(fn{1}) = GT.(fn{1}){idx};
        end
        neuron(iter).GT.somaloc = neuron(iter).GT.swcpixlocs(1,:);
    else
        neuron(iter).GTtrue = 0;
    end
    iter = iter+1;
end
%%
save(fullfile('/nrs/mouselight/seggui',['neuronaxon',tag]),'neuron')

%%
ran = max(subs);
branchcent=reshape([branches.cent]',3,[])';
tipsubs=[branches.tipsubs];tipsubs = permute(reshape(tipsubs,2,3,[]),[2 1 3]);tipsubs = reshape(tipsubs,3,[]);

% [sourcesGT,distGT] = knnsearch(tipsubs',GT.somalocs,'k',1);
somloc=zeros(length(neuron),3);for ii=1:length(neuron),somloc(ii,:)=neuron(ii).init.somaloc;end
[sourcesInit,distInit] = knnsearch(tipsubs',somloc,'k',1);
% sourcesInit = sourcesInit(index);
% distInit = distInit(index);
%%

viz = 0;
if viz
    % draw GT
    %     close all
    figure,
    %     myplot3(subs(1:10:end,:),{'k+','MarkerSize',2})
    hold on
    valinds=zeros(1,nsource);
    for ii=1:length(GT.swcpixlocs)
        if isempty(GT.swcpixlocs{ii})
            continue
        end
        valinds(ii)=1;
        % myplot3(swcpixlocs{ii},{'o','MarkerSize',2})
        gplot3(GT.connMatrix{ii},GT.swcpixlocs{ii},'-')
    end
    %     ran = max(subsin);
    xlim([0 ran(1)])
    ylim([0 ran(2)])
    zlim([0 ran(3)])
    legend(num2str([0 find(valinds)]'))
    leg= {'1','3','4','5','6','7','8','10','11'}
    legend(leg,'location','NorthEast')
    view([0 90])
    set(gca,'Ydir','reverse')
    set(gcf,'UserData',params)
    dcm_obj = datacursormode(gcf);
    set(dcm_obj,'UpdateFcn',@mycallback)
end
%%
% for each initialized neurons run a search
clear initassign dist_ii
sources = round(sourcesInit/2);
nsource =length(sources) ;
for ii=1:nsource
    [initassign{ii},dist_ii{ii}] = knnsearch(tipsubs',neuron(ii).init.swcpixlocs,'k',1);
end
%%
for ii=1:nsource
    tr=ceil(initassign{ii}/2);
    initlabelassign{ii} = unique(tr(dist_ii{ii}<10));
end
%%
% class labels of adjacent tiles
% D(1): Euclidean
% D(2): theta
% D(3): PCA
% D(4): KL
% D(5 [end]): hit

[w1,w2,w3]= find(distBr{1}); % based on euclidean tip distance
[w1,w2,w3]= find(connBr); % based on tip2shorthest distance
mm = w3<50;
mask = sparse(w1(mm),w2(mm),ones(sum(mm),1),size(distBr{1},1),size(distBr{1},1));
%
distin = connBr/1e3 .* distBr{3} .* mask;
nA = size(distin,1); % number of branches
% make is symetric
distin = max(distin,distin');

Yinit = zeros(nA,nsource);
for ii=1:nsource
    Yinit(initlabelassign{ii},ii)=1;
end

% uncertain labeling
Yinit(find(sum(Yinit,2)>1),:) = 0;

%%
%%%%%%%%%%%%%%%%%%%%%% QUERIES @@@@@@@@@@@@@@@@
% distin([38979],38826 ) = eps;
% distin(38826,38979) = eps;
% %%%%%%%%%%%%%%%%%%%%%% QUERIES - GLOBAL @@@@@@@@@@@@@@@@
% Yinit([132375, 71798,117809,118728,120008,155836,11342],2) = 1;
% Yinit([104562,153337],3) = 1;
% Yinit([132922,126964,131488,132585],1) = 1;
%%%%%%%%%%%%%%%%%%%%%% QUERIES - LOCAL @@@@@@@@@@@@@@@@
% Yinit([132375, 71798,117809,118728,120008,155836,11342],2) = 1;
% Yinit([104562,153337],3) = 1;
% Yinit([132922,126964,131488,132585],1) = 1;

% one pass
alpha = 1000;
conf = .5;
labY = Yinit;
sig = .01;

%%%%%%%%%%%%%%%%%%%%%% ------ @@@@@@@@@@@@@@@@@
L = getlaplacian(distin,sig);
[P,probs,clustL] = clusterGraph(L,alpha,labY,nsource);
clustL(isnan(P)) = 0;
%%
clustInit=clustL*0;
for ii=1:size(labY,1)
    ix = find(labY(ii,:));
    if ~isempty(ix)
        clustInit(ii) = ix;
    end
end
%%
% userdatin = {params,subs,nodeBrid,branchcent,opt}; % usefull for pointer callback functions
pointercallbackdata = {params,subs,nodeBrid,branchcent,opt}; % usefull for pointer callback functions
clustdata = {clustL,branches,neuron,ran,clustInit}

%%
userdata.pointercallbackdata.subs = subs;
userdata.pointercallbackdata.nodeBrid = nodeBrid;
userdata.pointercallbackdata.branchcent = branchcent;
userdata.pointercallbackdata.opt = opt;
userdata.clustdata.clustL = clustL;
userdata.clustdata.clustInit = clustInit;
userdata.clustdata.branches = branches;
userdata.clustdata.neuron = neuron;
userdata.clustdata.ran = ran;
%%
visualizeClusterGUI(userdata)
%%
fignum=1
vizClust(fignum,userdatin,clust,branches,GT,ran)



