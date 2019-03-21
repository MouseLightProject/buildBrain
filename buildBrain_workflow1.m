function buildBrain_workflow1(configfile)
%BUILDBRAIN Aggregrates skeletonization results and build connectivity
%graph of the brain
%
% [OUTPUTARGS] = BUILDBRAIN(INPUTARGS) Explain usage here
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

% $Author: base $	$Date: 2016/05/17 17:19:01 $	$Revision: 0.1 $
% Copyright: HHMI 2016
% if nargin<1
%     configfile = './config_files/config_buildBrain_20180801_prob0_repeat.cfg';
% %     configfile = './config_files/config_buildBrain_20180702_prob0_100.cfg';
% %        configfile = './config_files/config_buildBrain_20181001_prob0.cfg';
%     %     configfile = './config_files/config_buildBrain_20180309_prob1.cfg';
%     %     configfile = './config_files/config_buildBrain_20170911_prob0.cfg';
%     %     configfile = './config_files/config_buildBrain_20150619_octant12_prob0.cfg';
% end
% addpath(genpath('./common'))
% addpath(genpath('./functions'))
opt = configparser(configfile);
if ~isfield(opt,'sampling')
    opt.sampling = 'uni';
end
myh5 = opt.inputh5;
myh5prob = opt.h5prob;
%[~,name] = fileparts(myh5);
%h5infofile = fullfile('./h5infos',['h5inf_',name,'.mat']);
%mkdir(fileparts(h5infofile))

[brainSize,~,~,~] = h5parser(myh5,myh5prob);
%%
origin = h5read(opt.inputh5,[opt.h5prob,'_props/origin']);
spacing = h5read(opt.inputh5,[opt.h5prob,'_props/spacing']);
level = h5read(opt.inputh5,[opt.h5prob,'_props/level']);
params.outsiz = brainSize;
params.ox = origin(1);
params.oy = origin(2);
params.oz = origin(3);
params.sx = spacing(1);
params.sy = spacing(2);
params.sz = spacing(3);
params.level = level;

params.voxres = [params.sx params.sy params.sz]/2^(params.level)/1e3; % in um
opt.params = params;
%%
% create output folders
full_folder_path = fullfile(opt.outfolder,'full') ;
if ~exist(full_folder_path, 'file') ,    
    mkdir(full_folder_path) ;
end
frags_folder_path = fullfile(opt.outfolder,'frags') ;
if ~exist(frags_folder_path, 'file') ,    
    mkdir(frags_folder_path) ;
end
%%
[subs,~,A,~] = skel2graph(opt);
%subs_ori = subs;
%edges_ori = edges;
%weights_ori = weights;
%A_ori = A;
%%
%subs = subs_ori;
%edges = edges_ori;
%weights = weights_ori;
%A = A_ori;
%%
%maskids = []
% maskids = [56 672 1022 1031]; %ACB/Caudoputamen/Globus pallidus, external segment/Globus pallidus, internal segment
%maskids = [956 844 882 686 56 1022 1031 1021 1085 719 882 583 182305705 182305709 182305713]
%gt_swcfolder = '/nrs/mouselight/seggui/swcfiles/GT/2017-09-25_striatum_neurons_temp'
% if ~isempty(maskids)
%     %%
%     [hits_allen_brain] = maskWithAllenAtlas(params,subs_ori,maskids);
%     if exist('gt_swcfolder','var')
%         addpath(genpath('./scripts'))
%         [hits_gt,hits_delete,swcout,gtfile] = cropSectionBasedOnGT(params,gt_swcfolder,subs);
%         % keep allen_brain and gt then substract delete
%         keepthese = union(setdiff(hits_allen_brain,hits_delete),hits_gt);
%         figure, myplot3(subs(keepthese,:),'.')
%     else
%         keepthese = hits_allen_brain;
%     end
%     %%
%     subs = subs(keepthese,:);
%     A = A(keepthese,:);
%     A = A(:,keepthese);
% end
%
% delete junks due to ventricle
%[~,Comps] = graphconncomp(A,'DIRECTED',false);
%Y = histcounts(Comps,1:S+1);
% if 0
%     junk_locs = [[73336.7, 14360.0, 34718.8];
%     [71850.7, 14436.3, 34497.8];
%     [71870.7, 14308.7, 34508.5];
%     [73264.9, 14183.1, 34508.5];
%     [73328.9, 14801.0, 34508.5];
%     [71817.1, 14684.7, 34508.5];
%     [71875.3, 14339.0, 34448.2];
%     [73254.4, 14099.3, 34436.1];
%     [73366.1, 14453.1, 34753.0];
%     ];
% end
% if exist('junk_locs','var')
%     figure(100)
%     cla
%     ic_ids_all = zeros(1,length(Comps));
%     for ilocs = 1:size(junk_locs,1)
%         [aa,bb] = min(sum(abs(subs-um2pix(params,junk_locs(ilocs,:))),2));
%         ic_ids=Comps==Comps(bb);
%         ic_ids_all = ic_ids_all|ic_ids;
%         hold on
%         subs_ = subs(ic_ids,:);
%         A_ = A(ic_ids,:); A_ = A_(:,ic_ids);
%         gplot3(A_,subs_,'-');
%     end
%     keepthese = find(~ic_ids_all);
%     subs = subs(keepthese,:);
%     A = A(keepthese,:);
%     A = A(:,keepthese);
% end
%%
% junk_ids = junk_ids_070218sample();
% junk_ids = [];
% if exist('junk_ids','var') && ~isempty(junk_ids) ,
%     %%
%     if opt.viz
%         figure(100)
%         cla
%         for ilocs = 1:size(junk_ids,1)
%             ic_ids=Comps==junk_ids(ilocs);
%             hold on
%             subs_ = subs(ic_ids,:);
%             A_ = A(ic_ids,:); A_ = A_(:,ic_ids);
%             gplot3(A_,subs_,'-');
%         end
%     end
%     ic_ids_all = zeros(1,length(Comps));
%     for ilocs = 1:size(junk_ids,1)
%         ic_ids=Comps==junk_ids(ilocs);
%         ic_ids_all = ic_ids_all|ic_ids;
%     end
%     %%
%     keepthese = find(~ic_ids_all);
%     subs = subs(keepthese,:);
%     A = A(keepthese,:);
%     A = A(:,keepthese);
% end

%%
% [A,subs] = filterEdges(A,subs,params)

% %
% tstart = tic;
% affinityBuilder(opt,A,subs)
% sprintf('FINISHED IN: %d', round(toc(tstart)))
% return

%%
Gin = graph(max(A,A'));
workflow1(Gin,subs,opt)

% %% filter out results based on manual recons - this is for presentation only
% if 0
%     %swcfolder = '/groups/mouseprojects/home/base/CODE/MATLAB/recontree/malazgirt/singlespace'
%     %%
%     %     swcfolder = '/groups/mouseprojects/home/base/CODE/MATLAB/recontree/malazgirt/singlespace'
%     %     swcfiles = dir([swcfolder,'/*.swc']);
%     swcfolder = '/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/2017-09-11/cropSWC'
%     swcfiles = dir(fullfile(swcfolder,'*.swc'))
%     numSWCs = length(swcfiles);
%     %%
%     DATA = [];
%     upDATA = [];
%     for ii=1:numSWCs
%         %%
%         [swcData,offset,color, header] = loadSWC(fullfile(swcfolder,swcfiles(ii).name));
%         swcData(:,3:5) = swcData(:,3:5) + ones(size(swcData,1),1)*offset;
%         swcData(:,3:5) = swcData(:,3:5)*1000;
%         %%
%         % upsample points that are consequtive
%         upthese = diff(swcData(:,7))==1;
%         pd=round(sqrt(sum(diff(swcData(:,3:5)).^2,2)));
%         idxup = find(upthese);
%         updata_ = [];
%         for jj=idxup(:)'
%             st = swcData(jj,:);
%             en = swcData(jj+1,:);
%             sl = en-st;
%             sp = [[1:1e3:pd(jj)-1 pd(jj)]/pd(jj)];
%             %updata = [updata sl(:)*sp+st(:)*ones(1,length(sp))];
%             updata_{jj} = [sl(:)*sp+st(:)*ones(1,length(sp))];
%         end
%         upDATA = cat(1,upDATA,[updata_{:}]');
%         DATA=cat(1,DATA,swcData);
%     end
%     %%
%     pixsDATA = um2pix(params,DATA(:,3:5)/1000);
%     pixsupdata = um2pix(params,upDATA(:,3:5)/1000);
%     pixs = pixsupdata;
%     %     pixs = pixsDATA;
%     % pixs(pixs(:,2)<3000,:)=[];
%     %
%     [idxknn,dknn]=knnsearch(subs,pixs,'K',1);
%     validC = unique(Comps(idxknn(dknn<15)));
%     %%
%     if opt.viz % visualize clusters
%         %%
%         figure(22),
%         cla
%         plot3(pixs(:,1),pixs(:,2),pixs(:,3),'.k','MarkerSize',1)
%         hold on
%         axis equal tight
%         set(gca,'Ydir','reverse')
%         iter = 0;
%         for mC=ib(ia>100)%(1:100)%validC%1:size(Y,2)%ib(2:500)%
%             % for each cluster run reconstruction
%             if Y(mC)>opt.sizethreshold & ismember(mC,validC) %& ~(mC==ib(1))
%                 %             if Y(mC)>opt.sizethreshold & ismember(mC,validC)
%                 iter = iter+1;
%                 if ~rem(iter,round(length(validC)/100))
%                     iter
%                 end
%                 subidx = find(Comps==mC);
%                 plot3(subs(subidx,1),subs(subidx,2),subs(subidx,3),'.')
%             end
%             drawnow
%         end
%     end
% else
%     validC = 1:size(Y,2);
% end
% %%
% if 0
%     %%
%     figure(90),
%     plot3(subs(:,1),subs(:,2),subs(:,3),'.k','MarkerSize',.1)
%     hold on
%     axis equal tight
%     set(gca,'Ydir','reverse')
%     %%
%     figure(21),
%     cla
%     plot3(pixs(:,1),pixs(:,2),pixs(:,3),'.k','MarkerSize',.1)
%     hold on
%     axis equal tight
%     set(gca,'Ydir','reverse')
% end
% %%
% if 0
%     if 0
%         reconComp(opt,validC,Y,A_,Comps,subs,params)
%     else
%         Gin = graph(max(A,A'));
%         workflow1(Gin,subs,opt)
%     end
% else
%     tstart = tic;
%     affinityBuilder(opt,validC,Y,A_,Comps,subs)
%     sprintf('FINISHED IN: %d', round(toc(tstart)))
% end
% %%
% % reconComptest(opt,validC,Y,A_,Comps,subs,params)

end
