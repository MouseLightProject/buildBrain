function [hits_gt,hits_delete,swcout,gtfile] = cropSectionBasedOnGT(params,gt_swcfolder,subs)
hits_delete= [];
%% TODO upsample swc to make search better !!!!
[swcout,swcfiles] = loadSwcFolder(gt_swcfolder);
%%
% parse swcs into GT/init/mask
% make sure to add stuff for GT and delete stuff from mask
% PATCH
gtfile = zeros(1,length(swcout));
gtum = cell(1,length(swcout));
parfor iswc=1:length(swcout)
    swcfileparts = strsplit(swcfiles(iswc).name,'.');
    if ismember(swcfileparts{1}(1),{'R','G'})
        gtfile(iswc) = 1;
        % upsample swc to get better covarage
        e1=swcout{iswc}(2:end,1);
        e2=swcout{iswc}(2:end,7);
        E345=swcout{iswc}(1:end,3:5);
        
        G = sparse(e1,e2,1,e1(end),e1(end));
        L = getBranches(max(G,G'));
        %%
        updata=[];
        for il = 1:length(L)
            e345 = E345(L(il).set,:);
            pd=round(sqrt(sum(diff(e345,1,1).^2,2)));
            updata_ = [];
            for jj=1:length(pd)
                st = e345(jj,:);
                en = e345(jj+1,:);
                sl = en-st;
                sp = [[1:1e3:pd(jj)-1 pd(jj)]/pd(jj)];
                %updata = [updata sl(:)*sp+st(:)*ones(1,length(sp))];
                updata_{jj} = [sl(:)*sp+st(:)*ones(1,length(sp))];
            end
            if isempty(updata_)
                continue
            end
            updata_=[updata_{:}]';
            updata{il} = updata_;
        end
        gtum{iswc}=cat(1,updata{:});
    end
end
%%
gtum = cat(1,gtum{:});
% make sure to keep these indicies
hits_gt = maskWithInputLocations(subs,um2pix(params,gtum)+1);

%% stuf that enters and exists temporaryly
deleteum = [];
for iswc=1:length(swcout)
    swcfileparts = strsplit(swcfiles(iswc).name,'.');
    if ismember(swcfileparts{1}(1:4),{'mask'})
        deleteum{iswc} = swcout{iswc}(:,3:5);
    end
end
if ~isempty(deleteum)
    deleteum = cat(1,deleteum{:});
    hits_delete = maskWithInputLocations(subs,um2pix(params,deleteum));
end
% make sure to keep these indicies

