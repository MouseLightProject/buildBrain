function [outtree,selectedIdx] = sample_tree(intree, voxres, sampling_style, sampling_interval, lengthThr, largesampling)
%SAMPLETREE Summary of this function goes here
%
% [OUTPUTARGS] = SAMPLETREE(INPUTARGS) Explain usage here
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

% $Author: base $	$Date: 2016/07/18 16:22:30 $	$Revision: 0.1 $
% Copyright: HHMI 2016
[L,list] = get_branches(intree.dA);
XYZ = [intree.X intree.Y intree.Z];
R = intree.R;
D = intree.D;
%%
% down sample branches
selectedIdx = [];
for ii=1:length(L)
    %%
    set_ii = L(ii).set;
    % set_ii: +----x, includes child/excludes parent: indicies are
    % towards root
    % spacing is a function of branch length
    %%
    if isempty(set_ii)
        continue
    end
    %%
    if 1
        %%
        % flip indicies so that lower indicies are close to branching
        set_ii = [L(ii).parentnode set_ii(end:-1:1)];
        % get length of branch
        xyz = XYZ(set_ii,:);
        xyz = xyz.*(ones(size(xyz,1),1)*voxres);
        switch sampling_style
            case 'uni'
                p2 = sum((xyz(1:end-1,:)-xyz(2:end,:)).^2,2);
            case 'curv'
                p2 = (xyz(1:end-1,:)-xyz(2:end,:)).^2*[1;1;5];
            otherwise
                error('"%s" is not an allowed sampling style', sampling_style) ;
        end
        dists = sqrt(p2); % resolution is 0.33 um per pixel => multiply with 3 to make it unit um
        cdist = [0;cumsum(dists)];
        if cdist(end)>sampling_interval
            % make sure tip is in the list
            [aa,idx]=min(pdist2([0:sampling_interval:cdist(end)-sampling_interval]',cdist(:)),[],2);
            indicies = [set_ii(idx) set_ii(end)];
        else
            indicies = [set_ii(1) set_ii(end)];
        end
%         if length(indicies) > 3
%             selectedIdx{ii} = indicies(3:end)';
%         else
            selectedIdx{ii} = indicies(2:end)';
%         end
    else
        %%
        % flip indicies so that lower indicies are close to branching
        set_ii = set_ii(end:-1:1);
        % get length of branch
        xyz = XYZ(set_ii,:);
        xyz = xyz.*(ones(size(xyz,1),1)*voxres);
        switch sampling_style
            case 'uni'
                p2 = sum((xyz(1:end-1,:)-xyz(2:end,:)).^2,2);
            case 'curv'
                p2 = (xyz(1:end-1,:)-xyz(2:end,:)).^2*[1;1;50];
            otherwise
        end
        dists = sqrt(p2); % resolution is 0.33 um per pixel => multiply with 3 to make it unit um
        cdist = cumsum(dists);
        
        if isempty(cdist) | cdist(end)<=1.5*lengthThr %(um) : branch is a connactor and small
            % get the last index as an indicator
            indicies = set_ii(end); % keep the branch point (set flows from child towards parent)
        elseif cdist(end)>1.5*lengthThr & cdist(end)<largesampling
            %sample every 15(um): sampling interval for short branches
            sampling_style = min(1.5*lengthThr,15);
            [indicies] = sampleXum(cdist,sampling_style,set_ii);
        elseif cdist(end)>largesampling %(um)
            %sample every 50(um): sampling interval for long branches
            if 1.5*lengthThr>largesampling
                warning('pruning length is too big: %d compared to largesampling: %d',lengthThr,largesampling)
            end
            sampling_style = min(max(1.5*lengthThr,50),largesampling);
            [indicies] = sampleXum(cdist,sampling_style,set_ii);
        end
        selectedIdx{ii} = indicies;
    end
end
%%
clear edges
for ii=1:length(L)
    inds = [L(ii).parentnode;selectedIdx{ii}(:)];
    to = inds(1:end-1);
    from = inds(2:end);
    edges{ii} = [from(:) to(:)]';
end
E = [edges{:}]';
%%
[idx,ia,ic] = unique(E(:));
E_ = reshape(ic,[],2);
A_ = sparse(E_(:,1),E_(:,2),1,length(idx),length(idx));
%%
outtree.dA = A_;
outtree.X = XYZ(idx,1);
outtree.Y = XYZ(idx,2);
outtree.Z = XYZ(idx,3);
outtree.R = R(idx);
outtree.D = D(idx);

end