function branch2swc(params,branches,swcoutfolder)
%%
numbranches = length(branches);
colscomp = jet(numbranches);
colscomp = colscomp(randperm(numbranches),:);
%%
try parfor_progress(0);catch;end
parfor_progress(numbranches)
parfor iter = 1:numbranches
    parfor_progress;
    subs = branches(iter).subs;
    len = size(subs,1);
    % downsample for visualization on JW
    selected_inds = [1:5:round(len-5) len];
    if len<16
        selected_inds = round(len/2);
    else
        selected_inds = [6:5:round(len-10) len-5];
    end
    if isempty(selected_inds) | selected_inds<1
        continue
    end
    subs = subs(selected_inds,:);
    len = size(subs,1);
    
    e1 = 1:len;
    e7 = [1:len]-1;
    e7(1) = -1;
    e2 = ones(len,1);
    e6 = ones(len,1);
    e345 = pix2um(params,subs-1);
    swcData = cat(2,e1(:),e2,e345,e6,e7(:));
    
    pre = '';
    swcoutname = sprintf('auto%s',pre);
    fragname = sprintf('%s_br-%04d.swc',swcoutname,iter);
    outname = fullfile(swcoutfolder,fragname);
    if ~exist(fileparts(outname),'dir')
        mkdir(fileparts(outname))
    end
    % write swc
    fid = fopen(outname,'w');
    mytxt = sprintf('# Generated by workflow2\n');
    fprintf(fid,'%s',mytxt);
    mytxt = sprintf('# OFFSET %.6f %.6f %.6f\n',[0 0 0]);
    fprintf(fid,'%s',mytxt);
    mytxt = sprintf('# COLOR %f,%f,%f\n',colscomp(iter,1),colscomp(iter,2),colscomp(iter,3));
    fprintf(fid,'%s',mytxt);
    mytxt = sprintf('# NAME %s\n',fragname);
    fprintf(fid,'%s',mytxt);
    fprintf(fid,'%d %d %f %f %f %d %d\n',swcData');
    fclose(fid);
end
parfor_progress(0);

