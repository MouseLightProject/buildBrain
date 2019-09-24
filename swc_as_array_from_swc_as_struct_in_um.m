function tree_as_swc_array = swc_as_array_from_swc_as_struct_in_um(tree_as_struct)
    dA = tree_as_struct.dA ;
    node_count = size(dA,1) ;
    %XYZout = [outtree.X,outtree.Y,outtree.Z]-1;  % this is what
    %caused the horrble off-by-one issues in August 2019!
    xyz = [tree_as_struct.X,tree_as_struct.Y,tree_as_struct.Z] ;  % must be in um
    r = tree_as_struct.R;
    d = tree_as_struct.D;
    % transform location
    %XYZout = pix2um(params,XYZout); % center anisotropy to compansate imresize
    %

    At = dA+dA';
    [DISC,~,~] = graphtraverse(At,1,'Method','DFS');
    At(1:end,:) = At(DISC,:);
    At(:,1:end) = At(:,DISC);
    xyz = xyz(DISC,:);
    r = r(DISC,:);
    d = d(DISC,:);
    
    
    %%
%     if strcmp(version('-release'),'2015b')
%         [dist,pred] = graphalgs(algorithmkeys{algorithm},debug_level,directed,At,1);
%     else
    [~,~,pred] = graphshortestpath(At,1,'DIRECTED',false);
%     end
    tree_as_swc_array = [(1:node_count)' d xyz r pred(:)] ;
    tree_as_swc_array(1,7) = -1;
end
