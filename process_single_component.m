function process_single_component(output_folder_path, ...
                                  component_id, ...
                                  component, ...
                                  G, ...
                                  subs, ...
                                  size_threshold, ...
                                  length_threshold, ...
                                  do_visualize, ...
                                  origin_in_nm, ...
                                  spacing_in_nm, ...
                                  voxres, ...
                                  options)
    component_size= length(component) ;
    
    % Output some info
    fprintf('Processing component with id %d, size is %d nodes...\n', component_id, component_size) ;

    % Check if too small, return if so
    if component_size<=size_threshold ,
        return
    end                   

    % Check if the file already exists, return if so
    tree_mat_file_name = sprintf('auto-cc-%06d.mat', component_id) ;
    tree_mat_file_path = fullfile(output_folder_path, tree_mat_file_name);
    if exist(tree_mat_file_path, 'file') ,
        return
    end

    % Get the graph for this component
    subs_for_component = subs(component,:) ;
    G_for_component = G.subgraph(component) ;
    A_for_component = G_for_component.adjacency ;

    % Do something
    root_node_id = find(sum(A_for_component)==1, 1) ;
    % profile clear
    % profile -memory on
    eout = graphfuncs.buildgraph(A_for_component, root_node_id) ;
    % profile off
    
    % Package things up into an "SWC structure".
    % dA: directed adjacency graph.  A sparse matrix that represents the
    % directed edges that point from a node to the next node soma-ward (or,
    % more generally, rootward).
    % D: node descriptor.  This will become the 2nd column of the SWC, the
    % "structure identifier".  Used for different things in different
    % places.
    % R: The radius as the node, in um.  We basically don't use this, so it's
    % often just set to all ones for convenience.
    % X,Y,Z: The coordinates of each node, ideally in um, but sometimes we
    % have them in voxel indices.
    swc_struct.dA = sparse(eout(:,1),eout(:,2),1,component_size,component_size);
    swc_struct.D = ones(component_size,1);
    swc_struct.R = ones(component_size,1);
    swc_struct.X = subs_for_component(:,1);
    swc_struct.Y = subs_for_component(:,2);
    swc_struct.Z = subs_for_component(:,3);
    
    %%
    did_prune_some = true ;  % just to get prunTree() to be called at least once
    while did_prune_some ,
        [swc_struct, did_prune_some] = prunTree(swc_struct, length_threshold, voxres) ;
        if do_visualize
            hold on
            gplot3(swc_struct.dA,[swc_struct.X,swc_struct.Y,swc_struct.Z]);
            drawnow
        end
    end
    %%
    % shuffle root to one of the leafs for efficiency and not
    % splitting long stretches into half
    if size(swc_struct.dA,1)>1
        inupdate_A = max(swc_struct.dA, swc_struct.dA') ;  % make into an undirected graph adjacency matrix
        eoutprun = graphfuncs.buildgraph(inupdate_A, root_node_id) ;
        component_size = max(eoutprun(:));
        swc_struct.dA = sparse(eoutprun(:,1),eoutprun(:,2),1,component_size,component_size);
    else
    end
    if do_visualize
        hold on
        gplot3(swc_struct.dA,[swc_struct.X,swc_struct.Y,swc_struct.Z],'LineWidth',3);
        drawnow
    end
    %%
    if length(swc_struct.dA)<size_threshold
        %fprintf('Component with id %d fails to meet the size threshold after pruning, so discarding.\n', component_id) ;
        return
    end
    %%
    swc_struct_smoothed = smoothtree(swc_struct,options);
    %%
    % [L,list] = getBranches(inupdate_smoothed.dA);
%     if 0
%         outtree = inupdate_smoothed;
%     else
        % outtree_old = downSampleTree(inupdate_smoothed,opt);
    outtree_in_voxels = sampleTree(swc_struct_smoothed, options) ;
    outtree_in_voxels.units = 'voxels' ;
    if do_visualize
        cla
        gplot3(swc_struct_smoothed.dA,[swc_struct_smoothed.X,swc_struct_smoothed.Y,swc_struct_smoothed.Z],'LineWidth',3);
        hold on
        gplot3(outtree_in_voxels.dA,[outtree_in_voxels.X,outtree_in_voxels.Y,outtree_in_voxels.Z],'--','LineWidth',3);
        drawnow
    end        

    % Convert centerpoint from voxel coords to um
    outtree = convert_centerpoint_units_to_um(outtree_in_voxels, origin_in_nm, spacing_in_nm) ;

    % Write full tree as a .mat file
    save_tree_as_mat(tree_mat_file_path, component_id, outtree) ;
end    
