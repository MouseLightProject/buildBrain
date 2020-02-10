function named_tree = process_single_component_as_function_old(component_id, ...
                                                           component, ...
                                                           component_count, ...
                                                           A_for_component, ...
                                                           ijks_for_component, ...
                                                           size_threshold, ...
                                                           length_threshold, ...
                                                           do_visualize, ...
                                                           origin_in_nm, ...
                                                           spacing_in_nm, ...
                                                           voxres, ...
                                                           options, ...
                                                           params)
    component_size= length(component) ;
    
%     if ismember(component_id, [22361 26417 28631]) ,
%         do_visualize = true ;
%     end    
    
    % Output some info
    %fprintf('Processing component with id %d, size is %d nodes...\n', component_id, component_size) ;
    
%     % Get the graph for this component
%     A_for_component = G_for_component.adjacency ;

    % Do something
    % profile clear
    % profile -memory on
    eout = buildgraph_old(A_for_component) ;
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
    dA_struct.dA = sparse(eout(:,1),eout(:,2),1,component_size,component_size);
    dA_struct.D = ones(component_size,1);
    dA_struct.R = ones(component_size,1);
    dA_struct.X = ijks_for_component(:,1);
    dA_struct.Y = ijks_for_component(:,2);
    dA_struct.Z = ijks_for_component(:,3);

    if do_visualize
        figure() ;
        gplot3(dA_struct.dA,[dA_struct.X,dA_struct.Y,dA_struct.Z]);
        title(sprintf('Component %d', component_id)) ;
        drawnow
    end
    
    %%
    did_prune_some = true ;  % just to get prune_tree() to be called at least once
    while did_prune_some ,
        [dA_struct, did_prune_some] = prune_tree_old(dA_struct, length_threshold, voxres) ;
        if do_visualize
            hold on
            gplot3(dA_struct.dA,[dA_struct.X,dA_struct.Y,dA_struct.Z]);
            drawnow
        end
    end
    %%
    % shuffle root to one of the leafs for efficiency and not
    % splitting long stretches into half
    if size(dA_struct.dA,1)>1
        inupdate_A = max(dA_struct.dA, dA_struct.dA') ;  % make into an undirected graph adjacency matrix
        eoutprun = buildgraph_old(inupdate_A) ;
        component_size = max(eoutprun(:));
        dA_struct.dA = sparse(eoutprun(:,1),eoutprun(:,2),1,component_size,component_size);
    else
    end
    if do_visualize
        hold on
        gplot3(dA_struct.dA,[dA_struct.X,dA_struct.Y,dA_struct.Z],'LineWidth',3);
        drawnow
    end
    %%
    if length(dA_struct.X)<size_threshold
        fprintf('Component with id %d fails to meet the size threshold (%d) after pruning, so discarding.  (It contains %d nodes.)\n', ...
                component_id, ...
                size_threshold, ...
                length(dA_struct.dA)) ;
        return
    end
    %%
    swc_struct_smoothed = smoothtree_old(dA_struct,options);
    %%
    % [L,list] = getBranches(inupdate_smoothed.dA);
%     if 0
%         outtree = inupdate_smoothed;
%     else
        % outtree_old = downSampleTree(inupdate_smoothed,opt);
    outtree_in_voxels = sampleTree_old(swc_struct_smoothed, options, params) ;
    outtree_in_voxels.units = 'voxels' ;
    if do_visualize
        cla
        gplot3(swc_struct_smoothed.dA,[swc_struct_smoothed.X,swc_struct_smoothed.Y,swc_struct_smoothed.Z],'LineWidth',3);
        hold on
        gplot3(outtree_in_voxels.dA,[outtree_in_voxels.X,outtree_in_voxels.Y,outtree_in_voxels.Z],'--','LineWidth',3);
        drawnow
    end        

%     % Check for an off-by-one shift
%     tree_ijks = [outtree_in_voxels.X outtree_in_voxels.Y outtree_in_voxels.Z] ;
%     is_tree_point_a_skeleton_point_from_ijk = is_point_drawn_from_pool(tree_ijks, ijks_for_component, [1 1 1]) ;
%     fraction_of_tree_points_that_match_skeleton_points = mean(is_tree_point_a_skeleton_point_from_ijk) ;
    
    % Convert centerpoint from voxel coords to um
    outtree = convert_centerpoint_units_to_um(outtree_in_voxels, origin_in_nm, spacing_in_nm) ;

    % Convert to a named tree
    digits_needed_for_index = floor(log10(component_count)) + 1 ;
    tree_name_template = sprintf('auto-cc-%%0%dd', digits_needed_for_index) ;  % e.g. 'tree-%04d'
    tree_name = sprintf(tree_name_template, component_id) ;
    named_tree = named_tree_from_tree_as_dA_struct(outtree, tree_name) ;    
    
%     % Compute the output file path
%     tree_mat_file_name = sprintf('%s.mat', tree_name) ;
%     tree_mat_file_path = fullfile(output_folder_path, tree_mat_file_name);    
    
%     % Write full tree as a .mat file
%     save_named_tree_as_mat(tree_mat_file_path, named_tree) ;
end
