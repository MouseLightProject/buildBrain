function process_single_component(output_folder_path, ...
                                  component_id, ...
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
                                  sampling_style, ...
                                  sampling_interval, ...
                                  largesampling)
    component_size= length(component) ;
    
%     if ismember(component_id, [22361 26417 28631]) ,
%         do_visualize = true ;
%     end    
    
    % Output some info
    %fprintf('Processing component with id %d, size is %d nodes...\n', component_id, component_size) ;
    
%     % Get the graph for this component
%     A_for_component = G_for_component.adjacency ;

    if ~is_component_a_binary_tree(A_for_component) ,
        do_visualize = true ;
        keyboard
    end
    
    % Do something
    % profile clear
    % profile -memory on
    child_parent_node_id_pairs = rooted_spanning_tree_edges_from_graph(A_for_component) ;
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
    dA_struct.dA = sparse(child_parent_node_id_pairs(:,1),child_parent_node_id_pairs(:,2),1,component_size,component_size);
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
    dA_struct = prune_tree(dA_struct, length_threshold, voxres, do_visualize) ;    
    
    %%
    % shuffle root to one of the leafs for efficiency and not
    % splitting long stretches into half
    if size(dA_struct.dA,1)>1
        inupdate_A = max(dA_struct.dA, dA_struct.dA') ;  % make into an undirected graph adjacency matrix
        eoutprun = rooted_spanning_tree_edges_from_graph(inupdate_A) ;
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
    dA_struct_smoothed = smooth_tree(dA_struct, size_threshold);
    if do_visualize
        hold on
        gplot3(dA_struct_smoothed.dA,[dA_struct_smoothed.X,dA_struct_smoothed.Y,dA_struct_smoothed.Z],'LineWidth',3, 'Color', 'm');
        drawnow
    end        
    
    %%
    outtree_in_voxels = sample_tree(dA_struct_smoothed, voxres, sampling_style, sampling_interval, length_threshold, largesampling) ;
    outtree_in_voxels.units = 'voxels' ;
    if do_visualize
        cla
        gplot3(dA_struct_smoothed.dA,[dA_struct_smoothed.X,dA_struct_smoothed.Y,dA_struct_smoothed.Z],'LineWidth',3);
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
    
    % Compute the output file path
    tree_mat_file_name = sprintf('%s.mat', tree_name) ;
    tree_mat_file_path = fullfile(output_folder_path, tree_mat_file_name);    
    
    % Write full tree as a .mat file
    %save_named_tree_as_mat(tree_mat_file_path, named_tree) ;
end
