function named_tree = process_single_component_as_function(component_id, ...
                                                           component_count, ...
                                                           A_for_component, ...
                                                           xyz_for_component, ...
                                                           smoothing_filter_width, ...
                                                           length_threshold, ...
                                                           do_visualize, ...
                                                           sampling_interval)

    % Get a spanning tree
    dA_spanning = spanning_tree_adjacency_from_graph_adjacency(A_for_component) ;
    A_spanning = max(dA_spanning, dA_spanning') ;
    
    if do_visualize
        figure() ;
        gplot3(A_spanning, xyz_for_component) ;
        title(sprintf('Component %d', component_id)) ;
        drawnow
    end
    
    %%
    [A_pruned, xyz_pruned] = prune_tree(A_spanning, xyz_for_component, length_threshold, do_visualize) ;    
    if do_visualize
        hold on
        gplot3(A_pruned, xyz_pruned, 'LineWidth', 3) ;
        drawnow
    end
    
%     %%
%     pruned_node_count = size(xyz_pruned, 1) ;    
%     if pruned_node_count<size_threshold ,
%         fprintf('Component with id %d fails to meet the size threshold (%d) after pruning, so discarding.  (It contains %d nodes.)\n', ...
%                 component_id, ...
%                 size_threshold, ...
%                 pruned_node_count) ;
%         return
%     end

    % Decompose tree into chains
    pruned_node_ids_from_chain_id = chains_from_tree(A_pruned) ;

    % Smooth the z coords of chain nodes
    xyz_smoothed = smooth_chains(pruned_node_ids_from_chain_id, xyz_pruned, smoothing_filter_width) ;
    
    %%
    %[A_decimated, xyz_decimated] = decimate_tree(A_pruned, xyz_pruned, sampling_interval) ;    
        
    % Decimate chains
    decimated_node_ids_from_chain_id = decimate_chains(pruned_node_ids_from_chain_id, xyz_smoothed, sampling_interval) ;

    % Convert chains to edges
    edges_using_decimated_node_ids = edges_from_chains(decimated_node_ids_from_chain_id) ;

    % Defragment the node ids
    [edges_using_decimated_node_ids, xyz_decimated, pruned_node_id_from_decimated_node_id] = ...
        defragment_node_ids_in_edges(edges_using_decimated_node_ids, xyz_smoothed) ;

    % Convert edges to (sparse) adjacency
    decimated_node_count = length(pruned_node_id_from_decimated_node_id) ;
    A_decimated = undirected_adjacency_from_edges(edges_using_decimated_node_ids, decimated_node_count) ;    
    
    % Visualize the decimated tree
    if do_visualize
        hold on
        gplot3(A_decimated, xyz_decimated, '--', 'LineWidth', 3) ;
        drawnow
    end

    % Convert to a named tree
    digits_needed_for_index = floor(log10(component_count)) + 1 ;
    tree_name_template = sprintf('auto-cc-%%0%dd', digits_needed_for_index) ;  % e.g. 'tree-%04d'
    tree_name = sprintf(tree_name_template, component_id) ;
    named_tree = named_tree_from_undirected_graph(A_decimated, xyz_decimated, tree_name) ;    
    
%     % Compute the output file path
%     tree_mat_file_name = sprintf('%s.mat', tree_name) ;
%     tree_mat_file_path = fullfile(output_folder_path, tree_mat_file_name);    
    
    % Write full tree as a .mat file
    %save_named_tree_as_mat(tree_mat_file_path, named_tree) ;
end
