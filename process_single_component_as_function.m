function named_tree = process_single_component_as_function(component_id, ...
                                                           max_component_id, ...
                                                           A_for_component, ...
                                                           xyz_for_component, ...
                                                           smoothing_filter_width, ...
                                                           length_threshold, ...
                                                           do_visualize, ...
                                                           sampling_interval)

    % Plot the input                                                       
    if do_visualize ,
        f = figure('Color', 'w') ;
        ax = axes(f) ;
        l0 = gplot3(ax, A_for_component, xyz_for_component, 'b') ;  % blue
        hold(ax, 'on') ;
        drawnow
    end
                                                       
    % Get a spanning tree
    dA_spanning = spanning_tree_adjacency_from_graph_adjacency(A_for_component) ;
    A_spanning = max(dA_spanning, dA_spanning') ;
    
    if do_visualize
        l1 = gplot3(ax, A_spanning, xyz_for_component, 'Color', [0 0.5 0]) ;   % green
        drawnow
    end
    
    %%
    [A_pruned, xyz_pruned] = prune_tree(A_spanning, xyz_for_component, length_threshold, do_visualize) ;    
    if do_visualize
        l2 = gplot3(ax, A_pruned, xyz_pruned, 'r') ;  % red
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
    %xyz_smoothed = xyz_pruned ;

    % Visualize the smoothed tree
    if do_visualize        
        l3 = gplot3(ax, A_pruned, xyz_smoothed, 'Color', [1 2/3 0]) ;   % orange
        drawnow
    end  
    
    %%
    %[A_decimated, xyz_decimated] = decimate_tree(A_pruned, xyz_pruned, sampling_interval) ;    
        
    % Decimate chains
    decimated_node_ids_from_chain_id = decimate_chains(pruned_node_ids_from_chain_id, xyz_smoothed, sampling_interval) ;
    %decimated_node_ids_from_chain_id = pruned_node_ids_from_chain_id ;
    
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
        l4 = gplot3(ax, A_decimated, xyz_decimated, 'Color', [0.6 0 0.8]) ;   % purple
        legend([l0 l1 l2 l3 l4], {'original', 'tree', 'pruned', 'smoothed', 'decimated'}) ;        
        axis(ax, 'vis3d') ;
        title(ax, sprintf('Component %d', component_id)) ;
        xlabel(ax, 'x (um)') ;
        ylabel(ax, 'y (um)') ;
        zlabel(ax, 'z (um)') ;                
        drawnow
    end

    % Convert to a named tree
    digits_needed_for_index = floor(log10(max_component_id)) + 1 ;
    tree_name_template = sprintf('auto-cc-%%0%dd', digits_needed_for_index) ;  % e.g. 'tree-%04d'
    tree_name = sprintf(tree_name_template, component_id) ;
    named_tree = named_tree_from_undirected_graph(A_decimated, xyz_decimated, tree_name) ;    
    
%     % Compute the output file path
%     tree_mat_file_name = sprintf('%s.mat', tree_name) ;
%     tree_mat_file_path = fullfile(output_folder_path, tree_mat_file_name);    
    
    % Write full tree as a .mat file
    %save_named_tree_as_mat(tree_mat_file_path, named_tree) ;
end
