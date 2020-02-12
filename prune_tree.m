function [A_result, xyz_result] = prune_tree(A_tree, xyz_tree, length_threshold, do_visualize)
    A_result = A_tree ;
    xyz_result = xyz_tree ;
    did_prune_some = true ;  % just to get prune_tree_some() to be called at least once
    while did_prune_some ,
        [A_result, xyz_result, did_prune_some] = prune_tree_some(A_result, xyz_result, length_threshold) ;
        if do_visualize ,
            hold on
            gplot3(A_result, xyz_result);
            drawnow
        end
    end
end


function [A_tree_output, xyz_output, did_prune_some] = prune_tree_some(A_tree_input, xyz, length_threshold)
    % length_threshold should be in um.  Leaf chains shorter than that will be
    % trimmed off.
    
    % Find the leaf nodes
    degree = full(sum(A_tree_input)) ;
    is_leaf_node = (degree==1) ;
    leaf_node_ids = find(is_leaf_node) ;
    leaf_node_count = length(leaf_node_ids) ;
    
    % For each leaf node, trace inwards until you get to a nexus point.
    % If the chain of nodes before the nexus is too short, mark those nodes
    % as 'doomed'.
    node_count = size(xyz, 1) ;
    is_node_doomed = false(node_count, 1) ;        
    for leaf_node_index = 1 : leaf_node_count ,
        leaf_node_id = leaf_node_ids(leaf_node_index) ;
        [node_ids_in_chain, nexus_node_id] = trace_chain_in_undirected(A_tree_input, leaf_node_id) ;
        node_ids_in_chain_with_ends = [node_ids_in_chain nexus_node_id] ;
        xyz_chain = xyz(node_ids_in_chain_with_ends, :) ;
        dr_chain = diff(xyz_chain) ;
        ds_chain = sqrt(sum(dr_chain.^2, 2)) ;  % length of each edge in chain
        chain_length = sum(ds_chain) ;
        if chain_length < length_threshold ,
            is_node_doomed(node_ids_in_chain) = true ;
        end
    end             

    % Populate the outputs.  Basically the same as the inputs, but with
    % doomed nodes deleted.
    A_tree_output = A_tree_input ;
    A_tree_output(is_node_doomed,:) = [] ;
    A_tree_output(:,is_node_doomed) = [] ;
    xyz_output = xyz(~is_node_doomed,:) ;
    did_prune_some = any(is_node_doomed) ;
end
