function [A_result, xyz_result] = prune_tree(A_tree, xyz_tree, length_threshold, do_visualize)  %#ok<INUSD>
    A_result = A_tree ;
    xyz_result = xyz_tree ;
    did_prune_some = true ;  % just to get prune_tree_some() to be called at least once
    while did_prune_some ,
        [A_result, xyz_result, did_prune_some] = prune_tree_some(A_result, xyz_result, length_threshold) ;
%         if do_visualize ,
%             hold on
%             gplot3(A_result, xyz_result);
%             drawnow
%         end
    end
end


function [A_tree_output, xyz_output, did_prune_some] = prune_tree_some(A_tree_input, xyz, length_threshold)
    % length_threshold should be in um.  Leaf chains shorter than that will be
    % trimmed off.
    
    % Get the chains
    chains = chains_from_tree(A_tree_input) ;
    
    % Find the leaf nodes
    degree = full(sum(A_tree_input)) ;
    is_leaf_node = (degree==1) ;
    leaf_node_ids = find(is_leaf_node) ;
    leaf_node_count = length(leaf_node_ids) ;
    
    % If there are exactly two leaf nodes, the tree must be a chain, so we
    % don't bother to prune
    if leaf_node_count == 2 ,
        A_tree_output = A_tree_input ;
        xyz_output = xyz ;
        did_prune_some = false ;
    else
        % If get here, the tree is not a chain
        
        % For each leaf node, trace inwards until you get to a nexus point.
        % Record the length and the terminal node_id of the proximal branch
        % point.
        node_count = size(xyz, 1) ;
        length_from_leaf_node_index = zeros(leaf_node_count, 1) ;
        branch_node_id_from_leaf_node_index = zeros(leaf_node_count, 1) ;
        node_ids_in_chain_from_leaf_node_index  = cell(leaf_node_count, 1) ;
        for leaf_node_index = 1 : leaf_node_count ,
            leaf_node_id = leaf_node_ids(leaf_node_index) ;
            [node_ids_in_chain, branch_node_id] = trace_chain_in_undirected(A_tree_input, leaf_node_id) ;
            node_ids_in_chain_with_ends = [node_ids_in_chain branch_node_id] ;
            xyz_chain = xyz(node_ids_in_chain_with_ends, :) ;
            dr_chain = diff(xyz_chain) ;
            ds_chain = sqrt(sum(dr_chain.^2, 2)) ;  % length of each edge in chain
            leaf_chain_length = sum(ds_chain) ;
            length_from_leaf_node_index(leaf_node_index) = leaf_chain_length ;
            branch_node_id_from_leaf_node_index(leaf_node_index) = branch_node_id ;
            node_ids_in_chain_from_leaf_node_index{leaf_node_index} = node_ids_in_chain ;
        end

        % If there are multiple too-short leaf chains coming out of a single branch
        % node, we want to spare the longest of these, in hopes that it's a
        % valid on even though others are spurs.
        
        % Group the leaf chains by branch node id
        leaf_node_index_from_leaf_chain_index = (1:leaf_node_count) ;
        [leaf_node_indices_from_branch_index, branch_node_id_from_branch_index] = ...
            bag_items_by_id(branch_node_id_from_leaf_node_index, leaf_node_index_from_leaf_chain_index) ;
        branch_node_count = length(branch_node_id_from_branch_index) ;  % How many unique branch nodes the leaf chains go into
        
        % For each unique branch node, mark some or the associated leaf
        % chains as being doomed
        is_doomed_from_leaf_node_index = false(leaf_node_count, 1) ;        
        for nexus_node_index = 1 : branch_node_count ,
            leaf_node_indices = leaf_node_indices_from_branch_index{nexus_node_index} ;
            leaf_chain_lengths = length_from_leaf_node_index(leaf_node_indices) ;
            is_leaf_chain_stubby = (leaf_chain_lengths < length_threshold) ;
            if ~isscalar(leaf_node_indices) && all(is_leaf_chain_stubby) ,
                % don't want to delete them all---spare the longest one
                [~, i_max] = max(leaf_chain_lengths) ;
                is_leaf_chain_stubby(i_max) = false ;
            end
            % Delete the stubby leaf nodes
            stubby_leaf_node_indices = leaf_node_indices(is_leaf_chain_stubby) ;
            is_doomed_from_leaf_node_index(stubby_leaf_node_indices) = true ;            
        end
        
        % Mark all nodes of the doomed leaf node chains as doomed
        is_doomed_from_node_id = false(node_count, 1) ;
        for leaf_node_index = 1 : leaf_node_count ,
             if is_doomed_from_leaf_node_index(leaf_node_index) ,
                 node_ids_in_chain = node_ids_in_chain_from_leaf_node_index{leaf_node_index} ;
                 is_doomed_from_node_id(node_ids_in_chain) = true ;
             end
        end
        
        % Populate the outputs.  Basically the same as the inputs, but with
        % doomed nodes deleted.
        A_tree_output = A_tree_input ;
        A_tree_output(is_doomed_from_node_id,:) = [] ;
        A_tree_output(:,is_doomed_from_node_id) = [] ;
        xyz_output = xyz(~is_doomed_from_node_id,:) ;
        did_prune_some = any(is_doomed_from_node_id) ;
    end
end
