function [A_decimated, xyz_decimated, original_node_id_from_decimated_node_id] = decimate_tree(A_tree, xyz, desired_sampling_interval)
    % Decimates a tree represented as an undirected adjacency matrix
    % A_tree, with nodes at coords given by xyz, in um.  desired_sampling_interval
    % should be in um.
    %
    % Result is returned as an undirected tree graph, A_decimated, at roughly the desired
    % sampling interval.  The xyz coords of the nodes of A_decimated are in
    % the rows of xyz_decimated.  Each node of A_decimated is a node of
    % A_tree, and original_node_id_from_decimated_node_id allows the called
    % to map between them.
    
    % Special case for empty tree
    node_count = size(A_tree,1) ;
    if node_count==0 ,
        A_decimated = A_tree ;
        xyz_decimated = xyz ;
        original_node_id_from_decimated_node_id = [] ;
        return
    end
    
    % Decompose tree into chains
    node_ids_from_chain_id = chains_from_tree(A_tree) ;
    
    % Downsample each chain
    chain_count = length(node_ids_from_chain_id) ;
    decimated_node_ids_from_chain_id = cell(chain_count, 1) ;
    for chain_id = 1 : chain_count ,
        % Get the node ids for this chain
        node_ids_in_chain = node_ids_from_chain_id{chain_id} ;               
        decimated_node_ids_in_chain = decimate_chain(node_ids_in_chain, xyz, desired_sampling_interval) ;
        decimated_node_ids_from_chain_id{chain_id} = decimated_node_ids_in_chain ;
    end
    
    %
    % Now need to reassemble the decimated chains into a graph, and want to
    % get just the xyz coords used in the decimated graph
    %
    edges_from_chain_id = cell(chain_count, 1) ;
    for chain_id = 1 : chain_count ,
        decimated_node_ids = decimated_node_ids_from_chain_id{chain_id} ;
        to = decimated_node_ids(1:end-1) ;
        from = decimated_node_ids(2:end) ;
        edges_from_chain_id{chain_id} = [from(:) to(:)]' ;  % 2 x edge_count_for_this_chain
    end
    edges_using_original_node_ids = [edges_from_chain_id{:}]' ;  % edge_count x 2, all the edges, across all chains
    original_node_ids_in_edges = edges_using_original_node_ids(:) ;  % col vector
    
    edge_count = size(edges_using_original_node_ids, 1) ;  % edge count in the post-decimation graph
    [original_node_id_from_decimated_node_id,~,decimated_node_ids_in_edges] = unique(original_node_ids_in_edges) ;
    decimated_node_count = length(original_node_id_from_decimated_node_id) ;
    edges_using_decimated_node_ids = reshape(decimated_node_ids_in_edges, [edge_count 2]) ;

%     % Think this code would do the same thing as above, maybe is clearer    
%     original_node_id_from_decimated_node_id = unique(original_node_ids_in_edges) ;
%     max_original_node_id = original_node_id_from_decimated_node_id(end) ;  % works b/c unique sorts the result
%     decimate_node_id_from_original_node_id = invert_partial_map_array(original_node_id_from_decimated_node_id, max_original_node_id) ;
%     edges_using_decimated_node_ids = decimate_node_id_from_original_node_id(edges_using_original_node_ids) ;   
    
    dA_decimated = sparse(edges_using_decimated_node_ids(:,1), ...
                          edges_using_decimated_node_ids(:,2), ...
                          1, ...
                          decimated_node_count, ...
                          decimated_node_count) ;
    A_decimated = max(dA_decimated, dA_decimated') ;
    xyz_decimated = xyz(original_node_id_from_decimated_node_id,:) ;
end
